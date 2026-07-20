using Optim: Optim
using BlackBoxOptim: BlackBoxOptim
using Metaheuristics: Metaheuristics
using NLopt: NLopt
using NOMAD: NOMAD
using PolyChaos: PolyChaos
using LinearAlgebra
using Statistics
using Random

"""
    create_variant(io_settings, sim_params, project_config, sample_params)

Create a variant of the input file for the simulation run

# Arguments
- `io_settings::Dict{String,Any}`: IO settings
- `sim_params::Dict{String, Any}`: Simulation parameters
- `project_config::OrderedDict{String,Any}`: The project config read from the input file
- `sample_params::Dict{String,Any}`: The parameters that get changed for the simulation run
# Returns
- `OrderedDict{String,Any}`: Modified project_config
"""
function create_variant(io_settings::Dict{String,Any}, sim_params::Dict{String,Any},
                        project_config::OrderedDict{String,Any},
                        sample_params::Dict{String,Any})::OrderedDict{String,Any}
    # helper function to shorten parameter values for output file naming
    function compact_value(value)
        if abs(value) >= 100
            return string(round(Int, value))
        elseif abs(value) >= 10
            return string(round(value; digits=1))
        else
            return string(round(value; digits=2))
        end
    end

    # crete initial config for the current run
    cfg = deepcopy(project_config)

    # set up the parameters for this simulation variant
    run_name = ""
    for (key, value) in pairs(sample_params)
        uac, param_key = split(key, " ")
        if uac in keys(cfg["components"])
            cfg["components"][uac][param_key] = value
        else
            cfg[uac][param_key] = value
        end
        run_name *= uac * "_" * param_key * "_" * compact_value(value) * "_"
    end

    # rename outputs to clarify parameter values if outputfiles for each simulation should be generated
    for (key, value) in pairs(io_settings)
        if endswith(key, "file_path") && io_settings[key] != "nothing"
            dir, filename = splitdir(value)
            root, ext = splitext(filename)
            new_path = joinpath(dir, run_name * "_" * root * ext)
            cfg["io_settings"][key] = new_path
        end
    end

    #TODO postponed maybe this should be moved to profile processing to allow the profiles to be 
    # defined with "profiles" group, scale and addon without optimiser
    # if haskey(cfg, "profiles")
    #     profile_paths = Dict{String,String}()
    #     profile_scales = Dict{String,Float64}()
    #     profile_addons = Dict{String,Float64}()
    #     for (name, profile) in pairs(cfg["profiles"])
    #         profile_paths[name] = profile["path"]
    #         profile_scales[name] = profile["scale"]
    #         profile_addons[name] = profile["addon"]
    #     end
    # end

    return cfg
end

"""
    optim_func!(all_results, io_settings, sim_params, optim_results_path, project_config, 
                sample_values, run_lock, output_lock, results_lock)

Objective function called by algorithms. Wraps running of single simulations in a 
compatible format. Can be used as a batch function with Arrays for algorithms supporting it.

# Arguments
- `all_results::Vector{Any}`: Results of all runs
- `io_settings::Dict{String,Any}`: IO settings
- `sim_params::Dict{String,Any}`: Simulation parameters
- `optim_results_path::String`: Filepath for optim_results
- `project_config::OrderedDict{String,Any}`: The project config
- `sample_values::Union{Array{Float64},Float64}`: Values of the sample_params that get 
                                                  changed for the next simulation run or 
                                                  batch runs
- `run_lock::ReentrantLock`: Lock for writing to current_runs
- `output_lock::ReentrantLock`: Lock for file at optim_results_path
- `results_lock::ReentrantLock`: Lock for all_results
- `cancel_flag::Union{Nothing,Threads.Atomic{Bool}}`: Flag to pass STR+C down to all parallel runs
# Returns
- `Union{Array{Float64},Float64}`: Objective of the simulation run or batch runs
"""
function optim_func!(all_results::Vector{Any}, io_settings::Dict{String,Any},
                     sim_params::Dict{String,Any}, optim_results_path::String,
                     project_config::OrderedDict{String,Any},
                     sample_values::Union{Array{Float64},Float64}, run_lock::ReentrantLock,
                     output_lock::ReentrantLock,
                     results_lock::ReentrantLock;
                     cancel_flag::Union{Nothing,Threads.Atomic{Bool}}=nothing,
                     preparation_cache::Union{Nothing,PreparationCache}=nothing)::Union{Array{Float64},Float64}
    sample_params = Dict{String,Any}(zip(sim_params["optimisation"]["optim_params_keys"], sample_values))
    run_ID = uuid4()
    results = run_sample(io_settings, sim_params, optim_results_path, project_config,
                         sample_params, run_ID, run_lock, output_lock;
                         suppress_all_output=sim_params["optimisation"]["disable_all_simulation_outputs"],
                         cancel_flag=cancel_flag,
                         preparation_cache=preparation_cache)

    lock(results_lock) do
        push!(all_results, results)
    end

    return return objective_for_optimiser(results, sim_params["optimisation"])
end

function objective_for_optimiser(results::AbstractDict, optimiser::Dict{String,Any})
    objective = results["objective"]

    if optimiser["N_obj"] == 1
        return objective
    end

    # apply signs for each objective for multi-objective optimisation
    return Float64.(objective) .* optimiser["objective_signs"]
end

# transform normalized optim_params back to physical values
function to_physical_optim_values(sample_values, bounds::AbstractMatrix{<:Real})::Vector{Float64}
    values = sample_values isa Real ? [Float64(sample_values)] : vec(Float64.(sample_values))
    return bounds[:, 1] .+ values .* (bounds[:, 2] .- bounds[:, 1])
end

"""
    report_best_optimisation_result(all_results, optimiser)

Log the best successfully evaluated single-objective optimisation result.
Optimisation parameters are reported in their physical units.

# Arguments
- `all_results::Vector{Any}`: Results of all runs
- `optimiser::Dict{String,Any}`: The dict with the optimiser parameters
"""
function report_best_optimisation_result(all_results::Vector{Any}, optimiser::Dict{String,Any})
    if isempty(all_results)
        @globalInfo "No optimisation result is available."
        return
    end

    if optimiser["N_obj"] != 1
        # No single best solution can be reported for multi-objective optimisation.
        return
    end

    # Failed simulation runs use Inf as their objective and must not be selected.
    valid_results = filter(all_results) do result
        haskey(result, "objective") && result["objective"] isa Real && isfinite(result["objective"])
    end

    if isempty(valid_results)
        @globalInfo "No valid optimisation solution was found."
        return
    end

    objectives = Float64[result["objective"] for result in valid_results]
    best_objective, best_idx = findmin(objectives)
    best_result = valid_results[best_idx]

    format_value(value) = value isa Real ?
                          string(round(Float64(value); sigdigits=8)) :
                          string(value)

    parameter_lines = ["  $key = $(format_value(best_result[key]))" for key in optimiser["optim_params_keys"]]

    @globalInfo("Best optimisation result:\n" *
                "  Objective = $(format_value(best_objective))\n" *
                "  Physical parameter values:\n" *
                join(parameter_lines, "\n"))
end

"""
    monte_carlo_annealing!(all_results, obj, obj_lock, sim_params, optim_results_path, 
                           project_config, idx, run_ID, run_lock, output_lock, 
                           results_lock)

Combined monte carlo and simulated annealing algorithm. The temperature determines if a 
completely random or existing sample is used as starting point, determines the size of the 
neighborhood and the number of results (sorted by) global measure, from which a new sample 
is drawn.

# Arguments
- `all_results::Vector{Any}`: Results of all runs
- `io_settings::Dict{String,Any}`: IO settings used for simulation output and result writing.
- `obj::Array{Union{Float64,Nothing}}`: Objectives for optimisation
- `obj_lock::ReentrantLock`:: Lock for obj
- `sim_params::Dict{String,Any}`: Simulation parameters
- `optim_results_path::String`: Filepath for optim_results
- `project_config::OrderedDict{String,Any}`: The project config
- `idx::Int64`: The current run number 
- `run_ID::UUID`: The run ID used in the run registry
- `run_lock::ReentrantLock`: Lock for writing to current_runs
- `output_lock::ReentrantLock`: Lock for file at optim_results_path
- `results_lock::ReentrantLock`: Lock for all_results
- `cancel_flag::Union{Nothing,Threads.Atomic{Bool}}`: Flag to pass STR+C down to all parallel runs
"""
function monte_carlo_annealing!(all_results::Vector{Any}, io_settings::Dict{String,Any},
                                obj::Array{Union{Float64,Nothing}}, obj_lock::ReentrantLock,
                                sim_params::Dict{String,Any}, optim_results_path::String,
                                project_config::OrderedDict{String,Any}, idx::Int64,
                                run_ID::UUID, run_lock::ReentrantLock,
                                output_lock::ReentrantLock, results_lock::ReentrantLock;
                                cancel_flag::Union{Nothing,Threads.Atomic{Bool}}=nothing,
                                preparation_cache::Union{Nothing,PreparationCache}=nothing)
    optimiser = sim_params["optimisation"]
    # temperature schedule is simple inverse logistic curve
    temperature = 1.0 - 1.0 / (1.0 + exp(-8.0 * (idx / length(optimiser["iterator"]) - 0.5)))

    if length(all_results) == 0 || rand() < temperature
        # set parameters to equally distributed random values across whole parameter space
        sample_params = Dict{String,Any}()
        for (key, values) in zip(optimiser["optim_params_keys"], optimiser["optim_params_values"])
            sample_params[key] = rand(values)
        end
    else
        # set parameters to neighborhood of existing result, drawn from the top results
        # by global measure, where temperature determines the results pool and size of
        # neighborhood
        sample_idx = rand(1:max(1, Int(round(length(all_results) * temperature))))
        # sample = sample_idx >= 1 && sample_idx <= length(all_results) ? all_results[sample_idx] : all_results[1]
        sample = all_results[sample_idx]

        sample_params = Dict{String,Any}()
        for (key, values) in zip(optimiser["optim_params_keys"], optimiser["optim_params_values"])
            range = optimiser["nbh_scale"] * temperature * (maximum(values) - minimum(values))
            value = sample[key] + rand((-0.5 * range):(0.5 * range))
            sample_params[key] = clamp(value, minimum(values), maximum(values))
        end
    end

    # run sim and calculate objective results
    results = run_sample(io_settings, sim_params, optim_results_path, project_config,
                         sample_params, run_ID, run_lock, output_lock;
                         suppress_all_output=optimiser["disable_all_simulation_outputs"],
                         cancel_flag=cancel_flag,
                         preparation_cache=preparation_cache)

    # calculate minimum of results
    if any(!isnothing(obj))
        @lock obj_lock obj = results["objective"]
    else
        @lock obj_lock obj .= min.(obj, results["objective"])
    end

    # write output to all_results
    lock(results_lock) do
        push!(all_results, results)

        # calculate global measure and sort by it
        for res in all_results
            res["gm"] = norm(res[k] / m - 1 for (k, m) in zip(optimiser["objective_params_keys"], obj))
        end
        sort!(all_results; by=x -> x["gm"])
    end
end

"""
    calc_global_sensitivity!(model_function, bounds, all_results, sim_params)

Calculate the global sensitivity indices with polynomial chaos expansion (PCE). A 
surrogate 3rd degree polynomial model is fit to the existing data. If the existing data is 
doesn't produce a well enough fit more data is generated in batches until RMSE is < 0.1 or 
2x the max_runs is hit.

# Arguments
- `model_function::Function`: Function to run if more datapoints are needed
- `bounds::Array{Float64}`: Bounds in which to analyse parameters
- `all_results::Vector{Any}`: Results of all runs
- `optim_params_keys::Vector{String}`: Names of all variable parameters
- `sim_params::Dict{String,Any}`: Simulation parameters
# Returns
- `Float64`: Total-order Sobol sensitivity index
- `Float64`: First-order Sobol sensitivity index
- `Float64`: Relative root mean square error for the surrogate model
- `Float64`: R^2 for the surrogate model
"""
function calc_global_sensitivity!(model_function::Union{Nothing,Function},
                                  bounds::Array{Float64},
                                  all_results::Vector{Any},
                                  optim_params_keys::Array{String},
                                  sim_params::Dict{String,Any},
                                  cancel_optimisation::Threads.Atomic{Bool})::Tuple{Vector{Float64},Vector{Float64},
                                                                                    Float64,Float64}
    d = size(bounds, 1)
    deg = 3
    op = PolyChaos.Uniform01OrthoPoly(deg; Nrec=5 * deg)
    mop = PolyChaos.MultiOrthoPoly(fill(op, d), deg)

    # [-1,1] -> physical
    function to_phys(x, lo, hi)
        lo + (x + 1) * (hi - lo) / 2
    end

    # physical -> [-1,1]
    function to_std(x, lo, hi)
        (x - lo) / (hi - lo) * 2 - 1
    end

    # Fit PCE coefficients by least squares and estimate out-of-sample error via
    # leave-one-out CV using the hat-matrix shortcut:
    function fit_surrogate(X_std::Matrix{Float64}, y::Vector{Float64}, mop::PolyChaos.MultiOrthoPoly)
        Phi = zeros(length(y), mop.dim)
        for i in 1:length(y)
            Phi[i, :] = PolyChaos.evaluate(X_std[i, :], mop)
        end
        coeffs = Phi \ y
        resid = y .- Phi * coeffs
        H = Phi * pinv(Phi' * Phi) * Phi'     # hat matrix, only needed for its diagonal
        leverage = diag(H)
        loocv_resid = resid ./ (1 .- leverage)
        rel_rmse = sqrt(mean(loocv_resid .^ 2)) / std(y)
        r2 = 1 - sum(loocv_resid .^ 2) / sum((y .- mean(y)) .^ 2)

        return coeffs, rel_rmse, r2
    end

    # Take samples (standardized on [-1,1])
    n_existing = size(all_results, 1)
    if n_existing > 0
        keys = vcat("objective", sim_params["optimisation"]["optim_params_keys"]...)
        res_matrix = [d[k] for d in all_results, k in keys]
        X_phys = res_matrix[:, 2:end]
        y = res_matrix[:, 1]
        X_std = hcat([to_std.(X_phys[:, i], bounds[i, 1], bounds[i, 2]) for i in 1:d]...)
        coeffs, rel_rmse, r2 = fit_surrogate(X_std, y, mop)
    else
        rel_rmse = 1.0
        coeffs = 0.0
        r2 = 0.0
        y = []
        X_phys = Array{Float64}(undef, 0, d)
    end

    target_rel_rmse = 0.1
    target_r2 = 0.9

    if model_function !== nothing
        print_message = true
        maximum_sensitivity_runs = sim_params["optimisation"]["max_runs"] * 2
        status_lock = ReentrantLock()
        try
            while (rel_rmse > target_rel_rmse || r2 < target_r2) && length(y) < maximum_sensitivity_runs
                if print_message
                    @globalInfo("Performing additional runs for sensitivity analysis.\n" *
                                "Press Ctrl+C any time to stop creating additional runs and to " *
                                "calculate the sensitivity from the already completed runs.")
                    print_message = false
                end

                current_sample_count = length(y)
                if current_sample_count < mop.dim
                    n_new = max(mop.dim + 1 - current_sample_count, Threads.nthreads())
                else
                    n_new = Threads.nthreads()
                end

                # Do not exceed the configured sensitivity-run limit.
                remaining_runs = maximum_sensitivity_runs - current_sample_count
                n_new = min(n_new, remaining_runs)
                n_new > 0 || break

                X_std_new = rand(n_new, d) .* 2 .- 1
                X_phys_new = hcat([to_phys.(X_std_new[:, i], bounds[i, 1], bounds[i, 2]) for i in 1:d]...)

                y_new = zeros(Float64, n_new)
                completed_in_batch = Ref(0)
                Threads.@threads for i in 1:n_new
                    y_new[i] = model_function(X_phys_new[i, :])

                    lock(status_lock) do
                        completed_in_batch[] += 1
                        total_completed = current_sample_count + completed_in_batch[]
                        @globalInfo "Sensitivity status: $total_completed of up to $maximum_sensitivity_runs runs completed."
                    end
                end

                # Create candidate arrays first. They are committed only
                # after the new surrogate fit has completed successfully.
                X_phys_candidate = vcat(X_phys, X_phys_new)
                y_candidate = vcat(y, y_new)
                X_std_candidate = hcat([to_std.(X_phys_candidate[:, i], bounds[i, 1], bounds[i, 2]) for i in 1:d]...)

                coeffs_candidate, rel_rmse_candidate, r2_candidate = fit_surrogate(X_std_candidate, y_candidate, mop)

                # Commit only a completely simulated and fitted batch.
                X_phys = X_phys_candidate
                y = y_candidate
                coeffs = coeffs_candidate
                rel_rmse = rel_rmse_candidate
                r2 = r2_candidate

                @globalInfo("Current quality: relative RMSE=$(round(rel_rmse; digits=4)) (goal <= $target_rel_rmse), " *
                            "R²=$(round(r2; digits=4)) (goal >= $target_r2).")
            end
        catch e
            if e isa InterruptException
                cancel_optimisation[] = true
                @globalInfo("Additional sensitivity runs interrupted by Ctrl+C. " *
                            "Continuing sensitivity calculation with $(length(y)) completed runs.")
            else
                rethrow()
            end
        end
    end

    # Calculate Sobol indices from coefficients 
    total_var = sum(coeffs[2:end] .^ 2)

    S_first = zeros(d)
    S_total = zeros(d)
    for i in 1:d, k in 2:mop.dim
        degs = mop.ind[k, :]
        if degs[i] != 0
            S_total[i] += coeffs[k]^2
            if all(j == i || degs[j] == 0 for j in 1:d)
                S_first[i] += coeffs[k]^2
            end
        end
    end
    S_first ./= total_var
    S_total ./= total_var

    width = length.(optim_params_keys)
    @globalInfo "Global sensitivity results: \n" *
                "\t $(join(optim_params_keys, "\t")) \n" *
                "S_total\t $(join(rpad.(round.(S_total, digits=3), width), "\t")) \n" *
                "S_first\t $(join(rpad.(round.(S_first, digits=3), width), "\t")) \n" *
                "Surrogate RMSE: $(round(rel_rmse, digits=3)), R2: $(round(r2, digits=3)) \n" *
                "Important: Sobol indices depend on the selected parameter bounds!"

    return S_total, S_first, rel_rmse, r2
end

"""
    perform_optimisation(io_settings, sim_params, optim_results_path, project_config, all_results)

Run the configured optimisation or parameter variation and collect the resulting simulation
outputs.

This function initializes the locks required for parallel execution, clears the optimisation
results file, selects the optimisation backend from `sim_params["optimisation"]["type"]`,
and dispatches the corresponding optimisation workflow.

# Arguments
- `io_settings::Dict{String,Any}`: IO settings used for simulation output and result writing.
- `sim_params::Dict{String,Any}`: Simulation and optimisation parameters. 
- `project_config::OrderedDict{String,Any}`: Base project configuration used to generate
  simulation variants.

# Returns
- `Bool`: Flag is simulation was successful (true) or not (false)
- `Vector{Any}`: The updated `all_results` collection containing the results of all
  completed optimisation or parameter-variation runs.
"""
function perform_optimisation(io_settings::Dict{String,Any},
                              sim_params::Dict{String,Any},
                              project_config::OrderedDict{String,Any};
                              preparation_cache::Union{Nothing,PreparationCache}=nothing)::Tuple{Bool,Vector{Any}}
    # establish overarching locks for parallelization
    run_lock = ReentrantLock()
    output_lock = ReentrantLock()
    results_lock = ReentrantLock()

    optim_results_path = sim_params["run_path"](io_settings["optimisation_csv_file_path"])
    open(optim_results_path, "w") do f
    end

    optimiser = sim_params["optimisation"]

    # calculate inputs that should be cashed
    if preparation_cache !== nothing
        @globalInfo "Preparing reusable data for repeated simulation runs..."
        warmup_run_ID = uuid4()
        prepare_inputs(project_config, warmup_run_ID; preparation_cache=preparation_cache)
    end

    # warn if file outputs for multi-thread simulations should be created, as this can lead to troubles
    uses_threaded_sample_evaluation = length(optimiser["iterator"]) > 1 ||
                                      (optimiser["type"] == "Metaheuristics" && Threads.nthreads() > 1)
    if uses_threaded_sample_evaluation && !optimiser["disable_all_simulation_outputs"]
        @warn "Writing simulation outputs during multi-threaded optimisation may cause file-access conflicts or crashes. " *
              "To be safe, set `disable_all_simulation_outputs` to `true`, or run the optimisation with a single thread."
    end

    # handle interruption via STR+C for parallel runs and optimisation
    cancel_optimisation = Threads.Atomic{Bool}(false)

    # prepare result vector
    all_results = Vector{Any}()

    @globalInfo "Starting Simulations on $(Threads.nthreads()) Threads"
    main_start_time = now()

    if length(optimiser["iterator"]) > 1
        if optimiser["type"] == "monte_carlo_annealing"
            obj = Array{Union{Float64,Nothing}}(nothing)
            obj_lock = ReentrantLock()
        end
        try
            # handling status
            variation_start_time = now()
            completed_runs = Atomic{Int}(0)
            max_runs = length(optimiser["iterator"])
            worker_count = min(Threads.nthreads(), max_runs)

            @threads for sample_values in collect(optimiser["iterator"])
                if cancel_optimisation[]
                    continue
                end

                try
                    run_start_time = now()

                    # decide which algorithm to run based on type of optimiser
                    if optimiser["type"] == "parametervariation"
                        float_sample_values = sample_values isa Real ? Float64(sample_values) :
                                              Float64.(collect(sample_values))

                        optim_func!(all_results,
                                    io_settings,
                                    sim_params,
                                    optim_results_path,
                                    project_config,
                                    float_sample_values,
                                    run_lock,
                                    output_lock,
                                    results_lock;
                                    cancel_flag=cancel_optimisation,
                                    preparation_cache=preparation_cache)

                    elseif optimiser["type"] == "monte_carlo_annealing"
                        # TODO this is not working currently...
                        monte_carlo_annealing!(all_results,
                                               io_settings,
                                               obj,
                                               obj_lock,
                                               sim_params,
                                               optim_results_path,
                                               project_config,
                                               sample_values,
                                               sample_ID,
                                               run_lock,
                                               output_lock,
                                               results_lock;
                                               cancel_flag=cancel_optimisation,
                                               preparation_cache=preparation_cache)
                    end

                    runtime_seconds = round(Int, seconds(now() - run_start_time))
                    runtime_minutes, runtime_remaining_seconds = divrem(runtime_seconds, 60)
                    completed = atomic_add!(completed_runs, 1) + 1
                    elapsed_seconds = max(1, round(Int, seconds(now() - variation_start_time)))
                    if completed < worker_count
                        eta_text = "calculating..."
                    else
                        results_per_second = completed / elapsed_seconds
                        eta_seconds = round(Int, (max_runs - completed) / results_per_second)
                        eta_minutes, eta_remaining_seconds = divrem(max(eta_seconds, 0), 60)
                        eta_text = "$eta_minutes min $(lpad(eta_remaining_seconds, 2, '0')) s"
                    end

                    @globalInfo "[$completed/$max_runs] → completed in $runtime_minutes min " *
                                "$(lpad(runtime_remaining_seconds, 2, '0')) s. ETA: $eta_text"
                catch e
                    if e isa InterruptException
                        cancel_optimisation[] = true
                        continue
                    else
                        rethrow()
                    end
                end
            end
        catch e
            if e isa InterruptException
                # Handles Ctrl+C delivered to the task coordinating @threads.
                cancel_optimisation[] = true
            else
                rethrow()
            end
        end
    else
        # generic optimisation function
        # f_physical takes and returns the physical correct simulation parameter and results
        f_physical = function (sample_values)
            if cancel_optimisation[]
                throw(InterruptException())
            end

            return optim_func!(all_results, io_settings, sim_params, optim_results_path,
                               project_config, sample_values,
                               run_lock, output_lock, results_lock;
                               cancel_flag=cancel_optimisation,
                               preparation_cache=preparation_cache)
        end

        # f takes the normalised sample_values and converts them to physical values for the simulation run
        f = function (sample_values)
            physical_values = to_physical_optim_values(sample_values, optimiser["bounds"])
            return f_physical(physical_values)
        end

        # handle Logging for all algorithms
        progress_lock = ReentrantLock()
        progress_evaluations = Ref(0)
        progress_best = Ref(Inf)
        progress_every = get(optimiser, "progress_every", 1)

        # handles the progress logging of f
        f_progress = function (sample_values)
            cancel_optimisation[] && throw(InterruptException())
            result = f(sample_values)
            cancel_optimisation[] && throw(InterruptException())

            lock(progress_lock) do
                progress_evaluations[] += 1

                best_text = if optimiser["N_obj"] == 1
                    value = result isa Real ? Float64(result) : Float64(first(result))
                    progress_best[] = min(progress_best[], value)
                    string(round(progress_best[]; sigdigits=8))
                else
                    "multi-objective"
                end

                if progress_evaluations[] == 1 ||
                   progress_evaluations[] % progress_every == 0
                    @globalInfo("Optimisation progress: evaluations=$(progress_evaluations[]), best=$best_text",)
                end
            end

            return result
        end

        try
            if optimiser["type"] == "Optim"
                Optim.optimize(f_progress, optimiser["args"]...)
            elseif optimiser["type"] == "BlackBoxOptim"
                if optimiser["N_obj"] == 1
                    f_wrap = f_progress
                else
                    f_wrap(x) = Tuple(f_progress(x))
                end
                BlackBoxOptim.bboptimize(f_wrap, optimiser["args"]...; optimiser["kwargs"]...)
            elseif optimiser["type"] == "Metaheuristics"
                #TODO implement batch evaluation for other packages that need it        
                if Threads.nthreads() > 1
                    if optimiser["N_obj"] == 1
                        f_arr(x) = [f_progress(x)]
                    else
                        f_arr = f_progress
                    end
                    f_wrap = function (sample_values)
                        N_samples = size(sample_values, 1)
                        if N_samples > 1
                            objectives = zeros(N_samples, optimiser["N_obj"])
                            @threads for i in 1:N_samples
                                if cancel_optimisation[]
                                    continue
                                end

                                try
                                    objectives[i, :] = f_arr(sample_values[i, :])
                                catch e
                                    if e isa InterruptException
                                        cancel_optimisation[] = true
                                    else
                                        rethrow()
                                    end
                                end
                            end
                            if cancel_optimisation[]
                                throw(InterruptException())
                            end
                        else
                            objectives = f_arr(sample_values)
                        end
                        return objectives, zeros(N_samples, 1), zeros(N_samples, 1)
                    end
                else
                    f_wrap = f_progress
                end

                res = Metaheuristics.optimize(f_wrap, optimiser["args"]...)
                @globalInfo "Metaheuristics optimisation result:\n$res"

            elseif optimiser["type"] == "NLopt"
                f_nlopt = function (sample_values, gradient)
                    return f_progress(sample_values)
                end
                NLopt.min_objective!(optimiser["args"][1], f_nlopt)
                res = NLopt.optimize(optimiser["args"]...)
                @globalInfo "NLopt optimisation results: $res"
            elseif optimiser["type"] == "NOMAD"
                f_nomad = function (sample_values)
                    res = f_progress(sample_values)
                    success = res == Inf ? false : true
                    if length(res) == 1
                        res = [res]
                    end
                    return success, true, res
                end
                prob = NOMAD.NomadProblem(optimiser["args"][1:(end - 1)]..., f_nomad; optimiser["kwargs"]...)
                NOMAD.solve(prob, optimiser["args"][end])
            end
        catch e
            if e isa InterruptException
                cancel_optimisation[] = true
            else
                rethrow()
            end
        end
    end

    # handle info messages
    workflow_name = if optimiser["type"] == "parametervariation"
        "Parameter variation"
    elseif optimiser["type"] == "monte_carlo_annealing"
        "Monte Carlo annealing"
    else
        "Optimisation"
    end

    main_runtime_seconds = round(Int, seconds(now() - main_start_time))
    main_runtime_minutes, main_runtime_remaining_seconds = divrem(main_runtime_seconds, 60)
    main_run_count = length(all_results)

    if cancel_optimisation[]
        # use snapshot to avoid interference from threads that may still be running
        final_results = lock(results_lock) do
            copy(all_results)
        end
        main_run_count = length(final_results)
        @globalInfo "$workflow_name interrupted by user after $main_runtime_minutes min " *
                    "$(lpad(main_runtime_remaining_seconds, 2, '0')) s. $main_run_count runs completed. " *
                    "Recovering intermediate results..."

        # detect and output best simulation result, independent of any package-specific results
        report_best_optimisation_result(final_results, optimiser)
        return false, final_results
    end

    @globalInfo "$workflow_name completed in $main_runtime_minutes min " *
                "$(lpad(main_runtime_remaining_seconds, 2, '0')) s. $main_run_count runs completed."

    # detect and output best simulation result, independent of any package-specific results
    report_best_optimisation_result(all_results, optimiser)

    # Start sensitivity analysis based on the former and additional simulation runs.
    if optimiser["run_sensitivity"] && get(optimiser, "objective_function_name", "") == "multi-objective"
        @warn "Sensitivity analysis will not be performed because a multi-objective optimiser was selected."
        optimiser["run_sensitivity"] = false
    end

    if optimiser["run_sensitivity"]
        sensitivity_start_time = now()
        results_before_sensitivity = length(all_results)

        try
            if length(optimiser["iterator"]) > 1
                calc_global_sensitivity!(nothing, optimiser["bounds"][:, 1:2], all_results,
                                         optimiser["optim_params_keys"], sim_params, cancel_optimisation)

            else
                calc_global_sensitivity!(f_physical, optimiser["bounds"][:, 1:2], all_results,
                                         optimiser["optim_params_keys"], sim_params, cancel_optimisation)
            end
        catch e
            if e isa InterruptException
                cancel_optimisation[] = true

                sensitivity_runtime_seconds = round(Int, seconds(now() - sensitivity_start_time))
                sensitivity_runtime_minutes, sensitivity_runtime_remaining_seconds = divrem(sensitivity_runtime_seconds,
                                                                                            60)
                sensitivity_runs = length(all_results) - results_before_sensitivity

                @globalInfo "Global sensitivity analysis interrupted after $sensitivity_runtime_minutes min " *
                            "$(lpad(sensitivity_runtime_remaining_seconds, 2, '0')) s. $sensitivity_runs additional runs completed."

                return false, all_results
            else
                rethrow()
            end
        end

        sensitivity_runtime_seconds = round(Int, seconds(now() - sensitivity_start_time))
        sensitivity_runtime_minutes, sensitivity_runtime_remaining_seconds = divrem(sensitivity_runtime_seconds, 60)
        sensitivity_runs = length(all_results) - results_before_sensitivity
        @globalInfo "Global sensitivity analysis completed in $sensitivity_runtime_minutes min " *
                    "$(lpad(sensitivity_runtime_remaining_seconds, 2, '0')) s. $sensitivity_runs additional runs completed."

        overall_runtime_seconds = round(Int, seconds(now() - main_start_time))
        overall_runtime_minutes, overall_runtime_remaining_seconds = divrem(overall_runtime_seconds, 60)
        @globalInfo "Complete $workflow_name workflow finished in $overall_runtime_minutes min " *
                    "$(lpad(overall_runtime_remaining_seconds, 2, '0')) s. $(length(all_results)) total runs completed."
    end

    # write results to optimisation result file if not written continuously
    if !io_settings["write_optimisation_csv_continuously"] && !isempty(all_results)
        open(optim_results_path, "w") do file_handle
            # write header
            header = join(collect(keys(all_results[1])), ';') * "\n"
            write(file_handle, header)

            # write results
            for results in all_results
                row = join(collect(values(results)), ';') * "\n"
                row = replace(row, '.' => ',')
                write(file_handle, row)
            end
        end
    end

    return true, all_results
end
