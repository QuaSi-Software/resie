using Optim
using BlackBoxOptim
using Base.Threads


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
function create_variant(io_settings::Dict{String,Any}, sim_params::Dict{String, Any}, 
                        project_config::OrderedDict{String,Any}, 
                        sample_params::Dict{String,Any})::OrderedDict{String,Any}
    cfg = deepcopy(project_config)

    # set up the parameters for this simulation variant
    for (key, value) in pairs(sample_params)
        category, uac, param_key = split(key, " ")
        cfg[category][uac][param_key] = value
        # rename outputs to clarify parameter values if outputfiles for each simulation 
        # should be generated
        #TODO do same for all type of output files that get generated for each simulation
        # and add Info that this is happening with hint to set the files to nothing of if 
        # not needed for performance boost
        if io_settings["csv_output"] != "nothing"
            name, ext = rsplit(io_settings["csv_output_file"], '.'; limit=2)
            cfg["io_settings"]["csv_output_file"] = name * "_" * uac * "_" * param_key * 
                                                    "_" * string(value) * "." * ext
        end
        if io_settings["output_plot"] != "nothing"
            name, ext = rsplit(io_settings["output_plot_file"], '.'; limit=2)
            cfg["io_settings"]["output_plot_file"] = name * "_" * uac * "_" * param_key * 
                                                     "_" * string(value) * "." * ext
        end
        if io_settings["sankey_plot"] != "nothing"
            name, ext = rsplit(io_settings["sankey_plot_file"], '.'; limit=2)
            cfg["io_settings"]["sankey_plot_file"] = name * "_" * uac * "_" * param_key * 
                                                     "_" * string(value) * "." * ext
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
    
    #     # create correct profiles from price_profile_paths and add the to sim_output.
    #     # profiles will be overwritten for every run to make sure the profile_scales and 
    #     # profile_addons calculated correctly.
    #     # If multiple threads are used each thread gets their own profile.
    #     #TODO change directiory to something better
    #     profile_dir = sim_params["run_path"]("./profiles/parallel_runs")
    #     mkpath(profile_dir)

    #     profile_id = ifelse(Threads.nthreads() > 1, Threads.threadid(), 0)
    #     date_range = remove_leap_days(collect(sim_params["start_date"]:Second(sim_params["time_step_seconds"]):sim_params["end_date"]))
    #     new_paths = Dict{String,String}() 

    #     for (name, path) in pairs(profile_paths)
    #         if profile_scales[name] != 1 && profile_addons[name] != 0 && profile_id == 0
    #             new_paths[name] = path
    #         else
    #             profile = Profile(path, sim_params)
    #             values = [profile.data[dt] .* profile_scales[name] .+ profile_addons[name] for dt in date_range]
    #             new_path = profile_dir * "/" * split(path[1:end-4], '/')[end] * "_$profile_id.prf" 
    #             save_to_prf(collect(date_range), values, new_path)
    #             new_paths[name] = path
    #         end
    #     end

    #     # replace the profile names with the paths to the new profiles
    #     function replace_profiles!(cfg::AbstractDict, replacements::Dict{String,String})
    #         for (k, v) in cfg
    #             if v isa String && haskey(replacements, v)
    #                 cfg[k] = replacements[v]
    #             elseif v isa AbstractDict
    #                 replace_profiles!(v, replacements)
    #             end
    #         end

    #     end

    #     replace_profiles!(cfg, new_paths)
    # end

    return cfg
end

"""
    optim_func!(all_results, io_settings, sim_params, optim_results_path, project_config, 
                        sample_values, run_lock, output_lock, results_lock)

Objective function called by algorithms. Wraps running of single simulations in a 
compatible format. Can be used as a batch function with Arrays for algorithms supporting it.

# Arguments
- `all_results::Array{Union{Float64,Nothing}}`: Results of all runs
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
# Returns
- `Union{Array{Float64},Float64}`: Objective of the simulation run or batch runs
"""
function optim_func!(all_results::Array{Union{Float64,Nothing}}, 
                     io_settings::Dict{String,Any}, sim_params::Dict{String,Any}, 
                     optim_results_path::String, project_config::OrderedDict{String,Any}, 
                     sample_values::Union{Array{Float64},Float64}, run_lock::ReentrantLock, 
                     output_lock::ReentrantLock, 
                     results_lock::ReentrantLock)::Union{Array{Float64},Float64}
    #TODO implement batch evaluation for Metaheuristics e.g. if size(sample_values) > 1 then @threads return Vector
    sample_params = Dict{String, Any}(zip(sim_params["optimisation"]["optim_params_keys"], sample_values))
    run_ID = uuid4()
    results = run_sample(io_settings, sim_params, optim_results_path, project_config, 
                         sample_params, run_ID, run_lock, output_lock)

    lock(results_lock) do 
        push!(all_results, results)
    end

    return results["objective"]
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
- `all_results::Array{Union{Float64,Nothing}}`: Results of all runs
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
"""
function monte_carlo_annealing!(all_results::Array{Union{Float64,Nothing}}, 
                                obj::Array{Union{Float64,Nothing}}, obj_lock::ReentrantLock, 
                                sim_params::Dict{String,Any}, optim_results_path::String, 
                                project_config::OrderedDict{String,Any}, idx::Int64, 
                                run_ID::UUID, run_lock::ReentrantLock, 
                                output_lock::ReentrantLock, results_lock::ReentrantLock)
    optimiser = sim_params["optimisation"]
    # temperature schedule is simple inverse logistic curve
    temperature = 1.0 - 1.0 / (1.0 + exp(-8.0 * (idx / length(optimiser["iterator"]) - 0.5)))

    if length(all_results) == 0 || rand() < temperature
        # set parameters to equally distributed random values across whole parameter space
        sample_params = Dict{String, Any}()
        for (key, param) in pairs(optimiser["optim_params"])
            #TODO maybe define optim params also as ranges but use minimum(range) and maximum(range) as limits
            sample_params[key] = rand(range(start=param["min"], stop=param["max"], length=100))
        end
    else
        # set parameters to neighborhood of existing result, drawn from the top results
        # by global measure, where temperature determines the results pool and size of
        # neighborhood
        sample_idx = rand(1:max(1, Int(round(length(all_results) * temperature))))
        # sample = sample_idx >= 1 && sample_idx <= length(all_results) ? all_results[sample_idx] : all_results[1]
        sample = all_results[sample_idx]
        
        sample_params = Dict{String, Any}()
        for (key, param) in pairs(optimiser["optim_params"])
            range = optimiser["nbh_scale"] * temperature * (param["max"] - param["min"])
            value = sample[key] + rand((-0.5 * range):(0.5 * range))
            sample_params[key] = clamp(value, param["min"], param["max"])
        end
    end

    # run sim and calculate objective results
    results = run_sample(io_settings, sim_params, optim_results_path, project_config,
                         sample_params, run_ID, run_lock, output_lock)

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
            res["gm"] = norm(res[k] / m - 1 for (k, m) in zip(parse_outkeys(optimiser["objective_params"]), obj))
        end
        sort!(all_results; by=x -> x["gm"])
    end
end