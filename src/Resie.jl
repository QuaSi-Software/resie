module Resie

using Printf
using Dates: now, seconds
using UUIDs
using LinearAlgebra
using Base.Threads
using Optim
using BlackBoxOptim

"""
Contains the parameters, instantiated components and the order of operations for a simulation run.

In short, it bundles all the data required to perform a simulation from start to finish.
Through its fields it contains a lot of complexity in terms of hierarchical data structures.
However the struct should not be used as argument for complex calculations. Instead it is
intended to be used for the run registry only.

Because the types used to represent the field are defined in modules that depend upon the
definition of the struct in the first place, the fields are defined with the generic 'Any'
type.
"""
mutable struct SimulationRun
    parameters::Dict{String,Any}
    io_settings::Dict{String,Any}
    components::Dict{String,Any}
    order_of_operations::Vector{Any}
end

# this registry should be the only global state in the package and contains the state for
# ongoing or paused simulation runs
current_runs::Dict{UUID,SimulationRun} = Dict{UUID,SimulationRun}()

"""
    get_run(id)

Get the simulation run container for the given ID.

# Args
- `id::UUID`: The ID of the run
# Returns
- `SimulationRun`: The simulation run container
"""
function get_run(id::UUID)::SimulationRun
    return current_runs[id]
end

"""
    close_run(id)

Closes the given run, removing it from the run registry.

# Args
- `id:UUID`: The ID of the run
"""
function close_run(id::UUID)
    delete!(current_runs, id)
end

"""
Custom exception `InputError` used to signify that an input was not correctly set up,
outside the allowed range, etc.
Call with `throw(InputError("msg"))` or `throw(InputError())`.
"""
struct InputError <: Exception
    msg::Union{AbstractString,Nothing}
end
InputError() = InputError(nothing)

# note: includes that contain their own module, which have to be submodules of the Resie
# module, are included first, then can be accessed with the "using" keyword. files that
# contain code that is intended to be used in-place of their include statement (as part
# of the Resie module), are included after the "using" statements have been declared.
# this is done so the latter files can access the symbols of the submodules the same as
# if the code was inside this file.

include("profiles/solar_irradiance.jl")
using .SolarIrradiance
include("profiles/base.jl")
using .Profiles
include("profiles/weatherdata.jl")
using .Weatherdata

include("energy_systems/base.jl")
using .EnergySystems

include("project_loading.jl")
include("file_output.jl")
include("economy.jl")
include("emissions.jl")

include("resie_logger.jl")
using .Resie_Logger

using PlotlyJS
using ColorSchemes
using Colors
using Interpolations
using JSON
using Dates

"""
    run_simulation_loop(sim_params, io_settings, components, operations)

Performs the simulation as loop over time steps and records outputs.

# Arguments
-`sim_params::Dict{String,Any}`: Simulation parameters
-`io_settings::Dict{String,Any}`: IO settings
-`components::Grouping`: The energy system components
-`operations::OrderOfOperations`: Order of operations
"""
function run_simulation_loop(sim_params::Dict{String,Any},
                             io_settings::Dict{String,Any},
                             components::Grouping,
                             operations::OrderOfOperations)
    # get list of requested output keys for lineplot and csv export
    output_keys_lineplot,
    output_keys_to_CSV,
    output_keys_economic_emissions,
    output_keys_optimise = get_output_keys(io_settings,
                                           sim_params["economic_parameters"],
                                           sim_params["emissions_parameters"],
                                           sim_params["optimisation"],
                                           components)
    all_requested_output_keys = Vector{Resie.EnergySystems.OutputKey}(unique(vcat(something(output_keys_lineplot,
                                                                                            String[]),
                                                                                  something(output_keys_economic_emissions,
                                                                                            String[]),
                                                                                  something(output_keys_optimise,
                                                                                            String[]))))
    weather_data_keys = get_weather_data_keys(sim_params)
    do_create_plot_data = output_keys_lineplot !== nothing
    do_create_plot_weather = weather_data_keys !== nothing && io_settings["plot_weather_data"]
    do_write_CSV_weather = weather_data_keys !== nothing && io_settings["csv_output_weather"]
    weather_CSV_keys = do_write_CSV_weather ? weather_data_keys : nothing
    do_write_CSV = output_keys_to_CSV !== nothing || do_write_CSV_weather
    do_write_CSV_continuously = io_settings["write_csv_continuously"]
    do_write_summary_CSV = io_settings["write_summary_CSV"]
    csv_file_path = io_settings["csv_output_file"]
    csv_time_unit = io_settings["csv_time_unit"]
    do_calculate_economy = sim_params["economic_parameters"]["calculate_economy"]
    do_calculate_emissions = sim_params["emissions_parameters"]["calculate_emissions"]
    do_optimise = sim_params["optimisation"]["run_optimisation"]

    # Initialize the arrays for output
    output_weather_lineplot = do_create_plot_weather ?
                              zeros(Float64, sim_params["number_of_time_steps_output"], 1 + length(weather_data_keys)) :
                              nothing
    output_csv = do_write_CSV && !do_write_CSV_continuously ?
                 Matrix{String}(undef, sim_params["number_of_time_steps_output"],
                                1 + length(output_keys_to_CSV) + (do_write_CSV_weather ? length(weather_data_keys) : 0)) :
                 nothing
    output_data_all_requested = do_create_plot_data || do_calculate_economy || do_calculate_emissions ?
                                zeros(Float64, sim_params["number_of_time_steps_output"],
                                      1 + length(all_requested_output_keys)) :
                                nothing

    # write CSV file headers
    if do_write_CSV
        header = get_output_header(output_keys_to_CSV, weather_CSV_keys, csv_time_unit)
        # Reset the output file and add headers for the given outputs.
        open(sim_params["run_path"](csv_file_path), "w") do file_handle
            write(file_handle, join(header, ';') * "\n")
        end
    end

    # check if sankey should be plotted
    do_create_sankey = io_settings["sankey_plot"] !== "nothing"
    if do_create_sankey
        # get information about all interfaces for Sankey
        nr_of_interfaces,
        medium_of_interfaces,
        output_sourcenames_sankey,
        output_targetnames_sankey = get_interface_information(components)
        # preallocate for speed: Matrix with data of interfaces in every timestep
        output_interface_values = zeros(Float64, sim_params["number_of_time_steps_output"], nr_of_interfaces)
    end

    # export order of operation and other additional info like optional plots
    dump_auxiliary_outputs(io_settings, components, operations, sim_params)

    @info "-- Start time step loop"
    if sim_params["start_date_output"] == sim_params["start_date"]
        @info "-- No preheating activated. Starting output from beginning."
    else
        @info "-- Starting with preheating, no output will be written until completed."
    end
    start = now()
    for steps in 1:sim_params["number_of_time_steps"]
        # check if data should be output
        do_output = sim_params["current_date"] >= sim_params["start_date_output"]
        output_steps = Int(max(1,
                               Int(steps) - (Int(sim_params["number_of_time_steps"]) -
                                             Int(sim_params["number_of_time_steps_output"]))))
        if output_steps == 1 && do_output && sim_params["start_date_output"] != sim_params["start_date"]
            @info "-- Preheating completed. Starting output now."
        end

        operations_adjusted = reorder_operations_in_time_step(components, operations, sim_params)
        perform_operations(components, operations_adjusted, sim_params)

        if do_output
            # check if any component and/or interface was not balanced
            interface_warnings = check_balances_of_interfaces(components, sim_params["epsilon"])
            component_warnings = check_balances_of_components(components, sim_params["epsilon"])
            if length(interface_warnings) > 0
                for (key, balance) in interface_warnings
                    @balanceWarn "In timestep $(sim_params["current_date"]), the balance in interface " *
                                 "$key was not zero: $balance"
                end
            end
            if length(component_warnings) > 0
                for (key, balance) in component_warnings
                    @error "In timestep $(sim_params["current_date"]), the balance for component " *
                           "$key was not zero: $balance. This is probably caused by a bug in the component."
                end
            end

            # write requested output data to the CSV file if configured, or to output
            # storage if not
            if do_write_CSV
                row = get_output_row(output_keys_to_CSV, 
                                     weather_CSV_keys, 
                                     sim_params, 
                                     csv_time_unit)
                row[2:end] = replace.(row[2:end], '.' => ',')
                if do_write_CSV_continuously
                    # Write row to the output file
                    open(sim_params["run_path"](csv_file_path), "a") do file_handle
                        write(file_handle, join(row, ';') * "\n")
                    end
                else
                    output_csv[output_steps, :] = row
                end
            end

            # get the energy transported through each interface in every timestep for Sankey
            if do_create_sankey
                output_interface_values[output_steps, :] = collect_interface_energies(components, nr_of_interfaces)
            end
            # gather output data of weather
            if do_create_plot_weather
                output_weather_lineplot[output_steps, :] = gather_weather_data(weather_data_keys, sim_params)
            end
            # gather output data of each component for line plot, economy and emissions
            if do_create_plot_data || do_calculate_economy || do_calculate_emissions || do_optimise
                output_data_all_requested[output_steps, :] = gather_output_data(all_requested_output_keys,
                                                                                sim_params["time_since_output"])
            end

            # simulation update
            sim_params["time_since_output"] += Int(sim_params["time_step_seconds"])
        end

        # simulation update
        sim_params["time"] += Int(sim_params["time_step_seconds"])
        sim_params["current_date"] = add_ignoring_leap_days(sim_params["current_date"],
                                                            Second(sim_params["time_step_seconds"]))
        # progress report
        if sim_params["step_info_interval"] > 0 && steps % sim_params["step_info_interval"] == 0
            eta = ((sim_params["number_of_time_steps"] - steps)
                   *
                   seconds(now() - start)
                   /
                   sim_params["step_info_interval"])
            @info "Progress: $(steps)/$(sim_params["number_of_time_steps"]) in" *
                  " $(seconds(now() -  start)) s. ETA: $(@sprintf("%.4f", eta)) s"
            start = now()
        end
    end
    @info "-- Finished time step loop"

    # extract output data per category including the time step in the first column
    # TODO replace with parse_outkeys() or use this approach in do_optimise
    output_key_signature(k::EnergySystems.OutputKey) = (String(k.unit.uac), k.medium, k.value_key)
    function subset_cols(key_indexes, keys::Union{Nothing,Vector{EnergySystems.OutputKey}})
        keys === nothing && return Int[]
        return 1 .+ [key_indexes[output_key_signature(k)] for k in keys]
    end
    time_and_subset_cols(key_indexes, keys) = vcat(1, subset_cols(key_indexes, keys))
    key_indexes = Dict(output_key_signature(k) => i for (i, k) in pairs(all_requested_output_keys))

    if do_calculate_economy || do_calculate_emissions
        output_data_economic_emissions = Matrix(view(output_data_all_requested, :,
                                                     time_and_subset_cols(key_indexes, output_keys_economic_emissions)))
        economic_emissions_data = prepare_economic_emissions_data(components, output_keys_economic_emissions,
                                                                  output_data_economic_emissions)
        economic_result = do_calculate_economy ? calculate_economy(economic_emissions_data, sim_params) : nothing
        emissions_result = do_calculate_emissions ? calculate_emissions(economic_emissions_data, sim_params) : nothing
    end

    if do_optimise
        optim_results = Dict{String,Any}()
        # create keys consistent with csv_output for return data for return Dict
        output_data_header = get_output_header(all_requested_output_keys, nothing, csv_time_unit)
        output_data = OrderedDict{String, AbstractArray}(zip(output_data_header, eachcol(output_data_all_requested)))

        function write_optim_results!(params::Dict{String,Any}, 
                                      output_data::OrderedDict{String, AbstractArray}, 
                                      res::Dict{String,Any})::Dict{String,Any}
            for (func, spec) in pairs(params)
                if func == "sum"
                    keys = parse_outkeys(spec)
                    for key in keys
                        res["sum $key"] = sum(output_data[key])
                    end
                elseif func == "mean"
                    keys = parse_outkeys(spec)
                    for key in keys                    
                        res["mean $key"] = sum(output_data[key]) / length(output_data[key])
                    end
                elseif func == "economic"
                    for key in spec
                        res["$func $key"] = getfield(economic_result, Symbol(key))
                    end
                elseif func == "emissions"
                    for key in spec
                        res["$func $key"] = getfield(emissions_result, Symbol(key))
                    end
                end
            end
            return res
        end

        # write objective parameters in the global result dictionary that is returned by run_simulation_loop()
        write_optim_results!(sim_params["optimisation"]["objective_params"], output_data, 
                             optim_results)
        # calculate objective from the results 
        optim_results["objective"] = sim_params["optimisation"]["objective_function"](values(optim_results))
        
        # write necessary values for matrix_plot in the global result dictionary that is returned by run_simulation_loop()
        if io_settings["matrix_plot"] == "custom"
            write_optim_results!(io_settings["matrix_plot_spec"], output_data, 
                                 optim_results)
        end
        # TODO for validation/ test with Hafner
        # if !isnothing(findfirst(x -> x > 99.0, output_data["Hafner_Puffer_gross Load%"]))
        #     results["hours_to_full"] = findfirst(x -> x > 99.0, output_data["Hafner_Puffer_gross Load%"]) * sim_params["time_step_seconds"] / 3600
        # else 
        #     results["hours_to_full"] = NaN
        # end
        # if !isnothing(findfirst(isequal(0.0), output_data["Hafner_Puffer_gross Load%"]))
        #     results["hours_to_empty"] = findfirst(isequal(0.0), output_data["Hafner_Puffer_gross Load%"]) * sim_params["time_step_seconds"] / 3600
        # else
        #     results["hours_to_empty"] = NaN
        # end
        # results["Eigennutzungsgrad"] = sum(output_data["m_e_ac_230v EnergyFlow Hafner_PV_Freiflaeche->Hafner_WP"]) / sum(output_data["Hafner_PV_Freiflaeche Supply"])
        # results["Eigenversorgungsgrad"] = sum(output_data["m_e_ac_230v EnergyFlow Hafner_PV_Freiflaeche->Hafner_WP"]) / sum(output_data["Hafner_WP m_e_ac_230v IN"])
        # if haskey(output_data, "Grid_Price_IN Temperature")
        #     results["energy_cost_grid"] = sum(output_data["Hafner_Stromnetz_IN m_e_ac_230v OUT"] .* output_data["Grid_Price_IN Temperature"] ./ 10^6)
        #     results["energy_cost_pv"] = sum(output_data["m_e_ac_230v EnergyFlow Hafner_PV_Freiflaeche->Hafner_WP"]) .* 0.08 ./ 10^3
        #     results["energy_cost_wp_total"] = results["energy_cost_grid"] + results["energy_cost_pv"]
        #     results["objective"] = results["energy_cost_wp_total"]
        # else
        #     results["objective"] = results["Eigenversorgungsgrad"]
        # end
    end

    # write output to CSV if not done continuously
    if do_write_CSV 
        if do_write_CSV_continuously
            @info "CSV-file with outputs continuously written to $(csv_file_path)"
        else
            open(sim_params["run_path"](csv_file_path), "a") do file_handle
                for row_idx in 1:size(output_csv)[1]
                    write(file_handle, join(output_csv[row_idx, :], ";") * "\n")
                end
            end
            @info "CSV-file with outputs written to $(csv_file_path)"
        end
    end

    if do_write_summary_CSV
        output_path = replace(csv_file_path, r"\.csv$"i => "_aggregated.csv")

        success = aggregate_csv(csv_file_path,
                                output_path,
                                output_keys_to_CSV,
                                weather_CSV_keys,
                                ["Transfer", "Demand", "IN", "OUT", "Supply", "Losses", 
                                 "Gains", "EnergyFlow", "Balance", "Charge"], # energy terms
                                ["_sum"],                           # cumulative terms
                                ["COP", "Time_active", "Avg_PLR", "MixingTemperature_Input",
                                 "MixingTemperature_Output"],       # zero as missing terms
                                ["Time"],                           # time columns
                                ';',                                # separator
                                sim_params["epsilon"],              # threshold
                                io_settings["fixed_output_precision"])
        success && @info "Summary CSV-file with outputs created and written to $(output_path)"
    end

    # create profile line plot
    if do_create_plot_data || do_create_plot_weather
        output_data_lineplot = Matrix(view(output_data_all_requested, :,
                                           time_and_subset_cols(key_indexes, output_keys_lineplot)))
        create_profile_line_plots(output_data_lineplot,
                                  output_keys_lineplot,
                                  output_weather_lineplot,
                                  weather_data_keys,
                                  io_settings,
                                  sim_params)
        filepath = sim_params["run_path"](io_settings["output_plot_file"])
        @info "Line plot created and saved to $(sim_params["run_path"](filepath))"
    end

    # create Sankey diagram
    if do_create_sankey
        create_sankey(output_sourcenames_sankey,
                      output_targetnames_sankey,
                      output_interface_values,
                      medium_of_interfaces,
                      nr_of_interfaces,
                      io_settings,
                      sim_params)
        filepath = sim_params["run_path"](io_settings["sankey_plot_file"])
        @info "Sankey created and saved to $filepath"
    end

    # plot additional figures potentially available from components after simulation
    if io_settings["auxiliary_plots"]
        component_list = []
        output_path = sim_params["run_path"](io_settings["auxiliary_plots_path"])
        for component in components
            if plot_optional_figures_end(component[2], sim_params, output_path)
                push!(component_list, component[2].uac)
            end
        end
        if length(component_list) > 0
            @info "(Further) auxiliary plots are saved to folder $(output_path) for the " * 
                  "following components: $(join(component_list, ", "))"
        end
    end

    # output economic results
    if do_calculate_economy
        if io_settings["plot_economic_cashflows"]
            filepath = sim_params["run_path"](io_settings["economic_plot_cashflows_file_path"])
            success = plot_economic_results(economic_result, filepath, sim_params,
                                            io_settings["fixed_output_precision"], "cashflows")
            success && @info "Economy plot created and saved to $filepath"
        end

        if io_settings["plot_economic_present_values"]
            filepath = sim_params["run_path"](io_settings["economic_plot_present_values_file_path"])
            success = plot_economic_results(economic_result, filepath, sim_params,
                                            io_settings["fixed_output_precision"], "present_values")
            success && @info "Economy plot created and saved to $filepath"
        end

        # export economic results to CSV
        if io_settings["output_economic_CSV"]
            filepath = sim_params["run_path"](io_settings["economic_CSV_file_path"])
            success = write_economic_results_to_CSV(economic_result, filepath, sim_params)
            success && @info "Economic results exported to $filepath"
        end
    end

    # output emissions results
    if do_calculate_emissions
        # plot figure with yearly emissions
        if io_settings["plot_emission_results"]
            filepath = sim_params["run_path"](io_settings["emissions_plot_file_path"])
            success = plot_emissions_results(emissions_result, filepath, sim_params,
                                             io_settings["fixed_output_precision"])
            success && @info "Emissions plot created and saved to $filepath"
        end

        if io_settings["output_emissions_CSV"]
            # export emissions results to CSV
            filepath = sim_params["run_path"](io_settings["emissions_CSV_file_path"])
            success = write_emissions_results_to_CSV(emissions_result, filepath, sim_params)
            success && @info "Emissions results exported to $filepath"
        end
    end

    # plot utilized price and emission profiles
    if (do_calculate_economy || do_calculate_emissions) && io_settings["plot_price_and_emission_profiles"]
        filepath = sim_params["run_path"](io_settings["price_and_emission_profile_file_path"])
        success = plot_extended_price_and_emissions_profiles(economic_result, emissions_result, filepath, sim_params,
                                                             io_settings["fixed_output_precision"])
        success && @info "Utilized price and emission profiles exported as plot to $filepath"
    end

    return do_optimise ? optim_results : nothing
end

"""
    load_and_run(filepath, run_ID)

Load a project from the given file and run the simulation with it.

# Arguments
- `filepath::String`: Filepath to the project config file.
- `run_ID::UUID`: The run ID used in the run registry
# Returns
- `Bool`: `true` if the simulation was successful, `false` otherwise.
"""
function load_and_run(filepath::String, run_ID::UUID)::Bool
    start = now()
    @globalInfo "---- Simulation setup ----"
    @globalInfo "-- Starting simulation at $(start)"
    @globalInfo "-- Now reading project config"

    project_config = nothing

    try
        # we can't use run_path() here because that is only defined after loading the config
        # in the first place. so we hope we've been given a valid path and forbid upwards
        # path traversal for security reasons
        if occursin("..", filepath)
            @error "Project config filepath must not contain .. path traversal"
            return false
        end
        project_config = read_JSON(abspath(filepath))
    catch exc
        if isa(exc, MethodError)
            @error "Could not parse project config file at $(abspath(filepath))"
            return false
        end
    end

    if project_config === nothing
        @error "Could not find or parse project config file at $(abspath(filepath))"
        return false
    end

    @globalInfo "-- Now preparing inputs"

    # establish overarching locks for parallelization
    run_lock = ReentrantLock()
    output_lock = ReentrantLock()
    results_lock = ReentrantLock()

    io_settings = get_io_settings(project_config)
    sim_params = get_simulation_params(project_config, io_settings)

    if sim_params["optimisation"]["run_optimisation"]
        optim_results_path = sim_params["run_path"](io_settings["optimisation_csv_file_path"])
        open(optim_results_path, "w") do f end #TODO maybe remove
        optimiser = sim_params["optimisation"]

        all_results = []
        if optimiser["type"] == "monte_carlo_annealing"
            obj = Array{Union{Float64,Nothing}}(nothing)
            obj_lock = ReentrantLock()
        end
        nr_runs = Atomic{Int}(1)
        @globalInfo "Starting Simulations on $(Threads.nthreads()) Threads"
        #TODO find a way to cancel all the runs with STRG+C besides smashing the keys
        if length(optimiser["iterator"]) > 1
            @threads for sample_values in collect(optimiser["iterator"])
                run_nr = nr_runs[]
                atomic_add!(nr_runs, 1)
                sample_ID = uuid4()
                start_time = now()

                # decide which algorithm to run based on type of optimiser
                if optimiser["type"] == "parametervariation"
                    sample_params = Dict{String, Any}(zip(optimiser["optim_params_keys"], sample_values))
                    results = run_sample(io_settings, sim_params, optim_results_path, 
                                         project_config, sample_params, sample_ID, 
                                         run_lock, output_lock)
                    lock(results_lock) do 
                        push!(all_results, results) 
                    end

                elseif optimiser["type"] == "monte_carlo_annealing"
                    monte_carlo_annealing!(all_results, obj, obj_lock, 
                                           sim_params, optim_results_path, project_config, 
                                           sample_values, sample_ID, 
                                           run_lock, output_lock, results_lock)
                end
                runtime = round(Int, seconds(now() - start_time))
                max_runs = length(optimiser["iterator"])
                eta = round(Int, (max_runs - run_nr) * runtime / Threads.nthreads() / 60)
                @globalInfo "[$run_nr/$max_runs] → completed in $runtime s. ETA: $eta min"
            end
        else
            start_time = now()
            #TODO evtl. replace optimize and bboptimize with Optimization.jl interface
            if optimiser["type"] == "Optim"
                optimize(sample_values -> optim_func!(all_results, io_settings, sim_params, optim_results_path, 
                                                      project_config, sample_values, 
                                                      run_lock, output_lock, results_lock), 
                         optimiser["args"]...)
            elseif optimiser["type"] == "BlackBoxOptim"
                bboptimize(sample_values -> optim_func!(all_results, io_settings, sim_params, optim_results_path, 
                                                        project_config, sample_values, 
                                                        run_lock, output_lock, results_lock),
                           optimiser["args"]...; optimiser["kwargs"]...)
            end
            runtime = round(Int, seconds(now() - start_time))
            @globalInfo "[$(length(total_results)) runs → completed in $runtime s."
        end

        if !io_settings["write_optimisation_csv_continuously"]
            open(optim_results_path, "w") do file_handle
                # write header
                header = join(collect(keys(all_results[1])), ';') * "\n"
                write(file_handle, header)

                # write rows TODO can probably speed up by collecting data and writing once
                for results in all_results            
                    row = join(collect(values(results)), ';') * "\n"
                    row = replace(row, '.' => ',')
                    write(file_handle, row)
                end
            end
        end

        if io_settings["matrix_plot"] != "nothing"
            create_matrix_plot(all_results, io_settings, sim_params)
        end
    else
        _ = run_sample(io_settings, sim_params, nothing, project_config, 
                       nothing, run_ID, run_lock, output_lock)
    end

    return true
end

#TODO move everything connected to optimisation to new file
function monte_carlo_annealing!(all_results, obj, obj_lock, sim_params, 
                                optim_results_path, project_config, idx,
                                run_ID, run_lock, output_lock, results_lock)
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

function optim_func!(all_results, io_settings, sim_params, optim_results_path, project_config, 
                     sample_values, run_lock, output_lock, results_lock)
    #TODO implement batch evaluation for Metaheuristics e.g. if size(sample_values) then @threads?
    sample_params = Dict{String, Any}(zip(sim_params["optimisation"]["optim_params_keys"], sample_values))
    run_ID = uuid4()
    results = run_sample(io_settings, sim_params, optim_results_path, project_config, 
                         sample_params, run_ID, run_lock, output_lock)

    lock(results_lock) do 
        push!(all_results, results)
    end

    return results["objective"]
end

function run_sample(io_settings, sim_params, optim_results_path, project_config, 
                    sample_params, run_ID, run_lock, output_lock)
    start = now()
    if sample_params !== nothing
        project_config = create_variant(io_settings, sim_params, project_config, sample_params)
    end

    results = OrderedDict()

    try
        sim_params, io_settings, components, operations = prepare_inputs(project_config, run_ID)
        @info "-- Simulation setup complete in $(seconds(now() - start)) s"

        lock(run_lock) do 
            current_runs[run_ID] = SimulationRun(sim_params, io_settings, components, operations)
        end

        start = now()
        @info "---- Simulation loop ----"
        
        if !isnothing(sample_params)
            for (key, value) in pairs(sample_params)
                results[key] = value
            end 
        end

        sim_output = run_simulation_loop(sim_params, io_settings, components, operations)
        if !isnothing(sim_output)
            for (key, value) in pairs(sim_output)
                results[key] = value
            end
        end
 
        results["error"] = ""

    catch e
        if !isnothing(optim_results_path) && filesize(optim_results_path) == 0
            throw(e)
        end
        # save excact error message to output file
        error_message = sprint(showerror, e)
        full_error_message = error_message * "\n" * sprint(Base.show_backtrace, catch_backtrace())
        @globalInfo full_error_message
        results["error"] =  "\"" * replace(full_error_message, "\"" => "\"\"") * "\"\n"

        if sim_params["optimisation"]["run_optimisation"]
            for key in parse_outkeys(sim_params["optimisation"]["output_keys_sum"])
                results[key] = NaN
            end
            for key in parse_outkeys(sim_params["optimisation"]["output_keys_mean"])
                results[key] = NaN
            end 
            if haskey(sim_params["optimisation"], "objective_function")
                results["objective"] = NaN
            end
        end
    end

    if sim_params["optimisation"]["run_optimisation"] && io_settings["write_optimisation_csv_continuously"]
        # Write results to seperate file after all simulations are finished.
        row = join(collect(values(results)), ';') * "\n"
        row = replace(row, '.' => ',')
        # Lock the file writing
        lock(output_lock) do 
            # create header if file is empty
            if filesize(optim_results_path) == 0
                header = join(collect(keys(results)), ';') * "\n"
                open(optim_results_path, "w") do file_handle
                    write(file_handle, header)
                end
            end
            open(optim_results_path, "a") do file_handle
                write(file_handle, row)
            end
        end
    end
        
    @info "-- Simulation loop complete in $(seconds(now() - start)) s"
    lock(run_lock) do 
        close_run(run_ID)
    end

    return results
end

function create_variant(io_settings::Dict{String,Any}, sim_params::Dict{String, Any}, 
                        project_config::OrderedDict{String,Any}, 
                        sample_params::Dict{String,Any})
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

    #TODO maybe this should be moved to profile processing to allow the profiles to be 
    # defined with "profiles" group, scale and addon without optimiser
    if haskey(cfg, "profiles")
        profile_paths = Dict{String,String}()
        profile_scales = Dict{String,Float64}()
        profile_addons = Dict{String,Float64}()
        for (name, profile) in pairs(cfg["profiles"])
            profile_paths[name] = profile["path"]
            profile_scales[name] = profile["scale"]
            profile_addons[name] = profile["addon"]
        end
    
        # create correct profiles from price_profile_paths and add the to sim_output.
        # profiles will be overwritten for every run to make sure the profile_scales and 
        # profile_addons calculated correctly.
        # If multiple threads are used each thread gets their own profile.
        #TODO change directiory to something better
        profile_dir = sim_params["run_path"]("./profiles/parallel_runs")
        mkpath(profile_dir)

        profile_id = ifelse(Threads.nthreads() > 1, Threads.threadid(), 0)
        date_range = remove_leap_days(collect(sim_params["start_date"]:Second(sim_params["time_step_seconds"]):sim_params["end_date"]))
        new_paths = Dict{String,String}() 

        for (name, path) in pairs(profile_paths)
            if profile_scales[name] != 1 && profile_addons[name] != 0 && profile_id == 0
                new_paths[name] = path
            else
                profile = Profile(path, sim_params)
                values = [profile.data[dt] .* profile_scales[name] .+ profile_addons[name] for dt in date_range]
                new_path = profile_dir * "/" * split(path[1:end-4], '/')[end] * "_$profile_id.prf" 
                save_to_prf(collect(date_range), values, new_path)
                new_paths[name] = path
            end
        end

        # replace the profile names with the paths to the new profiles
        function replace_profiles!(cfg::AbstractDict, replacements::Dict{String,String})
            for (k, v) in cfg
                if v isa String && haskey(replacements, v)
                    cfg[k] = replacements[v]
                elseif v isa AbstractDict
                    replace_profiles!(v, replacements)
                end
            end

        end

        replace_profiles!(cfg, new_paths)
    end

    return cfg
end

end # module

