module Resie

using Printf
using Dates: now, seconds
using UUIDs
using Base.Threads
using Logging

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

include("resie_logger.jl")
using .Resie_Logger

include("project_loading.jl")
include("simulation_output.jl")
include("parameter_study_output.jl")
include("economy.jl")
include("emissions.jl")
include("parameter_study.jl")

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
- `cancel_flag::Union{Nothing,Threads.Atomic{Bool}}`: Flag to pass STR+C down to all parallel runs
"""
function run_simulation_loop(sim_params::Dict{String,Any},
                             io_settings::Dict{String,Any},
                             components::Grouping,
                             operations::OrderOfOperations;
                             suppress_all_output::Bool=false,
                             cancel_flag::Union{Nothing,Threads.Atomic{Bool}}=nothing)
    collect_objective_results = sim_params["parameter_study"]["runtime"]["enabled"]
    # get list of requested output keys for lineplot and csv export
    output_keys_lineplot,
    output_keys_to_CSV,
    output_keys_economic_emissions,
    output_keys_parameter_study = get_output_keys(io_settings,
                                                  sim_params["economic_parameters"],
                                                  sim_params["emissions_parameters"],
                                                  sim_params["parameter_study"],
                                                  components,
                                                  suppress_all_output)
    all_requested_output_keys = Vector{Resie.EnergySystems.OutputKey}(unique(vcat(something(output_keys_lineplot,
                                                                                            String[]),
                                                                                  something(output_keys_economic_emissions,
                                                                                            String[]),
                                                                                  something(output_keys_parameter_study,
                                                                                            String[]))))
    weather_data_keys = get_weather_data_keys(sim_params, suppress_all_output)
    do_create_plot_data = output_keys_lineplot !== nothing
    do_create_plot_weather = weather_data_keys !== nothing && io_settings["plot_weather_data"]
    do_write_CSV_weather = weather_data_keys !== nothing && io_settings["csv_output_weather"]
    weather_CSV_keys = do_write_CSV_weather ? weather_data_keys : nothing
    do_write_CSV = output_keys_to_CSV !== nothing || do_write_CSV_weather
    do_write_CSV_continuously = io_settings["write_csv_continuously"]
    do_write_summary_CSV = !suppress_all_output && io_settings["write_summary_csv"]
    csv_file_path = io_settings["csv_output_file_path"]
    csv_time_unit = io_settings["csv_time_unit"]
    do_calculate_economy = sim_params["economic_parameters"]["calculate_economy"]
    do_calculate_emissions = sim_params["emissions_parameters"]["calculate_emissions"]
    do_create_sankey = !suppress_all_output && io_settings["sankey_plot"] !== "nothing"

    # Initialize the arrays for output
    output_weather_lineplot = do_create_plot_weather ?
                              zeros(Float64, sim_params["number_of_time_steps_output"], 1 + length(weather_data_keys)) :
                              nothing
    output_csv = do_write_CSV && !do_write_CSV_continuously ?
                 Matrix{String}(undef, sim_params["number_of_time_steps_output"],
                                1 + length(output_keys_to_CSV) + (do_write_CSV_weather ? length(weather_data_keys) : 0)) :
                 nothing
    output_data_all_requested = do_create_plot_data || do_calculate_economy || do_calculate_emissions ||
                                collect_objective_results ?
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
    dump_auxiliary_outputs(io_settings, components, operations, sim_params, suppress_all_output)

    @info "-- Start time step loop"
    if sim_params["start_date_output"] == sim_params["start_date"]
        @info "-- No preheating activated. Starting output from beginning."
    else
        @info "-- Starting with preheating, no output will be written until completed."
    end
    start = now()
    for steps in 1:sim_params["number_of_time_steps"]
        # handle interruption from STRG+C for multi-thread simulations
        if cancel_flag !== nothing && cancel_flag[]
            throw(InterruptException())
        end

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
                                     csv_time_unit,
                                     io_settings)
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
            if do_create_plot_data || do_calculate_economy || do_calculate_emissions || collect_objective_results
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

    # Map each requested OutputKey to its column in output_data_all_requested.
    key_indexes = output_key_indexes(all_requested_output_keys)

    if do_calculate_economy || do_calculate_emissions
        economic_emissions_columns = output_data_columns(key_indexes, output_keys_economic_emissions; include_time=true)
        output_data_economic_emissions = Matrix(view(output_data_all_requested, :, economic_emissions_columns))
        economic_emissions_data = prepare_economic_emissions_data(components, output_keys_economic_emissions,
                                                                  output_data_economic_emissions)
        economic_result = do_calculate_economy ? calculate_economy(economic_emissions_data, sim_params) : nothing
        emissions_result = do_calculate_emissions ? calculate_emissions(economic_emissions_data, sim_params) : nothing
    end

    if collect_objective_results
        objective_results = Dict{String,Union{Array{Float64},Float64}}()
        parameter_study_columns = output_data_column_map(key_indexes, output_keys_parameter_study)

        function write_objective_results!(params::Array{String},
                                          output_data::AbstractMatrix{<:Real},
                                          output_columns::AbstractDict{String,Int},
                                          res::Dict{String,Union{Array{Float64},Float64}})
            for key in params
                func, spec = split(key, " "; limit=2)
                if func == "sum" || func == "mean"
                    values = view(output_data, :, output_columns[spec])
                    res[key] = func == "sum" ? sum(values) : sum(values) / length(values)
                elseif func == "economic"
                    res[key] = getfield(economic_result, Symbol(spec))
                elseif func == "emissions"
                    res[key] = getfield(emissions_result, Symbol(spec))
                end
            end
        end
        # write objective parameters in the global result dictionary that is returned by run_simulation_loop()
        write_objective_results!(sim_params["parameter_study"]["runtime"]["objective_params_keys"],
                                 output_data_all_requested,
                                 parameter_study_columns,
                                 objective_results)
        # calculate objective from the results 
        objective_values = [objective_results[key]
                            for key in sim_params["parameter_study"]["runtime"]["objective_params_keys"]]
        objective_results["objective"] = sim_params["parameter_study"]["runtime"]["objective_function"](objective_values)
    end

    # write output to CSV if not done continuously
    if do_write_CSV
        csv_file_path_abs = sim_params["run_path"](csv_file_path)
        if do_write_CSV_continuously
            @info "CSV-file with outputs continuously written to $(csv_file_path_abs)"
        else
            open(csv_file_path_abs, "a") do file_handle
                for row_idx in 1:size(output_csv)[1]
                    write(file_handle, join(output_csv[row_idx, :], ";") * "\n")
                end
            end
            @info "CSV-file with outputs written to $(csv_file_path_abs)"
        end
    end

    if do_write_CSV && do_write_summary_CSV
        output_path = sim_params["run_path"](replace(csv_file_path, r"\.csv$"i => "_aggregated.csv"))

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
        lineplot_columns = output_data_columns(key_indexes, output_keys_lineplot; include_time=true)
        output_data_lineplot = Matrix(view(output_data_all_requested, :, lineplot_columns))
        create_profile_line_plots(output_data_lineplot,
                                  output_keys_lineplot,
                                  output_weather_lineplot,
                                  weather_data_keys,
                                  io_settings,
                                  sim_params)
        filepath = sim_params["run_path"](io_settings["output_plot_file_path"])
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
        filepath = sim_params["run_path"](io_settings["sankey_plot_file_path"])
        @info "Sankey created and saved to $filepath"
    end

    # plot additional figures potentially available from components after simulation
    if io_settings["auxiliary_plots"] && !suppress_all_output
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
    if do_calculate_economy && !suppress_all_output
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
        if io_settings["output_economic_csv"]
            filepath = sim_params["run_path"](io_settings["economic_csv_file_path"])
            success = write_economic_results_to_CSV(economic_result, filepath, sim_params)
            success && @info "Economic results exported to $filepath"
        end
    end

    # output emissions results
    if do_calculate_emissions && !suppress_all_output
        # plot figure with yearly emissions
        if io_settings["plot_emission_results"]
            filepath = sim_params["run_path"](io_settings["emissions_plot_file_path"])
            success = plot_emissions_results(emissions_result, filepath, sim_params,
                                             io_settings["fixed_output_precision"])
            success && @info "Emissions plot created and saved to $filepath"
        end

        if io_settings["output_emissions_csv"]
            # export emissions results to CSV
            filepath = sim_params["run_path"](io_settings["emissions_csv_file_path"])
            success = write_emissions_results_to_CSV(emissions_result, filepath, sim_params)
            success && @info "Emissions results exported to $filepath"
        end
    end

    # plot utilized price and emission profiles
    if (do_calculate_economy || do_calculate_emissions) && io_settings["plot_price_and_emission_profiles"] &&
       !suppress_all_output
        filepath = sim_params["run_path"](io_settings["price_and_emission_profile_file_path"])
        success = plot_extended_price_and_emissions_profiles(economic_result, emissions_result, filepath, sim_params,
                                                             io_settings["fixed_output_precision"])
        success && @info "Utilized price and emission profiles exported as plot to $filepath"
    end

    return collect_objective_results ? objective_results : nothing
end

"""
    load_and_run(filepath, run_ID)

Load a project from the given file and run the simulation with it.

# Arguments
- `filepath::String`: Filepath to the project config file.
- `run_ID::UUID`: The run ID used in the run registry
- `logger::Union{Nothing,Resie_Logger.CustomLogger}`: Logger used for ReSiE
# Returns
- `Bool`: `true` if the simulation was successful, `false` otherwise.
"""
function load_and_run(filepath::String, run_ID::UUID; logger::Union{Nothing,Resie_Logger.CustomLogger}=nothing)::Bool
    start = now()
    success = true
    @globalInfo "---- Simulation setup ----"
    @globalInfo "Starting simulation at $(start)"
    @globalInfo "Now reading project config"

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

    @globalInfo "Now preparing inputs"
    preparation_cache = PreparationCache()

    # set log level by operating mode
    if logger !== nothing
        Resie_Logger.set_min_log_level!(logger, get_min_log_level(project_config, logger))
    end
    io_settings = get_io_settings(project_config)
    sim_params = get_simulation_params(project_config, io_settings; preparation_cache=preparation_cache)

    if sim_params["parameter_study"]["runtime"]["enabled"]
        # perform the configured parameter study
        success,
        evaluated_parameter_sets,
        local_evaluated_parameter_sets = perform_parameter_study(io_settings,
                                                                 sim_params,
                                                                 project_config;
                                                                 preparation_cache=preparation_cache)

        figure_results = if sim_params["parameter_study"]["sensitivity_analysis"]["include_local_sensitivity_results_in_figures"]
            vcat(evaluated_parameter_sets, local_evaluated_parameter_sets)
        else
            evaluated_parameter_sets
        end

        if !isempty(figure_results)
            create_parameter_study_diagnostic_plots(figure_results, io_settings, sim_params)
        end
    else
        # perform single simulation run
        # establish overarching locks for parallelization
        run_lock = ReentrantLock()
        output_lock = ReentrantLock()
        try
            _ = run_simulation_sample(io_settings, sim_params, nothing, project_config,
                                      nothing, run_ID, run_lock, output_lock;
                                      suppress_all_output=false,
                                      preparation_cache=preparation_cache)
        catch e
            if e isa InterruptException
                @globalInfo "Simulation interrupted by user."
                success = false
            else
                rethrow()
            end
        end
    end

    return success
end

"""
    run_simulation_sample(io_settings, sim_params, parameter_study_results_path, project_config,
                          parameter_set, run_ID, run_lock, output_lock; suppress_all_output=false)

Run a single simulation sample with given parameters.

# Arguments
- `io_settings::Dict{String,Any}`: IO settings
- `sim_params::Dict{String,Any}`: Simulation parameters
- `parameter_study_results_path::Union{String,Nothing}`: File path for parameter-study results
- `project_config::OrderedDict{String,Any}`: The project config
- `parameter_set::Dict{String, Any}`: Values and names of the sample parameters that get 
                                      used for the next simulation run
- `run_ID::UUID`: The run ID used in the run registry
- `run_lock::ReentrantLock`: Lock for writing to current_runs
- `output_lock::ReentrantLock`: Lock for file at parameter_study_results_path
- `suppress_all_output::Bool=false`: Bool that can be set to suppress the generation of all outputs written to hard drive.
- `cancel_flag::Union{Nothing,Threads.Atomic{Bool}}`: Flag to pass STR+C down to all parallel runs
# Returns
- `OrderedDict{String,Union{Float64, Int64, String}}`: Results of the simulation run
"""
function run_simulation_sample(io_settings::Dict{String,Any}, sim_params::Dict{String,Any},
                               parameter_study_results_path::Union{String,Nothing},
                               project_config::OrderedDict{String,Any},
                               parameter_set::Union{Dict{String,Any},Nothing}, run_ID::UUID,
                               run_lock::ReentrantLock, output_lock::ReentrantLock;
                               suppress_all_output::Bool=false,
                               cancel_flag::Union{Nothing,Threads.Atomic{Bool}}=nothing,
                               preparation_cache::Union{Nothing,PreparationCache}=nothing)::OrderedDict{String,Any}
    start = now()
    if !isnothing(parameter_set)
        project_config = create_parameter_variant(io_settings, sim_params, project_config, parameter_set)
    end

    results = OrderedDict{String,Any}()

    if !isnothing(parameter_set)
        for (key, value) in pairs(parameter_set)
            results[key] = value
        end
    end

    try
        sim_params, io_settings, components, operations = prepare_inputs(project_config, run_ID;
                                                                         preparation_cache=preparation_cache)
        @info "-- Simulation setup complete in $(seconds(now() - start)) s"

        lock(run_lock) do
            current_runs[run_ID] = SimulationRun(sim_params, io_settings, components, operations)
        end

        start = now()
        @info "---- Simulation loop ----"

        sim_output = run_simulation_loop(sim_params, io_settings, components, operations;
                                         suppress_all_output=suppress_all_output,
                                         cancel_flag=cancel_flag)
        if !isnothing(sim_output)
            for (key, value) in pairs(sim_output)
                results[key] = value
            end
        end

        results["error"] = ""

    catch e
        if e isa InterruptException
            if cancel_flag !== nothing
                cancel_flag[] = true
            end
            rethrow()
        end

        if !isnothing(parameter_study_results_path) && filesize(parameter_study_results_path) == 0
            throw(e)
        end

        if sim_params["parameter_study"]["runtime"]["enabled"]
            for key in sim_params["parameter_study"]["runtime"]["objective_params_keys"]
                results[key] = NaN
            end

            results["objective"] = Inf
        end

        # save exact error message to output file
        error_message = sprint(showerror, e)
        full_error_message = error_message * "\n" * sprint(Base.show_backtrace, catch_backtrace())
        @globalInfo full_error_message
        results["error"] = "\"" * replace(full_error_message, "\"" => "\"\"") * "\"\n"
    finally
        lock(run_lock) do
            close_run(run_ID)
        end
    end

    if sim_params["parameter_study"]["runtime"]["enabled"] &&
       io_settings["write_parameter_study_csv_continuously"] &&
       !isnothing(parameter_study_results_path)
        # Write results to file after the single simulation has finished.
        row = join(collect(values(results)), ';') * "\n"
        row = replace(row, ',' => ' ')
        row = replace(row, '.' => ',')
        # Lock the file writing
        lock(output_lock) do
            # create header if file is empty
            if filesize(parameter_study_results_path) == 0
                header = join(collect(keys(results)), ';') * "\n"
                open(parameter_study_results_path, "w") do file_handle
                    write(file_handle, header)
                end
            end
            open(parameter_study_results_path, "a") do file_handle
                write(file_handle, row)
            end
        end
    end

    @info "-- Simulation loop complete in $(seconds(now() - start)) s"
    return results
end

end # module
