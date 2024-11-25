module Resie

# note: includes that contain their own module, which have to be submodules of the Resie
# module, are included first, then can be accessed with the "using" keyword. files that
# contain code that is intended to be used in-place of their include statement (as part
# of the Resie module), are included after the "using" statements have been declared.
# this is done so the latter files can access the symbols of the submodules the same as
# if the code was inside this file.

include("profiles/base.jl")
using .Profiles
include("profiles/weatherdata.jl")
using .Weatherdata

include("energy_systems/base.jl")
using .EnergySystems

include("project_loading.jl")
include("file_output.jl")

include("resie_logger.jl")
using .Resie_Logger

using PlotlyJS
using ColorSchemes
using Colors
using Interpolations
using JSON
using Dates
using UUIDs

const HOURS_PER_SECOND::Float64 = 1.0 / 3600.0
const SECONDS_PER_HOUR::Float64 = 3600.0

"""
Contains the parameters, instantiated components and the order of operations for a simulation run.

In short, it bundles all the data required to perform a simulation from start to finish.
Through its fields it contains a lot of complexity in terms of hierarchical data structures.
However the struct should not be used as argument for complex calculations. Instead it is
mainly used for indexing the various sub data structures.
"""
mutable struct SimulationRun
    parameters::Dict{String,Any}
    components::Grouping
    order_of_operations::StepInstructions
end

# these variables should be the only global state in the package and contain the state for
# ongoing or paused simulation runs
current_runs::Dict{UUID,SimulationRun} = Dict{UUID,SimulationRun}()

"""
    get_simulation_params(project_config)

Constructs the dictionary of simulation parameters.

# Arguments
-`project_config::Dict{AbstractString,Any}`: The project config
# Returns
-`Dict{String,Any}`: The simulation parameter dictionary
"""
function get_simulation_params(project_config::Dict{AbstractString,Any})::Dict{String,Any}
    time_step, start_date, end_date, nr_of_steps = get_timesteps(project_config["simulation_parameters"])

    sim_params = Dict{String,Any}(
        "time" => 0,
        "current_date" => start_date,
        "time_step_seconds" => time_step,
        "number_of_time_steps" => nr_of_steps,
        "start_date" => start_date,
        "end_date" => end_date,
        "epsilon" => 1e-9,
        "latitude" => default(project_config["simulation_parameters"], "latitude", nothing),
        "longitude" => default(project_config["simulation_parameters"], "longitude", nothing),
    )

    # add helper functions to convert power to work and vice-versa. this uses the time step
    # of the simulation as the duration required for the conversion.
    sim_params["watt_to_wh"] = function (watts::Float64)
        return watts * time_step * HOURS_PER_SECOND
    end
    sim_params["wh_to_watts"] = function (wh::Float64)
        return wh * SECONDS_PER_HOUR / time_step
    end

    weather_file_path = default(project_config["simulation_parameters"],
                                "weather_file_path",
                                nothing)

    if weather_file_path !== nothing
        sim_params["weather_data"], lat, long = WeatherData(weather_file_path, sim_params)

        if sim_params["latitude"] === nothing || sim_params["longitude"] === nothing
            sim_params["latitude"] = lat
            sim_params["longitude"] = long
        else
            @info "The coordinates given in the weather file where overwritten by the " *
                  "ones given in the input file:\n" *
                  "Latidude: $(sim_params["latitude"]); Longitude: $(sim_params["longitude"])"
        end
    end

    return sim_params
end

"""
    prepare_inputs(project_config)

Construct and prepare parameters, energy system components and the order of operation.

# Arguments
-`project_config::Dict{AbstractString,Any}`: The project config
# Returns
-`Dict{String,Any}`: Simulation parameters
-`Grouping`: The constructed energy system components
-`StepInstructions`: Order of operations
"""
function prepare_inputs(project_config::Dict{AbstractString,Any})
    sim_params = get_simulation_params(project_config)

    components = load_components(project_config["components"], sim_params)

    if haskey(project_config, "order_of_operation") && length(project_config["order_of_operation"]) > 0
        step_order = load_order_of_operations(project_config["order_of_operation"], components)
        @info "The order of operations was successfully imported from the input file.\n" *
              "Note that the order of operations has a major impact on the simulation " *
              "result and should only be changed by experienced users!"
    else
        step_order = calculate_order_of_operations(components)
    end

    return sim_params, components, step_order
end

"""
    run_simulation_loop()

Performs the simulation as loop over time steps and records outputs.

# Arguments
-`project_config::Dict{AbstractString,Any}`: The project config
-`sim_params::Dict{String,Any}`: Simulation parameters
-`components::Grouping`: The energy system components
-`step_order::StepInstructions`:: Order of operations
"""
function run_simulation_loop(project_config::Dict{AbstractString,Any},
                             sim_params::Dict{String,Any},
                             components::Grouping,
                             step_order::StepInstructions)
    # get list of requested output keys for lineplot and csv export
    output_keys_lineplot, output_keys_to_csv = get_output_keys(project_config["io_settings"], components)
    do_create_plot = !(output_keys_lineplot === nothing)
    do_write_CSV = !(output_keys_to_csv === nothing)
    csv_output_file_path = default(project_config["io_settings"],
                                   "csv_output_file",
                                   "./output/out.csv")
    csv_time_unit = default(project_config["io_settings"], "csv_time_unit", "seconds")
    if !(csv_time_unit in ["seconds", "minutes", "hours", "date"])
        @error "The `csv_time_unit` has to be one of: seconds, minutes, hours, date!"
        throw(IntputError)
    end

    # Initialize the array for output plots
    if do_create_plot
        output_data_lineplot = zeros(Float64, sim_params["number_of_time_steps"], 1 + length(output_keys_lineplot))
    end
    # reset CSV file
    if do_write_CSV
        reset_file(csv_output_file_path, output_keys_to_csv, csv_time_unit)
    end

    # check if sankey should be plotted
    do_create_sankey = haskey(project_config["io_settings"], "sankey_plot") &&
                       project_config["io_settings"]["sankey_plot"] !== "nothing"
    if do_create_sankey
        # get infomration about all interfaces for Sankey
        nr_of_interfaces,
        medium_of_interfaces,
        output_sourcenames_sankey,
        output_targetnames_sankey = get_interface_information(components)
        # preallocate for speed: Matrix with data of interfaces in every timestep
        output_interface_values = zeros(Float64, sim_params["number_of_time_steps"], nr_of_interfaces)
    end

    # export order of operation and other additional info like optional plots
    dump_auxiliary_outputs(project_config, components, step_order, sim_params)

    for steps in 1:sim_params["number_of_time_steps"]
        # perform the simulation
        perform_steps(components, step_order, sim_params)

        # check if any component was not balanced
        warnings = check_balances(components, sim_params["epsilon"])
        if length(warnings) > 0
            for (key, balance) in warnings
                @balanceWarn "Balance for component $key was not zero in timestep $(sim_params["time"]): $balance"
            end
        end

        # write requested output data of the components to CSV-file
        # This is currently done in every time step to keep data even if 
        # an error occurs.
        if do_write_CSV
            write_to_file(csv_output_file_path, output_keys_to_csv, sim_params, csv_time_unit)
        end

        # get the energy transported through each interface in every timestep for Sankey
        if do_create_sankey
            output_interface_values[steps, :] = collect_interface_energies(components, nr_of_interfaces)
        end

        # gather output data of each component for line plot
        if do_create_plot
            output_data_lineplot[steps, :] = geather_output_data(output_keys_lineplot, sim_params["time"])
        end

        # simulation update
        sim_params["time"] += Int(sim_params["time_step_seconds"])
        sim_params["current_date"] = add_ignoring_leap_days(sim_params["current_date"],
                                                            Second(sim_params["time_step_seconds"]))
    end

    ### create profile line plot
    if do_create_plot
        create_profile_line_plots(output_data_lineplot, output_keys_lineplot, project_config, sim_params)
        @info "Line plot created and saved to .output/output_plot.html"
    end

    ### create Sankey diagram
    if do_create_sankey
        create_sankey(output_sourcenames_sankey,
                      output_targetnames_sankey,
                      output_interface_values,
                      medium_of_interfaces,
                      nr_of_interfaces,
                      project_config["io_settings"])
        @info "Sankey created and saved to .output/output_sankey.html"
    end

    if do_write_CSV
        @info "CSV-file with outputs written to $(csv_output_file_path)"
    end
end

"""
    load_and_run(filepath)

Load a project from the given file and run the simulation with it.

# Arguments
- `filepath::String`: Filepath to the project config file.
"""
function load_and_run(filepath::String)
    project_config = nothing

    try
        project_config = read_JSON(abspath(filepath))
    catch exc
        if isa(exc, MethodError)
            @error "Could not parse project config file at $(filepath)"
            return
        end
    end

    if project_config === nothing
        @error "Could not find or parse project config file at $(filepath)"
        return
    end

    sim_params, components, step_order = prepare_inputs(project_config)
    run_ID = uuid1()
    sim_params["run_ID"] = run_ID
    current_runs[run_ID] = SimulationRun(sim_params, components, step_order)

    run_simulation_loop(project_config, sim_params, components, step_order)
end

end # module
