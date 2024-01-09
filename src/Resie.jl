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

"""
    run_simulation()

Read inputs, perform the simulation calculation and write outputs.

Due to the complexity of required inputs of a simulation and how the outputs are persisted
(to file), this function takes only one argument, namely the project config, and returns
nothing.
"""
function run_simulation(project_config::Dict{AbstractString,Any})
    # get time steps from input file
    time_step, start_timestamp, end_timestamp, start_date = get_timesteps(project_config["simulation_parameters"])  
    nr_of_steps = UInt(max(1, (end_timestamp - start_timestamp) / time_step))

    sim_params = Dict{String,Any}(
        "time" => start_timestamp,
        "time_step_seconds" => time_step,
        "number_of_time_steps" => nr_of_steps,
        "start_date" => start_date,
        "epsilon" => 1e-9
    )
    EnergySystems.set_timestep(sim_params["time_step_seconds"])

    # load weather data 
    if haskey(project_config["io_settings"], "weather_file_path") && length(project_config["io_settings"]["weather_file_path"]) > 0
        sim_params["weatherdata"] = WeatherData(project_config["io_settings"]["weather_file_path"], sim_params)
    end

    components = load_components(project_config["components"], sim_params)

    if haskey(project_config, "order_of_operation") && length(project_config["order_of_operation"]) > 0
        step_order = load_order_of_operations(project_config["order_of_operation"], components)
        @info "The order of operations was successfully imported from the input file.\n" *
              "Note that the order of operations has a major impact on the simulation result and should only be changed by experienced users!"
    else
        step_order = calculate_order_of_operations(components)
    end

    # get list of requested output keys for lineplot and csv export
    output_keys_lineplot, output_keys_to_csv = get_output_keys(project_config["io_settings"], components)
    do_create_plot = !(output_keys_lineplot === nothing)
    do_write_CSV = !(output_keys_to_csv === nothing)

    # Initialize the array for output plots
    if do_create_plot
        output_data_lineplot = zeros(Float64, nr_of_steps, 1 + length(output_keys_lineplot))
    end
    # reset CSV file
    if do_write_CSV
        reset_file(project_config["io_settings"]["output_file"], output_keys_to_csv)
    end
   
    # get infomration about all interfaces for Sankey
    nr_of_interfaces, medium_of_interfaces, output_sourcenames_sankey, output_targetnames_sankey = get_interface_information(components)
   
    # preallocate for speed: Matrix with data of interfaces in every timestep
    output_interface_values = zeros(Float64, nr_of_steps, nr_of_interfaces)

    # export order or operatin (OoO)
    if project_config["io_settings"]["dump_info"]
        dump_info(
            project_config["io_settings"]["dump_info_file"],
            components, step_order, sim_params
        )
    end

    for steps = 1:nr_of_steps
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
            write_to_file(project_config["io_settings"]["output_file"], output_keys_to_csv, sim_params["time"])
        end

        # get the energy transported through each interface in every timestep for Sankey
        output_interface_values[steps, :] = collect_interface_energies(components, nr_of_interfaces)
      
        # gather output data of each component for line plot
        if do_create_plot
            output_data_lineplot[steps, :] = geather_output_data(output_keys_lineplot, sim_params["time"])
        end

        # simulation update
        sim_params["time"] += Int(sim_params["time_step_seconds"])
    end

    ### create profile line plot
    if do_create_plot
        create_profile_line_plots(output_data_lineplot, output_keys_lineplot, project_config, sim_params)
        @info "Line plot created and saved to .output/output_plot.html"

    end

    ### create Sankey diagram
    create_sankey(output_sourcenames_sankey, output_targetnames_sankey, output_interface_values, medium_of_interfaces, nr_of_interfaces)
    @info "Sankey created and saved to .output/output_sankey.html"
    
    if do_write_CSV
        @info "CSV-file with outputs written to $(project_config["io_settings"]["output_file"])"
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

    run_simulation(project_config)
end

end # module
