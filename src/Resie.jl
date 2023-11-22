module Resie

# note: includes that contain their own module, which have to be submodules of the Resie
# module, are included first, then can be accessed with the "using" keyword. files that
# contain code that is intended to be used in-place of their include statement (as part
# of the Resie module), are included after the "using" statements have been declared.
# this is done so the latter files can access the symbols of the submodules the same as
# if the code was inside this file.

include("profiles/base.jl")
using .Profiles

include("energy_systems/base.jl")
using .EnergySystems

include("project_loading.jl")
include("file_output.jl")

using PlotlyJS
using ColorSchemes

"""
    run_simulation()

Read inputs, perform the simulation calculation and write outputs.

Due to the complexity of required inputs of a simulation and how the outputs are persisted
(to file), this function takes only one argument, namely the project config, and returns
nothing.
"""
function run_simulation(project_config::Dict{AbstractString,Any})
    components = load_components(project_config["components"])

    if haskey(project_config, "order_of_operation") && length(project_config["order_of_operation"]) > 0
        step_order = load_order_of_operations(project_config["order_of_operation"], components)
        println("The order of operations was successfully imported from the input file.\nNote that the order of operations has a major impact on the simulation result and should only be changed by experienced users!")
    else
        step_order = calculate_order_of_operations(components)
    end

    time_step = 900
    if "time_step_seconds" in keys(project_config["simulation_parameters"])
        time_step = UInt(project_config["simulation_parameters"]["time_step_seconds"])
    end

    start_timestamp = 0
    if "start" in keys(project_config["simulation_parameters"])
        start_timestamp = Integer(project_config["simulation_parameters"]["start"])
    end

    end_timestamp = 900
    if "end" in keys(project_config["simulation_parameters"])
        end_timestamp = Integer(project_config["simulation_parameters"]["end"])
    end

    nr_of_steps = UInt(max(1, (end_timestamp - start_timestamp) / time_step))

    parameters = Dict{String,Any}(
        "time" => start_timestamp,
        "time_step_seconds" => time_step,
        "epsilon" => 1e-9
    )
    EnergySystems.set_timestep(parameters["time_step_seconds"])

    # read in flags, if all outputs should be printed and/or plotted
    do_plot_all_outputs = false
    do_create_plot = false
    if haskey(project_config["io_settings"], "output_plot")
        do_create_plot = true  # if the key exists, set do_create_plot to true by default
        if project_config["io_settings"]["output_plot"] == "all"
            do_plot_all_outputs = true
        elseif project_config["io_settings"]["output_plot"] == "nothing"
            do_create_plot = false
        end
    end

    do_write_all_CSV_outputs = false
    do_write_CSV = false
    if haskey(project_config["io_settings"], "output_keys")
        do_write_CSV = true  # if the key exists, set do_create_plot to true by default
        if project_config["io_settings"]["output_keys"] == "all"
            do_write_all_CSV_outputs = true
        elseif project_config["io_settings"]["output_keys"] == "nothing"
            do_write_CSV = false
        end
    end

    # collect all possible outputs of all units if needed
    if do_plot_all_outputs || do_write_all_CSV_outputs
        all_output_keys = Vector{EnergySystems.OutputKey}()
        for unit in components
            temp_dict = Dict{String, Any}(unit[2].uac => output_values(unit[2]))
            append!(all_output_keys, output_keys(components, temp_dict))
        end
    end

    # collect output keys for lineplot and csv output
    if do_create_plot  # line plot
        output_keys_lineplot = do_plot_all_outputs ? all_output_keys : Vector{EnergySystems.OutputKey}()
        # Collect output keys if not plotting all outputs
        if !do_plot_all_outputs
            for plot in project_config["io_settings"]["output_plot"]
                key = plot[2]["key"]
                append!(output_keys_lineplot, output_keys(components, key))
            end
        end
        # Initialize the array for output plots
        output_data_lineplot = zeros(Float64, nr_of_steps, 1 + length(output_keys_lineplot))
    end

    if do_write_CSV  # csv export
        if do_write_all_CSV_outputs # gather all outputs
            output_keys_to_csv = all_output_keys
        else  # get only requested output keys from input file for CSV-export
            output_keys_to_csv = output_keys(components, project_config["io_settings"]["output_keys"])
        end
        reset_file(project_config["io_settings"]["output_file"], output_keys_to_csv)
    end

    ### prepare array for output of all energy flow of all system interfaces for Sankey
    # get number of system interfaces for preallocation and medium, source and target of
    # each interface for sankey diagram
    nr_of_interfaces = 0
    medium_of_interfaces = []
    output_sourcenames_sankey = []
    output_targetnames_sankey = []
    for each_component in components
        for each_outputinterface in each_component[2].output_interfaces
            medium = nothing  # reset medium
            if isa(each_outputinterface, Pair) # some output_interfaces are wrapped in a Touple
                medium = each_outputinterface[1] # then, the medium is stored separately
                each_outputinterface = each_outputinterface[2]
            end

            if isdefined(each_outputinterface, :target)
                # count interface
                nr_of_interfaces += 1

                #get name of source and sink
                push!(output_sourcenames_sankey, each_outputinterface.source.uac)
                push!(output_targetnames_sankey, each_outputinterface.target.uac)

                # get name of medium
                if isdefined(each_outputinterface.target, :medium)
                    push!(medium_of_interfaces, each_outputinterface.target.medium)
                elseif isdefined(each_outputinterface.source, :medium)
                    push!(medium_of_interfaces, each_outputinterface.source.medium)
                elseif !(medium === nothing)
                    push!(medium_of_interfaces, medium)
                else
                    println("Warning: The name of the medium was not detected. This may lead to wrong colouring in Sankey plot.")
                end
            end
        end
    end
    println(
        length(medium_of_interfaces) !== nr_of_interfaces
        ? "Warning: error in extracting information from input file for sankey plot."
        : ""
    )
    # preallocate for speed: Matrix with data of interfaces in every timestep
    output_interface_values = zeros(Float64, nr_of_steps, nr_of_interfaces)

    # export order or operatin (OoO)
    if project_config["io_settings"]["dump_info"]
        dump_info(
            project_config["io_settings"]["dump_info_file"],
            components, step_order, parameters
        )
    end

    for steps = 1:nr_of_steps
        # perform the simulation
        perform_steps(components, step_order, parameters)

        # check if any component was not balanced
        warnings = check_balances(components, parameters["epsilon"])
        if length(warnings) > 0
            print("Time is $(parameters["time"])\n")
            for (key, balance) in warnings
                print("Warning: Balance for component $key was not zero: $balance\n")
            end
        end

        # write requested output data of the components to CSV-file
        # This is currently done in every time step to keep data even if 
        # an error occurs.
        if do_write_CSV
            write_to_file(
                            project_config["io_settings"]["output_file"],
                            output_keys_to_csv,
                            parameters["time"]
                          )
        end

        # get the energy transported through each interface in every timestep for Sankey
        # if the balance of an interface was not zero, the actual energy that was flowing
        # is written to the outputs.
        # Attention: This can lead to overfilling of demands which is currenlty not visible
        # in the sankey diagram!
        n = 1
        for each_component in components
            for each_outputinterface in each_component[2].output_interfaces
                if isa(each_outputinterface, Pair) # some output_interfaces are wrapped in a Touple
                    if isdefined(each_outputinterface[2], :target)
                        output_interface_values[steps, n] = calculate_energy_flow(each_outputinterface[2])  
                        n += 1
                    end
                elseif isdefined(each_outputinterface, :target)
                    output_interface_values[steps, n] = calculate_energy_flow(each_outputinterface) 
                    n += 1
                end
            end
        end

        # gather output data of each component for line plot
        if do_create_plot
            output_data_lineplot[steps, :] = geather_output_data(output_keys_lineplot, parameters["time"])
        end

        # simulation update
        parameters["time"] += Int(parameters["time_step_seconds"])
    end

    ### create profile line plot
    if do_create_plot
        create_profile_line_plots(  output_data_lineplot,
                                    output_keys_lineplot,
                                    do_plot_all_outputs,
                                    do_plot_all_outputs ? nothing : project_config["io_settings"]["output_plot"]
                                 )
    end

    ### create Sankey diagram
    create_sankey(output_sourcenames_sankey, output_targetnames_sankey, output_interface_values, medium_of_interfaces, nr_of_interfaces)

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
            println("Could not parse project config file")
            return
        end
    end

    if project_config === nothing
        println("Could not find or parse project config file")
        return
    end

    run_simulation(project_config)
end

end # module
