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
    _ = load_medien(("user_defined_media" in keys(project_config)
                     ? project_config["user_defined_media"]
                     : Array{Any}(undef, 0)))
    systems = load_systems(project_config["energy_systems"])
    step_order = order_of_operations(systems)

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

    outputs = output_keys(systems, project_config["io_settings"]["output_keys"])
    reset_file(project_config["io_settings"]["output_file"], outputs)

    ### set data for profile line plot
    # check if profile line plot should be created (checks if output_plot is present or if
    # it is "nothing")
    plot_bool = !(
        !(haskey(project_config["io_settings"], "output_plot"))
        ||
        project_config["io_settings"]["output_plot"] === "nothing"
    )

    if plot_bool
        # set keys for output plots    
        outputs_plot_keys = Vector{EnergySystems.OutputKey}()
        for plot in project_config["io_settings"]["output_plot"]
            append!(outputs_plot_keys, output_keys(systems, plot[2]["key"]))
        end

        # prepare array for output plots
        outputs_plot_data = zeros(Float64, nr_of_steps, 1 + length(outputs_plot_keys))
    end

    ### prepare array for output of all energy flow of all system interfaces
    # get number of system interfaces for preallocation and medium, source and target of
    # each interface for sankey diagram
    nr_of_interfaces = 0
    medium_of_interfaces = []
    output_all_sourcenames = []
    output_all_targetnames = []
    for each_system in systems
        for each_outputinterface in each_system[2].output_interfaces
            if isa(each_outputinterface, Pair) # some output_interfaces are wrapped in a Touple
                if isdefined(each_outputinterface[2], :target)
                    # count interface
                    nr_of_interfaces += 1

                    #get name of source and sink
                    push!(output_all_sourcenames, each_outputinterface[2].source.uac)
                    push!(output_all_targetnames, each_outputinterface[2].target.uac)

                    # get name of medium
                    if isdefined(each_outputinterface[2].target, :medium)
                        push!(medium_of_interfaces, each_outputinterface[2].target.medium)
                    elseif isdefined(each_outputinterface[2].source, :medium)
                        push!(medium_of_interfaces, each_outputinterface[2].source.medium)
                    else
                        push!(medium_of_interfaces, each_outputinterface[1])
                    end
                end
            elseif isdefined(each_outputinterface, :target)
                # count interface
                nr_of_interfaces += 1

                #get name of source and sink
                push!(output_all_sourcenames, each_outputinterface.source.uac)
                push!(output_all_targetnames, each_outputinterface.target.uac)

                # get name of medium
                if isdefined(each_outputinterface.target, :medium)
                    push!(medium_of_interfaces, each_outputinterface.target.medium)
                elseif isdefined(each_outputinterface.source, :medium)
                    push!(medium_of_interfaces, each_outputinterface.source.medium)
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
    output_all_values = zeros(Float64, nr_of_steps, nr_of_interfaces)

    if project_config["io_settings"]["dump_info"]
        dump_info(
            project_config["io_settings"]["dump_info_file"],
            systems, step_order, parameters
        )
    end

    for steps = 1:nr_of_steps
        # perform the simulation
        perform_steps(systems, step_order, parameters)

        # check if any energy system was not balanced
        warnings = check_balances(systems, parameters["epsilon"])
        if length(warnings) > 0
            print("Time is $(parameters["time"])\n")
            for (key, balance) in warnings
                print("Warning: Balance for system $key was not zero: $balance\n")
            end
        end

        # output to file
        write_to_file(
            project_config["io_settings"]["output_file"],
            outputs,
            parameters["time"]
        )

        # get all data of all interfaces in every timestep for Sankey
        # if the balance of an interface was not zero, the actual energy that was flowing
        # is written to the outputs.
        # Attention: This can lead to overfilling of demands which is currenlty not visible
        # in the sankey diagram!
        n = 1
        for each_system in systems
            for each_outputinterface in each_system[2].output_interfaces
                if isa(each_outputinterface, Pair) # some output_interfaces are wrapped in a Touple
                    if isdefined(each_outputinterface[2], :target)
                        output_all_values[steps, n] = (each_outputinterface[2].sum_abs_change + each_outputinterface[2].balance) / 2 #@ToDo: correct???
                        n += 1
                    end
                elseif isdefined(each_outputinterface, :target)
                    output_all_values[steps, n] = (each_outputinterface.sum_abs_change + each_outputinterface.balance) / 2 #@ToDo: correct???
                    n += 1
                end
            end
        end

        # gather data for profile line plot (@ToDo: may extract from all data in post processing)
        if plot_bool
            outputs_plot_data[steps, :] = geather_output_data(
                outputs_plot_keys,
                parameters["time"]
            )
        end

        # simulation update
        parameters["time"] += Int(parameters["time_step_seconds"])
    end

    ### create profile line plot
    if plot_bool
        create_profile_line_plots(outputs_plot_data, outputs_plot_keys, project_config["io_settings"]["output_plot"])
    end

    ### create Sankey diagram
    create_sankey(output_all_sourcenames, output_all_targetnames, output_all_values, medium_of_interfaces, nr_of_interfaces)

end

"""
    main()

Entry point into the simulation engine. The simulation is controlled and configured by a
config file.

# Command line arguments
## Positional arguments
- `String`: Filepath to the project config file (see documentation on file format). Can be
a path relative to the CWD of the caller.
"""
function main()
    if length(ARGS) > 0
        filepath = ARGS[1]
        if filepath !== nothing && filepath != ""
            project_config = nothing

            try
                project_config = read_JSON(abspath(filepath))
            catch exc
                if isa(exc, MethodError)
                    println("Could not parse project config file")
                    return
                end
            end

            run_simulation(project_config)
            return
        end

        println("Could not find or access project config file")
        return
    end

    println("No project config file argument given")
end
main()

end # module
