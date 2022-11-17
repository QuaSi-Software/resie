module Bran

# note: includes that contain their own module, which have to be submodules of the Bran
# module, are included first, then can be accessed with the "using" keyword. files that
# contain code that is intended to be used in-place of their include statement (as part
# of the Bran module), are included after the "using" statements have been declared.
# this is done so the latter files can access the symbols of the submodules the same as
# if the code was inside this file.

include("profiles/base.jl")
using .Profiles

include("energy_systems/base.jl")
using .EnergySystems

include("project_loading.jl")
include("file_output.jl")

"""
    run_simulation()

Read inputs, perform the simulation calculation and write outputs.

Due to the complexity of required inputs of a simulation and how the outputs are persisted
(to file), this function takes only one argument, namely the project config, and returns
nothing.
"""
function run_simulation(project_config :: Dict{AbstractString, Any})
    systems = load_systems(project_config["energy_systems"])
    step_order = order_of_steps(systems, project_config["energy_systems"])

    time_step = 900
    if "time_step_seconds" in keys(project_config["simulation_parameters"])
        time_step = UInt(project_config["simulation_parameters"]["time_step_seconds"])
    end

    parameters = Dict{String, Any}(
        "time" => 0,
        "time_step_seconds" => time_step,
        "epsilon" => 1e-9
    )

    outputs = output_keys(systems, project_config["io_settings"]["output_keys"])
    reset_file(project_config["io_settings"]["output_file"], outputs)

    if project_config["io_settings"]["dump_info"]
        dump_info(
            project_config["io_settings"]["dump_info_file"],
            systems, step_order, parameters
        )
    end

    for i = 1:(96*7)
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

        # output
        write_to_file(
            project_config["io_settings"]["output_file"],
            outputs,
            parameters["time"]
        )

        # simulation update
        parameters["time"] += Int(parameters["time_step_seconds"])
    end
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
