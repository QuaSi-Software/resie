module Bran
"""
The time step, in seconds, used by the simulation.

@TODO: Move this into the input parameters to make it customizable at runtime.
"""
const TIME_STEP = UInt(900)

include("energy_systems/base.jl")

using .EnergySystems
import JSON

"""
    output_keys(from_config)

Transform the output keys definition in the project config file into a list of OutputKey
items. This is done to speed up selection of values for the output in each time step,
as this transformation has to be done only once at the beginning.
"""
function output_keys(
    systems :: Grouping,
    from_config :: Dict{String, Any}
) :: Vector{EnergySystems.OutputKey}
    outputs = Vector{EnergySystems.OutputKey}()

    for unit_key in keys(from_config)
        unit = systems[unit_key]

        for entry in from_config[unit_key]
            splitted = split(String(entry))
            if length(splitted) > 1
                medium_key = splitted[1]
                medium = getproperty(EnergySystems, Symbol(String(medium_key)))
                value_key = splitted[2]
            else
                medium = nothing
                value_key = splitted[1]
            end

            push!(outputs, EnergySystems.OutputKey(
                unit=unit,
                medium=medium,
                value_key=value_key
            ))
        end
    end

    return outputs
end

"""
    read_JSON(filepath)

Read and parse the JSON-encoded Dict in the given file.
"""
function read_JSON(filepath :: String) :: Dict{AbstractString, Any}
    open(filepath, "r") do file_handle
        content = read(file_handle, String)
        return JSON.parse(content)
    end
end

"""
    reset_file(filepath, output_keys)

Reset the output file and add headers for the given outputs.
"""
function reset_file(
    filepath :: String,
    output_keys :: Vector{EnergySystems.OutputKey}
)
    open(abspath(filepath), "w") do file_handle
        write(file_handle, "Time [s]")

        for outkey in output_keys
            if outkey.medium === nothing
                header = "$(outkey.unit.uac) $(outkey.value_key)"
            else
                header = "$(outkey.unit.uac) $(outkey.medium) $(outkey.value_key)"
            end
            write(file_handle, ";$header")
        end

        write(file_handle, "\n")
    end
end

"""
    write_to_file(filepath, output_keys, time)

Write the given outputs for the given time to file.
"""
function write_to_file(
    filepath :: String,
    output_keys :: Vector{EnergySystems.OutputKey},
    time :: Int
)
    open(abspath(filepath), "a") do file_handle
        write(file_handle, "$time")

        for outkey in output_keys
            value = output_value(outkey.unit, outkey)
            value = replace("$value", "." => ",")
            write(file_handle, ";$value")
        end

        write(file_handle, "\n")
    end
end

"""
    load_systems(config)

Construct instances of energy systems from the given config.

The config must have the structure:
```
{
    "UAC key": {
        "type": "PVPlant",
        ...
    },
    ...
}
```

The required parameters to construct an energy system from one entry in the config must
match what is required for the particular system. The `type` parameter must be present and
must match the symbol of the energy system class exactly. The structure is described in
more detail in the accompanying documentation on the project file.
"""
function load_systems(config :: Dict{String, Any}) :: Grouping
    systems = Grouping()

    for (unit_key, entry) in pairs(config)
        symbol = Symbol(String(entry["type"]))
        unit_class = getproperty(EnergySystems, symbol)
        if unit_class <: EnergySystems.EnergySystem
            instance = unit_class(unit_key, entry)
            systems[unit_key] = instance
        end
    end

    for (unit_key, entry) in pairs(config)
        if length(entry["control_refs"]) > 0
            others = Grouping(key => systems[key] for key in entry["control_refs"])
            link_control_with(systems[unit_key], others)
        end

        if length(entry["production_refs"]) > 0
            others = Grouping(key => systems[key] for key in entry["production_refs"])
            link_production_with(systems[unit_key], others)
        end
    end

    return systems
end

"""
    dump_info(file_path, systems, order_of_steps, parameters)

Dump a bunch of information to file that might be useful to explain the result of a run.

This is mostly used for debugging and development purposes, but might prove useful in
general to find out why the systems behave in the simulation as they do.
"""
function dump_info(
    file_path :: String,
    systems :: Grouping,
    order_of_steps :: StepInstructions,
    parameters :: Dict{String, Any}
)
    open(abspath(file_path), "w") do file_handle
        write(file_handle, "# Simulation step order\n")

        for entry in order_of_steps
            for step in entry[2:lastindex(entry)]
                write(file_handle, "1. `$(entry[1]) $(entry[2])`\n")
            end
        end
    end
end

"""
    order_of_steps(systems)

Calculate the order of steps that need to be performed to simulate the given systems.

This function works by an algorithm described in more detail in the accompanying
documentation. The result of this are step-by-step instructions telling the simulation
engine in which order the system functions are performed for each unit. This algorithm is
not trivial and might not work for each possible grouping of systems.

# Args
- `system::Grouping`: The systems for which an order is required
# Returns
- `StepInstructions`: The order in the structure:
    ```
    [
        ["UAC Key", s_step],
        ...
    ]
    ```
"""
function order_of_steps(systems :: Grouping) :: StepInstructions
    systems_by_function = [
        [unit for unit in each(systems)
            if unit.sys_function == EnergySystems.sf_fixed_source],
        [unit for unit in each(systems)
            if unit.sys_function == EnergySystems.sf_fixed_sink],
        [unit for unit in each(systems)
            if unit.sys_function == EnergySystems.sf_bus],
        [unit for unit in each(systems)
            if unit.sys_function == EnergySystems.sf_transformer],
        [unit for unit in each(systems)
            if unit.sys_function == EnergySystems.sf_storage],
        [unit for unit in each(systems)
            if unit.sys_function == EnergySystems.sf_dispatchable_source],
        [unit for unit in each(systems)
            if unit.sys_function == EnergySystems.sf_dispatchable_sink],
    ]

    simulation_order = []
    initial_nr = length(systems) * 7

    # reset all systems, order doesn't matter
    for sf_order = 1:7
        for unit in values(systems_by_function[sf_order])
            push!(simulation_order, (initial_nr, (unit.uac, EnergySystems.s_reset)))
            initial_nr -= 1
        end
    end

    # control/produce all systems except dispatchable sources/sinks. the order corresponds
    # to the general order of system functions
    for sf_order = 1:5
        for unit in values(systems_by_function[sf_order])
            push!(simulation_order, (initial_nr, (unit.uac, EnergySystems.s_control)))
            initial_nr -= 1
            push!(simulation_order, (initial_nr, (unit.uac, EnergySystems.s_produce)))
            initial_nr -= 1
        end
    end

    # sandwich loading of storage systems between the other systems and dispatchable ones
    for unit in values(systems_by_function[5])
        push!(simulation_order, (initial_nr, (unit.uac, EnergySystems.s_load)))
        initial_nr -= 1
    end

    # handle dispatchable sources/sinks
    for sf_order = 6:7
        for unit in values(systems_by_function[sf_order])
            push!(simulation_order, (initial_nr, (unit.uac, EnergySystems.s_control)))
            initial_nr -= 1
            push!(simulation_order, (initial_nr, (unit.uac, EnergySystems.s_produce)))
            initial_nr -= 1
        end
    end

    # finally, distribute bus systems
    for unit in values(systems_by_function[3])
        push!(simulation_order, (initial_nr, (unit.uac, EnergySystems.s_distribute)))
        initial_nr -= 1
    end

    fn_first = function(entry) return entry[1] end
    return [u[2] for u in sort(simulation_order, by=fn_first, rev=true)]
end

"""
    run_simulation()

Read inputs, perform the simulation calculation and write outputs.

Due to the complexity of required inputs of a simulation and how the outputs are persisted
(to file), this function takes only one argument, namely the project config, and returns
nothing.
"""
function run_simulation(project_config :: Dict{AbstractString, Any})
    systems = load_systems(project_config["energy_systems"])
    step_order = order_of_steps(systems)

    parameters = Dict{String, Any}(
        "time" => 0,
        "time_step_seconds" => TIME_STEP,
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
        parameters["time"] += Int(TIME_STEP)
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
