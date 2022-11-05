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
    end

    link_production_with(
        systems["TST_01_HZG_01_GRI"],
        Grouping("TST_01_HZG_01_CHP" => systems["TST_01_HZG_01_CHP"])
    )
    link_production_with(
        systems["TST_01_ELT_01_GRI"],
        Grouping("TST_01_ELT_01_BUS" => systems["TST_01_ELT_01_BUS"])
    )
    link_production_with(
        systems["TST_01_HZG_01_BFT"],
        Grouping("TST_01_HZG_01_BUS" => systems["TST_01_HZG_01_BUS"])
    )
    link_production_with(
        systems["TST_01_ELT_01_BAT"],
        Grouping("TST_01_ELT_01_BUS" => systems["TST_01_ELT_01_BUS"])
    )
    link_production_with(
        systems["TST_01_HZG_01_CHP"],
        Grouping(
            "TST_01_HZG_01_BUS" => systems["TST_01_HZG_01_BUS"],
            "TST_01_ELT_01_BUS" => systems["TST_01_ELT_01_BUS"]
        )
    )
    link_production_with(
        systems["TST_01_HZG_01_HTP"],
        Grouping("TST_01_HZG_01_BUS" => systems["TST_01_HZG_01_BUS"])
    )
    link_production_with(
        systems["TST_01_ELT_01_PVP"],
        Grouping("TST_01_ELT_01_BUS" => systems["TST_01_ELT_01_BUS"])
    )
    link_production_with(
        systems["TST_01_HZG_01_BUS"],
        Grouping(
            "TST_01_HZG_01_DEM" => systems["TST_01_HZG_01_DEM"],
            "TST_01_HZG_01_BFT" => systems["TST_01_HZG_01_BFT"]
        )
    )
    link_production_with(
        systems["TST_01_ELT_01_BUS"],
        Grouping(
            "TST_01_ELT_01_DEM" => systems["TST_01_ELT_01_DEM"],
            "TST_01_HZG_01_HTP" => systems["TST_01_HZG_01_HTP"],
            "TST_01_ELT_01_BAT" => systems["TST_01_ELT_01_BAT"],
            "TST_01_ELT_01_GRO" => systems["TST_01_ELT_01_GRO"]
        )
    )

    return systems
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

    simulation_order = [
        ["TST_01_ELT_01_PVP", EnergySystems.s_reset], # limited_source
        ["TST_01_HZG_01_DEM", EnergySystems.s_reset], # limited_sink
        ["TST_01_ELT_01_DEM", EnergySystems.s_reset], # limited_sink
        ["TST_01_HZG_01_BUS", EnergySystems.s_reset], # bus
        ["TST_01_ELT_01_BUS", EnergySystems.s_reset], # bus
        ["TST_01_HZG_01_CHP", EnergySystems.s_reset], # transformer
        ["TST_01_HZG_01_HTP", EnergySystems.s_reset], # transformer
        ["TST_01_HZG_01_BFT", EnergySystems.s_reset], # storage
        ["TST_01_ELT_01_BAT", EnergySystems.s_reset], # storage
        ["TST_01_HZG_01_GRI", EnergySystems.s_reset], # infinite_source
        ["TST_01_ELT_01_GRI", EnergySystems.s_reset], # infinite_source
        ["TST_01_ELT_01_GRO", EnergySystems.s_reset], # infinite_sink
        ["TST_01_ELT_01_PVP", EnergySystems.s_control, EnergySystems.s_produce], # limited_source
        ["TST_01_HZG_01_DEM", EnergySystems.s_control, EnergySystems.s_produce], # limited_sink
        ["TST_01_ELT_01_DEM", EnergySystems.s_control, EnergySystems.s_produce], # limited_sink
        ["TST_01_HZG_01_BUS", EnergySystems.s_control, EnergySystems.s_produce], # bus
        ["TST_01_ELT_01_BUS", EnergySystems.s_control, EnergySystems.s_produce], # bus
        ["TST_01_HZG_01_CHP", EnergySystems.s_control, EnergySystems.s_produce], # transformer
        ["TST_01_HZG_01_HTP", EnergySystems.s_control, EnergySystems.s_produce], # transformer
        ["TST_01_HZG_01_BFT", EnergySystems.s_control, EnergySystems.s_produce], # storage
        ["TST_01_ELT_01_BAT", EnergySystems.s_control, EnergySystems.s_produce], # storage
        ["TST_01_HZG_01_BFT", EnergySystems.s_load], # storage
        ["TST_01_ELT_01_BAT", EnergySystems.s_load], # storage
        ["TST_01_HZG_01_GRI", EnergySystems.s_control, EnergySystems.s_produce], # infinite_source
        ["TST_01_ELT_01_GRI", EnergySystems.s_control, EnergySystems.s_produce], # infinite_source
        ["TST_01_ELT_01_GRO", EnergySystems.s_control, EnergySystems.s_produce], # infinite_sink
        ["TST_01_HZG_01_BUS", EnergySystems.s_distribute], # bus
        ["TST_01_ELT_01_BUS", EnergySystems.s_distribute], # bus
    ]

    parameters = Dict{String, Any}(
        "time" => 0,
        "time_step_seconds" => TIME_STEP,
        "epsilon" => 1e-9
    )

    outputs = output_keys(systems, project_config["io_settings"]["output_keys"])
    reset_file(project_config["io_settings"]["output_file"], outputs)

    for i = 1:(96*7)
        # perform the simulation
        perform_steps(systems, simulation_order, parameters)

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
