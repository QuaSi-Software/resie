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
Holds the options which output values should be recorded.

This is a specific data structure intended to speed up recording output by avoiding the
need to parse the user-submitted config options for every time step.
"""
Base.@kwdef struct OutputKey
    unit :: EnergySystem
    medium :: Union{Nothing, MediumCategory}
    value_key :: String
end

"""
    output_keys(from_config)

Transform the output keys definition in the project config file into a list of OutputKey
items. This is done to speed up selection of values for the output in each time step,
as this transformation has to be done only once at the beginning.
"""
function output_keys(
    systems :: Grouping,
    from_config :: Dict{String, Any}
) :: Vector{OutputKey}
    outputs = Vector{OutputKey}()

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

            push!(outputs, OutputKey(
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
    reset_file(systems)

Reset the output file and add headers for the given systems
"""
function reset_file(systems :: Grouping)
    open("./out.csv", "w") do file_handle
        write(file_handle, "Time [s]")

        for (key, unit) in pairs(systems)
            for val in specific_values(unit, Int(0))
                write(file_handle, ";$key $(val[1])")
            end

            if isa(unit, Bus) continue end

            for (medium, inface) in pairs(unit.input_interfaces)
                if inface === nothing continue end
                write(file_handle, ";$key $medium IN")
            end

            for (medium, outface) in pairs(unit.output_interfaces)
                if outface === nothing continue end
                write(file_handle, ";$key $medium OUT")
            end
        end

        write(file_handle, "\n")
    end
end

"""
    write_to_file(systems, time)

Write the energy transfer values and additional state to the output file.
"""
function write_to_file(systems :: Grouping, time :: Int)
    open("./out.csv", "a") do file_handle
        write(file_handle, "$time")

        for unit in each(systems)

            for val in specific_values(unit, time)
                write(file_handle, replace(
                    replace(";$(val[2])", "/" => ";"),
                    "." => ","
                ))
            end

            if isa(unit, Bus) continue end

            for inface in values(unit.input_interfaces)
                if inface === nothing continue end
                write(file_handle, replace(
                    replace(";$(inface.sum_abs_change * 0.5)", "/" => ";"),
                    "." => ","
                ))
            end

            for outface in values(unit.output_interfaces)
                if outface === nothing continue end
                write(file_handle, replace(
                    replace(";$(outface.sum_abs_change * 0.5)", "/" => ";"),
                    "." => ","
                ))
            end
        end

        write(file_handle, "\n")
    end
end

"""
    run_simulation()

Read inputs, perform the simulation calculation and write outputs.

This is the entry point to the simulation engine. Due to the complexity of required inputs
and how the outputs are written (to file), this function doesn't take any arguments and
returns nothing.
"""
function run_simulation(project_config :: Dict{AbstractString, Any})
    systems = Grouping(
        "TST_01_HZG_01_GRI" => make_GridConnection(EnergySystems.m_c_g_natgas, true),
        "TST_01_ELT_01_GRI" => make_GridConnection(EnergySystems.m_e_ac_230v, true),
        "TST_01_ELT_01_GRO" => make_GridConnection(EnergySystems.m_e_ac_230v, false),
        "TST_01_HZG_01_BFT" => make_BufferTank(40000.0, 20000.0),
        "TST_01_ELT_01_BAT" => make_Battery("Economical discharge", 10000.0, 5000.0),
        "TST_01_HZG_01_CHP" => make_CHPP("Ensure storage", 12500.0),
        "TST_01_HZG_01_HTP" => make_HeatPump("Ensure storage", 20000.0, 3.0),
        "TST_01_ELT_01_PVP" => make_PVPlant(15000.0),
        "TST_01_ELT_01_BUS" => make_Bus(EnergySystems.m_e_ac_230v),
        "TST_01_HZG_01_BUS" => make_Bus(EnergySystems.m_h_w_60c),
        "TST_01_HZG_01_DEM" => make_Demand(EnergySystems.m_h_w_60c, 10000.0),
        "TST_01_ELT_01_DEM" => make_Demand(EnergySystems.m_e_ac_230v, 15000.0),
    )

    for (key, unit) in pairs(systems)
        unit.uac = key
    end

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

    link_control_with(
        systems["TST_01_HZG_01_CHP"],
        Grouping("TST_01_HZG_01_BFT" => systems["TST_01_HZG_01_BFT"])
    )
    link_control_with(
        systems["TST_01_HZG_01_HTP"],
        Grouping("TST_01_HZG_01_BFT" => systems["TST_01_HZG_01_BFT"])
    )
    link_control_with(
        systems["TST_01_ELT_01_BAT"],
        Grouping("TST_01_ELT_01_PVP" => systems["TST_01_ELT_01_PVP"])
    )

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

    parameters = Dict{String, Any}(
        "time" => 0,
        "time_step_seconds" => TIME_STEP,
        "epsilon" => 1e-9
    )

    reset_file(systems)
    outputs = output_keys(systems, project_config["output_keys"])

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
        write_to_file(systems, parameters["time"])

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
