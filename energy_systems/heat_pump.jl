Base.@kwdef mutable struct HeatPump <: ControlledSystem
    controller :: StateMachine

    input_interfaces :: InterfaceMap
    output_interfaces :: InterfaceMap

    last_consumed_e :: Float64 = 0.0
    last_produced_h :: Float64 = 0.0

    power :: Float64
    min_power_fraction :: Float64 = 0.2
    cop :: Float64
end

function make_HeatPump(strategy :: String, power :: Float64, cop :: Float64) :: HeatPump
    if strategy == "Ensure storage"
        return HeatPump(
            StateMachine( # HeatPump.controller
                state=UInt(1),
                state_names=Dict{UInt, String}(
                    1 => "Off",
                    2 => "Load"
                ),
                time_in_state=UInt(0),
                transitions=Dict{UInt, TruthTable}(
                    1 => TruthTable( # State: Off
                        conditions=[
                            Condition(
                                "Buffer < X%",
                                Dict{String, Any}(
                                    "percentage" => 0.1
                                )
                            ),
                        ],
                        table_data=Dict{Tuple, UInt}(
                            (true,) => 2,
                            (false,) => 1
                        )
                    ),

                    2 => TruthTable( # State: Load
                        conditions=[
                            Condition(
                                "Buffer >= X%",
                                Dict{String, Any}(
                                    "percentage" => 0.5
                                )
                            ),
                            Condition(
                                "Would overfill thermal buffer",
                                Dict{String, Any}()
                            ),
                        ],
                        table_data=Dict{Tuple, UInt}(
                            (false, false) => 2,
                            (false, true) => 1,
                            (true, false) => 1,
                            (true, true) => 1,
                        )
                    ),
                )
            ),
            InterfaceMap( # input_interfaces
                m_e_ac_230v => nothing
            ),
            InterfaceMap( # output_interfaces
                m_h_w_60c => nothing
            ),
            0.0, # last_consumed_e
            0.0, # last_produced_h
            power, # power
            0.2, # min_power_fraction
            cop, # electricity_fraction
        )
    else
        return HeatPump(controller=StateMachine(), power=power, cop=cop)
    end
end

function specific_values(unit :: HeatPump, time :: Int) :: Vector{Tuple}
    return [
        ("Consumption E", "$(unit.last_consumed_e)"),
        ("Production H", "$(unit.last_produced_h)")
    ]
end

export HeatPump, make_HeatPump, specific_values