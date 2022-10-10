Base.@kwdef mutable struct HeatPump <: ControlledSystem
    controller :: StateMachine
    is_storage :: Bool

    input_interfaces :: InterfaceMap
    output_interfaces :: InterfaceMap

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
            false, # is_storage
            InterfaceMap( # input_interfaces
                m_e_ac_230v => nothing
            ),
            InterfaceMap( # output_interfaces
                m_h_w_60c => nothing
            ),
            power, # power
            0.2, # min_power_fraction
            cop, # electricity_fraction
        )
    else
        return HeatPump(controller=StateMachine(), power=power, cop=cop)
    end
end

function produce(unit :: HeatPump, parameters :: Dict{String, Any}, watt_to_wh :: Function)
    if unit.controller.state == 2
        max_produce_h = watt_to_wh(unit.power)

        usage_fraction = 1.0 # @TODO: implement partial load depending on space in buffer

        add!(unit.output_interfaces[m_h_w_60c], max_produce_h * usage_fraction)
        sub!(unit.input_interfaces[m_e_ac_230v], max_produce_h * usage_fraction / unit.cop)
    end
end

function specific_values(unit :: HeatPump, time :: Int) :: Vector{Tuple}
    return []
end

export HeatPump, make_HeatPump, specific_values