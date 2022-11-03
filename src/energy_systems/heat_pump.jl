"""
Implementation of a heat pump energy system.

For the moment this remains a simple implementation that requires no heat input and
produces heat of medium m_h_w_60c from electricity. Has a fixed coefficient of performance
(COP) of 3 and a minimum power fraction of 20%. The power parameters is considered the
maximum power of heat output the heat pump can produce.

The only currently implemented operation strategy involves checking the load of a linked
buffer tank and en-/disabling the heat pump when a threshold is reached, in addition to an
overfill shutoff condition.
"""
mutable struct HeatPump <: ControlledSystem
    uac :: String
    controller :: StateMachine
    sys_function :: SystemFunction

    input_interfaces :: InterfaceMap
    output_interfaces :: InterfaceMap

    power :: Float64
    min_power_fraction :: Float64
    cop :: Float64

    function HeatPump(uac :: String, config :: Dict{String, Any})
        if config["strategy"] == "Ensure storage"
            controller = StateMachine(
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
        )
        else
            controller = StateMachine()
        end

        return new(
            uac, # uac
            controller, # controller
            transformer, # sys_function
            InterfaceMap( # input_interfaces
                m_e_ac_230v => nothing
            ),
            InterfaceMap( # output_interfaces
                m_h_w_60c => nothing
            ),
            config["power"], # power
            0.2, # min_power_fraction
            3.0 # cop
        )
    end
end

function produce(unit :: HeatPump, parameters :: Dict{String, Any}, watt_to_wh :: Function)
    if unit.controller.state == 2
        max_produce_h = watt_to_wh(unit.power)

        balance, potential = balance_on(
            unit.output_interfaces[m_h_w_60c],
            unit.output_interfaces[m_h_w_60c].target
        )
        if balance + potential >= 0.0
            return # don't add to a surplus of energy
        end

        usage_fraction = min(1.0, abs(balance + potential) / max_produce_h)
        if usage_fraction < unit.min_power_fraction
            return
        end

        add!(unit.output_interfaces[m_h_w_60c], max_produce_h * usage_fraction)
        sub!(unit.input_interfaces[m_e_ac_230v], max_produce_h * usage_fraction / unit.cop)
    end
end

export HeatPump, output_values, output_value