"""
Implementation of a gas boiler producing heat from chemical energy in gaseous form.

The only currently implemented operation strategy involves checking the load of a linked
buffer tank and en-/disabling the boiler when a threshold is reached, in addition to an
overfill shutoff condition.
"""
mutable struct GasBoiler <: ControlledSystem
    uac :: String
    controller :: Controller
    sys_function :: SystemFunction

    input_interfaces :: InterfaceMap
    output_interfaces :: InterfaceMap

    power :: Float64
    min_power_fraction :: Float64

    function GasBoiler(uac :: String, config :: Dict{String, Any})
        if config["strategy"]["name"] == "Storage-driven"
            strategy = config["strategy"]["name"]

            machine = StateMachine(
                state=UInt(1),
                state_names=Dict{UInt, String}(
                    1 => "Off",
                    2 => "Load",
                ),
                time_in_state=UInt(0),
                transitions=Dict{UInt, TruthTable}(
                    1 => TruthTable( # State: Off
                        conditions=[
                            Condition(
                                "Buffer < X%",
                                Dict{String, Any}(
                                    "percentage" => config["strategy"]["low_threshold"]
                                )
                            ),
                        ],
                        table_data=Dict{Tuple, UInt}(
                            (false,) => 1,
                            (true,) => 2,
                        )
                    ),

                    2 => TruthTable( # State: Load
                        conditions=[
                            Condition(
                                "Buffer >= X%",
                                Dict{String, Any}(
                                    "percentage" => config["strategy"]["high_threshold"]
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
            strategy = "Default"
            machine = StateMachine()
        end

        return new(
            uac, # uac
            Controller(strategy, machine), # controller
            sf_transformer, # sys_function
            InterfaceMap( # input_interfaces
                m_c_g_natgas => nothing
            ),
            InterfaceMap( # output_interfaces
                m_h_w_ht1 => nothing
            ),
            config["power"], # power
            0.1, # min_power_fraction
        )
    end
end

function produce(unit :: GasBoiler, parameters :: Dict{String, Any}, watt_to_wh :: Function)
    if unit.controller.state_machine.state == 2
        max_produce_h = watt_to_wh(unit.power)

        balance, potential = balance_on(
            unit.output_interfaces[m_h_w_ht1],
            unit.output_interfaces[m_h_w_ht1].target
        )
        if balance + potential >= 0.0
            return # don't add to a surplus of energy
        end

        usage_fraction = min(1.0, abs(balance + potential) / max_produce_h)
        if usage_fraction < unit.min_power_fraction
            return
        end

        add!(unit.output_interfaces[m_h_w_ht1], max_produce_h * usage_fraction)
        sub!(unit.input_interfaces[m_c_g_natgas], watt_to_wh(unit.power * usage_fraction))
    end
end

export GasBoiler