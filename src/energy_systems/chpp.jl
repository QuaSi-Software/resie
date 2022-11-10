"""
Implementation of a combined-heat-power-plant (CHPP) energy system.

For the moment this remains a simple implementation that converts natural gas into
electricity and heat (as medium m_h_w_ht1) at a defined ratio of 1:0.4:0.6. Has a minimum
run time of 1800s taken into consideration in its control behaviour and a minimum power
fraction of 20%. The power is considered the maximum amount of both heat and electricity
that the CHPP can produce.

The only currently implemented operation strategy involves checking the load of a linked
buffer tank and en-/disabling the CHPP when a threshold is reached, in addition to an
overfill shutoff condition.
"""
mutable struct CHPP <: ControlledSystem
    uac :: String
    controller :: StateMachine
    sys_function :: SystemFunction

    input_interfaces :: InterfaceMap
    output_interfaces :: InterfaceMap

    power :: Float64
    electricity_fraction :: Float64
    min_power_fraction :: Float64
    min_run_time :: UInt

    function CHPP(uac :: String, config :: Dict{String, Any})
        if config["strategy"]["name"] == "Ensure storage"
            controller = StateMachine(
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
                                "Min run time",
                                Dict{String, Any}()
                            ),
                            Condition(
                                "Would overfill thermal buffer",
                                Dict{String, Any}()
                            ),
                        ],
                        table_data=Dict{Tuple, UInt}(
                            (false, false, false) => 2,
                            (false, true, false) => 2,
                            (true, false, false) => 2,
                            (true, true, false) => 1,
                            (false, false, true) => 1,
                            (false, true, true) => 1,
                            (true, false, true) => 1,
                            (true, true, true) => 1,
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
            sf_transformer, # sys_function
            InterfaceMap( # input_interfaces
                m_c_g_natgas => nothing
            ),
            InterfaceMap( # output_interfaces
                m_h_w_ht1 => nothing,
                m_e_ac_230v => nothing
            ),
            config["power"], # power
            0.4, # electricity_fraction
            0.2, # min_power_fraction
            1800 # min_run_time
        )
    end
end

function produce(unit :: CHPP, parameters :: Dict{String, Any}, watt_to_wh :: Function)
    if unit.controller.state == 2
        max_produce_h = watt_to_wh(unit.power * (1.0 - unit.electricity_fraction))
        max_produce_e = watt_to_wh(unit.power * unit.electricity_fraction)

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

        add!(unit.output_interfaces[m_e_ac_230v], max_produce_e * usage_fraction)
        add!(unit.output_interfaces[m_h_w_ht1], max_produce_h * usage_fraction)
        sub!(unit.input_interfaces[m_c_g_natgas], watt_to_wh(unit.power * usage_fraction))
    end
end

export CHPP