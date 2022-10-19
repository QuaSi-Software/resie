Base.@kwdef mutable struct CHPP <: ControlledSystem
    controller :: StateMachine
    sys_function :: SystemFunction

    input_interfaces :: InterfaceMap
    output_interfaces :: InterfaceMap

    power :: Float64
    electricity_fraction :: Float64 = 0.4
    min_power_fraction :: Float64
    min_run_time :: UInt = 1800
end

function make_CHPP(strategy :: String, power :: Float64) :: CHPP
    if strategy == "Ensure storage"
        return CHPP(
            StateMachine( # CHPP.controller
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
                                    "percentage" => 0.2
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
                                    "percentage" => 0.9
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
            ),
            transformer, # sys_function
            InterfaceMap( # input_interfaces
                m_c_g_natgas => nothing
            ),
            InterfaceMap( # output_interfaces
                m_h_w_60c => nothing,
                m_e_ac_230v => nothing
            ),
            power, # power
            0.4, # electricity_fraction
            0.2, # min_power_fraction
            1800 # min_run_time
        )
    else
        return CHPP(controller=StateMachine(), power=power)
    end
end

function produce(unit :: CHPP, parameters :: Dict{String, Any}, watt_to_wh :: Function)
    if unit.controller.state == 2
        max_produce_h = watt_to_wh(unit.power * (1.0 - unit.electricity_fraction))
        max_produce_e = watt_to_wh(unit.power * unit.electricity_fraction)

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

        add!(unit.output_interfaces[m_e_ac_230v], max_produce_e * usage_fraction)
        add!(unit.output_interfaces[m_h_w_60c], max_produce_h * usage_fraction)
        sub!(unit.input_interfaces[m_c_g_natgas], watt_to_wh(unit.power * usage_fraction))
    end
end

function specific_values(unit :: CHPP, time :: Int) :: Vector{Tuple}
    return []
end

export CHPP, make_CHPP, specific_values