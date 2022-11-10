"""
Implementation of an electrolyser, turning electricity and water into H2, O2 and heat.

For the moment this remains a simple implementation that converts electricity into
the gases and heat (as medium m_h_w_60c) at a defined ratio of 1:0.6:0.4. Has a minimum
run time of 3600s taken into consideration in its control behaviour and a minimum power
fraction of 20%. The power is considered the maximum amount of electricity that the
electrolyser can consume.

At the moment there is no operation strategy is implemented and the production of the
electrolyser is controlled by the demand it is linked to requires.
"""
mutable struct Electrolyser <: ControlledSystem
    uac :: String
    controller :: StateMachine
    sys_function :: SystemFunction

    input_interfaces :: InterfaceMap
    output_interfaces :: InterfaceMap

    power :: Float64
    heat_fraction :: Float64
    min_power_fraction :: Float64
    min_run_time :: UInt

    function Electrolyser(uac :: String, config :: Dict{String, Any})
        return new(
            uac, # uac
            StateMachine(), # controller
            sf_transformer, # sys_function
            InterfaceMap( # input_interfaces
                m_e_ac_230v => nothing
            ),
            InterfaceMap( # output_interfaces
                m_h_w_60c => nothing,
                m_c_g_h2 => nothing,
                m_c_g_o2 => nothing
            ),
            config["power"], # power
            0.4, # heat_fraction
            0.2, # min_power_fraction
            3600 # min_run_time
        )
    end
end

function produce(unit :: Electrolyser, parameters :: Dict{String, Any}, watt_to_wh :: Function)
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

export Electrolyser