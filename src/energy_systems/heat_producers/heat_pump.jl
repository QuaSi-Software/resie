"""
Implementation of a heat pump energy system.

For the moment this remains a simple implementation that requires a low temperature heat
and electricity input and produces high temperature heat. Has a fixed coefficient of
performance (COP) of 3 and a minimum power fraction of 20%. The power parameter is
considered the maximum power of heat output the heat pump can produce.

The only currently implemented operation strategy involves checking the load of a linked
buffer tank and en-/disabling the heat pump when a threshold is reached, in addition to an
overfill shutoff condition.
"""
mutable struct HeatPump <: ControlledSystem
    uac :: String
    controller :: Controller
    sys_function :: SystemFunction

    input_interfaces :: InterfaceMap
    output_interfaces :: InterfaceMap

    power :: Float64
    min_power_fraction :: Float64
    min_run_time :: UInt
    cop :: Float64

    function HeatPump(uac :: String, config :: Dict{String, Any})
        return new(
            uac, # uac
            controller_for_strategy( # controller
                config["strategy"]["name"], config["strategy"]
            ),
            sf_transformer, # sys_function
            InterfaceMap( # input_interfaces
                m_h_w_lt1 => nothing,
                m_e_ac_230v => nothing
            ),
            InterfaceMap( # output_interfaces
                m_h_w_ht1 => nothing
            ),
            config["power"], # power
            0.2, # min_power_fraction
            0, # min_run_time
            3.0 # cop
        )
    end
end

function produce(unit :: HeatPump, parameters :: Dict{String, Any}, watt_to_wh :: Function)
    if unit.controller.strategy == "storage_driven" && unit.controller.state_machine.state == 2
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
        sub!(unit.input_interfaces[m_e_ac_230v], max_produce_h * usage_fraction / unit.cop)
        sub!(
            unit.input_interfaces[m_h_w_lt1],
            max_produce_h * usage_fraction * (1.0 - 1.0 / unit.cop)
        )

    elseif unit.controller.strategy == "supply_driven"
        balance, _ = balance_on(unit.input_interfaces[m_h_w_lt1], unit)

        if balance < parameters["epsilon"]
            return # do nothing if there is no heat to consume
        end

        max_consume_h = min(unit.power * (1.0 - 1.0 / unit.cop), balance)
        consume_e = max_consume_h / (unit.cop - 1.0)
        produce_h = max_consume_h + consume_e

        add!(unit.output_interfaces[m_h_w_ht1], produce_h)
        sub!(unit.input_interfaces[m_h_w_lt1], max_consume_h)
        sub!(unit.input_interfaces[m_e_ac_230v], consume_e)

    elseif unit.controller.strategy == "demand_driven"
        max_produce_h = watt_to_wh(unit.power)

        balance, _ = balance_on(
            unit.output_interfaces[m_h_w_ht1],
            unit.output_interfaces[m_h_w_ht1].target
        )
        if balance >= 0.0
            return # don't add to a surplus of energy
        end

        usage_fraction = min(1.0, abs(balance) / max_produce_h)
        if usage_fraction < unit.min_power_fraction
            return
        end

        add!(unit.output_interfaces[m_h_w_ht1], max_produce_h * usage_fraction)
        sub!(unit.input_interfaces[m_e_ac_230v], max_produce_h * usage_fraction / unit.cop)
        sub!(
            unit.input_interfaces[m_h_w_lt1],
            max_produce_h * usage_fraction * (1.0 - 1.0 / unit.cop)
        )
    end
end

export HeatPump