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
    controller :: Controller
    sys_function :: SystemFunction

    input_interfaces :: InterfaceMap
    output_interfaces :: InterfaceMap

    power :: Float64
    electricity_fraction :: Float64
    min_power_fraction :: Float64
    min_run_time :: UInt

    function CHPP(uac :: String, config :: Dict{String, Any})
        return new(
            uac, # uac
            controller_for_strategy( # controller
                config["strategy"]["name"], config["strategy"]
            ),
            sf_transformer, # sys_function
            InterfaceMap( # input_interfaces
                m_c_g_natgas => nothing
            ),
            InterfaceMap( # output_interfaces
                m_h_w_ht1 => nothing,
                m_e_ac_230v => nothing
            ),
            config["power"], # power
            "electricity_fraction" in keys(config) # electricity_fraction
                ? config["electricity_fraction"]
                : 0.4,
            "min_power_fraction" in keys(config) # min_power_fraction
                ? config["min_power_fraction"]
                : 0.2,
            "min_run_time" in keys(config) # min_run_time
                ? config["min_run_time"]
                : 1800,
        )
    end
end

function produce(unit :: CHPP, parameters :: Dict{String, Any}, watt_to_wh :: Function)
    strategy = unit.controller.strategy
    if strategy == "storage_driven" && unit.controller.state_machine.state != 2
        return
    end

    max_produce_h = watt_to_wh(unit.power * (1.0 - unit.electricity_fraction))
    max_produce_e = watt_to_wh(unit.power * unit.electricity_fraction)

    balance, potential, _ = balance_on(
        unit.output_interfaces[m_h_w_ht1],
        unit.output_interfaces[m_h_w_ht1].target
    )

    demand_to_meet = (
        strategy == "storage_driven"
        ? balance + potential
        : balance
    )

    if demand_to_meet >= 0.0
        return
    end

    usage_fraction = min(1.0, abs(balance + potential) / max_produce_h)
    if usage_fraction < unit.min_power_fraction
        return
    end

    add!(unit.output_interfaces[m_e_ac_230v], max_produce_e * usage_fraction)
    add!(unit.output_interfaces[m_h_w_ht1], max_produce_h * usage_fraction)
    sub!(unit.input_interfaces[m_c_g_natgas], watt_to_wh(unit.power * usage_fraction))
end

export CHPP