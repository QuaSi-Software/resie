"""
Implementation of a gas boiler producing heat from chemical energy in gaseous form.

The only currently implemented operation strategy involves checking the load of a linked
buffer tank and en-/disabling the boiler when a threshold is reached, in addition to an
overfill shutoff condition.
"""
mutable struct GasBoiler <: ControlledSystem
    uac::String
    controller::Controller
    sys_function::SystemFunction

    input_interfaces::InterfaceMap
    output_interfaces::InterfaceMap

    power::Float64
    min_power_fraction::Float64
    min_run_time::UInt

    function GasBoiler(uac::String, config::Dict{String,Any})
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
                m_h_w_ht1 => nothing
            ),
            config["power"], # power
            "min_power_fraction" in keys(config) # min_power_fraction
                ? config["min_power_fraction"]
                : 0.1,
            "min_run_time" in keys(config) # min_run_time
                ? config["min_run_time"]
                : 0,
        )
    end
end

function produce(unit::GasBoiler, parameters::Dict{String,Any}, watt_to_wh::Function)
    strategy = unit.controller.strategy
    if strategy == "storage_driven" && unit.controller.state_machine.state != 2
        return
    end

    max_produce_h = watt_to_wh(unit.power)

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

    usage_fraction = min(1.0, abs(demand_to_meet) / max_produce_h)
    if usage_fraction < unit.min_power_fraction
        return
    end

    add!(unit.output_interfaces[m_h_w_ht1], max_produce_h * usage_fraction)
    sub!(unit.input_interfaces[m_c_g_natgas], watt_to_wh(unit.power * usage_fraction))
end

export GasBoiler