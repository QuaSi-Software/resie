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
    uac::String
    controller::Controller
    sys_function::SystemFunction

    input_interfaces::InterfaceMap
    output_interfaces::InterfaceMap

    m_el_in::Symbol
    m_heat_out::Symbol
    m_heat_in::Symbol

    power::Float64
    min_power_fraction::Float64
    min_run_time::UInt
    fixed_cop::Float64
    cop::Float64

    function HeatPump(uac::String, config::Dict{String,Any})
        m_el_in = Symbol(default(config, "m_el_in", "m_e_ac_230v"))
        m_heat_out = Symbol(default(config, "m_heat_out", "m_h_w_ht1"))
        m_heat_in = Symbol(default(config, "m_heat_in", "m_h_w_lt1"))
        register_media([m_el_in, m_heat_out, m_heat_in])

        return new(
            uac, # uac
            controller_for_strategy( # controller
                config["strategy"]["name"], config["strategy"]
            ),
            sf_transformer, # sys_function
            InterfaceMap( # input_interfaces
                m_heat_in => nothing,
                m_el_in => nothing
            ),
            InterfaceMap( # output_interfaces
                m_heat_out => nothing
            ),
            m_el_in,
            m_heat_out,
            m_heat_in,
            config["power"], # power
            default(config, "min_power_fraction", 0.2),
            default(config, "min_run_time", 0),
            default(config, "fixed_cop", 3.0),
            0.0, # cop
        )
    end
end

function output_values(unit::HeatPump)::Vector{String}
    return ["OUT", "Max_Power", "Temperature"]
end

function output_value(unit::HeatPump, key::OutputKey)::Float64
    if key.value_key == "IN"
        return unit.input_interfaces[key.medium].sum_abs_change * 0.5
    elseif key.value_key == "OUT"
        return unit.output_interfaces[key.medium].sum_abs_change * 0.5
    elseif key.value_key == "COP"
        return unit.cop
    end
    throw(KeyError(key.value_key))
end

function dynamic_cop(in_temp::Temperature, out_temp::Temperature)::Union{Nothing,Float64}
    if (in_temp === nothing || out_temp === nothing)
        return nothing
    end

    delta_t = out_temp - in_temp
    return 8.0 * exp(-0.08 * delta_t) + 1
end

function produce(unit::HeatPump, parameters::Dict{String,Any}, watt_to_wh::Function)
    in_blnc, _, in_temp = balance_on(
        unit.input_interfaces[unit.m_heat_in],
        unit
    )

    out_blnc, out_pot, out_temp = balance_on(
        unit.output_interfaces[unit.m_heat_out],
        unit.output_interfaces[unit.m_heat_out].target
    )

    cop = dynamic_cop(in_temp, out_temp)
    unit.cop = cop === nothing ? unit.fixed_cop : cop

    if unit.controller.strategy == "storage_driven" && unit.controller.state_machine.state == 2
        max_produce_h = watt_to_wh(unit.power)

        if out_blnc + out_pot >= 0.0
            return # don't add to a surplus of energy
        end

        usage_fraction = min(1.0, abs(out_blnc + out_pot) / max_produce_h)
        if usage_fraction < unit.min_power_fraction
            return
        end

        add!(unit.output_interfaces[unit.m_heat_out], max_produce_h * usage_fraction)
        sub!(unit.input_interfaces[unit.m_el_in], max_produce_h * usage_fraction / unit.cop)
        sub!(
            unit.input_interfaces[unit.m_heat_in],
            max_produce_h * usage_fraction * (1.0 - 1.0 / unit.cop)
        )

    elseif unit.controller.strategy == "supply_driven"
        if in_blnc < parameters["epsilon"]
            return # do nothing if there is no heat to consume
        end

        max_consume_h = min(unit.power * (1.0 - 1.0 / unit.cop), in_blnc)
        consume_e = max_consume_h / (unit.cop - 1.0)
        produce_h = max_consume_h + consume_e

        add!(unit.output_interfaces[unit.m_heat_out], produce_h)
        sub!(unit.input_interfaces[unit.m_heat_in], max_consume_h)
        sub!(unit.input_interfaces[unit.m_el_in], consume_e)

    elseif unit.controller.strategy == "demand_driven"
        max_produce_h = watt_to_wh(unit.power)

        if out_blnc >= 0.0
            return # don't add to a surplus of energy
        end

        usage_fraction = min(1.0, abs(out_blnc) / max_produce_h)
        if usage_fraction < unit.min_power_fraction
            return
        end

        add!(
            unit.output_interfaces[unit.m_heat_out],
            max_produce_h * usage_fraction,
            out_temp
        )
        sub!(unit.input_interfaces[unit.m_el_in], max_produce_h * usage_fraction / unit.cop)
        sub!(
            unit.input_interfaces[unit.m_heat_in],
            max_produce_h * usage_fraction * (1.0 - 1.0 / unit.cop)
        )
    end
end

export HeatPump