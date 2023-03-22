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
    return ["IN", "OUT", "COP"]
end

function output_value(unit::HeatPump, key::OutputKey)::Float64
    if key.value_key == "IN"
        return calculate_energy_flow(unit.input_interfaces[key.medium])
    elseif key.value_key == "OUT"
        return calculate_energy_flow(unit.output_interfaces[key.medium])
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
    # get balance on in- and outputs, but only if they act as limitations (default: all are limiting, equals true)
    # electricity 
    if unit.controller.parameter["m_el_in"] == true 
        InterfaceInfo = balance_on(
            unit.input_interfaces[unit.m_el_in],
            unit.input_interfaces[unit.m_el_in].source
        )
        potential_energy_el = InterfaceInfo.balance + InterfaceInfo.energy_potential
        potential_storage_el = InterfaceInfo.storage_potential
        if (unit.controller.parameter["unload_storages"] ? potential_energy_el + potential_storage_el : potential_energy_el) <= parameters["epsilon"]
            return # do nothing if there is no electricity to consume
        end
    else # unlimited demand in interface is assumed
        potential_energy_el = Inf
        potential_storage_el = Inf
    end
   
    # heat in 
    if unit.controller.parameter["m_heat_in"] == true 
        InterfaceInfo = balance_on(
            unit.input_interfaces[unit.m_heat_in],
            unit.input_interfaces[unit.m_heat_in].source
        )
        potential_energy_heat_in = InterfaceInfo.balance + InterfaceInfo.energy_potential
        potential_storage_heat_in = InterfaceInfo.storage_potential
        in_temp = InterfaceInfo.temperature
        if (unit.controller.parameter["unload_storages"] ? potential_energy_heat_in + potential_storage_heat_in : potential_energy_heat_in) <= parameters["epsilon"]
            return # do nothing if there is no heat to consume
        end
    else # unlimited demand in interface is assumed
        potential_energy_heat_in = Inf
        potential_storage_heat_in = Inf
    end

    # heat out
    if unit.controller.parameter["m_heat_out"] == true 
        InterfaceInfo = balance_on(
            unit.output_interfaces[unit.m_heat_out],
            unit.output_interfaces[unit.m_heat_out].target
        )
        potential_energy_heat_out = InterfaceInfo.balance + InterfaceInfo.energy_potential
        potential_storage_heat_out = InterfaceInfo.storage_potential
        out_temp = InterfaceInfo.temperature
        if (unit.controller.parameter["load_storages"] ? potential_energy_heat_out + potential_storage_heat_out : potential_energy_heat_out) >= -parameters["epsilon"]
            return # don't add to a surplus of heat
        end
    else # unlimited demand in interface is assumed
        potential_energy_heat_out = Inf
        potential_storage_heat_out = Inf
    end
   
    # check if temperature has already been read from input and output interface
    if !(@isdefined in_temp)
        InterfaceInfo = balance_on(
            unit.input_interfaces[unit.m_heat_in],
            unit
        )
        in_temp = InterfaceInfo.temperature
    end

    if !(@isdefined out_temp)
        InterfaceInfo = balance_on(
            unit.output_interfaces[unit.m_heat_out],
            unit.output_interfaces[unit.m_heat_out].target
        )
        out_temp = InterfaceInfo.temperature
    end

    # calculate COP
    cop = dynamic_cop(in_temp, out_temp)
    unit.cop = cop === nothing ? unit.fixed_cop : cop

    # maximum possible in and outputs of heat pump, not regarding any external limits!
    max_produce_heat = watt_to_wh(unit.power)
    max_consume_heat = max_produce_heat * (1.0 - 1.0 / unit.cop)
    max_consume_el = max_produce_heat - max_consume_heat

    if unit.controller.strategy == "storage_driven" && unit.controller.state_machine.state == 2

        usage_fraction_heat_out = -((unit.controller.parameter["load_storages"] ? potential_energy_heat_out + potential_storage_heat_out : potential_energy_heat_out) / max_produce_heat)
        usage_fraction_heat_in = +((unit.controller.parameter["unload_storages"] ? potential_energy_heat_in + potential_storage_heat_in : potential_energy_heat_in) / max_consume_heat)
        usage_fraction_el = +((unit.controller.parameter["unload_storages"] ? potential_energy_el + potential_storage_el : potential_energy_el) / max_consume_el)
        other_limitations = 1.0

    elseif unit.controller.strategy == "storage_driven" 
        return # do not start due to statemachine!

    elseif unit.controller.strategy == "supply_driven"

        usage_fraction_heat_out = -((unit.controller.parameter["load_storages"] ? potential_energy_heat_out + potential_storage_heat_out : potential_energy_heat_out) / max_produce_heat)
        usage_fraction_heat_in = +((unit.controller.parameter["unload_storages"] ? potential_energy_heat_in + potential_storage_heat_in : potential_energy_heat_in) / max_consume_heat)
        usage_fraction_el = +((unit.controller.parameter["unload_storages"] ? potential_energy_el + potential_storage_el : potential_energy_el) / max_consume_el)
        other_limitations = 1.0

    elseif unit.controller.strategy == "demand_driven"

        usage_fraction_heat_out = -((unit.controller.parameter["load_storages"] ? potential_energy_heat_out + potential_storage_heat_out : potential_energy_heat_out) / max_produce_heat)
        usage_fraction_heat_in = +((unit.controller.parameter["unload_storages"] ? potential_energy_heat_in + potential_storage_heat_in : potential_energy_heat_in) / max_consume_heat)
        usage_fraction_el = +((unit.controller.parameter["unload_storages"] ? potential_energy_el + potential_storage_el : potential_energy_el) / max_consume_el)
        other_limitations = 1.0

    end

    # get smallest usage fraction
    usage_fraction = min(
        1.0, 
        usage_fraction_heat_out,
        usage_fraction_heat_in, 
        usage_fraction_el,
        other_limitations
        )

    # exit if usage_fraction is below min_power_fraciton 
    if usage_fraction < unit.min_power_fraction
        return
    end

    # write consumed and produced energy in interfaces
    add!(
        unit.output_interfaces[unit.m_heat_out],
        max_produce_heat * usage_fraction,
        out_temp
    )
    sub!(unit.input_interfaces[unit.m_el_in], max_produce_heat * usage_fraction / unit.cop)
    sub!(
        unit.input_interfaces[unit.m_heat_in],
        max_produce_heat * usage_fraction * (1.0 - 1.0 / unit.cop)
    )
end

export HeatPump