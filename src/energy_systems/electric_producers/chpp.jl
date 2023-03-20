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
    uac::String
    controller::Controller
    sys_function::SystemFunction

    input_interfaces::InterfaceMap
    output_interfaces::InterfaceMap

    m_gas_in::Symbol
    m_heat_out::Symbol
    m_el_out::Symbol

    power::Float64
    electricity_fraction::Float64
    min_power_fraction::Float64
    min_run_time::UInt

    function CHPP(uac::String, config::Dict{String,Any})
        m_gas_in = Symbol(default(config, "m_gas_in", "m_c_g_natgas"))
        m_heat_out = Symbol(default(config, "m_heat_out", "m_h_w_ht1"))
        m_el_out = Symbol(default(config, "m_el_out", "m_e_ac_230v"))
        register_media([m_gas_in, m_heat_out, m_el_out])

        return new(
            uac, # uac
            controller_for_strategy( # controller
                config["strategy"]["name"], config["strategy"]
            ),
            sf_transformer, # sys_function
            InterfaceMap( # input_interfaces
                m_gas_in => nothing
            ),
            InterfaceMap( # output_interfaces
                m_heat_out => nothing,
                m_el_out => nothing
            ),
            m_gas_in,
            m_heat_out,
            m_el_out,
            config["power"], # power
            default(config, "electricity_fraction", 0.4),
            default(config, "min_power_fraction", 0.2),
            default(config, "min_run_time", 1800),
        )
    end
end

function produce(unit::CHPP, parameters::Dict{String,Any}, watt_to_wh::Function)
    strategy = unit.controller.strategy
    if strategy == "storage_driven" && unit.controller.state_machine.state != 2
        return
    end

    max_produce_heat = watt_to_wh(unit.power * (1.0 - unit.electricity_fraction))
    max_produce_el = watt_to_wh(unit.power * unit.electricity_fraction)
    max_consume_gas = max_produce_heat + max_produce_el

    # get balance on in- and outputs, but only if they act as limitations (default: all are limiting, equals true)
    # Gas input 
    if unit.controller.parameter["m_gas_in"] == true 
        InterfaceInfo = balance_on(
            unit.input_interfaces[unit.m_gas_in],
            unit.input_interfaces[unit.m_gas_in].source
        )
        potential_energy_gas_in = InterfaceInfo.balance + InterfaceInfo.energy_potential
        potential_storage_gas_in = InterfaceInfo.storage_potential
        if (unit.controller.parameter["unload_storages"] ? potential_energy_gas_in + potential_storage_gas_in : potential_energy_gas_in) <= parameters["epsilon"]
            return # do nothing if there is no gas to consume
        end
    else # unlimited demand in interface is assumed
        potential_energy_gas_in = -Inf
        potential_storage_gas_in = -Inf
    end
   
    # electricity output 
    if unit.controller.parameter["m_el_out"] == true 
        InterfaceInfo = balance_on(
            unit.output_interfaces[unit.m_el_out],
            unit.output_interfaces[unit.m_el_out].target
        )
        potential_energy_el_out = InterfaceInfo.balance + InterfaceInfo.energy_potential
        potential_storage_el_out = InterfaceInfo.storage_potential
        if (unit.controller.parameter["load_storages"] ? potential_energy_el_out + potential_storage_el_out : potential_energy_el_out) >= -parameters["epsilon"]
            return # don't add to a surplus of electricity
        end
    else # unlimited demand in interface is assumed
        potential_energy_el_out = Inf
        potential_storage_el_out = Inf
    end

    # heat output
    if unit.controller.parameter["m_heat_out"] == true 
        InterfaceInfo = balance_on(
            unit.output_interfaces[unit.m_heat_out],
            unit.output_interfaces[unit.m_heat_out].target
        )
        potential_energy_heat_out = InterfaceInfo.balance + InterfaceInfo.energy_potential
        potential_storage_heat_out = InterfaceInfo.storage_potential
        if (unit.controller.parameter["load_storages"] ? potential_energy_heat_out + potential_storage_heat_out : potential_energy_heat_out) >= -parameters["epsilon"]
            return # don't add to a surplus of heat
        end
    else # unlimited demand in interface is assumed
        potential_energy_heat_out = Inf
        potential_storage_heat_out = Inf
    end
   
    # calculate usage fractions
    usage_fraction_heat_out = -((unit.controller.parameter["load_storages"] ? potential_energy_heat_out + potential_storage_heat_out : potential_energy_heat_out) / max_produce_heat)
    usage_fraction_el_out = -((unit.controller.parameter["unload_storages"] ? potential_energy_el_out + potential_storage_el_out : potential_energy_el_out) / max_produce_el)
    usage_fraction_gas_in = +((unit.controller.parameter["unload_storages"] ? potential_energy_gas_in + potential_storage_gas_in : potential_energy_gas_in) / max_consume_gas)
    other_limitations = 1.0

    # get smallest usage fraction
    usage_fraction = min(
        1.0, 
        usage_fraction_heat_out,
        usage_fraction_el_out, 
        usage_fraction_gas_in,
        other_limitations
        )

    # exit if usage_fraction is below min_power_fraciton 
    if usage_fraction < unit.min_power_fraction
        return
    end

    add!(unit.output_interfaces[unit.m_el_out], max_produce_el * usage_fraction)
    add!(unit.output_interfaces[unit.m_heat_out], max_produce_heat * usage_fraction)
    sub!(unit.input_interfaces[unit.m_gas_in], max_consume_gas * usage_fraction)
end

export CHPP