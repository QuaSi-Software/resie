"""
Implementation of an electrolyser, turning electricity and water into H2, O2 and heat.

For the moment this remains a simple implementation that converts electricity into
the gases and heat (as medium m_h_w_ht1) at a defined ratio of 1:0.6:0.4. Has a minimum
run time of 3600s taken into consideration in its control behaviour and a minimum power
fraction of 20%. The power is considered the maximum amount of electricity that the
electrolyser can consume.

At the moment there is no operation strategy is implemented and the production of the
electrolyser is controlled by the demand it is linked to requires.
"""
mutable struct Electrolyser <: ControlledSystem
    uac::String
    controller::Controller
    sys_function::SystemFunction

    input_interfaces::InterfaceMap
    output_interfaces::InterfaceMap

    m_el_in::Symbol
    m_heat_out::Symbol
    m_h2_out::Symbol
    m_o2_out::Symbol

    power::Float64
    heat_fraction::Float64
    min_power_fraction::Float64
    min_run_time::UInt
    output_temperature::Temperature

    function Electrolyser(uac::String, config::Dict{String,Any})
        m_el_in = Symbol(default(config, "m_el_in", "m_e_ac_230v"))
        m_heat_out = Symbol(default(config, "m_heat_out", "m_h_w_lt1"))
        m_h2_out = Symbol(default(config, "m_h2_out", "m_c_g_h2"))
        m_o2_out = Symbol(default(config, "m_o2_out", "m_c_g_o2"))
        register_media([m_el_in, m_heat_out, m_h2_out, m_o2_out])

        return new(
            uac, # uac
            controller_for_strategy( # controller
                config["strategy"]["name"], config["strategy"]
            ),
            sf_transformer, # sys_function
            InterfaceMap( # input_interfaces
                m_el_in => nothing
            ),
            InterfaceMap( # output_interfaces
                m_heat_out => nothing,
                m_h2_out => nothing,
                m_o2_out => nothing,
            ),
            m_el_in,
            m_heat_out,
            m_h2_out,
            m_o2_out,
            config["power"], # power
            default(config, "heat_fraction", 0.4),
            default(config, "min_power_fraction", 0.2),
            default(config, "min_run_time", 3600),
            default(config, "output_temperature", 55.0),
        )
    end
end

function produce(unit::Electrolyser, parameters::Dict{String,Any}, watt_to_wh::Function)
    # get maximum energy demand and supply of electrolyser, not regarding external bounds:
    max_produce_heat = watt_to_wh(unit.power * unit.heat_fraction)
    max_produce_h2 = watt_to_wh(unit.power * (1.0 - unit.heat_fraction))
    max_produce_o2 = 0.5 * max_produce_h2  # @TODO: handle O2 calculation if it ever becomes relevant. for now use molar ratio
    max_consume_el = watt_to_wh(unit.power)

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

    # hydrogen
    if unit.controller.parameter["m_h2_out"] == true   
        InterfaceInfo = balance_on(
            unit.output_interfaces[unit.m_h2_out],
            unit.output_interfaces[unit.m_h2_out].target
        )
        potential_energy_h2 = InterfaceInfo.balance + InterfaceInfo.energy_potential
        potential_storage_h2 = InterfaceInfo.storage_potential
        if (unit.controller.parameter["load_storages"] ? potential_energy_h2 + potential_storage_h2 : potential_energy_h2) >= -parameters["epsilon"]
            return  # don't add to a surplus of hydrogen
        end
    else # unlimited demand in interface is assumed
        potential_energy_h2 = -Inf
        potential_storage_h2 = -Inf
    end

    # oxygen
    if unit.controller.parameter["m_o2_out"] == true 
        InterfaceInfo = balance_on(
            unit.output_interfaces[unit.m_o2_out],
            unit.output_interfaces[unit.m_o2_out].target
        )
        potential_energy_o2 = InterfaceInfo.balance + InterfaceInfo.energy_potential
        potential_storage_o2 = InterfaceInfo.storage_potential
        if (unit.controller.parameter["load_storages"] ? potential_energy_o2 + potential_storage_o2 : potential_energy_o2) >= -parameters["epsilon"]
            return  # don't add to a surplus of oxygen
        end
    else # unlimited demand in interface is assumed
        potential_energy_o2 = -Inf
        potential_storage_o2 = -Inf
    end

    # heat
    if unit.controller.parameter["m_heat_out"] == true  
        InterfaceInfo = balance_on(
            unit.output_interfaces[unit.m_heat_out],
            unit.output_interfaces[unit.m_heat_out].target
        )
        potential_energy_heat = InterfaceInfo.balance + InterfaceInfo.energy_potential
        potential_storage_heat = InterfaceInfo.storage_potential
        if (unit.controller.parameter["load_storages"] ? potential_energy_heat + potential_storage_heat : potential_energy_heat) >= 0
            return  # don't add to a surplus of heat
        end
    else # unlimited supply in interface is assumed
        potential_energy_heat = -Inf
        potential_storage_heat = -Inf
    end

    # get usage fraction of external profile (normalized from 0 to 1)
    usage_fraction_operation_profile = unit.controller.parameter["operation_profile_path"] === nothing ? 1.0 : value_at_time(unit.controller.parameter["operation_profile"], parameters["time"])
    if usage_fraction_operation_profile <= 0.0
        return # no operation allowed from external profile
    end

    # get usage_factions depending on control strategy
    if unit.controller.strategy == "storage_driven" && unit.controller.state_machine.state == 2
        usage_fraction_el = +(unit.controller.parameter["unload_storages"] ? potential_energy_el + potential_storage_el : potential_energy_el) / max_consume_el 
        usage_fraction_h2 = -(unit.controller.parameter["load_storages"] ? potential_energy_h2 + potential_storage_h2 : potential_energy_h2) / max_produce_h2 
        usage_fraction_o2 = -(unit.controller.parameter["load_storages"] ? potential_energy_o2 + potential_storage_o2 : potential_energy_o2) / max_produce_o2 
        usage_fraction_heat = -(unit.controller.parameter["load_storages"] ? potential_energy_heat + potential_storage_heat : potential_energy_heat) / max_produce_heat

    elseif unit.controller.strategy == "storage_driven" 
        return # do not start due to statemachine!
    
    elseif unit.controller.strategy == "supply_driven"
        usage_fraction_el = +(unit.controller.parameter["unload_storages"] ? potential_energy_el + potential_storage_el : potential_energy_el) / max_consume_el 
        usage_fraction_h2 = -(unit.controller.parameter["load_storages"] ? potential_energy_h2 + potential_storage_h2 : potential_energy_h2) / max_produce_h2 
        usage_fraction_o2 = -(unit.controller.parameter["load_storages"] ? potential_energy_o2 + potential_storage_o2 : potential_energy_o2) / max_produce_o2 
        usage_fraction_heat = -(unit.controller.parameter["load_storages"] ? potential_energy_heat + potential_storage_heat : potential_energy_heat) / max_produce_heat

    elseif unit.controller.strategy == "demand_driven"
        usage_fraction_el = +(unit.controller.parameter["unload_storages"] ? potential_energy_el + potential_storage_el : potential_energy_el) / max_consume_el 
        usage_fraction_h2 = -(unit.controller.parameter["load_storages"] ? potential_energy_h2 + potential_storage_h2 : potential_energy_h2) / max_produce_h2 
        usage_fraction_o2 = -(unit.controller.parameter["load_storages"] ? potential_energy_o2 + potential_storage_o2 : potential_energy_o2) / max_produce_o2 
        usage_fraction_heat = -(unit.controller.parameter["load_storages"] ? potential_energy_heat + potential_storage_heat : potential_energy_heat) / max_produce_heat

    else
        throw(ArgumentError("Error: No valid control strategy chosen for electrolyser. Must be one of storage_driven, supply_driven, demand_driven."))
    end

    # get smallest usage fraction
    usage_fraction = min(
        1.0, 
        usage_fraction_el,
        usage_fraction_h2, 
        usage_fraction_o2,
        usage_fraction_heat,
        usage_fraction_operation_profile
        )

    # exit if usage_fraction is below min_power_fraciton 
    if usage_fraction < unit.min_power_fraction
        return
    end

    # write production and demand in interfaces
    add!(unit.output_interfaces[unit.m_h2_out], max_produce_h2 * usage_fraction)
    add!(unit.output_interfaces[unit.m_o2_out], max_produce_o2 * usage_fraction)
    add!(
        unit.output_interfaces[unit.m_heat_out],
        max_produce_heat * usage_fraction,
        unit.output_temperature
    )
    sub!(unit.input_interfaces[unit.m_el_in], max_consume_el * usage_fraction)

end

export Electrolyser