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

    m_gas_in::Symbol
    m_heat_out::Symbol

    power::Float64
    min_power_fraction::Float64
    min_run_time::UInt

    function GasBoiler(uac::String, config::Dict{String,Any})
        m_gas_in = Symbol(default(config, "m_gas_in", "m_c_g_natgas"))
        m_heat_out = Symbol(default(config, "m_heat_out", "m_h_w_ht1"))
        register_media([m_gas_in, m_heat_out])

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
                m_heat_out => nothing
            ),
            m_gas_in,
            m_heat_out,
            config["power"], # power
            default(config, "min_power_fraction", 0.1),
            default(config, "min_run_time", 0),
        )
    end
end

function produce(unit::GasBoiler, parameters::Dict{String,Any}, watt_to_wh::Function)
    strategy = unit.controller.strategy
    if strategy == "storage_driven" && unit.controller.state_machine.state != 2
        return
    end
    
    max_produce_heat = watt_to_wh(unit.power)
    max_consume_gas = max_produce_heat 

    # get balance on in- and outputs, but only if they act as limitations (default: all are limiting, equals true)
    # Gas input 
    if unit.controller.parameter["m_gas_in"] == true 
        exchange = balance_on(
            unit.input_interfaces[unit.m_gas_in],
            unit.input_interfaces[unit.m_gas_in].source
        )
        potential_energy_gas_in = exchange.balance + exchange.energy_potential
        potential_storage_gas_in = exchange.storage_potential
        if (unit.controller.parameter["unload_storages"] ? potential_energy_gas_in + potential_storage_gas_in : potential_energy_gas_in) <= parameters["epsilon"]
            return # do nothing if there is no gas to consume
        end
    else # unlimited demand in interface is assumed
        potential_energy_gas_in = -Inf
        potential_storage_gas_in = -Inf
    end

    # heat output
    if unit.controller.parameter["m_heat_out"] == true 
        exchange = balance_on(
            unit.output_interfaces[unit.m_heat_out],
            unit.output_interfaces[unit.m_heat_out].target
        )
        potential_energy_heat_out = exchange.balance + exchange.energy_potential
        potential_storage_heat_out = exchange.storage_potential
        if (unit.controller.parameter["load_storages"] ? potential_energy_heat_out + potential_storage_heat_out : potential_energy_heat_out) >= -parameters["epsilon"]
            return # don't add to a surplus of heat
        end
    else # unlimited demand in interface is assumed
        potential_energy_heat_out = Inf
        potential_storage_heat_out = Inf
    end
   
    # get usage fraction of external profile (normalized from 0 to 1)
    usage_fraction_operation_profile = unit.controller.parameter["operation_profile_path"] === nothing ? 1.0 : value_at_time(unit.controller.parameter["operation_profile"], parameters["time"])
    if usage_fraction_operation_profile <= 0.0
        return # no operation allowed from external profile
    end
    
    # calculate usage fractions
    usage_fraction_heat_out = -((unit.controller.parameter["load_storages"] ? potential_energy_heat_out + potential_storage_heat_out : potential_energy_heat_out) / max_produce_heat)
    usage_fraction_gas_in = +((unit.controller.parameter["unload_storages"] ? potential_energy_gas_in + potential_storage_gas_in : potential_energy_gas_in) / max_consume_gas)

    # get smallest usage fraction
    usage_fraction = min(
        1.0, 
        usage_fraction_heat_out,
        usage_fraction_gas_in,
        usage_fraction_operation_profile
        )

    if usage_fraction < unit.min_power_fraction
        return
    end

    add!(unit.output_interfaces[unit.m_heat_out], max_produce_heat * usage_fraction)
    sub!(unit.input_interfaces[unit.m_gas_in], max_consume_gas * usage_fraction)
end

export GasBoiler