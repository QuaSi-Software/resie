"""
Implementation of a seasonal thermal storage system.

This is a simplified model, which mostly deals with amounts of energy and considers
temperatures only for the available temperature as the tank is depleted.
"""
mutable struct SeasonalThermalStorage <: ControlledSystem
    uac::String
    controller::Controller
    sys_function::SystemFunction

    input_interfaces::InterfaceMap
    output_interfaces::InterfaceMap

    m_heat_in::Symbol
    m_heat_out::Symbol

    capacity::Float64
    load::Float64

    use_adaptive_temperature::Bool
    switch_point::Float64
    high_temperature::Float64
    low_temperature::Float64

    function SeasonalThermalStorage(uac::String, config::Dict{String,Any})
        m_heat_in = Symbol(default(config, "m_heat_in", "m_h_w_ht1"))
        m_heat_out = Symbol(default(config, "m_heat_out", "m_h_w_lt1"))
        register_media([m_heat_in, m_heat_out])

        return new(
            uac, # uac
            controller_for_strategy( # controller
                config["strategy"]["name"], config["strategy"]
            ),
            sf_storage, # sys_function
            InterfaceMap( # input_interfaces
                m_heat_in => nothing
            ),
            InterfaceMap( # output_interfaces
                m_heat_out => nothing
            ),
            m_heat_in,
            m_heat_out,
            config["capacity"], # capacity
            config["load"], # load
            default(config, "use_adaptive_temperature", false),
            default(config, "switch_point", 0.25),
            default(config, "high_temperature", 90.0),
            default(config, "low_temperature", 15.0),
        )
    end
end

function temperature_at_load(unit::SeasonalThermalStorage)::Temperature
    if unit.use_adaptive_temperature
        partial_load = min(1.0, unit.load / (unit.capacity * unit.switch_point))
        return (unit.high_temperature - unit.low_temperature) * partial_load + unit.low_temperature
    else
        return unit.high_temperature
    end
end

function balance_on(
    interface::SystemInterface,
    unit::SeasonalThermalStorage
)::NamedTuple{}

    caller_is_input = false   # ==true if interface is input of unit (caller puts energy in unit); 
                              # ==false if interface is output of unit (caller gets energy from unit)
    
    # check if caller is input or output of unit
    for (_, input_uac) in pairs(unit.input_interfaces)
        if input_uac == interface.source.uac
            caller_is_input = true
            break
        end
        if input_uac.source.uac == interface.source.uac
            caller_is_input = true
            break
        end
    end

    return (
            balance = interface.balance,
            storage_potential = caller_is_input ? -(unit.capacity-unit.load) : unit.load,
            energy_potential = 0.0,
            temperature = interface.temperature
            )
end

function produce(unit::SeasonalThermalStorage, parameters::Dict{String,Any}, watt_to_wh::Function)
    outface = unit.output_interfaces[unit.m_heat_out]
    InterfaceInfo = balance_on(outface, outface.target)
    demand_temp = InterfaceInfo.temperature

    if InterfaceInfo.balance >= 0.0
        return # produce is only concerned with moving energy to the target
    end

    if demand_temp !== nothing && demand_temp > temperature_at_load(unit)
        return # we can only supply energy if it's at a higher temperature,
        # effectively reducing the tank's capacity for any demand at
        # a temperature higher than the lower limit of the tank
    end

    if unit.load > abs(InterfaceInfo.balance)
        unit.load += InterfaceInfo.balance
        add!(outface, abs(InterfaceInfo.balance), demand_temp)
    else
        add!(outface, unit.load, demand_temp)
        unit.load = 0.0
    end
end

function load(unit::SeasonalThermalStorage, parameters::Dict{String,Any}, watt_to_wh::Function)
    inface = unit.input_interfaces[unit.m_heat_in]
    InterfaceInfo = balance_on(inface, inface.source)
    supply_temp = InterfaceInfo.temperature

    if InterfaceInfo.balance <= 0.0
        return # load is only concerned with receiving energy from the target
    end

    if supply_temp !== nothing && supply_temp < unit.low_temperature
        return # we can only take in energy if it's at a higher temperature than the
        # tank's lower limit
    end

    diff = unit.capacity - unit.load
    if diff > InterfaceInfo.balance
        unit.load += InterfaceInfo.balance
        sub!(inface, InterfaceInfo.balance, supply_temp)
    else
        unit.load = unit.capacity
        sub!(inface, diff, supply_temp)
    end
end

function output_values(unit::SeasonalThermalStorage)::Vector{String}
    return ["IN", "OUT", "Load", "Capacity"]
end

function output_value(unit::SeasonalThermalStorage, key::OutputKey)::Float64
    if key.value_key == "IN"
        return unit.input_interfaces[key.medium].sum_abs_change * 0.5
    elseif key.value_key == "OUT"
        return unit.output_interfaces[key.medium].sum_abs_change * 0.5
    elseif key.value_key == "Load"
        return unit.load
    elseif key.value_key == "Capacity"
        return unit.capacity
    end
    throw(KeyError(key.value_key))
end

export SeasonalThermalStorage