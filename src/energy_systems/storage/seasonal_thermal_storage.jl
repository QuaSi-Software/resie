"""
Implementation of a seasonal thermal storage component.

This is a simplified model, which mostly deals with amounts of energy and considers
temperatures only for the available temperature as the tank is depleted.
"""
mutable struct SeasonalThermalStorage <: Component
    uac::String
    controller::Controller
    sys_function::SystemFunction

    input_interfaces::InterfaceMap
    output_interfaces::InterfaceMap

    m_heat_in::Symbol
    m_heat_out::Symbol

    capacity::Float64
    load::Float64
    losses::Float64

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
            0.0, # losses
            default(config, "use_adaptive_temperature", false),
            default(config, "switch_point", 0.25),
            default(config, "high_temperature", 90.0),
            default(config, "low_temperature", 15.0),
        )
    end
end

function control(
    unit::SeasonalThermalStorage,
    components::Grouping,
    parameters::Dict{String,Any}
)
    move_state(unit, components, parameters)
    unit.output_interfaces[unit.m_heat_out].temperature = highest_temperature(temperature_at_load(unit), unit.output_interfaces[unit.m_heat_out].temperature)
    unit.input_interfaces[unit.m_heat_in].temperature = highest_temperature(unit.high_temperature, unit.input_interfaces[unit.m_heat_in].temperature)

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
    
    caller_is_input = unit.uac == interface.target.uac ? true : false
    # ==true if interface is input of unit (caller puts energy in unit); 
    # ==false if interface is output of unit (caller gets energy from unit)
    
    return (
            balance = interface.balance,
            storage_potential = caller_is_input ? -(unit.capacity-unit.load) : unit.load,
            energy_potential = 0.0,
            temperature = interface.temperature
            )
end

function process(unit::SeasonalThermalStorage, parameters::Dict{String,Any})
    outface = unit.output_interfaces[unit.m_heat_out]
    exchange = balance_on(outface, outface.target)
    demand_temp = exchange.temperature

    if unit.controller.parameter["name"] == "default"
        energy_demand = exchange.balance
    elseif unit.controller.parameter["name"] == "extended_storage_control"
        if unit.controller.parameter["load_any_storage"]
            energy_demand = exchange.balance + exchange.storage_potential
        else
            energy_demand = exchange.balance
        end
    else
        energy_demand = exchange.balance
    end

    if energy_demand >= 0.0
        return # process is only concerned with moving energy to the target
    end

    if demand_temp !== nothing && demand_temp > temperature_at_load(unit)
        return # we can only supply energy if it's at a higher temperature,
        # effectively reducing the tank's capacity for any demand at
        # a temperature higher than the lower limit of the tank
    end

    if unit.load > abs(energy_demand)
        unit.load += energy_demand
        add!(outface, abs(energy_demand), demand_temp)
    else
        add!(outface, unit.load, demand_temp)
        unit.load = 0.0
    end
end

function load(unit::SeasonalThermalStorage, parameters::Dict{String,Any})
    inface = unit.input_interfaces[unit.m_heat_in]
    exchange = balance_on(inface, inface.source)
    supply_temp = exchange.temperature
    energy_available = exchange.balance

    if energy_available <= 0.0
        return # load is only concerned with receiving energy from the target
    end

    if supply_temp !== nothing && supply_temp < unit.low_temperature
        return # we can only take in energy if it's at a higher temperature than the
        # tank's lower limit
    end

    diff = unit.capacity - unit.load
    if diff > energy_available
        unit.load += energy_available
        sub!(inface, energy_available, supply_temp)
    else
        unit.load = unit.capacity
        sub!(inface, diff, supply_temp)
    end
end

function output_values(unit::SeasonalThermalStorage)::Vector{String}
    return [string(unit.m_heat_in)*" IN",
            string(unit.m_heat_out)*" OUT",
            "Load",
            "Load%",
            "Capacity",
            "Losses"]
end

function output_value(unit::SeasonalThermalStorage, key::OutputKey)::Float64
    if key.value_key == "IN"
        return calculate_energy_flow(unit.input_interfaces[key.medium])
    elseif key.value_key == "OUT"
        return calculate_energy_flow(unit.output_interfaces[key.medium])
    elseif key.value_key == "Load"
        return unit.load
    elseif key.value_key == "Load%"
        return 100 * unit.load / unit.capacity
    elseif key.value_key == "Capacity"
        return unit.capacity
    elseif key.value_key == "Losses"
        return unit.losses
    end
    throw(KeyError(key.value_key))
end

export SeasonalThermalStorage