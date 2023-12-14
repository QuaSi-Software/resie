"""
Implementation of a buffer tank holding hot water for heating or DHW purposes.

This is a simplified model, which mostly deals with amounts of energy and considers
temperatures only for the available temperature as the tank is depleted.
"""
mutable struct BufferTank <: Component
    uac::String
    controller::Controller
    sys_function::SystemFunction

    input_interfaces::InterfaceMap
    output_interfaces::InterfaceMap

    medium::Symbol

    capacity::Float64
    load::Float64

    use_adaptive_temperature::Bool
    switch_point::Float64
    high_temperature::Float64
    low_temperature::Float64

    function BufferTank(uac::String, config::Dict{String,Any})
        medium = Symbol(default(config, "medium", "m_h_w_ht1"))
        register_media([medium])

        return new(
            uac, # uac
            controller_for_strategy( # controller
                config["strategy"]["name"], config["strategy"]
            ),
            sf_storage, # sys_function
            InterfaceMap( # input_interfaces
                medium => nothing
            ),
            InterfaceMap( # output_interfaces
                medium => nothing
            ),
            medium,
            config["capacity"], # capacity
            config["load"], # load
            default(config, "use_adaptive_temperature", false),
            default(config, "switch_point", 0.15),
            default(config, "high_temperature", 75.0),
            default(config, "low_temperature", 20),
        )
    end
end

function control(
    unit::BufferTank,
    components::Grouping,
    parameters::Dict{String,Any}
)
    move_state(unit, components, parameters)

    # revise! This is not correct anymore!
    unit.output_interfaces[unit.medium].temperature = highest(temperature_at_load(unit), unit.output_interfaces[unit.medium].temperature)
    unit.input_interfaces[unit.medium].temperature = highest(unit.high_temperature, unit.input_interfaces[unit.medium].temperature)

end

function temperature_at_load(unit::BufferTank)::Temperature
    if unit.use_adaptive_temperature
        partial_load = min(1.0, unit.load / (unit.capacity * unit.switch_point))
        return (unit.high_temperature - unit.low_temperature) * partial_load + unit.low_temperature
    else
        return unit.high_temperature
    end
end

function balance_on(
    interface::SystemInterface,
    unit::BufferTank
)::Vector{EnergyExchange}
    caller_is_input = unit.uac == interface.target.uac

    return [EnEx(
        balance=interface.balance,
        uac=unit.uac,
        energy_potential=0.0,
        storage_potential=caller_is_input ? -(unit.capacity - unit.load) : unit.load,
        temperature=interface.temperature,
        pressure=nothing,
        voltage=nothing,
    )]
end

function process(unit::BufferTank, parameters::Dict{String,Any})
    outface = unit.output_interfaces[unit.medium]
    exchanges = balance_on(outface, outface.target)
    demand_temp = temperature_first(exchanges)

    if unit.controller.parameter["name"] == "default"
        energy_demand = balance(exchanges)
    elseif unit.controller.parameter["name"] == "extended_storage_control"
        if unit.controller.parameter["load_any_storage"]
            energy_demand = balance(exchanges) + storage_potential(exchanges)
        else
            energy_demand = balance(exchanges)
        end
    else
        energy_demand = balance(exchanges)
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

function load(unit::BufferTank, parameters::Dict{String,Any})
    inface = unit.input_interfaces[unit.medium]
    exchanges = balance_on(inface, inface.source)
    supply_temp = temperature_first(exchanges)
    energy_available = balance(exchanges)

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

function output_values(unit::BufferTank)::Vector{String}
    return ["IN", "OUT", "Load", "Load%", "Capacity"]
end

function output_value(unit::BufferTank, key::OutputKey)::Float64
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
    end
    throw(KeyError(key.value_key))
end

export BufferTank