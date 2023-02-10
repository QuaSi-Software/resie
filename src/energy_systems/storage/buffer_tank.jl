"""
Implementation of a buffer tank holding hot water for heating or DHW purposes.

This is a simplified model, which mostly deals with amounts of energy and considers
temperatures only for the available temperature as the tank is depleted.
"""
mutable struct BufferTank <: ControlledSystem
    uac :: String
    controller :: Controller
    sys_function :: SystemFunction

    input_interfaces :: InterfaceMap
    output_interfaces :: InterfaceMap

    capacity :: Float64
    load :: Float64

    high_temperature :: Float64
    low_temperature :: Float64

    function BufferTank(uac :: String, config :: Dict{String, Any})
        return new(
            uac, # uac
            controller_for_strategy( # controller
                config["strategy"]["name"], config["strategy"]
            ),
            sf_storage, # sys_function
            InterfaceMap( # input_interfaces
                m_h_w_ht1 => nothing
            ),
            InterfaceMap( # output_interfaces
                m_h_w_ht1 => nothing
            ),
            config["capacity"], # capacity
            config["load"], # load
            "high_temperature" in keys(config) # high_temperature
                ? config["high_temperature"]
                : 75.0,
            "low_temperature" in keys(config) # low_temperature
                ? config["low_temperature"]
                : 20.0
        )
    end
end

function temperature_at_load(unit :: BufferTank) :: Temperature
    switch_point = 0.15
    partial_load = min(1.0, unit.load / (unit.capacity * switch_point))
    return (unit.high_temperature - unit.low_temperature) * partial_load + unit.low_temperature
end

function balance_on(
    interface :: SystemInterface,
    unit :: BufferTank
) :: Tuple{Float64, Float64, Temperature}
    return interface.balance, unit.capacity - unit.load, interface.temperature
end

function produce(unit :: BufferTank, parameters :: Dict{String, Any}, watt_to_wh :: Function)
    outface = unit.output_interfaces[m_h_w_ht1]
    balance, _, demand_temp = balance_on(outface, outface.target)

    if balance >= 0.0
        return # produce is only concerned with moving energy to the target
    end

    if demand_temp !== nothing && demand_temp > temperature_at_load(unit)
        return # we can only supply energy if it's at a higher temperature,
               # effectively reducing the tank's capacity for any demand at
               # a temperature higher than the lower limit of the tank
    end

    if unit.load > abs(balance)
        unit.load += balance
        add!(outface, abs(balance), demand_temp)
    else
        add!(outface, unit.load, demand_temp)
        unit.load = 0.0
    end
end

function load(unit :: BufferTank, parameters :: Dict{String, Any}, watt_to_wh :: Function)
    inface = unit.input_interfaces[m_h_w_ht1]
    balance, _, supply_temp = balance_on(inface, inface.source)

    if balance <= 0.0
        return # load is only concerned with receiving energy from the target
    end

    if supply_temp !== nothing && supply_temp < unit.low_temperature
        return # we can only take in energy if it's at a higher temperature than the
               # tank's lower limit
    end

    diff = unit.capacity - unit.load
    if diff > balance
        unit.load += balance
        sub!(inface, balance, supply_temp)
    else
        unit.load = unit.capacity
        sub!(inface, diff, supply_temp)
    end
end

function output_values(unit :: BufferTank) :: Vector{String}
    return ["IN", "OUT", "Load", "Capacity"]
end

function output_value(unit :: BufferTank, key :: OutputKey) :: Float64
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

export BufferTank