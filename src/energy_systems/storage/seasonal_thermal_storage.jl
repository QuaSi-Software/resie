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

    function SeasonalThermalStorage(uac::String, config::Dict{String,Any}, sim_params::Dict{String,Any})
        m_heat_in = Symbol(default(config, "m_heat_in", "m_h_w_ht1"))
        m_heat_out = Symbol(default(config, "m_heat_out", "m_h_w_lt1"))
        register_media([m_heat_in, m_heat_out])

        return new(
            uac, # uac
            controller_for_strategy( # controller
                config["strategy"]["name"], config["strategy"], sim_params
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
    sim_params::Dict{String,Any}
)
    move_state(unit, components, sim_params)

    if unit.output_interfaces[unit.m_heat_out].temperature === nothing
        set_temperature!(unit.output_interfaces[unit.m_heat_out], temperature_at_load(unit))
    end
    if unit.input_interfaces[unit.m_heat_in].temperature === nothing
        set_temperature!(unit.input_interfaces[unit.m_heat_in], unit.high_temperature)
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

function process(unit::SeasonalThermalStorage, sim_params::Dict{String,Any})
    outface = unit.output_interfaces[unit.m_heat_out]
    exchanges = balance_on(outface, outface.target)
    energy_demanded = balance(exchanges)

    if (
        unit.controller.parameter["name"] == "extended_storage_control"
        && unit.controller.parameter["load_any_storage"]
    )
        energy_demanded += storage_potential(exchanges)
    end

    # shortcut if there is no energy demanded
    if energy_demanded >= -sim_params["epsilon"]
        return
    end

    for exchange in exchanges
        demanded_on_interface = exchange.balance
        if (
            unit.controller.parameter["name"] == "extended_storage_control"
            && unit.controller.parameter["load_any_storage"]
        )
            demanded_on_interface += exchange.storage_potential
        end

        if demanded_on_interface >= -sim_params["epsilon"]
            continue
        end

        demand_temp = exchange.temperature
        if demand_temp !== nothing && demand_temp > temperature_at_load(unit)
            # we can only supply energy at a temperature at or below the tank's current
            # output temperature
            continue
        end

        used_heat = min(abs(energy_demanded), abs(demanded_on_interface))

        if unit.load > used_heat
            unit.load -= used_heat
            add!(outface, used_heat, temperature_at_load(unit))
            energy_demanded += used_heat
        else
            add!(outface, unit.load, temperature_at_load(unit))
            energy_demanded += unit.load
            unit.load = 0.0
        end
    end
end

function load(unit::SeasonalThermalStorage, sim_params::Dict{String,Any})
    inface = unit.input_interfaces[unit.m_heat_in]
    exchanges = balance_on(inface, inface.source)
    energy_available = balance(exchanges)

    # shortcut if there is no energy to be used
    if energy_available <= sim_params["epsilon"]
        return
    end

    for exchange in exchanges
        if exchange.balance < sim_params["epsilon"]
            continue
        end

        supply_temp = exchange.temperature
        if supply_temp !== nothing && supply_temp < unit.high_temperature
            # we can only take in energy if it's at a higher/equal temperature than the
            # storage's upper limit for temperatures
            continue
        end

        used_heat = min(energy_available, exchange.balance)
        diff = unit.capacity - unit.load

        if diff > used_heat
            unit.load += used_heat
            sub!(inface, used_heat, unit.high_temperature)
            energy_available -= used_heat
        else
            unit.load = unit.capacity
            sub!(inface, diff, unit.high_temperature)
            energy_available -= diff
        end
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