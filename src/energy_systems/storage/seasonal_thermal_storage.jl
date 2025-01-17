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
    load_end_of_last_timestep::Float64
    losses::Float64

    use_adaptive_temperature::Bool
    switch_point::Float64
    high_temperature::Float64
    low_temperature::Float64

    function SeasonalThermalStorage(uac::String, config::Dict{String,Any}, sim_params::Dict{String,Any})
        m_heat_in = Symbol(default(config, "m_heat_in", "m_h_w_ht1"))
        m_heat_out = Symbol(default(config, "m_heat_out", "m_h_w_lt1"))
        register_media([m_heat_in, m_heat_out])

        return new(uac, # uac
                   Controller(default(config, "control_parameters", nothing)),
                   sf_storage, # sys_function
                   InterfaceMap(m_heat_in => nothing),  # input_interfaces
                   InterfaceMap(m_heat_out => nothing), # output_interfaces
                   m_heat_in,
                   m_heat_out,
                   config["capacity"], # capacity
                   config["load"],     # load
                   0.0,                # load_end_of_last_timestep
                   0.0,                # losses
                   default(config, "use_adaptive_temperature", false),
                   default(config, "switch_point", 0.25),
                   default(config, "high_temperature", 90.0),
                   default(config, "low_temperature", 15.0))
    end
end

function initialise!(unit::SeasonalThermalStorage, sim_params::Dict{String,Any})
    set_storage_transfer!(unit.input_interfaces[unit.m_heat_in],
                          unload_storages(unit.controller, unit.m_heat_in))
    set_storage_transfer!(unit.output_interfaces[unit.m_heat_out],
                          load_storages(unit.controller, unit.m_heat_out))

    unit.load_end_of_last_timestep = copy(unit.load)
end

function control(unit::SeasonalThermalStorage,
                 components::Grouping,
                 sim_params::Dict{String,Any})
    update(unit.controller)

    set_temperature!(unit.output_interfaces[unit.m_heat_out],
                     nothing,
                     temperature_at_load(unit))
    set_temperature!(unit.input_interfaces[unit.m_heat_in],
                     unit.high_temperature,
                     unit.high_temperature)

    set_max_energy!(unit.input_interfaces[unit.m_heat_in], unit.capacity - unit.load)
    set_max_energy!(unit.output_interfaces[unit.m_heat_out], unit.load)
end

function temperature_at_load(unit::SeasonalThermalStorage)::Temperature
    if unit.use_adaptive_temperature
        partial_load = min(1.0, unit.load / (unit.capacity * unit.switch_point))
        return (unit.high_temperature - unit.low_temperature) * partial_load + unit.low_temperature
    else
        return unit.high_temperature
    end
end

function balance_on(interface::SystemInterface,
                    unit::SeasonalThermalStorage)::Vector{EnergyExchange}
    caller_is_input = unit.uac == interface.target.uac
    balance_written = interface.max_energy.max_energy[1] === nothing || interface.sum_abs_change > 0.0
    purpose_uac = unit.uac == interface.target.uac ? interface.target.uac : interface.source.uac

    return [EnEx(;
                 balance=interface.balance,
                 energy_potential=balance_written ? 0.0 : (caller_is_input ? -(unit.capacity - unit.load) : unit.load),
                 purpose_uac=purpose_uac,
                 temperature_min=interface.temperature_min,
                 temperature_max=interface.temperature_max,
                 pressure=nothing,
                 voltage=nothing)]
end

function process(unit::SeasonalThermalStorage, sim_params::Dict{String,Any})
    outface = unit.output_interfaces[unit.m_heat_out]
    exchanges = balance_on(outface, outface.target)
    energy_demanded = balance(exchanges) + energy_potential(exchanges)

    # shortcut if there is no energy demanded
    if energy_demanded >= -sim_params["epsilon"]
        set_max_energy!(unit.output_interfaces[unit.m_heat_out], 0.0)
        return
    end

    for exchange in exchanges
        demanded_on_interface = exchange.balance + exchange.energy_potential

        if demanded_on_interface >= -sim_params["epsilon"]
            continue
        end

        tank_temp = temperature_at_load(unit)
        if (exchange.temperature_min !== nothing &&
            exchange.temperature_min > tank_temp)
            # we can only supply energy at a temperature at or below the tank's current
            # output temperature
            continue
        end

        used_heat = min(abs(energy_demanded), abs(demanded_on_interface))

        if unit.load > used_heat
            unit.load -= used_heat
            add!(outface, used_heat, tank_temp)
            energy_demanded += used_heat
        else
            add!(outface, unit.load, tank_temp)
            energy_demanded += unit.load
            unit.load = 0.0
        end
    end
end

function load(unit::SeasonalThermalStorage, sim_params::Dict{String,Any})
    inface = unit.input_interfaces[unit.m_heat_in]
    exchanges = balance_on(inface, inface.source)
    energy_available = balance(exchanges) + energy_potential(exchanges)

    # shortcut if there is no energy to be used
    if energy_available <= sim_params["epsilon"]
        set_max_energy!(unit.input_interfaces[unit.m_heat_in], 0.0)
        unit.load_end_of_last_timestep = copy(unit.load)
        return
    end

    for exchange in exchanges
        exchange_energy_available = exchange.balance + exchange.energy_potential

        if exchange_energy_available < sim_params["epsilon"]
            continue
        end

        if (exchange.temperature_min !== nothing &&
            exchange.temperature_min > unit.high_temperature
            ||
            exchange.temperature_max !== nothing &&
            exchange.temperature_max < unit.high_temperature)
            # we can only take in energy if it's at a higher/equal temperature than the
            # storage's upper limit for temperatures
            continue
        end

        used_heat = min(energy_available, exchange_energy_available)
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

    unit.load_end_of_last_timestep = copy(unit.load)
end

function output_values(unit::SeasonalThermalStorage)::Vector{String}
    return [string(unit.m_heat_in) * " IN",
            string(unit.m_heat_out) * " OUT",
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
