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

    model_type::Symbol

    capacity::Floathing
    volume::Floathing
    rho_medium::Float64
    cp_medium::Float64
    high_temperature::Float64
    low_temperature::Float64
    max_load_rate::Floathing
    max_unload_rate::Floathing
    max_input_energy::Floathing
    max_output_energy::Floathing
    consider_losses::Bool

    # for losses
    h_to_r::Float64
    surface_lid_bottom::Float64
    surface_barrel::Float64
    thermal_transmission_lid::Float64
    thermal_transmission_barrel::Float64
    thermal_transmission_bottom::Float64
    ambient_temperature_profile::Union{Profile,Nothing}
    ambient_temperature::Temperature
    ground_temperauture::Temperature

    # for model type "balanced"
    switch_point::Float64

    current_max_output_temperature::Float64
    load::Float64
    losses::Float64

    function BufferTank(uac::String, config::Dict{String,Any}, sim_params::Dict{String,Any})
        medium = Symbol(default(config, "medium", "m_h_w_ht1"))
        register_media([medium])

        input_model_type = default(config, "model_type", "ideally_stratified")
        if input_model_type == "ideally_stratified"
            model_type = :stratified
        elseif input_model_type == "balanced"
            model_type = :balanced
        elseif input_model_type = "ideally_mixed"
            model_type = :ideally_mixed
        else
            @error "For the buffer tank $uac, the model_type could not be detected. Is has to be one of: " *
                   "'ideally_stratified', 'balanced', 'ideally_mixed'."
            throw(InputError)
        end

        constant_ambient_temperature = default(config, "constant_ambient_temperature", nothing)
        if constant_ambient_temperature !== nothing
            ambient_temperature_profile = get_ambient_temperature_profile_from_config(config, sim_params, uac)
        else
            ambient_temperature_profile = nothing
        end

        return new(uac, # uac
                   Controller(default(config, "control_parameters", nothing)),
                   sf_storage,                                 # sys_function
                   InterfaceMap(medium => nothing),            # input_interfaces
                   InterfaceMap(medium => nothing),            # output_interfaces
                   medium,                                     # medium in the buffer tank
                   model_type,                                 # can be one of :ideally_stratified, :balanced, :ideally_mixed
                   default(config, "capacity", nothing),       # capacity of the buffer tank [Wh] (either capacity or volume has to be given)
                   default(config, "volume", nothing),         # volume of the buffer tank [m^3] (either capacity or volume has to be given)
                   default(config, "rho_medium", 1000.0),      # density of the medium [kg/m^3]
                   default(config, "cp_medium", 4.18),         # specific thermal capacity of medium [kJ/kgK]
                   default(config, "high_temperature", 75.0),  # upper temperature of the buffer tank [°C]
                   default(config, "low_temperature", 20),     # lower temperature of the buffer tank [°C]
                   default(config, "max_load_rate", nothing),  # maximum load rate given in 1/h
                   default(config, "max_unload_rate", nothing),# maximum unload rate given in 1/h
                   nothing,                                    # maximum input energy per time step [Wh]
                   nothing,                                    # maximum output energy per time step [Wh]
                   default(config, "consider_losses", false),  # consider_losses [Bool]
                   default(config, "h_to_r", 2.0),             # ratio of height to radius of the cylinder [-]
                   0.0,               # surface_lid_bottom, surface of the lid and the bottom of the cylindder [m^2]
                   0.0,               # surface_barrel, surface of the barrel of the cylindder [m^2]
                   default(config, "thermal_transmission_lid", 0.2),     # [W/mK]
                   default(config, "thermal_transmission_barrel", 0.2),  # [W/mK]
                   default(config, "thermal_transmission_bottom", 0.2),  # [W/mK]
                   ambient_temperature_profile,                          # [°C]
                   constant_ambient_temperature,                         # ambient_temperature [°C]
                   default(config, "ground_temperauture", 12),           # [°C]
                   default(config, "switch_point", 0.15),                # [%/100] load where to change from stratified to mixed mode, only for model_type :balanced
                   0.0,                                        # current_max_output_temperature at the beginning of the time step
                   default(config, "initial_load", 0.0),       # current load, set to inital_load at the beginning [%]
                   0.0)                                        # losses in current time step [Wh]
    end
end

function initialise!(unit::BufferTank, sim_params::Dict{String,Any})
    set_storage_transfer!(unit.input_interfaces[unit.medium],
                          unload_storages(unit.controller, unit.medium))
    set_storage_transfer!(unit.output_interfaces[unit.medium],
                          load_storages(unit.controller, unit.medium))

    # calculate volume and capacity
    if unit.capacity === nothing && unit.volume === nothing
        @error "For the buffer tank $(unit.uac), either a volume or a capacity has to be given, but none of them is given."
        throw(InputError)
    elseif unit.capacity === nothing && unit.volume !== nothing
        unit.capacity = unit.volume * unit.rho_medium / 3.6 * unit.cp_medium *
                        (unit.high_temperature - unit.low_temperature)             # [Wh]
    elseif unit.capacity !== nothing && unit.volume === nothing
        unit.volume = unit.capacity /
                      (unit.rho_medium / 3.6 * unit.cp_medium * (unit.high_temperature - unit.low_temperature))  # [m^3]
    else
        @error "For the buffer tank $(unit.uac), either a volume or a capacity has to be given, but both are given."
        throw(InputError)
    end

    # calculate maximum input and output energy
    if unit.max_load_rate !== nothing
        unit.max_input_energy = unit.max_load_rate * unit.capacity / (sim_params["time_step_seconds"] / 60 / 60)  # [Wh per timestep]
    end
    if unit.max_unload_rate !== nothing
        unit.max_output_energy = unit.max_unload_rate * unit.capacity / (sim_params["time_step_seconds"] / 60 / 60)  # [Wh per timestep]
    end

    # calculate surfaces of the butter tank cylinder
    radius = cbrt(unit.volume / (unit.h_to_r * pi))  # [m]
    height = radius * unit.h_to_r                    # [m]
    unit.surface_barrel = 2 * pi * radius * height   # [m^2]
    unit.surface_lid_bottom = radius^2 * pi          # [m^2]
end

function control(unit::BufferTank,
                 components::Grouping,
                 sim_params::Dict{String,Any})
    update(unit.controller)

    unit.current_max_output_temperature = temperature_at_load(unit)

    set_temperature!(unit.output_interfaces[unit.medium],
                     nothing,
                     unit.current_max_output_temperature)
    set_temperature!(unit.input_interfaces[unit.medium],
                     unit.high_temperature,
                     unit.high_temperature)

    set_max_energy!(unit.input_interfaces[unit.medium], unit.capacity - unit.load)
    set_max_energy!(unit.output_interfaces[unit.medium], unit.load)

    if unit.constant_ambient_temperature !== nothing
        unit.ambient_temperature = unit.constant_ambient_temperature
    else
        unit.ambient_temperature = Profiles.value_at_time(unit.ambient_temperature_profile, sim_params)
    end
end

function temperature_at_load(unit::BufferTank)::Temperature
    if unit.use_adaptive_temperature
        partial_load = min(1.0, unit.load / (unit.capacity * unit.switch_point))
        return (unit.high_temperature - unit.low_temperature) * partial_load + unit.low_temperature
    else
        return unit.high_temperature
    end
end

function balance_on(interface::SystemInterface,
                    unit::BufferTank)::Vector{EnergyExchange}
    caller_is_input = unit.uac == interface.target.uac
    purpose_uac = unit.uac == interface.target.uac ? interface.target.uac : interface.source.uac

    return [EnEx(;
                 balance=interface.balance,
                 energy_potential=caller_is_input ? -(unit.capacity - unit.load) : unit.load,
                 purpose_uac=purpose_uac,
                 temperature_min=interface.temperature_min,
                 temperature_max=interface.temperature_max,
                 pressure=nothing,
                 voltage=nothing)]
end

function process(unit::BufferTank, sim_params::Dict{String,Any})
    outface = unit.output_interfaces[unit.medium]
    exchanges = balance_on(outface, outface.target)
    energy_demanded = balance(exchanges) + energy_potential(exchanges)

    # shortcut if there is no energy demanded
    if energy_demanded >= -sim_params["epsilon"]
        set_max_energy!(unit.output_interfaces[unit.medium], 0.0)
        return
    end

    for exchange in exchanges
        demanded_on_interface = exchange.balance + exchange.energy_potential

        if demanded_on_interface >= -sim_params["epsilon"]
            continue
        end

        tank_temp = temperature_at_load(unit)
        if (exchange.temperature_min !== nothing
            &&
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

function load(unit::BufferTank, sim_params::Dict{String,Any})
    inface = unit.input_interfaces[unit.medium]
    exchanges = balance_on(inface, inface.source)
    energy_available = balance(exchanges) + energy_potential(exchanges)

    # shortcut if there is no energy to be used
    if energy_available <= sim_params["epsilon"]
        set_max_energy!(unit.input_interfaces[unit.medium], 0.0)
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
            # tank's upper limit for temperatures
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
end

function output_values(unit::BufferTank)::Vector{String}
    return [string(unit.medium) * " IN",
            string(unit.medium) * " OUT",
            "Load",
            "Load%",
            "Capacity",
            "Losses",
            "CurrentMaxOutTemp"]
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
    elseif key.value_key == "Losses"
        return unit.losses
    elseif key.value_key == "CurrentMaxOutTemp"
        return unit.current_max_output_temperature
    end
    throw(KeyError(key.value_key))
end

export BufferTank
