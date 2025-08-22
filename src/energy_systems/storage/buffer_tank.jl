"""
Implementation of a buffer tank holding hot water for heating or DHW purposes.

This is a simplified model of a buffer tank with three different model types:
* ideally_stratified: An ideally stratified cylindrical tank providing thermal energy always 
                      at the tank's upper temperature. Losses, if activated, are only 
                      energy losses reducing the amount of available energy (no exergy losses). 
* ideally_mixed:      An ideally mixed cylindrical tank providing thermal energy at a temperature
                      between the tank's upper and lower temperature. Note that the upper 
                      temperature is only supplied at 100% load. Losses, if activated, are
                      considered as energy and exergy losses reducing the energy and the 
                      current supply temperature.
* balanced:           The balanced model of the cylindrical buffer tank is a mix of the 
                      ideally stratified and ideally mixed model. At a load higher than the 
                      user-defined switch point, the ideally stratified model is used. At a 
                      load less than the switch point, the model switches to the ideally mixed model, 
                      representing a more realistic temperature behaviour of the energy supply.
Note that all three models can only be loaded with energy at a temperature of at least the 
upper temperature of the tank.

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
    ground_temperature::Temperature

    # for model type "balanced"
    switch_point::Float64

    current_max_output_temperature::Float64
    initial_load::Float64
    load::Float64
    load_end_of_last_timestep::Float64
    losses::Float64
    process_done::Bool
    load_done::Bool

    function BufferTank(uac::String, config::Dict{String,Any}, sim_params::Dict{String,Any})
        medium = Symbol(default(config, "medium", "m_h_w_ht1"))
        register_media([medium])

        input_model_type = default(config, "model_type", "ideally_stratified")
        if input_model_type == "ideally_stratified"
            model_type = :ideally_stratified
        elseif input_model_type == "balanced"
            model_type = :balanced
        elseif input_model_type == "ideally_mixed"
            model_type = :ideally_mixed
        else
            @error "For the buffer tank $uac, the model_type could not be detected. Is has to be one of: " *
                   "'ideally_stratified', 'balanced', 'ideally_mixed'."
            throw(InputError)
        end

        consider_losses = default(config, "consider_losses", false)
        if consider_losses
            constant_ambient_temperature,
            ambient_temperature_profile = get_parameter_profile_from_config(config,
                                                                            sim_params,
                                                                            "ambient_temperature",
                                                                            "ambient_temperature_profile_file_path",
                                                                            "ambient_temperature_from_global_file",
                                                                            "constant_ambient_temperature",
                                                                            uac;
                                                                            required=true)

        else
            constant_ambient_temperature = nothing
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
                   consider_losses,                            # consider_losses [Bool]
                   default(config, "h_to_r", 2.0),             # ratio of height to radius of the cylinder [-]
                   0.0,                                        # surface_lid_bottom, surface of the lid and the bottom of the cylinder [m^2]
                   0.0,                                        # surface_barrel, surface of the barrel of the cylinder [m^2]
                   default(config, "thermal_transmission_lid", 1.2),     # [W/(m^2K)]
                   default(config, "thermal_transmission_barrel", 1.2),  # [W/(m^2K)]
                   default(config, "thermal_transmission_bottom", 1.2),  # [W/(m^2K)]
                   ambient_temperature_profile,                          # [°C]
                   constant_ambient_temperature,                         # ambient_temperature [°C]
                   default(config, "ground_temperature", 12),            # [°C]
                   default(config, "switch_point", 0.15),                # [%/100] load where to change from stratified to mixed mode, only for model_type :balanced
                   0.0,                                        # current_max_output_temperature at the beginning of the time step
                   default(config, "initial_load", 0.0),       # initial_load [%/100]
                   0.0,                                        # load, set to inital_load at the beginning [Wh]
                   0.0,                                        # load_end_of_last_timestep, stores the load of the previous time step without losses
                   0.0,                                        # losses in current time step [Wh]
                   false,                                      # process_done, bool indicating if the process step has already been performed in the current time step
                   false)                                      # load_done, bool indicating if the load step has already been performed in the current time step
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
        if unit.consider_losses
            unit.volume = unit.capacity /
                          (unit.rho_medium / 3.6 * unit.cp_medium * (unit.high_temperature - unit.low_temperature))  # [m^3]
        end
    else
        @error "For the buffer tank $(unit.uac), either a volume or a capacity has to be given, but both are given."
        throw(InputError)
    end

    # calculate maximum input and output energy
    if unit.max_load_rate === nothing
        unit.max_input_energy = Inf
    else
        unit.max_input_energy = unit.max_load_rate * unit.capacity * (sim_params["time_step_seconds"] / 60 / 60)  # [Wh per timestep]
    end
    if unit.max_unload_rate === nothing
        unit.max_output_energy = Inf
    else
        unit.max_output_energy = unit.max_unload_rate * unit.capacity * (sim_params["time_step_seconds"] / 60 / 60)  # [Wh per timestep]
    end

    # calculate surfaces of the buffer tank cylinder
    if unit.consider_losses
        radius = cbrt(unit.volume / (unit.h_to_r * pi))  # [m]
        height = radius * unit.h_to_r                    # [m]
        unit.surface_barrel = 2 * pi * radius * height   # [m^2]
        unit.surface_lid_bottom = radius^2 * pi          # [m^2]
    end

    # set initial state
    unit.load = unit.initial_load * unit.capacity
    unit.load_end_of_last_timestep = copy(unit.load)
end

function control(unit::BufferTank,
                 components::Grouping,
                 sim_params::Dict{String,Any})
    update(unit.controller)

    unit.current_max_output_temperature = temperature_at_load(unit)

    set_max_energy!(unit.input_interfaces[unit.medium], min(unit.capacity - unit.load, unit.max_input_energy),
                    unit.high_temperature, nothing)
    set_max_energy!(unit.output_interfaces[unit.medium], min(unit.load, unit.max_output_energy), nothing,
                    unit.current_max_output_temperature)

    if unit.ambient_temperature_profile !== nothing && unit.consider_losses
        unit.ambient_temperature = Profiles.value_at_time(unit.ambient_temperature_profile, sim_params)
    end
end

function temperature_at_load(unit::BufferTank)::Temperature
    if unit.model_type == :ideally_stratified
        # always high temperature is available
        return unit.high_temperature
    elseif unit.model_type == :balanced
        # When the storage is loaded above the switch_point, the high temperature is available. 
        # At loads below the switch_point, a linear course between high and low temperature can be supplied.
        partial_load = min(1.0, unit.load / (unit.capacity * unit.switch_point))
        return (unit.high_temperature - unit.low_temperature) * partial_load + unit.low_temperature
    elseif unit.model_type == :ideally_mixed
        # A linear course between high and low temperature is supplied, depending on the current load.
        return unit.low_temperature + (unit.load / unit.capacity) * (unit.high_temperature - unit.low_temperature)
    end
end

function calculate_losses!(unit::BufferTank, sim_params)
    if !unit.consider_losses
        unit.load_end_of_last_timestep = copy(unit.load)
        return
    end

    function calculate_energy_loss(unit, sim_params)
        if unit.load <= sim_params["epsilon"]
            # no gains or losses if storage is completely empty for ideally_stratified model
            return 0.0
        else
            barrel_surface = unit.surface_barrel * unit.load / unit.capacity
            temperature_difference = unit.high_temperature - unit.ambient_temperature

            return (unit.thermal_transmission_lid * unit.surface_lid_bottom * temperature_difference +
                    unit.thermal_transmission_barrel * barrel_surface * temperature_difference) *
                   sim_params["time_step_seconds"] / 60 / 60
        end
    end

    if unit.model_type == :ideally_stratified
        # losses only through layer with high temperature (energy losses)
        unit.losses = calculate_energy_loss(unit, sim_params)
    elseif unit.model_type == :balanced
        # losses above switch_point result in energy decrease, below switch_point in exergy losses
        if unit.load / unit.capacity > unit.switch_point
            unit.losses = calculate_energy_loss(unit, sim_params)
        else
            current_tank_temperature = unit.low_temperature +
                                       (unit.load / (unit.switch_point * unit.capacity)) *
                                       (unit.high_temperature - unit.low_temperature)
            temperature_difference_air = current_tank_temperature - unit.ambient_temperature
            barrel_surface = unit.surface_barrel * unit.switch_point

            unit.losses = (unit.thermal_transmission_lid * unit.surface_lid_bottom * temperature_difference_air +
                           unit.thermal_transmission_barrel * barrel_surface * temperature_difference_air) *
                          sim_params["time_step_seconds"] / 60 / 60
        end
    elseif unit.model_type == :ideally_mixed
        # losses are exergy losses through the whole storage
        current_tank_temperature = unit.low_temperature +
                                   (unit.load / unit.capacity) * (unit.high_temperature - unit.low_temperature)
        temperature_difference_air = current_tank_temperature - unit.ambient_temperature
        temperature_difference_ground = current_tank_temperature - unit.ground_temperature

        unit.losses = (unit.thermal_transmission_lid * unit.surface_lid_bottom * temperature_difference_air +
                       unit.thermal_transmission_bottom * unit.surface_lid_bottom * temperature_difference_ground +
                       unit.thermal_transmission_barrel * unit.surface_barrel * temperature_difference_air) *
                      sim_params["time_step_seconds"] / 60 / 60
    end

    # save load at the end of the current time step before applying the losses for the control modules
    unit.load_end_of_last_timestep = copy(unit.load)

    # update load of storage and limit losses to current load
    unit.losses = min(unit.losses, unit.load)
    unit.load = unit.load - unit.losses
end

function process(unit::BufferTank, sim_params::Dict{String,Any})
    outface = unit.output_interfaces[unit.medium]
    exchanges = balance_on(outface, outface.target)
    energy_demanded = balance(exchanges) + energy_potential(exchanges)

    # shortcut if there is no energy demanded
    if energy_demanded >= -sim_params["epsilon"]
        set_max_energy!(unit.output_interfaces[unit.medium], 0.0)
        handle_component_update!(unit, "process", sim_params)
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
            add!(outface, used_heat, nothing, tank_temp)
            energy_demanded += used_heat
        else
            add!(outface, unit.load, nothing, tank_temp)
            energy_demanded += unit.load
            unit.load = 0.0
        end
    end

    handle_component_update!(unit, "process", sim_params)
end

function handle_component_update!(unit::BufferTank, step::String, sim_params::Dict{String,Any})
    if step == "process"
        unit.process_done = true
    elseif step == "load"
        unit.load_done = true
    end
    if unit.process_done && unit.load_done
        # update component
        calculate_losses!(unit, sim_params)
        # reset 
        unit.process_done = false
        unit.load_done = false
    end
end

function load(unit::BufferTank, sim_params::Dict{String,Any})
    inface = unit.input_interfaces[unit.medium]
    exchanges = balance_on(inface, inface.source)
    energy_available = balance(exchanges) + energy_potential(exchanges)

    # shortcut if there is no energy to be used
    if energy_available <= sim_params["epsilon"]
        handle_component_update!(unit, "load", sim_params)
        set_max_energy!(unit.input_interfaces[unit.medium], 0.0)
        return
    end

    for exchange in exchanges
        exchange_energy_available = exchange.balance + exchange.energy_potential

        if exchange_energy_available < sim_params["epsilon"]
            continue
        end

        if exchange.temperature_max !== nothing &&
           exchange.temperature_max < unit.high_temperature
            # we can only take in energy if it's at a higher/equal temperature than the
            # tank's upper limit for temperatures
            continue
        end

        used_heat = min(energy_available, exchange_energy_available)
        diff = unit.capacity - unit.load

        if diff > used_heat
            unit.load += used_heat
            sub!(inface, used_heat, unit.high_temperature, nothing)
            energy_available -= used_heat
        else
            unit.load = unit.capacity
            sub!(inface, diff, unit.high_temperature, nothing)
            energy_available -= diff
        end
    end

    handle_component_update!(unit, "load", sim_params)
end

function output_values(unit::BufferTank)::Vector{String}
    return [string(unit.medium) * " IN",
            string(unit.medium) * " OUT",
            "Load",
            "Load%",
            "Capacity",
            "LossesGains",
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
    elseif key.value_key == "LossesGains"
        return -unit.losses
    elseif key.value_key == "CurrentMaxOutTemp"
        return unit.current_max_output_temperature
    end
    throw(KeyError(key.value_key))
end

export BufferTank
