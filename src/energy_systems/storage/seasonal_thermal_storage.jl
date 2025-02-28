"""
Implementation of a seasonal thermal storage component.

This is a simplified model, which mostly deals with amounts of energy and considers
temperatures only for the available temperature as the tank is depleted.
"""
mutable struct SeasonalThermalStorage <: Component
    # general
    uac::String
    controller::Controller
    sys_function::SystemFunction

    input_interfaces::InterfaceMap
    output_interfaces::InterfaceMap
    m_heat_in::Symbol
    m_heat_out::Symbol

    # geometry and physical properties
    capacity::Floathing
    volume::Floathing
    hr_ratio::Float64
    sidewall_angle::Float64
    rho_medium::Float64
    cp_medium::Float64
    diffussion_coefficient::Float64
    surface_area_lid::Float64
    surface_area_bottom::Float64
    surface_area_barrel_segments::Vector{Float64}
    volume_segments::Vector{Float64}
    height::Float64
    number_of_layer_total::Int64
    number_of_layer_above_ground::Int64

    # loading and unloading
    high_temperature::Float64
    low_temperature::Float64
    max_load_rate::Floathing
    max_unload_rate::Floathing
    max_input_energy::Floathing
    max_output_energy::Floathing

    # losses
    thermal_transmission_lid::Float64
    thermal_transmission_barrel::Float64
    thermal_transmission_bottom::Float64
    ambient_temperature_profile::Union{Profile,Nothing}
    ambient_temperature::Temperature
    ground_temperature::Temperature

    # other
    current_max_output_temperature::Float64
    initial_load::Float64
    load::Float64
    load_end_of_last_timestep::Float64
    losses::Float64
    temperature_segments::Vector{Float64}

    function SeasonalThermalStorage(uac::String, config::Dict{String,Any}, sim_params::Dict{String,Any})
        m_heat_in = Symbol(default(config, "m_heat_in", "m_h_w_ht1"))
        m_heat_out = Symbol(default(config, "m_heat_out", "m_h_w_lt1"))
        register_media([m_heat_in, m_heat_out])

        constant_ambient_temperature,
        ambient_temperature_profile = get_parameter_profile_from_config(config,
                                                                        sim_params,
                                                                        "ambient_temperature",
                                                                        "ambient_temperature_profile_file_path",
                                                                        "ambient_temperature_from_global_file",
                                                                        "constant_ambient_temperature",
                                                                        uac;
                                                                        required=true)

        return new(uac,                                            # uac
                   Controller(default(config, "control_parameters", nothing)),
                   sf_storage,                                     # sys_function
                   InterfaceMap(m_heat_in => nothing),             # input_interfaces
                   InterfaceMap(m_heat_out => nothing),            # output_interfaces
                   m_heat_in,                                      # medium in the STES   
                   m_heat_out,                                     # medium out of the STES
                   # geometry and physical properties
                   0.0,                                            # capacity of the STES [Wh]
                   default(config, "volume", nothing),             # volume of the STES [m^3]
                   default(config, "hr_ratio", 0.5),               # ratio of the height to the mean radius of the STES
                   default(config, "sidewall_angle", 0.0),         # angle of the sidewall of the STES with respect to the horizont [°]
                   default(config, "rho_medium", 1000.0),          # density of the medium [kg/m^3]
                   default(config, "cp_medium", 4.18),             # specific thermal capacity of medium [kJ/kgK]
                   default(config, "diffussion_coefficient", 0.0), # diffussion coefficient of the medium [m^2/s]
                   0.0,                                            # surface_area_lid, surface of the lid of the STES [m^2]
                   0.0,                                            # surface_area_bottom, surface of the bottom of the STES [m^2]
                   [],                                             # surface_area_barrel_segments, surface of the barrel segments of the STES [m^2]
                   [],                                             # volume_segments, volume of the segments of the STES [m^3]
                   0.0,                                            # height of the STES [m]
                   default(config, "number_of_layer_total", 25),   # number of layers in the STES
                   default(config, "number_of_layer_above_ground", 5), # number of layers above ground in the STES
                   # loading and unloading
                   default(config, "high_temperature", 75.0),      # upper temperature of the STES [°C]
                   default(config, "low_temperature", 20),         # lower temperature of the STES [°C]
                   default(config, "max_load_rate", nothing),      # maximum load rate given in 1/h
                   default(config, "max_unload_rate", nothing),    # maximum unload rate given in 1/h
                   nothing,                                        # max_input_energy, maximum input energy per time step [Wh]
                   nothing,                                        # max_output_energy, maximum output energy per time step [Wh]
                   # Losses
                   default(config, "thermal_transmission_lid", 1.2),     # [W/(m^2K)]
                   default(config, "thermal_transmission_barrel", 1.2),  # [W/(m^2K)]
                   default(config, "thermal_transmission_bottom", 1.2),  # [W/(m^2K)]
                   ambient_temperature_profile,                          # [°C]
                   constant_ambient_temperature,                         # ambient_temperature [°C]
                   default(config, "ground_temperature", 12),            # [°C]
                   # other
                   0.0,                                            # current_max_output_temperature at the beginning of the time step
                   default(config, "initial_load", 0.0),           # initial_load [%/100]
                   0.0,                                            # load, set to inital_load at the beginning [Wh]
                   0.0,                                            # load_end_of_last_timestep, stores the load of the previous time step without losses
                   0.0,                                            # losses in current time step [Wh]
                   [])                                             # temperature_segments: temperatures of the segments
    end
end

function initialise!(unit::SeasonalThermalStorage, sim_params::Dict{String,Any})
    set_storage_transfer!(unit.input_interfaces[unit.m_heat_in],
                          unload_storages(unit.controller, unit.m_heat_in))
    set_storage_transfer!(unit.output_interfaces[unit.m_heat_out],
                          load_storages(unit.controller, unit.m_heat_out))

    # set temperature vector: assuming a uniform temperature profile
    mean_temperature = unit.initial_load * (unit.high_temperature - unit.low_temperature) + unit.low_temperature
    unit.temperature_segments = [mean_temperature for _ in 1:(unit.number_of_layer_total)]

    # calculate geometry of the STES
    unit.surface_area_lid,
    unit.surface_area_barrel_segments,
    unit.surface_area_bottom,
    unit.volume_segments,
    unit.height = calc_cylinder_geometry(unit.uac,
                                         unit.volume,
                                         unit.sidewall_angle,
                                         unit.hr_ratio,
                                         unit.number_of_layer_total)

    unit.capacity = unit.volume * unit.rho_medium * unit.cp_medium * (unit.high_temperature - unit.low_temperature)  # [Wh]

    # calculate maximum input and output energy
    if unit.max_load_rate === nothing
        unit.max_input_energy = Inf
    else
        unit.max_input_energy = unit.max_load_rate * unit.capacity / (sim_params["time_step_seconds"] / 60 / 60)  # [Wh per timestep]
    end
    if unit.max_unload_rate === nothing
        unit.max_output_energy = Inf
    else
        unit.max_output_energy = unit.max_unload_rate * unit.capacity / (sim_params["time_step_seconds"] / 60 / 60)  # [Wh per timestep]
    end

    # set initial state
    unit.load = unit.initial_load * unit.capacity
    unit.load_end_of_last_timestep = copy(unit.load)
end

"""
    calc_cylinder_geometry(uac::String, volume::Float64, alpha::Float64, hr::Float64, segments::Int64)

Calculate the surface areas and volumes of sections of a cylindrical or truncated cone-shaped storage tank.

# Arguments
- `uac::String`: Unique identifier for the storage tank.
- `volume::Float64`: Total volume of the storage tank [m^3].
- `alpha::Float64`: Slope angle of the truncated cone in degrees. If `alpha` is 90, the tank is a cylinder.
- `hr::Float64`: Height-to-radius ratio of the tank.
- `segments::Int64`: Number of segments to divide the tank into for calculation.

# Returns
- `a_lid::Float64`: Surface area of the top lid of the tank [m^2].
- `a_barrel::Vector{Float64}`: Surface areas of the barrel sections of the tank [m^2].
- `a_bottom::Float64`: Surface area of the bottom of the tank [m^2].
- `v_section::Vector{Float64}`: Volumes of the sections of the tank [m^3].
- `height::Float64`: Calculated height of the tank [m].

"""
function calc_cylinder_geometry(uac::String, volume::Float64, alpha::Float64, hr::Float64, segments::Int64)
    a_barrel = zeros(segments)
    v_section = zeros(segments)

    if alpha == 90  # cylinder
        # calculate radius and height of cylinder
        radius = cbrt(volume / (pi * hr))
        height = hr * radius

        # calculate surfaces and volumes of the sections
        a_lid = pi * radius^2
        for n in 1:segments
            a_barrel[n] = 2 * pi * radius * height / segments
            v_section[n] = pi * radius^2 * height / segments
        end
        a_bottom = a_lid
    else  # truncated cone
        # helper functions
        deg2rad(deg) = deg * (π / 180)
        rad2deg(rad) = rad * (180 / π)

        alpha_rad = deg2rad(alpha)
        alpha_tan = tan(alpha_rad)

        radius_large = cbrt((3 * volume) / (pi * alpha_tan * (1 - ((2 * alpha_tan - hr) / (2 * alpha_tan + hr))^3)))
        radius_small = radius_large * (2 * alpha_tan - hr) / (2 * alpha_tan + hr)

        if radius_small <= 0
            alpha_min = rad2deg(atan(hr / 2))
            alpha_max = 180 - alpha_min
            hr_max = 2 * alpha_tan

            @error "For the STES $(uac), reduce the h/r ratio for the truncated cone or increase the slope angle.  " *
                   "For a given h/r ratio of $(hr), the slope angle must be greater than $(round(alpha_min, 2))° and " *
                   "less than $(round(alpha_max, 2))°. For a given slope angle of $(alpha)°, the h/r ratio must be " *
                   "less than $(round(hr_max, 2))°."
            throw(InputError)
        end

        height = hr * (radius_large + radius_small) / 2
        a_lid = pi * radius_large^2
        a_bottom = pi * radius_small^2

        for n in 1:segments
            radius_seg_large = radius_large - (segments - n) * (radius_large - radius_small) / segments
            radius_seg_small = radius_small + (n - 1) * (radius_large - radius_small) / segments
            a_barrel[n] = (radius_seg_large + radius_seg_small) * pi *
                          sqrt((radius_seg_large - radius_seg_small)^2 + (height / segments)^2)
            v_section[n] = (height / segments) * pi / 3 *
                           (radius_seg_large^2 + radius_seg_large * radius_seg_small + radius_seg_small^2)
        end
    end

    return a_lid, a_barrel, a_bottom, v_section, height
end

function control(unit::SeasonalThermalStorage,
                 components::Grouping,
                 sim_params::Dict{String,Any})
    update(unit.controller)

    unit.current_max_output_temperature = temperature_at_load(unit)

    set_temperature!(unit.output_interfaces[unit.m_heat_out],
                     nothing,
                     unit.current_max_output_temperature)
    set_temperature!(unit.input_interfaces[unit.m_heat_in],
                     unit.low_temperature,
                     unit.high_temperature)

    set_max_energy!(unit.input_interfaces[unit.m_heat_in], min(unit.capacity - unit.load, unit.max_input_energy))
    set_max_energy!(unit.output_interfaces[unit.m_heat_out], min(unit.load, unit.max_output_energy))

    if unit.ambient_temperature_profile !== nothing
        unit.ambient_temperature = Profiles.value_at_time(unit.ambient_temperature_profile, sim_params)
    end
end

function temperature_at_load(unit::SeasonalThermalStorage)::Temperature
    # the current maximal output temperature
    return unit.high_temperature   # TODO
end

function update_STES!(unit::SeasonalThermalStorage, sim_params)
    # calculate losses and update temperatures of the slices
    unit.losses = nothing  # TODO

    # save load at the end of the current time step before applying the losses for the control modules
    unit.load_end_of_last_timestep = copy(unit.load)

    # update load of storage and limit losses to current load
    unit.losses = min(unit.losses, unit.load)
    unit.load = unit.load - unit.losses
end

function balance_on(interface::SystemInterface,
                    unit::SeasonalThermalStorage)::Vector{EnergyExchange}
    caller_is_input = unit.uac == interface.target.uac
    balance_written = interface.max_energy.max_energy[1] === nothing || interface.sum_abs_change > 0.0
    purpose_uac = unit.uac == interface.target.uac ? interface.target.uac : interface.source.uac

    # The losses should be applied to the load at the beginning of the next timestep, but they
    # might already be included in the load here. To avoid other components calling the balance_on() 
    # in between, what then would include the current losses that are intended to be applied 
    # in the next timestep, the losses must be added here again.
    # In the first time step and ahead of their calculation, losses are zero.
    current_load = unit.load + unit.losses # TODO

    return [EnEx(;
                 balance=interface.balance,
                 energy_potential=balance_written ? 0.0 :
                                  (caller_is_input ? -(unit.capacity - current_load) : current_load),
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
        update_STES!(unit, sim_params)
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

    update_STES!(unit, sim_params)
end

function output_values(unit::SeasonalThermalStorage)::Vector{String}
    return [string(unit.m_heat_in) * " IN",
            string(unit.m_heat_out) * " OUT",
            "Load",
            "Load%",
            "Capacity",
            "Losses",
            "CurrentMaxOutTemp"]
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
    elseif key.value_key == "CurrentMaxOutTemp"
        return unit.current_max_output_temperature
    end
    throw(KeyError(key.value_key))
end

export SeasonalThermalStorage
