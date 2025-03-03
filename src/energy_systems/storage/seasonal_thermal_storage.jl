"""
Implementation of a seasonal thermal storage component.

This is a simplified model, which mostly deals with amounts of energy and considers
temperatures only for the available temperature as the tank is depleted.
"""

using Plots: plot, scatter, savefig
using Plots: Plots

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
    radius_small::Float64
    radius_large::Float64
    number_of_layer_total::Int64
    number_of_layer_above_ground::Int64
    output_layer_from_top::Int64
    dz::Vector{Float64}
    dz_normalized::Vector{Float64}
    sigma::Vector{Float64}
    lambda::Vector{Float64}
    phi::Vector{Float64}
    theta::Vector{Float64}

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
    effective_ambient_temperature::Temperature

    # state variables
    current_max_output_temperature::Float64
    initial_load::Float64
    load::Float64
    load_end_of_last_timestep::Float64
    losses::Float64
    temperature_segments::Vector{Float64}
    current_energy_input_output::Vector{Float64}
    current_temperature_input_output::Vector{Temperature}

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

        # Note: layer numbering starts at the bottom with index 1 and ends at the top with index number_of_layer_total
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
                   0.0,                                            # radius_small, small (lower) radius of the STES [m]
                   0.0,                                            # radius_large, large (upper) radius of the STES [m]
                   default(config, "number_of_layer_total", 25),   # number of layers in the STES
                   default(config, "number_of_layer_above_ground", 5), # number of layers above ground in the STES
                   default(config, "output_layer_from_top", 1),    # layer number of the output layer, counted from the top
                   [],                                             # dz, thickness of the layers of the STES [m]
                   [],                                             # dz_normalized, normalized dz with respect to to the volume of each section
                   [],                                             # [1/h]  factor for losses to ambiente: area_of_losses * U[kJ/m^2K] / (roh * cp * volume_segment)
                   [],                                             # [K/kJ] factor for input/output energy:  1 / (roh * cp * volume_segment ) 
                   [],                                             # [1/kg] factor for input/output mass flow:  1 / (roh  * volume_segment )
                   [],                                             # volume-ratios of sections: V_section[n-1] / (V_section[n] + V_section[n-1])

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
                   [],                                                   # effective_ambient_temperature corresponding to each layer [°C]          
                   # other
                   0.0,                                            # current_max_output_temperature at the beginning of the time step
                   default(config, "initial_load", 0.0),           # initial_load [%/100] assuming perfectly mixed storage at the begin
                   0.0,                                            # load, set to initial_load at the beginning [Wh]
                   0.0,                                            # load_end_of_last_timestep, stores the load of the previous time step without losses
                   0.0,                                            # losses in current time step [Wh]
                   [],                                             # temperature_segments: temperatures of the segments
                   [],                                             # current_energy_input_output, energy input (+) or output(-) in current time step [Wh]
                   [])                                             # current_temperature_input_output, temperature of the input in current time step [°C]
    end
end

function initialise!(unit::SeasonalThermalStorage, sim_params::Dict{String,Any})
    set_storage_transfer!(unit.input_interfaces[unit.m_heat_in],
                          unload_storages(unit.controller, unit.m_heat_in))
    set_storage_transfer!(unit.output_interfaces[unit.m_heat_out],
                          load_storages(unit.controller, unit.m_heat_out))

    # set temperature vector: assuming a uniform temperature profile (mixed storage)
    mean_temperature = unit.initial_load * (unit.high_temperature - unit.low_temperature) + unit.low_temperature
    unit.temperature_segments = [mean_temperature for _ in 1:(unit.number_of_layer_total)]

    # calculate thermal transmission coefficients # TODO other transmission below ground?
    unit.thermal_transmission_barrels = [unit.thermal_transmission_barrel for _ in 1:(unit.number_of_layer_total)]

    # calculate geometry of the STES
    unit.surface_area_lid,
    unit.surface_area_barrel_segments,
    unit.surface_area_bottom,
    unit.volume_segments,
    unit.height,
    unit.radius_small,
    unit.radius_large = calc_cylinder_geometry(unit.uac,
                                               unit.volume,
                                               unit.sidewall_angle,
                                               unit.hr_ratio,
                                               unit.number_of_layer_total)

    # calculate (equally spaced) thickness of the layers of the STES [m]
    unit.dz = fill(unit.height / unit.number_of_layer_total, unit.number_of_layer_total)
    # calculate the normalized dz with respect to to the volume of each section
    unit.dz_normalized = unit.dz .* sqrt.((unit.volume_segments .* unit.number_of_layer_total) ./ unit.volume)

    # calculate coefficient for losses to ambiente
    unit.sigma = zeros(unit.number_of_layer_total)
    unit.sigma[1] = unit.surface_area_bottom * unit.thermal_transmission_bottom * 3.6 /
                    (unit.rho_medium * unit.cp_medium * unit.volume_segments[1])          # [1/h] losses to ambiente through bottom
    unit.sigma[end] = unit.surface_area_lid * unit.thermal_transmission_lid * 3.6 /
                      (unit.rho_medium * unit.cp_medium * unit.volume_segments[end])      # [1/h] losses to ambiente through lid
    unit.sigma = unit.sigma .+
                 unit.surface_area_barrel_segments .* unit.thermal_transmission_barrel .* 3.6 ./
                 (unit.rho_medium * unit.cp_medium * unit.volume_segments)                # [1/h]  losses to ambiente though barrel

    # calculate coefficient for input/output energy --> needed? TODO
    unit.lambda = 1 / (unit.rho_medium * unit.cp_medium * unit.volume_segments)          # [K/kJ] 

    # coefficient for input/output mass flow, assuming water as fluid
    cp_water = 4.18                                                                      # [kJ/kgK]
    unit.phi = cp_water / (unit.cp_medium * unit.rho_medium * unit.volume_segments)      # [1/kg]

    # coefficient for buoyancy effects
    unit.theta = [unit.volume_segments[n - 1] / (unit.volume_segments[n] + unit.volume_segments[n - 1])
                  for n in 2:(unit.number_of_layer_total)]
    pushfirst!(unit.theta, 0.0)  # Set first element to 0

    unit.capacity = unit.volume * unit.rho_medium * unit.cp_medium * (unit.high_temperature - unit.low_temperature)  # [Wh]
    unit.load = unit.initial_load * unit.capacity

    # TODO move
    unit.load = (weighted_mean(unit.temperature_segments, unit.volume_segments) - unit.low_temperature) /
                (unit.high_temperature - unit.low_temperature)

    # calculate maximum input and output energy
    if unit.max_load_rate === nothing
        unit.max_input_energy = Inf
    else
        unit.max_input_energy = unit.max_load_rate * unit.capacity / (sim_params["time_step_seconds"] / 60 / 60)     # [Wh per timestep]
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
    weighted_mean(values::Vector{Any}, weights::Vector{Any}) -> Float64

Calculate the weighted mean of a set of values.

# Arguments
- `values::Vector{Any}`: A vector containing the values.
- `weights::Vector{Any}`: A vector containing the weights corresponding to each value.

# Returns
- `Float64`: The weighted mean of the values.
"""
function weighted_mean(values::Vector{Any}, weights::Vector{Any})
    return sum(values .* weights) / sum(weights)
end

"""
    calc_cylinder_geometry(uac::String, volume::Float64, alpha::Float64, hr::Float64, n_segments::Int64)

Calculate the geometry of a seasonal thermal energy storage (STES) system, either as a cylinder or a truncated cone.

# Arguments
- `uac::String`: Unique identifier for the STES
- `volume::Float64`: Total volume of the storage [m^3]
- `alpha::Float64`: Slope angle of the truncated cone in degrees with respect to the horizontal. 
                    If `alpha` is 90, the storage is a cylinder.
- `hr::Float64`: Height-to-radius ratio
- `n_segments::Int64`: Number of segments to divide the storage into for calculation

# Returns
- `a_lid::Float64`: Surface area of the top lid [m^2]
- `a_barrel::Vector{Float64}`: Surface areas of the barrel sections [m^2]
- `a_bottom::Float64`: Surface area of the bottom lid [m^2]
- `v_section::Vector{Float64}`: Volumes of the sections [m^3]
- `height::Float64`: Height of the storage [m]
- `radius_small::Float64`: Radius of the smaller base (for truncated cone) [m]
- `radius_large::Float64`: Radius of the larger base (for truncated cone) [m]

"""
function calc_cylinder_geometry(uac::String, volume::Float64, alpha::Float64, hr::Float64, n_segments::Int64)
    a_barrel = zeros(n_segments)
    v_section = zeros(n_segments)

    if alpha == 90  # cylinder
        # calculate radius and height of cylinder
        radius = cbrt(volume / (pi * hr))
        height = hr * radius

        # calculate surfaces and volumes of the sections
        a_lid = pi * radius^2
        for n in 1:n_segments
            a_barrel[n] = 2 * pi * radius * height / n_segments
            v_section[n] = pi * radius^2 * height / n_segments
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

        for n in 1:n_segments
            radius_seg_large = radius_large - (n_segments - n) * (radius_large - radius_small) / n_segments
            radius_seg_small = radius_small + (n - 1) * (radius_large - radius_small) / n_segments
            a_barrel[n] = (radius_seg_large + radius_seg_small) * pi *
                          sqrt((radius_seg_large - radius_seg_small)^2 + (height / n_segments)^2)
            v_section[n] = (height / n_segments) * pi / 3 *
                           (radius_seg_large^2 + radius_seg_large * radius_seg_small + radius_seg_small^2)
        end
    end

    return a_lid, a_barrel, a_bottom, v_section, height, radius_small, radius_large
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

    # TODO: check also if energy is more than the bottom-secton can hold!
    # Or change algorithm to consider more than one layer if energy exceeds the maximum?
    set_max_energy!(unit.input_interfaces[unit.m_heat_in], min(unit.capacity - unit.load, unit.max_input_energy))
    set_max_energy!(unit.output_interfaces[unit.m_heat_out], min(unit.load, unit.max_output_energy))

    if unit.ambient_temperature_profile !== nothing
        unit.ambient_temperature = Profiles.value_at_time(unit.ambient_temperature_profile, sim_params)
    end

    # effective ambient temperature for each layer of the STES
    unit.effective_ambient_temperature = vcat(fill(unit.ambient_temperature, unit.number_of_layer_above_ground),
                                              fill(unit.ground_temperature,
                                                   unit.number_of_layer_total - unit.number_of_layer_above_ground))
end

function temperature_at_load(unit::SeasonalThermalStorage)::Temperature
    # the current maximal output temperature
    return unit.temperature_segments[end]   # TODO
end

"""
    update_STES!(unit::SeasonalThermalStorage, sim_params)

This function updates the temperature segments of the STES unit by calculating the new temperatures 
for each layer based on thermal diffusion, losses, thermal input/output, and mass input/output
by solving a 1D partial differential equation using an explicit Euler method.
The function also accounts for mixing due to buoyancy effects if a temperature gradient is present.

The method is based on:
Lago, J. et al. (2019): A 1-dimensional continuous and smooth model for thermally stratified storage tanks including 
                        mixing and buoyancy, Applied Energy 248, S. 640-655: doi: 10.1016/j.apenergy.2019.04.139                   
Steinacker, H. (2022): Entwicklung eines dynamischen Simulationsmodells zur Optimierung von wärmegekoppelten 
                       Wasserstoffkonzepten für die klimaneutrale Quartiersversorgung, unpublished master thesis, 
                       University of Stuttgart.

# Arguments
- `unit::SeasonalThermalStorage`: The STES unit to be updated
- `sim_params`: simulation parameters

"""
function update_STES!(unit::SeasonalThermalStorage, sim_params)
    # set alias for better readability
    t_old = unit.temperature_segments
    dt = sim_params["time_step_seconds"] / 60 / 60  # [h]
    t_new = zeros(unit.number_of_layer_total)

    # set lower node always to 1 for charging and discharging to always use the whole storage
    lower_node = 1

    # find index of highest layer that is below the temperature of the charging input temperature
    # This represents a ideal charging system (lance e.g.)
    upper_node_charging = unit.number_of_layer_total
    for i in length(t_old):-1:1
        if t_old[i] < temp_charging_in
            upper_node_charging = i
            break
        end
    end

    # index of layer for discharging, currently only one is supported
    upper_node_discharging = unit.number_of_layer_total - unit.output_layer_from_top + 1

    # mass_in and mass_out are both positive here! TODO
    # temp_charging_in = 
    # mass_in = 
    # temp_discharging_out =
    # mass_out =

    # mass flow and temperatures for charging
    mass_in_temp = vcat(zeros(lower_node - 1),
                        t_old[(lower_node + 1):upper_node_charging],
                        [temp_charging_in],
                        zeros(unit.number_of_layer_total - upper_node_charging))
    mass_in_vec = [((lower_node <= n <= upper_node_charging) ? mass_in : 0.0) for n in 1:(unit.number_of_layer_total)]
    # mass flow and temperatures for discharging
    mass_out_temp = vcat(zeros(lower_node - 1),
                         [temp_discharging_out],
                         t_old[lower_node:(upper_node_discharging - 1)],
                         zeros(unit.number_of_layer_total - upper_node_discharging))
    mass_out_vec = [((lower_node <= n <= upper_node_discharging) ? mass_out : 0.0)
                    for n in 1:(unit.number_of_layer_total)]

    for n in unit.number_of_layer_total
        if n == 0                # bottom layer, single-side
            t_new[n] = t_old[n] +
                       (60 * 60 * unit.diffussion_coefficient * (t_old[n + 1] - t_old[n]) / unit.dz_normalized[n]^2 +    # thermal diffusion
                        unit.beta[n] * (unit.effective_ambient_temperature[n] - t_old[n]) +                              # losses through bottom and side walls
                        # unit.lamda[n] * (Q_in_out)[n] +                                                                # thermal input and output
                        unit.phi[n] * mass_in_vec[n] * (mass_in_temp[n] - t_old[n]) +                                    # mass input
                        unit.phi[n] * mass_out_vec[n] * (mass_out_temp[n] - t_old[n])) * dt                              # mass output
        elseif n == n_sec - 1    # top layer, single-side
            t_new[n] = t_old[n] +
                       (60 * 60 * unit.diffussion_coefficient * (t_old[n - 1] - t_old[n]) / unit.dz_normalized[n]^2 +    # thermal diffusion
                        unit.beta[n] * (unit.effective_ambient_temperature[n] - t_old[n]) +                              # losses through lid and side walls
                        # unit.lamda[n] * Q_in_out[n] +                                                                  # thermal input and output
                        unit.phi[n] * mass_in_vec[n] * (mass_in_temp[n] - t_old[n]) +                                    # mass input
                        unit.phi[n] * mass_out_vec[n] * (mass_out_temp[n] - t_old[n])) * dt                              # mass output
        else                     # mid layer
            t_new[n] = t_old[n] +
                       (60 * 60 * unit.diffussion_coefficient * (t_old[n + 1] + t_old[n - 1] - 2 * t_old[n]) /
                        unit.dz_normalized[n]^2 +                                                                        # thermal diffusion
                        unit.beta[n] * (unit.effective_ambient_temperature[n] - t_old[n]) +                              # losses through side walls
                        # unit.lamda[n] * Q_in_out[n] +                                                                  # thermal input and output
                        unit.phi[n] * mass_in_vec[n] * (mass_in_temp[n] - t_old[n]) +                                    # mass input
                        unit.phi[n] * mass_out_vec[n] * (mass_out_temp[n] - t_old[n])) * dt                              # mass output
        end

        if n > 0   # mixing due to buoancy effecs, if temperature gradient is present
            adjust_temp = max(0, t_new[n - 1] - t_new[n])
            t_new[n] = t_new[n] + theta[n] * adjust_temp
            t_new[n - 1] = t_new[n - 1] - (1 - theta[n]) * adjust_temp
        end
    end

    load_new = (weighted_mean(t_new, unit.volume_segments) - unit.low_temperature) /
               (unit.high_temperature - unit.low_temperature)
    unit.losses = unit.load + sum(Q_in_out) + sum(Mass_in_out .* unit.cp_medium .* (Mass_in_out_temp .- t_old)) / 3.6 -
                  load_new

    # save load at the end of the current time step before applying the losses for the control modules
    unit.load_end_of_last_timestep = copy(unit.load)
    unit.load = load_new
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
