import Interpolations as ip
using Roots
using Resie.SolarIrradiance

"""
Implementation of a solarthermal collector.
Works for flat plate collectors, vacuum tube collectors and PVT modules.
The calculation is based on EN ISO 9806:2017 for quasi-dynamic models.
"""

mutable struct SolarthermalCollector <: Component
    # general
    uac::String
    controller::Controller
    sys_function::SystemFunction

    input_interfaces::InterfaceMap
    output_interfaces::InterfaceMap
    m_heat_out::Symbol
    # collector installation
    collector_gross_area::Float64
    tilt_angle::Float64
    azimuth_angle::Float64
    ground_reflectance::Float64
    # collector parameters
    eta_0_b::Float64
    K_b_array::Array{Union{Missing,Float64},2}
    K_b_itp::Tuple{AbstractArray, AbstractArray}
    K_b::Float64
    K_d::Float64
    a_params::Array{Float64,1}
    sigma::Float64
    vol_heat_cap::Float64
    # weather profiles
    ambient_temperature_profile::Union{Profile,Nothing}
    beam_solar_hor_irradiance_profile::Union{Profile,Nothing}
    diffuse_solar_hor_irradiance_profile::Union{Profile,Nothing}
    long_wave_irradiance_profile::Union{Profile,Nothing}
    wind_speed_profile::Union{Profile,Nothing}
    # weather parameters
    wind_speed_reduction::Float64 
    ambient_temperature::Temperature
    beam_solar_irradiance_in_plane::Float64
    diffuse_solar_irradiance_in_plane::Float64
    reduced_wind_speed::Floathing
    beam_solar_hor_irradiance::Floathing
    diffuse_solar_hor_irradiance::Floathing
    long_wave_irradiance::Floathing
    direct_normal_irradiance::Float64
    # operation parameters
    delta_T::Union{Float64,Nothing}
    spec_flow_rate::Floathing
    delta_T_min::Union{Float64,Nothing}
    spec_flow_rate_min::Floathing
    # results
    used_energy::Float64
    output_temperature::Temperature
    average_temperature::Temperature
    last_average_temperature::Temperature
    spec_flow_rate_actual::Floathing
    # internal parameters for temperature layers
    available_energies::Array{Floathing}
    average_temperatures::Array{Floathing}
    output_temperatures::Array{Floathing}
    runtimes::Array{Floathing}
    temperature_energy_pairs::Dict{Temperature,Tuple}
    calc_uacs::Array{Stringing}
    mean_ambient_temperature::Float64

    function SolarthermalCollector(uac::String, config::Dict{String,Any}, sim_params::Dict{String,Any})
        const_ambient_temperature, 
        ambient_temperature_profile = get_parameter_profile_from_config(config, sim_params, "ambient_temperature",
                                                                        "ambient_temperature_profile_file_path",
                                                                        "ambient_temperature_from_global_file",
                                                                        "constant_ambient_temperature", uac; 
                                                                        required=true) # °C
        const_beam_solar_hor_irradiance, 
        beam_solar_hor_irradiance_profile = get_parameter_profile_from_config(config, sim_params,
                                                                              "beam_solar_radiation",
                                                                              "beam_solar_radiation_profile_file_path",
                                                                              "beam_solar_radiation_from_global_file",
                                                                              "constant_beam_solar_radiation", uac;
                                                                              required=true) # Wh/m²
        const_long_wave_irradiance, 
        long_wave_irradiance_profile = get_parameter_profile_from_config(config, sim_params,
                                                                         "infrared_sky_radiation",
                                                                         "infrared_sky_radiation_profile_file_path",
                                                                         "infrared_sky_radiation_from_global_file",
                                                                         "constant_infrared_sky_radiation", uac;
                                                                         required=true) # Wh/m²
        const_diffuse_hor_irradiance, 
        diffuse_solar_hor_irradiance_profile = get_parameter_profile_from_config(config, sim_params,
                                                                                 "diffuse_solar_radiation",
                                                                                 "diffuse_solar_radiation_profile_file_path",
                                                                                 "diffuse_solar_radiation_from_global_file",
                                                                                 "constant_diffuse_solar_radiation",
                                                                                 uac; required=true) # Wh/m²
        const_wind_speed, 
        wind_speed_profile = get_parameter_profile_from_config(config, sim_params, "wind_speed",
                                                               "wind_speed_profile_file_path", "wind_speed_from_global_file",
                                                               "constant_wind_speed", uac; required=true) # m/s
        if const_wind_speed === nothing
            const_wind_speed = 0
        else
            const_wind_speed = max(const_wind_speed * default(config, "wind_speed_reduction", 0.5) - 3, 0)
        end

        m_heat_out = Symbol(default(config, "m_heat_out", "m_h_w_ht1"))
        register_media([m_heat_out])

        return new(uac, # uac
                   Controller(default(config, "control_parameters", nothing)),
                   sf_bounded_source, # sys_function
                   InterfaceMap(), # input_interfaces
                   InterfaceMap(m_heat_out => nothing), # output_interfaces
                   m_heat_out, # medium name of output interface
                   # collector installation
                   config["collector_gross_area"], # gross area of collector
                   config["tilt_angle"], # tilt angle
                   config["azimuth_angle"], # azimuth angle or orientation angle
                   default(config, "ground_reflectance", 0.2), # reflectance of the ground around the collector
                   # collector parameters
                   config["eta_0_b"],
                   vcat([10, 20, 30, 40, 50, 60, 70, 80, 90]', config["K_b_t_array"]', config["K_b_l_array"]'),
                        # row1: angle 0° to 90°/vertical/west, 
                        # row2: transversal/tilt, 
                        # row3: longitudinal/azimuth/orientation
                   ([], []), # initialise Interpolation functions for K_b
                   1, # Calculate K_b from K_b_array in control()
                   config["K_d"], # collector parameter for diffuse irradiance 
                   config["a_params"], # collector sim_params a1 to a8 corresponding to EN ISO 9806:2017
                   5.670374419 * 10^-8, # Stefan Boltzmann Constant
                   default(config, "vol_heat_capacity", 4.2e6), # volumetric heat capacity of the fluid in the collector in [J/(m³*K)]
                   # weather profiles
                   ambient_temperature_profile,
                   beam_solar_hor_irradiance_profile,
                   diffuse_solar_hor_irradiance_profile,
                   long_wave_irradiance_profile,
                   wind_speed_profile, 
                   # weather parameters
                   default(config, "wind_speed_reduction", 0.5), # adjust the wind speed by this factor to account for different wind conditions compared to measured wind speed at 10m height
                   const_ambient_temperature, # ambient temperature [°C]
                   0.0, # beam_solar_irradiance in collector plane [W/m²]
                   0.0, # diffuse_solar_irradiance_in_plane [W/m²]
                   const_wind_speed, # wind_speed multiplied with wind speed reduction [m/s]
                   const_beam_solar_hor_irradiance, # beam horizontal irradiance [W/m²]
                   const_diffuse_hor_irradiance, # diffuse horizontal irradiance [W/m²]
                   const_long_wave_irradiance, # long wave irradiance [W/m²]
                   0.0, # direct_normal_irradiance calculated from sun position [W/m²]
                   # operation parameters
                   default(config, "delta_T", nothing), # delta_T between input and output temperature in [K]
                   default(config, "spec_flow_rate", nothing), # nominal specific volume flow of the thermal collector in [m³/(m²*s)]
                   default(config, "delta_T_min", 2), # minimal delta_T between input and output temperature for the collector to start producing energy in K; used together with spec_flow_rate
                   default(config, "spec_flow_rate_min", 0.000002), # minimal specific volume flow of the thermal collector to start producing energy in m³/(m²*s); used together with delta_T
                   # results
                   0.0, # energy that was removed from the system in current time step
                   0.0, # output temperature in current time step, calculated in control()
                   0.0, # average temperature in current time step, calculated in control()
                   0.0, # average temperature from last time step, calculated in control()
                   0.0, # actual flow rate resulting from simulation and used energy
                   # internal parameters for temperature layers
                   [], # list of available_energies in for each component connected to the collector
                   [], # list of average_temperatures in for each component connected to the collector
                   [], # list of output_temperatures in for each component connected to the collector
                   [], # list of runtimes at different temperature levels for each component connected to the collector
                   Dict{Temperature,Tuple}(), # list of all calculated output temperature and output energy pairs
                   [], # uacs for which energy was calculated
                   10.0)
    end
end

function initialise!(unit::SolarthermalCollector, sim_params::Dict{String,Any})
    set_storage_transfer!(unit.output_interfaces[unit.m_heat_out],
                          load_storages(unit.controller, unit.m_heat_out))

    if sim_params["longitude"] === nothing || sim_params["latitude"] === nothing
        @error "Longitude and latitude must be provided through a weather file or in the 
        simulation parameters to calculate the sun position."
        throw(InputError)
    end

    unit.output_temperature = Profiles.value_at_time(unit.ambient_temperature_profile, sim_params)
    unit.average_temperature = unit.output_temperature
    unit.last_average_temperature = unit.output_temperature
    if unit.ambient_temperature_profile !== nothing
        unit.mean_ambient_temperature = _mean(collect(values(unit.ambient_temperature_profile.data)))
    else
        unit.mean_ambient_temperature = unit.ambient_temperature
    end

    unit.K_b_itp = init_K_b(unit.K_b_array)
end

function control(unit::SolarthermalCollector,
                 components::Grouping,
                 sim_params::Dict{String,Any})
    update(unit.controller)

    # get values from profiles if available and convert irradiances from energy to power
    if unit.beam_solar_hor_irradiance_profile !== nothing
        unit.beam_solar_hor_irradiance = Profiles.power_at_time(unit.beam_solar_hor_irradiance_profile,
                                                                sim_params)
    end

    if unit.diffuse_solar_hor_irradiance_profile !== nothing
        unit.diffuse_solar_hor_irradiance = Profiles.power_at_time(unit.diffuse_solar_hor_irradiance_profile,
                                                                   sim_params)
    end

    if unit.ambient_temperature_profile !== nothing
        unit.ambient_temperature = Profiles.value_at_time(unit.ambient_temperature_profile, sim_params)
    end

    if unit.wind_speed_profile !== nothing
        unit.reduced_wind_speed = max(Profiles.value_at_time(unit.wind_speed_profile, sim_params) * unit.wind_speed_reduction - 3, 0)
    end
    
    if unit.long_wave_irradiance_profile !== nothing
        unit.long_wave_irradiance = Profiles.power_at_time(unit.long_wave_irradiance_profile,
                                                           sim_params)                                                                      
    end

    unit.beam_solar_irradiance_in_plane, unit.diffuse_solar_irradiance_in_plane,
    unit.direct_normal_irradiance, angle_of_incidence, longitudinal_angle, transversal_angle,
    solar_zenith, solar_azimuth = irr_in_plane(sim_params, unit.tilt_angle, unit.azimuth_angle,
                                               unit.beam_solar_hor_irradiance, unit.diffuse_solar_hor_irradiance, 
                                               nothing, 1.0, unit.mean_ambient_temperature, unit.ground_reflectance)

    # Calculate K_b based on a table with provided values for longitudinal and transversal angles.
    # Table with provided values is mirrored for negative angles.
    if abs(transversal_angle) < 90 && abs(longitudinal_angle) < 90
        unit.K_b = unit.K_b_itp[1](abs(transversal_angle)) * unit.K_b_itp[2](abs(longitudinal_angle))
    else
        unit.K_b = 0
    end

    # reset list of all calculated output temperature and output energy pairs; 
    # used as lookup table to speed up calculation
    unit.temperature_energy_pairs = Dict{Temperature,Tuple}()

    unit.available_energies, has_calculated_all_maxima, unit.average_temperatures,
    unit.output_temperatures, unit.runtimes = calculate_energies(unit, components, sim_params)

    unit.output_temperature = highest(unit.output_temperatures)

    set_max_energy!(unit.output_interfaces[unit.m_heat_out],
                    unit.available_energies[1:end-1],
                    [nothing for _ in 1:length(unit.available_energies)-1],
                    unit.output_temperatures[1:end-1],
                    unit.calc_uacs[1:end-1],
                    has_calculated_all_maxima,
                    false)
end

function process(unit::SolarthermalCollector, sim_params::Dict{String,Any})
    exchanges = balance_on(unit.output_interfaces[unit.m_heat_out], unit.output_interfaces[unit.m_heat_out].target)
    energy_demands = [abs(e.balance + e.energy_potential) for e in exchanges]

    output_temperature_interface = lowest(unit.output_temperature, temp_min_highest(exchanges))
    unit.used_energy = 0.0

    if sum(energy_demands) > 0 && _sum(unit.available_energies) > 0
        uac_idx = 0
        for (uac_idx, e) in enumerate(exchanges)
            comp_energy_demand = abs(e.balance + e.energy_potential)
            used_energy_comp = min(comp_energy_demand, unit.available_energies[uac_idx])

            # recalculate the average_temperature if not all energy is used to keep a constant 
            # flow rate
            if unit.delta_T === nothing && unit.spec_flow_rate !== nothing &&
               used_energy_comp < unit.available_energies[uac_idx]
               # end of condition
               unit.average_temperatures[uac_idx] = e.temperature_min -
                                                    sim_params["wh_to_watts"](used_energy_comp) /
                                                    unit.collector_gross_area /
                                                    (unit.spec_flow_rate * 2 * unit.vol_heat_cap)
            end

            unit.used_energy += used_energy_comp
            unit.runtimes[uac_idx] = comp_energy_demand / unit.available_energies[uac_idx] * unit.runtimes[uac_idx]
        end

        unit.runtimes[end] = 0
        add!(unit.output_interfaces[unit.m_heat_out],
             unit.used_energy,
             output_temperature_interface)

    else
        unit.used_energy = 0.0
        set_max_energy!(unit.output_interfaces[unit.m_heat_out], 0.0)
    end

    for (uac_idx, uac) in enumerate(unit.calc_uacs)
        if !(uac in [e.purpose_uac for e in exchanges] || uac == unit.uac) && !isnothing(unit.average_temperatures[uac_idx])
            unit.runtimes[uac_idx] = 0
        elseif isnothing(unit.average_temperatures[uac_idx])
            unit.runtimes[uac_idx] = nothing
        end
    end

    if _sum(unit.runtimes) == 0
        stagnation_temp = get_stagnation_temperature(unit, sim_params)
        unit.output_temperature = stagnation_temp
        unit.average_temperature = stagnation_temp
    else
        unit.output_temperature = _weighted_mean(unit.output_temperatures, unit.runtimes)
        unit.average_temperature = _weighted_mean(unit.average_temperatures, unit.runtimes)
    end

    unit.last_average_temperature = unit.average_temperature

    if unit.output_temperature == unit.average_temperature
        unit.spec_flow_rate_actual = 0
    else
        unit.spec_flow_rate_actual = (unit.used_energy * 3600.0 / sim_params["time_step_seconds"]) /
                                     unit.collector_gross_area /
                                     ((unit.output_temperature - unit.average_temperature) * 2 * unit.vol_heat_cap)
    end
end

"""
    calculate_energies(unit::SolarthermalCollector, components, sim_params::Dict{String,Any})

Function used in the control() function of the solarthermal collector to calculate the 
available energies, average temperatures, output temperatures and partial runtime for each 
connected component. Results are written as dictionaries with the uac of the component as key.

# Arguments
- `unit::SolarthermalCollector`: The solarthermal collector for which the calculation is 
                                 performed.
- `components::Grouping`: All components of the energy system.
- `sim_params::Dict{String,Any}`: Simulation parameters.

# Returns
- `available_energies::Dict{String,Floathing}`: The available energy to each connected 
                                                component with uac as keys.
- `calc_max::Bool`: A Bool indicating if max energy was calculated or not.
- `average_temperatures::Dict{String,Temperature}`: The average collector temperatures for 
                                                     each connected component with uac as keys.
- `output_temperatures::Dict{String,Temperature}`: The output collector temperatures for each  
                                                   connected component with uac as keys.
- `runtimes::Dict{String,Floathing}`: The partial time as part of 1 each connected component  
                                      gets delivered energy with uac as keys.
"""
function calculate_energies(unit::SolarthermalCollector, components, sim_params::Dict{String,Any})
    left_energy_factor = 1.0
    component_idx = 1

    exchanges = balance_on(unit.output_interfaces[unit.m_heat_out], unit.output_interfaces[unit.m_heat_out].target)

    available_energies = Array{Floathing}(nothing, length(exchanges) + 1)
    average_temperatures = Array{Floathing}(nothing, length(exchanges) + 1)
    output_temperatures = Array{Floathing}(nothing, length(exchanges) + 1)
    runtimes = Array{Floathing}(nothing, length(exchanges) + 1)

    potential_energies = Float64[]
    unit.calc_uacs = Stringing[]
    target_temperatures = Temperature[]
    for exchange in exchanges
        # no check for possibly written max energy --> is done by set_max_energy!() for 
        # 1-to-1 connections and for a connection to a bus, the interface is only used by 
        # the current component
        success, target_temperature, max_energy = determine_temperature_and_energy(unit.controller,
                                                                                   components,
                                                                                   unit.uac,
                                                                                   exchange.purpose_uac,
                                                                                   sim_params)
        if !success
            # no control module is provided between source and target
            # set temperature to exchange.temperature_min --> can also be nothing, but in 
            # case it is given we will use it if no control module is used
            target_temperature, max_energy = check_temperature_and_get_max_energy(unit,
                                                                                  sim_params,
                                                                                  exchange.temperature_min,
                                                                                  false)
        end

        push!(unit.calc_uacs, exchange.purpose_uac)
        push!(target_temperatures, target_temperature)
        push!(potential_energies, max_energy)
    end

    while component_idx <= length(exchanges) && left_energy_factor > 0
        if potential_energies[component_idx] == 0
            t_avg = nothing
            t_target = nothing
            produced_energy = 0.0
        else
            produced_energy = calculate_output_energy_from_output_temperature(unit, target_temperatures[component_idx],
                                                                              sim_params) * left_energy_factor
            _, t_avg, t_target = unit.temperature_energy_pairs[target_temperatures[component_idx]]
        end

        # check if more energy can be produced for a certain component than is demanded
        if any(isinf, potential_energies)
            available_energies[component_idx] = produced_energy
            runtimes[component_idx] = 1

        elseif produced_energy == 0
            available_energies[component_idx] = produced_energy
            runtimes[component_idx] = nothing

        elseif produced_energy > potential_energies[component_idx]
            available_energies[component_idx] = potential_energies[component_idx]
            left_energy_factor *= (produced_energy - potential_energies[component_idx]) / produced_energy
            runtimes[component_idx] = (produced_energy - potential_energies[component_idx]) / produced_energy

        else
            available_energies[component_idx] = produced_energy
            runtimes[component_idx] = left_energy_factor
            left_energy_factor = 0
        end

        average_temperatures[component_idx] = t_avg
        output_temperatures[component_idx] = t_target
        component_idx += 1
    end

    push!(unit.calc_uacs, unit.uac)

    # if no energy at expected temperature levels is available then use max_output_temperature 
    # which equals average temperature with a used energy of 0
    if _sum(available_energies) == 0
        stagnation_temp = get_stagnation_temperature(unit, sim_params)
        calc_max = false
        average_temperatures[end] = stagnation_temp
        output_temperatures[end] = stagnation_temp
        runtimes[end] = 1

    elseif any(isinf, potential_energies)
        stagnation_temp = get_stagnation_temperature(unit, sim_params)
        calc_max = true
        average_temperatures[end] = stagnation_temp
        output_temperatures[end] = stagnation_temp
        runtimes[end] = 1
    else
        calc_max = false
        average_temperatures[end] = nothing
        output_temperatures[end] = nothing
        runtimes[end] = nothing
    end

    return (available_energies,
            calc_max,
            average_temperatures,
            output_temperatures,
            runtimes)
end

"""
    get_output_temperature_bounds(unit::SolarthermalCollector, 
                                 sim_params::Dict{String,Any})

Function used by the control module negotiate_temperature to get the current minimum 
and maximum temperatures that can be supplied.

# Arguments
- `unit::SolarthermalCollector`: The solarthermal collector for which the calculation is 
                                 performed.
- `sim_params::Dict{String,Any}`: simulation parameters

# Returns
- `output_min_temperature::Temperature`: The current minimal temperature that can be given 
                                         in the output
- `output_max_temperature::Temperature`: The current maximum temperature that can be given 
                                         in the output

"""
function get_output_temperature_bounds(unit::SolarthermalCollector, 
                                       sim_params::Dict{String,Any})::Tuple{Temperature,Temperature}
    if unit.delta_T_min == 0.0 || unit.spec_flow_rate_min == 0
        current_max_temperature = get_stagnation_temperature(unit, sim_params)
        energy_at_max_temp = 0

    elseif unit.delta_T !== nothing && unit.spec_flow_rate === nothing
        avg_temp = find_zero((t_avg -> spec_thermal_power_func(t_avg, t_avg + unit.delta_T / 2, 
                                                               unit.last_average_temperature,
                                                               unit.spec_flow_rate_min, unit, 
                                                               sim_params["time_step_seconds"]),
                              t_avg -> derivate_spec_thermal_power_func(t_avg, unit.spec_flow_rate_min, unit,
                                                                        sim_params["time_step_seconds"])),
                             unit.last_average_temperature,
                             Roots.Newton())
        current_max_temperature = avg_temp + unit.delta_T / 2
        energy_at_max_temp = unit.spec_flow_rate_min * unit.delta_T * unit.vol_heat_cap

    elseif unit.delta_T === nothing && unit.spec_flow_rate !== nothing
        avg_temp = find_zero((t_avg -> spec_thermal_power_func(t_avg, t_avg + unit.delta_T_min / 2, 
                                                               unit.last_average_temperature,
                                                               unit.spec_flow_rate, unit, 
                                                               sim_params["time_step_seconds"]),
                              t_avg -> derivate_spec_thermal_power_func(t_avg, unit.spec_flow_rate, unit,
                                                                        sim_params["time_step_seconds"])),
                             unit.last_average_temperature,
                             Roots.Newton())
        current_max_temperature = avg_temp + unit.delta_T_min / 2
        energy_at_max_temp = unit.spec_flow_rate * unit.delta_T_min * unit.vol_heat_cap

    else
        @error "Error in config file: Exclusively delta_T OR spec_flow_rate must have a value"
        throw(InputError)
    end
    
    unit.temperature_energy_pairs[current_max_temperature] = (energy_at_max_temp, 
                                                              current_max_temperature, 
                                                              current_max_temperature)

    return nothing, current_max_temperature
end

"""
    function get_stagnation_temperature(unit::SolarthermalCollector, sim_params::Dict{String,Any})

Convenience function to calculation stagnation temperature of the solarthermal collector. 
This is the temperature when no energy is taken from the collector by setting the flow rate 
to 0.

# Arguments
- `unit::SolarthermalCollector`: The solarthermal collector for which the calculation is 
                                 performed.
- `sim_params::Dict{String,Any}`: simulation parameters

# Returns
- `stagnation_temp::Temperature`: Stagnation temperature with no energy taken from the collector
"""
function get_stagnation_temperature(unit::SolarthermalCollector, sim_params::Dict{String,Any})::Temperature
    stagnation_temp = find_zero((t_avg -> spec_thermal_power_func(t_avg, t_avg, unit.last_average_temperature,
                                                                        0, unit, sim_params["time_step_seconds"]),
                                        t_avg -> derivate_spec_thermal_power_func(t_avg, 0, unit,
                                                                                sim_params["time_step_seconds"])),
                                    unit.last_average_temperature,
                                    Roots.Newton())
    return stagnation_temp
end

"""
    calculate_output_energy_from_output_temperature(unit::SolarthermalCollector,
                                                    minimum_output_temperature::Temperature, 
                                                    sim_params::Dict{String,Any})

Calculates the maximum energy that can be taken out of the collector field while not  
undercutting the given minimum_output_temperature.
If no converging result can be found, the maximum energy is set to the maximum energy of the 
collector field. The maximum energy is limited to the user-input of the max_energy.

# Arguments
- `unit::SolarthermalCollector`: the unit struct of the solarthermal collector with all its 
                                 parameters
- `minimum_output_temperature::Temperature`: The desired minimum output temperature of the 
                                             fluid temperature
- `sim_params::Dict{String,Any}`: simulation parameters
# Returns
- `current_max_output_energy::Float64`: the maximum energy that can be taken out of the  
                                        solarthermal collector while not undercutting the 
                                        given minimum_output_temperature.     
"""
function calculate_output_energy_from_output_temperature(unit::SolarthermalCollector,
                                                         minimum_output_temperature::Temperature,
                                                         sim_params::Dict{String,Any})::Float64

    # calculate specific thermal power depending on control mode; reuse already calculated values if possible
    if minimum_output_temperature in keys(unit.temperature_energy_pairs)
        p_spec_th, t_avg, t_target = unit.temperature_energy_pairs[minimum_output_temperature]

    elseif unit.delta_T !== nothing && unit.spec_flow_rate === nothing
        p_spec_th, t_avg, t_target = calc_thermal_power_fixed_delta_T(unit, sim_params, minimum_output_temperature)

    elseif unit.delta_T === nothing && unit.spec_flow_rate !== nothing
        p_spec_th, t_avg, t_target = calc_thermal_power_fixed_flow_rate(unit, sim_params, minimum_output_temperature)

    else
        @error "Error in config file: Exclusively delta_T OR spec_flow_rate must have a value"
        throw(InputError)
    end

    unit.temperature_energy_pairs[minimum_output_temperature] = (p_spec_th, t_avg, t_target)
    max_energy = isnothing(p_spec_th) ? 0.0 : sim_params["watt_to_wh"](p_spec_th * unit.collector_gross_area)

    return max_energy
end

"""
    check_temperature_and_get_max_energy(unit::SolarthermalCollector,
                                         sim_params::Dict{String,Any},
                                         temperature_output::Temperature)

Checks if a requested output temperature is within the allowed operational range for the 
solarthermal collector. If no temperature is given, the function sets the temperature to the 
mean temperature that will be reached with maximum energy draw during the whole time step.
The function also determines the maximum extractable energy for the given temperature.

# Arguments
- `unit::SolarthermalCollector`: The solarthermal collector unit for which the calculation 
                                 is performed.
- `sim_params::Dict{String,Any}`: Simulation parameters.
- `temperature_output::Temperature`: The requested output temperature of the solarthermal 
                                     collector.
- `limit_max_output_energy_to_avoid_pulsing::Bool`: Not used in solarthermal collector

# Returns
Returns a tuple containing:
- The temperature that can be used for energy extraction.
- The maximum energy that can be extracted under the given temperature.
"""
function check_temperature_and_get_max_energy(unit::SolarthermalCollector,
                                              sim_params::Dict{String,Any},
                                              temperature_output::Temperature,
                                              limit_max_output_energy_to_avoid_pulsing::Bool)::Tuple{Temperature,
                                                                                                     Float64}
    # get max output temperature of solarthermal collector
    _, source_max_out_temperature = get_output_temperature_bounds(unit, sim_params)

    if temperature_output > source_max_out_temperature || temperature_output === nothing
        # the requested temperature is higher than the current output temperature or
        # no temperature information is given, the collector doesn't run and produces no energy
        temperature_output = nothing
        max_energy = 0.0
    else
        # calculate the maximum energy that can be delivered while not changing the output temperature for the next time step
        max_energy = calculate_output_energy_from_output_temperature(unit, temperature_output, sim_params)
    end

    return temperature_output, max_energy
end

"""
    init_K_b(K_b_array::Array{Union{Missing,Float64},2})

Calculate K_b based on a table with provided values for longitudinal and transversal angles 
and the angle of irradiance on the plane of the solarthermal collector.
Table with provided values is mirrored for negative angles.

# Arguments
- `K_b_array::Array{Union{Missing,Float64},2}`: Simulation parameters.

# Returns
Returns a tuple containing two interpolation objects:
- interpolation object for transversal K_b value.
- interpolation object for longitudinal K_b value.
"""
function init_K_b(K_b_array::Array{Union{Missing,Float64},2})
    K_b_array = hcat([0, 1, 1], K_b_array)
    K_b_array[:, 10] = [90, 0, 0]
    angle_range = K_b_array[1, 1]:10:K_b_array[1, end]

    if any(ismissing, K_b_array)
        K_b_filtered = []
        missing_idx = []

        for col_idx in 1:length(K_b_array[1, :])
            if any(ismissing, K_b_array[:, col_idx])
                append!(missing_idx, col_idx)
            else
                append!(K_b_filtered, K_b_array[:, col_idx])
            end
        end
        K_b_filtered = reshape(K_b_filtered, (3, :))

        itp_t = ip.interpolate((K_b_filtered[1, :],), K_b_filtered[2, :], ip.Gridded(ip.Linear()))
        itp_l = ip.interpolate((K_b_filtered[1, :],), K_b_filtered[3, :], ip.Gridded(ip.Linear()))

        for idx in missing_idx
            K_b_array[2, idx] = itp_t(angle_range[idx])
            K_b_array[3, idx] = itp_l(angle_range[idx])
        end
    end
    interp_transversal = ip.scale(ip.interpolate(K_b_array[2, :], ip.BSpline(ip.Cubic(ip.Flat(ip.OnGrid())))),
                                  angle_range)
    interp_longitudinal = ip.scale(ip.interpolate(K_b_array[3, :], ip.BSpline(ip.Cubic(ip.Flat(ip.OnGrid())))),
                                   angle_range)

    return interp_transversal, interp_longitudinal
end

"""
    calc_thermal_power_fixed_delta_T(unit::SolarthermalCollector, 
                                     sim_params::Dict{String,Any}, 
                                     target_temperature::Temperature)

Calculate the specific thermal power of the solarthermal collector under the assumption that 
the collector has a fixed delta_T and the flow rate is variable 

# Arguments
- `unit::SolarthermalCollector`: The solarthermal collector unit for which the calculation 
                                 is performed.
- `sim_params::Dict{String,Any}`: Simulation parameters.
- `target_temperature::Temperature`: The targeted output temperature of the solarthermal 
                                     collector.

# Returns
- `p_spec_th::Floathing`: The specific thermal power of the collector if target temperature 
                          was reached otherwise nothing.
- `average_temperature::Temperature`: The average temperature of the collector if target   
                                      temperature was reached otherwise nothing.
- `target_temperature::Temperature`: The targeted output temperature of the collector if it 
                                     was reached otherwise nothing.
"""
function calc_thermal_power_fixed_delta_T(unit::SolarthermalCollector, 
                                          sim_params::Dict{String,Any}, 
                                          target_temperature::Temperature)
    average_temperature = target_temperature - unit.delta_T / 2
    p_spec_th = spec_thermal_power_func(average_temperature, target_temperature, unit.last_average_temperature, 0, unit,
                                        sim_params["time_step_seconds"])
    spec_flow_rate = p_spec_th / (unit.delta_T * unit.vol_heat_cap)

    if spec_flow_rate >= unit.spec_flow_rate_min
        return (p_spec_th,
                average_temperature,
                target_temperature)
    else
        return (nothing,
                nothing,
                nothing)
    end
end

"""
    calc_thermal_power_fixed_delta_T(unit::SolarthermalCollector, 
                                     sim_params::Dict{String,Any}, 
                                     target_temperature::Temperature)

Calculate the specific thermal power of the solarthermal collector under the assumption that 
the collector has a fixed flow rate and the delta_T is variable 

# Arguments
- `unit::SolarthermalCollector`: The solarthermal collector unit for which the calculation 
                                 is performed.
- `sim_params::Dict{String,Any}`: Simulation parameters.
- `target_temperature::Temperature`: The targeted output temperature of the solarthermal 
                                     collector.

# Returns
- `p_spec_th::Floathing`: The specific thermal power of the collector if target temperature 
                          was reached otherwise nothing.        
- `average_temperature::Temperature`: The average temperature of the collector if target   
                                      temperature was reached otherwise nothing.
- `target_temperature::Temperature`: The targeted output temperature of the collector if it  
                                     was reached otherwise nothing.
"""
function calc_thermal_power_fixed_flow_rate(unit::SolarthermalCollector, 
                                            sim_params::Dict{String,Any},
                                            target_temperature::Temperature)
    # go through the expected temperature levels and from high to low and check if the 
    # solarthermal collector can deliver energy at this level
    try
        average_temperature = find_zero((t_avg -> spec_thermal_power_func(t_avg, target_temperature,
                                                                          unit.last_average_temperature,
                                                                          unit.spec_flow_rate, unit,
                                                                          sim_params["time_step_seconds"]),
                                         t_avg -> derivate_spec_thermal_power_func(t_avg, unit.spec_flow_rate, unit,
                                                                                   sim_params["time_step_seconds"])),
                                        target_temperature - unit.delta_T_min,
                                        Roots.Newton())
        delta_T = (target_temperature - average_temperature) * 2
        if delta_T >= unit.delta_T_min
            p_spec_th = unit.spec_flow_rate * (target_temperature - average_temperature) * 2 * unit.vol_heat_cap
            return (p_spec_th,
                    average_temperature,
                    target_temperature)
        else
            return (nothing,
                    nothing,
                    nothing)
        end
    catch
        return (nothing,
                nothing,
                nothing)
    end
end

"""
    spec_thermal_power_func(t_avg::Temperature, t_target::Temperature, 
                            t_avg_last::Temperature, spec_flow_rate::Float64, 
                            unit::SolarthermalCollector, time_step:UInt64)

Function for calculating the thermal power output of a solarthermal collector.
Can be used to solve after t_avg and find average_temperature for specific used energy

# Arguments
- `t_avg::Temperature`: Average temperature of the collector. Value to solve after.
- `t_target::Temperature`: Target temperature of the collector 
- `t_avg_last::Temperature`: Average temperature of the collector from the last timestep
- `spec_flow_rate::Float64`: Specific flow rate through the collector in m³/(s*m²)
- `unit::SolarthermalCollector`: The solarthermal collector unit for which the calculation 
                                 is performed.
- `time_step::UInt64`: Time_step of the simulation in seconds

# Returns
- 'Float64': Difference between generated energy, losses and power output
"""
function spec_thermal_power_func(t_avg, t_target, t_avg_last, spec_flow_rate, unit::SolarthermalCollector, time_step)
    unit.eta_0_b * unit.K_b * unit.beam_solar_irradiance_in_plane +
    unit.eta_0_b * unit.K_d * unit.diffuse_solar_irradiance_in_plane -
    unit.a_params[1] * (t_avg - unit.ambient_temperature) -
    unit.a_params[2] * (t_avg - unit.ambient_temperature)^2 -
    unit.a_params[3] * unit.reduced_wind_speed * (t_avg - unit.ambient_temperature) +
    unit.a_params[4] * (unit.long_wave_irradiance - unit.sigma * (unit.ambient_temperature + 273.15)^4) -
    unit.a_params[5] * ((t_avg - t_avg_last) / time_step) -
    unit.a_params[6] * unit.reduced_wind_speed *
    (unit.beam_solar_irradiance_in_plane + unit.diffuse_solar_irradiance_in_plane) -
    unit.a_params[7] * unit.reduced_wind_speed *
    (unit.long_wave_irradiance - unit.sigma * (unit.ambient_temperature + 273.15)^4) -
    unit.a_params[8] * (t_avg - unit.ambient_temperature)^4 -
    spec_flow_rate * (t_target - t_avg) * 2 * unit.vol_heat_cap
end

"""
    spec_thermal_power_func(t_avg::Temperature, spec_flow_rate::Float64, 
                            unit::SolarthermalCollector, time_step:Number)

Derivative of the thermal power output function to speed up solve function for 
spec_thermal_power_func

# Arguments
- `t_avg::Temperature`: Average temperature of the collector. Value to solve after.
- `spec_flow_rate::Float64`: Specific flow rate through the collector in m³/(s*m²)
- `unit::SolarthermalCollector`: The solarthermal collector unit for which the calculation 
                                 is performed.
- `time_step::Number`: time_step of the simulation in seconds

# Returns
- 'Float64': Derivative of difference between generated energy, losses and power output
"""
function derivate_spec_thermal_power_func(t_avg, spec_flow_rate, unit::SolarthermalCollector, time_step)
    unit.a_params[1] * (-1) -
    unit.a_params[2] * 2 * (t_avg - unit.ambient_temperature) -
    unit.a_params[3] * unit.reduced_wind_speed -
    unit.a_params[5] * 1 / time_step -
    unit.a_params[8] * 4 * (t_avg - unit.ambient_temperature)^3 +
    spec_flow_rate * 2 * unit.vol_heat_cap
end

function output_values(unit::SolarthermalCollector)::Vector{String}
    return [string(unit.m_heat_out) * " OUT",
            "Temperature_Output",
            "Temperature_Mean_Collector",
            "direct_normal_irradiance",
            "beam_solar_irradiance_in_plane",
            "diffuse_solar_irradiance_in_plane",
            "delta_T",
            "spec_flow_rate"]
end

function output_value(unit::SolarthermalCollector, key::OutputKey)::Float64
    if key.value_key == "OUT"
        return calculate_energy_flow(unit.output_interfaces[key.medium])
    elseif key.value_key == "Temperature_Output"
        return unit.output_temperature
    elseif key.value_key == "Temperature_Mean_Collector"
        return unit.average_temperature
    elseif key.value_key == "direct_normal_irradiance"
        return unit.direct_normal_irradiance
    elseif key.value_key == "beam_solar_irradiance_in_plane"
        return unit.beam_solar_irradiance_in_plane
    elseif key.value_key == "diffuse_solar_irradiance_in_plane"
        return unit.diffuse_solar_irradiance_in_plane
    elseif key.value_key == "delta_T"
        return (unit.output_temperature - unit.average_temperature) * 2
    elseif key.value_key == "spec_flow_rate"
        return unit.spec_flow_rate_actual
    end
    throw(KeyError(key.value_key))
end

export SolarthermalCollector