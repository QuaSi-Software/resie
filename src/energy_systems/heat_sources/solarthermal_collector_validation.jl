import Interpolations as ip
using Roots
using Resie.SolarIrradiance

"""
Implementation of a solarthermal collector.
Works for flat plate collectors, vacuum tube collectors and PVT modules.
The calculation is based on EN ISO 9806:2017 for quasi-dynamic models.

Stagnation is ignored under the assumption that the system has either measures to prevent 
stagnation harming the collectors or the designed size is small enough in comparision to 
the global system for stagnation not to become a problem.

## ATTENTION: solarthermal heat collector is currently work in progress and not completed!!
"""

mutable struct SolarthermalCollectorVal <: Component
    # general
    uac::String
    controller::Controller
    sys_function::SystemFunction

    input_interfaces::InterfaceMap
    output_interfaces::InterfaceMap
    m_heat_out::Symbol
    # collector installation
    collector_gross_area::Float64
    tilt_angle::Float32
    azimuth_angle::Float32
    ground_reflectance::Float64
    # collector parameters
    eta_0_b::Float32
    K_b_array::Array{Union{Missing,Float32},2}
    K_b_itp::Any
    K_b::Float32
    K_d::Float32
    a_params::Array{Float32,1}
    sigma::Float64
    vol_heat_cap::Float64
    # weather profiles
    ambient_temperature_profile::Union{Profile,Nothing}
    beam_solar_hor_irradiance_profile::Union{Profile,Nothing}
    diffuse_solar_hor_irradiance_profile::Union{Profile,Nothing}
    long_wave_irradiance_profile::Union{Profile,Nothing}
    wind_speed_profile::Union{Profile,Nothing}
    # weather parameters
    wind_speed_reduction::Float32
    ambient_temperature::Temperature
    beam_solar_irradiance_in_plane::Float32
    diffuse_solar_irradiance_in_plane::Float32
    reduced_wind_speed::Floathing
    beam_solar_hor_irradiance::Floathing
    diffuse_solar_hor_irradiance::Floathing
    long_wave_irradiance::Floathing
    direct_normal_irradiance::Float32
    # operation parameters
    delta_T::Union{Float32,Nothing}
    spec_flow_rate::Union{Float64,Nothing}
    delta_T_min::Union{Float32,Nothing}
    spec_flow_rate_min::Union{Float64,Nothing}
    # results
    spec_thermal_power::Float64
    max_energy::Union{Float64,Nothing}
    used_energy::Float64
    output_temperature::Temperature
    average_temperature::Temperature
    last_average_temperature::Temperature
    spec_flow_rate_actual::Union{Float64,Nothing}
    # internal parameters for temperature layers
    available_energies::Dict{String,Floathing}
    average_temperatures::Dict{String,Temperature}
    output_temperatures::Dict{String,Temperature}
    runtimes::Dict{String,Floathing}
    temperature_energy_pairs::Dict{Temperature,Tuple}


    # for validation
    flow_rate_profile::Union{Profile,Nothing}
    input_temp_profile::Union{Profile,Nothing}
    input_temp::Float32
    zenith_angle::Float32
    angle_of_incidence::Float32
    beam_irr_profile::Union{Profile,Nothing}
    diff_irr_profile::Union{Profile,Nothing}
    dni_profile::Union{Profile,Nothing}
    spec_flow_rate_profile_value::Float64

    function SolarthermalCollectorVal(uac::String, config::Dict{String,Any}, sim_params::Dict{String,Any})
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
            (), # initialise Interpolation functions for K_b
            1, # Calculate K_b from K_b_array in control()
            config["K_d"], # collector parameter for diffuse irradiance 
            config["a_params"], # collector sim_params a1 to a8 corresponding to EN ISO 9806:2017
            5.670374419 * 10^-8, # Stefan Bolzmann Constant
            default(config, "vol_heat_capacity", 4.2e6), # volumnetric heat capacity of the fluid in the collector in [J/(m³*K)]
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
            # max(const_wind_speed * default(config, "wind_speed_reduction", 0.5) - 3, 0), # wind_speed multiplied with wind speed reduction [m/s]
            max(const_wind_speed * default(config, "wind_speed_reduction", 0.5), 0), # wind_speed multiplied with wind speed reduction [m/s]
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
            0.0, # specific thermal power of the solarthermal collector
            0.0, # maximum available energy
            0.0, # energy that was removed from the system in current time step
            0.0, # output temperature in current time step, calculated in control()
            0.0, # average temperature in current time step, calculated in control()
            0.0, # average temperature from last time step, calculated in control()
            0.0, # actual flow rate resulting from simulation and used energy
            # internal parameters for temperature layers
            Dict{String,Float64}(), # list of available_energies in for each component connected to the collector
            Dict{String,Temperature}(), # list of average_temperatures in for each component connected to the collector
            Dict{String,Temperature}(), # list of output_temperatures in for each component connected to the collector
            Dict{String,Float64}(), # list of runtimes at different temperature levels for each component connected to the collector
            Dict{Temperature,Tuple}(), # list of all calculated output temperature and output energy pairs
                
            # for validation
            nothing,
            nothing,
            -1.5,
            0.0,
            0.0,
            nothing,
            nothing,
            nothing,
            0.0
            )
    end
end

function initialise!(unit::SolarthermalCollectorVal, sim_params::Dict{String,Any})
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

    unit.K_b_itp = init_K_b(unit.K_b_array)

    # for validation
    # unit.flow_rate_profile = Profile("./profiles/validation/TRNSYS_flow_rate_l_h.prf", sim_params)
    # unit.flow_rate_profile = Profile("./profiles/validation/Dannenmannstr/absorber_flow_rate_calc_l_h.prf", sim_params)
    # unit.flow_rate_profile = Profile("./profiles/validation/RSTV/WMZ_flow_rate_l_h.prf", sim_params)
    unit.flow_rate_profile = Profile("./profiles/validation/RSTV/calc_flow_rate_l_h.prf", sim_params)
    # unit.input_temp_profile = Profile("./profiles/validation/TRNSYS_temperature_in.prf", sim_params)
    # unit.input_temp_profile = Profile("./profiles/validation/Dannenmannstr/absorber_temperature_in.prf", sim_params)
    unit.input_temp_profile = Profile("./profiles/validation/RSTV/all_temperature_in.prf", sim_params)
    # unit.input_temp_profile = Profile("./profiles/validation/RSTV/calc_temperature_in_haus4a.prf", sim_params)
    # unit.beam_irr_profile = Profile("./profiles/validation/TRNSYS_beam_irr.prf", sim_params, shift=Dates.Second(3600))
    # unit.diff_irr_profile = Profile("./profiles/validation/TRNSYS_diff_irr.prf", sim_params, shift=Dates.Second(3600))
    # unit.dni_profile = Profile("./profiles/validation/TRNSYS_dni.prf", sim_params)

end

function control(unit::SolarthermalCollectorVal,
                 components::Grouping,
                 sim_params::Dict{String,Any})
    update(unit.controller)

    # for validation
    unit.spec_flow_rate = Profiles.value_at_time(unit.flow_rate_profile, sim_params) / 
                          3600 / 1000 / unit.collector_gross_area
    # unit.spec_flow_rate = ifelse(unit.spec_flow_rate < 2.78e-5 * 0.8, 2.78e-5, unit.spec_flow_rate)
    runtime_percent = 1
    # runtime_percent = (Profiles.value_at_time(unit.flow_rate_profile, sim_params) / 3600 / unit.collector_gross_area
    #     ) / unit.spec_flow_rate

    unit.input_temp = Profiles.value_at_time(unit.input_temp_profile, sim_params)

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
        mean_ambient_temperature = _mean(collect(values(unit.ambient_temperature_profile.data)))
    else
        mean_ambient_temperature = unit.ambient_temperature
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
    solar_zenith, solar_azimuth = 
    irr_in_plane(sim_params, unit.tilt_angle, unit.azimuth_angle, unit.beam_solar_hor_irradiance, 
                 unit.diffuse_solar_hor_irradiance, nothing, 1.0, mean_ambient_temperature, 
                 unit.ground_reflectance
                 )

    # for validation
    # unit.beam_solar_irradiance_in_plane = (Profiles.value_at_time(unit.beam_irr_profile, sim_params))
    # unit.diffuse_solar_irradiance_in_plane = (Profiles.value_at_time(unit.diff_irr_profile, sim_params))
    # unit.direct_normal_irradiance = (Profiles.value_at_time(unit.dni_profile, sim_params))
    unit.spec_flow_rate_profile_value = unit.spec_flow_rate

    unit.zenith_angle = solar_zenith
    unit.angle_of_incidence = angle_of_incidence
    
    # Calculate K_b based on a table with provided values for longitudinal and transversal angles.
    # Table with provided values is mirrored for negative angles.
    if abs(transversal_angle) < 90 && abs(longitudinal_angle) < 90
        unit.K_b = unit.K_b_itp[1](abs(transversal_angle)) * unit.K_b_itp[2](abs(longitudinal_angle))
    else
        unit.K_b = 0
    end

    # reset list of all calculated output temperature and output energy pairs; used as lookup table to speed up calculation
    # calculation temperature when no energy is taken from collector
    unit.temperature_energy_pairs = Dict{Temperature,Tuple}()
    current_max_temperature = find_zero((t_avg -> spec_thermal_power_func(t_avg, t_avg, unit.last_average_temperature,
                                                                          0, unit, sim_params["time_step_seconds"]),
                                         t_avg -> derivate_spec_thermal_power_func(t_avg, 0, unit,
                                                                                   sim_params["time_step_seconds"])),
                                        unit.last_average_temperature, 
                                        Roots.Newton())
    unit.temperature_energy_pairs[current_max_temperature] = (0, current_max_temperature, current_max_temperature)
 
    unit.available_energies, has_calculated_all_maxima, unit.average_temperatures, 
    unit.output_temperatures, unit.runtimes = calculate_energies(unit, components, sim_params)

    unit.max_energy = ifelse(all(isnothing.(values(unit.available_energies))), 0, 
                             _sum(collect(values(unit.available_energies)))) * runtime_percent
    unit.output_temperature = highest(collect(values(unit.output_temperatures)))
    
    set_max_energy!(unit.output_interfaces[unit.m_heat_out], 
                    collect(values(sort(unit.available_energies))) * runtime_percent, 
                    [nothing for _ in 1:length(unit.available_energies)],
                    collect(values(sort(unit.output_temperatures))), 
                    collect(keys(sort(unit.available_energies))), 
                    has_calculated_all_maxima, 
                    false)
end

function process(unit::SolarthermalCollectorVal, sim_params::Dict{String,Any})
    exchanges = balance_on(unit.output_interfaces[unit.m_heat_out], unit.output_interfaces[unit.m_heat_out].target)
    energy_demands = [abs(e.balance + e.energy_potential) for e in exchanges]
    
    output_temperature_interface = lowest(unit.output_temperature, temp_min_highest(exchanges))
    unit.used_energy = 0.0

    if sum(energy_demands) > 0 && _sum(collect(values(unit.available_energies))) > 0
        for e in exchanges
            comp_uac = e.purpose_uac
            comp_energy_demand = abs(e.balance + e.energy_potential)
            used_energy_comp = min(comp_energy_demand, unit.available_energies[comp_uac])

            # recalculate the average_temperature if not all energy is used to keep a constant 
            # flow rate
            if unit.delta_T === nothing && unit.spec_flow_rate !== nothing &&
               used_energy_comp < unit.available_energies[comp_uac]
               # end of condition
                unit.average_temperatures[comp_uac] = e.temperature_min -
                                                      sim_params["wh_to_watts"](used_energy_comp) /
                                                      unit.collector_gross_area /
                                                      (unit.spec_flow_rate * 2 * unit.vol_heat_cap)
            end

            unit.used_energy += used_energy_comp
            unit.runtimes[comp_uac] = comp_energy_demand / unit.available_energies[comp_uac] * unit.runtimes[comp_uac]
        end

        unit.runtimes[unit.uac] = 0
        add!(unit.output_interfaces[unit.m_heat_out],
            unit.used_energy,
             output_temperature_interface)

    else
        unit.used_energy = 0.0
        set_max_energy!(unit.output_interfaces[unit.m_heat_out], 0.0)
    end
    
    for uac in keys(unit.runtimes)
        if !(uac in [e.purpose_uac for e in exchanges] || uac == unit.uac) && !isnothing(unit.average_temperatures[uac])
            unit.runtimes[uac] = 0
        elseif isnothing(unit.average_temperatures[uac])
            unit.runtimes[uac] = nothing
        end
    end

    unit.output_temperature = _weighted_mean(collect(values(sort(unit.output_temperatures))),
                                             collect(values(sort(unit.runtimes))))
    unit.average_temperature = _weighted_mean(collect(values(sort(unit.average_temperatures))),
                                              collect(values(sort(unit.runtimes))))
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
    calculate_energies(unit::SolarthermalCollectorVal, components, sim_params::Dict{String,Any})

Function used in the control() function of the solarthermal collector to calculate the 
available energeries, average temperatures, output temperatures and partial runtime for each 
connected component. Results are written as dictionares with the uac of the component as key.

# Arguments
- `unit::SolarthermalCollectorVal`: The solarthermal collector for which the calculation is 
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
function calculate_energies(unit::SolarthermalCollectorVal, components, sim_params::Dict{String,Any})
    available_energies = Dict{String,Floathing}()
    average_temperatures = Dict{String,Temperature}()
    output_temperatures = Dict{String,Temperature}()
    runtimes = Dict{String,Floathing}()
    left_energy_factor = 1.0
    component_idx = 1

    exchanges = balance_on(unit.output_interfaces[unit.m_heat_out], unit.output_interfaces[unit.m_heat_out].target)

    potential_energies = Float64[]
    uacs = Stringing[]
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
            #TODO for validation
            target_temperature = exchange.temperature_min
            max_energy = abs(exchange.balance + exchange.energy_potential)
        end

        push!(uacs, exchange.purpose_uac)
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
            available_energies[uacs[component_idx]] = produced_energy
            runtimes[uacs[component_idx]] = 1

        elseif produced_energy == 0
            available_energies[uacs[component_idx]] = produced_energy
            runtimes[uacs[component_idx]] = nothing

        elseif produced_energy > potential_energies[component_idx]
            available_energies[uacs[component_idx]] = potential_energies[component_idx]
            left_energy_factor *= (produced_energy - potential_energies[component_idx]) / produced_energy
            runtimes[uacs[component_idx]] = (produced_energy - potential_energies[component_idx]) / produced_energy

        else
            available_energies[uacs[component_idx]] = produced_energy
            runtimes[uacs[component_idx]] = left_energy_factor
            left_energy_factor = 0
        end

        average_temperatures[uacs[component_idx]] = t_avg
        output_temperatures[uacs[component_idx]] = t_target
        component_idx += 1
    end

    # if no energy at expected temperature levels is available then use max_output_temperature 
    # which equals average temperature with a used energy of 0
    if _sum(collect(values(available_energies))) == 0
        _, average_temperature = get_output_temperature_bounds(unit)

        calc_max = false
        average_temperatures[unit.uac] = average_temperature
        map!(v -> nothing, values(output_temperatures))
        output_temperatures[unit.uac] = average_temperature
        map!(v -> nothing, values(runtimes))
        runtimes[unit.uac] = 1

    elseif any(isinf, potential_energies)
        _, average_temperature = get_output_temperature_bounds(unit)

        calc_max = true
        average_temperatures[unit.uac] = average_temperature
        output_temperatures[unit.uac] = average_temperature
        runtimes[unit.uac] = 1
    else
        calc_max = false
        average_temperatures[unit.uac] = nothing
        output_temperatures[unit.uac] = nothing
        runtimes[unit.uac] = nothing
    end

    return (available_energies,
        calc_max,
        average_temperatures,
        output_temperatures,
            runtimes)
end

"""
    get_input_temperature_bounds(unit::SolarthermalCollectorVal)

Function used by the control module negotiate_temperature to get the current minimum 
and maximum temperatures that can be supplied.

# Arguments
- `unit::SolarthermalCollectorVal`: The solarthermal collector for which the calculation is performed.

# Returns
- `output_min_temperature::Temperature`: The current minimal temperature that can be given 
                                         in the output
- `output_max_temperature::Temperature`: The current maximum temperature that can be given 
                                         in the output

"""
function get_output_temperature_bounds(unit::SolarthermalCollectorVal)::Tuple{Temperature,Temperature}
    temps = collect(keys(unit.temperature_energy_pairs))
    zero_energy_mask = [v[1] == 0 for v in values(unit.temperature_energy_pairs)]
    return nothing, lowest(temps[zero_energy_mask])
end

"""
    calculate_output_energy_from_output_temperature(unit::SolarthermalCollectorVal,
                                                    minimum_output_temperature::Temperature, 
                                                    sim_params::Dict{String,Any})

Calculates the maximum energy that can be taken out of the collector field while not undercutting the given minimum_output_temperature.
If no converging result can be found, the maximum energy is set to the maximum energy of the collector field.
The maximum energy is limited to the user-input of the max_energy.

# Arguments
- `unit::SolarthermalCollectorVal`: the unit struct of the solarthermal collector with all its 
                                 parameters
- `minimum_output_temperature::Temperature`: The desired minimum output temperature of the 
                                             fluid temperature
- `sim_params::Dict{String,Any}`: simulation parameters
# Returns
- `current_max_output_energy::Float64`: the maximum energy that can be taken out of the  
                                        solarthermal collector while not undercutting the 
                                        given minimum_output_temperature.     
"""
function calculate_output_energy_from_output_temperature(unit::SolarthermalCollectorVal,
                                                         minimum_output_temperature::Temperature,
                                                         sim_params::Dict{String,Any})::Float64

    # calculate specific thermal power depending on control mode; reuse already caluclated values if possible
    if minimum_output_temperature in keys(unit.temperature_energy_pairs)
        p_spec_th, t_avg, t_target = unit.temperature_energy_pairs[minimum_output_temperature]

    elseif unit.delta_T !== nothing && unit.spec_flow_rate === nothing
        p_spec_th, t_avg, t_target = calc_thermal_power_fixed_delta_T(unit, sim_params, minimum_output_temperature)

    elseif unit.delta_T === nothing && unit.spec_flow_rate !== nothing
        # p_spec_th, t_avg, t_target = calc_thermal_power_fixed_flow_rate(unit, sim_params, minimum_output_temperature)
        p_spec_th, t_avg, t_target = calc_thermal_power_fixed_flow_input_temp!(unit, sim_params, minimum_output_temperature)

    else
        @error "Error in config file: Exclusively delta_T OR spec_flow_rate must have a value"
        throw(InputError)
    end

    unit.temperature_energy_pairs[minimum_output_temperature] = (p_spec_th, t_avg, t_target)
    max_energy = isnothing(p_spec_th) ? 0.0 : sim_params["watt_to_wh"](p_spec_th * unit.collector_gross_area)
                                                 
    return max_energy
end

"""
    check_temperature_and_get_max_energy(unit::SolarthermalCollectorVal,
                                         sim_params::Dict{String,Any},
                                         temperature_output::Temperature)

Checks if a requested output temperature is within the allowed operational range for the 
solarthermal collector. If no temperautre is given, the function sets the temperature to the 
mean temperature that will be reached with maximum energy draw during the whole time step.
The function also determines the maximum extractable energy for the given temperature.

# Arguments
- `unit::SolarthermalCollectorVal`: The solarthermal collector unit for which the calculation 
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
function check_temperature_and_get_max_energy(unit::SolarthermalCollectorVal,
                                              sim_params::Dict{String,Any},
                                              temperature_output::Temperature,
                                              limit_max_output_energy_to_avoid_pulsing::Bool)::Tuple{Temperature,Float64}                               
    # get max output temperature of solarthermal collector
    _, source_max_out_temperature = get_output_temperature_bounds(unit)

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
    init_K_b(K_b_array::Array{Union{Missing,Float32},2})

Calculate K_b based on a table with provided values for longitudinal and transversal angles 
and the angle of irradiance on the plane of the solarthermal collector.
Table with provided values is mirrored for negative angles.

# Arguments
- `K_b_array::Array{Union{Missing,Float32},2}`: Simulation parameters.

# Returns
Returns a tuple containing two interpolation objects:
- interpolation object for transversal K_b value.
- interpolation object for longitudinal K_b value.
"""
function init_K_b(K_b_array::Array{Union{Missing,Float32},2})
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
    calc_thermal_power_fixed_delta_T(unit::SolarthermalCollectorVal, 
                                     sim_params::Dict{String,Any}, 
                                     target_temperature::Temperature)

Calculate the specific thermal power of the solarthermal collector under the assumption that 
the collector has a fixed delta_T and the flow rate is variable 

# Arguments
- `unit::SolarthermalCollectorVal`: The solarthermal collector unit for which the calculation 
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
function calc_thermal_power_fixed_delta_T(unit::SolarthermalCollectorVal, 
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
    calc_thermal_power_fixed_delta_T(unit::SolarthermalCollectorVal, 
                                     sim_params::Dict{String,Any}, 
                                     target_temperature::Temperature)

Calculate the specific thermal power of the solarthermal collector under the assumption that 
the collector has a fixed flow rate and the delta_T is variable 

# Arguments
- `unit::SolarthermalCollectorVal`: The solarthermal collector unit for which the calculation 
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
function calc_thermal_power_fixed_flow_rate(unit::SolarthermalCollectorVal, sim_params::Dict{String,Any}, target_temperature)
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
Calculate the specific thermal power of the solarthermal collector under the assumption that 
the collector has a fixed flow rate and the delta_T is variable.
Simulates smaller timesteps to decrease influence of thermal capacity.
"""
function calc_thermal_power_fixed_flow_rate_sub!(unit::SolarthermalCollectorVal, sim_params::Dict{String,Any})
    t_in_min = 0

    # go through the expected temperature levels and from high to low and check if the 
    # solarthermal collector can deliver energy at this level
    t_avg_last = unit.last_average_temperature
    sub_time_step = 900
    average_temperature = 0.0
    try
        for time_step in sub_time_step:sub_time_step:sim_params["time_step_seconds"]
            average_temperature = find_zero(
                (t_avg->spec_thermal_power_func(t_avg, target_temperature, t_avg, unit.spec_flow_rate, unit, sub_time_step), 
                t_avg->derivate_spec_thermal_power_func(t_avg, unit.spec_flow_rate, unit, sub_time_step)
                ), 
                target_temperature - unit.delta_T_min, 
                Roots.Newton()
                )

            delta_T = (target_temperature - average_temperature) * 2
            if delta_T >= unit.delta_T_min && (target_temperature - delta_T) >= t_in_min
                t_avg_last = average_temperature
            else
                average_temperature = find_zero(
                    (t_avg->spec_thermal_power_func(t_avg, t_avg, t_avg, 0, unit, sub_time_step), 
                        t_avg->derivate_spec_thermal_power_func(t_avg, 0, unit, sub_time_step)
                    ), 
                    t_avg_last, 
                    Roots.Newton()
                    )
                t_avg_last = average_temperature
            end
        end
        delta_T = (target_temperature - average_temperature) * 2
        if delta_T >= unit.delta_T_min && (target_temperature - delta_T) >= t_in_min
            p_spec_th = unit.spec_flow_rate * (target_temperature - average_temperature) * 2 * unit.vol_heat_cap
            unit.average_temperature = average_temperature
            unit.output_temperature = target_temperature
            return (
                p_spec_th,
                average_temperature
            )
        else
            return (
                0,
                nothing,
                nothing
            )
        end
    catch
        return (
            0,
            nothing,
            nothing
        )
    end
end

# for validation
"""
Calculate the specific thermal power of the solarthermal collector under the assumption that 
the collector has a fixed flow rate and the input temperature is set additionaly 
"""
function calc_thermal_power_fixed_flow_input_temp!(unit, sim_params, target_temperature)
    # go through the expected temperature levels and from high to low and check if the 
    # solarthermal collector can deliver energy at this level

    try
        average_temperature = find_zero(
            (t_avg->spec_thermal_power_func_input_temp(t_avg, unit.input_temp, unit.last_average_temperature, unit.spec_flow_rate, unit, sim_params["time_step_seconds"]), 
            t_avg->derivate_spec_thermal_power_func_input_temp(t_avg, unit.spec_flow_rate, unit, sim_params["time_step_seconds"])
            ), 
            unit.input_temp + unit.delta_T_min, 
            Roots.Newton()
            )

        delta_T = (average_temperature - unit.input_temp) * 2
        # if true # TODO for debugging
        if unit.spec_flow_rate > 0
            p_spec_th = unit.spec_flow_rate * delta_T * unit.vol_heat_cap
            return (
                p_spec_th,
                average_temperature,
                unit.input_temp + delta_T
            )
        else
            return (
                nothing,
                nothing,
                nothing
            )
        end
    catch
        try
            average_temperature = find_zero(
                t_avg -> spec_thermal_power_func_input_temp(t_avg, unit.input_temp, unit.last_average_temperature, unit.spec_flow_rate, unit, sim_params["time_step_seconds"]), 
                (unit.input_temp-500, unit.input_temp+500),
                Bisection())

            delta_T = (average_temperature - unit.input_temp) * 2
            # if true # TODO for debugging
            if unit.spec_flow_rate > 0
                p_spec_th = unit.spec_flow_rate * delta_T * unit.vol_heat_cap
                return (
                    p_spec_th,
                    average_temperature,
                    unit.input_temp + delta_T
                )
            else
                return (
                    nothing,
                    nothing,
                    nothing
                )
            end
        catch
            @error "Solar thermal collector could not be calculated check inputs ->
            time_step:$(sim_params["current_date"]); 
            beam_solar_hor_irradiance:$(sim_params["wh_to_watts"](Profiles.value_at_time(
                unit.beam_solar_hor_irradiance_profile, sim_params)));
            beam_solar_irradiance_in_plane:$(unit.beam_solar_irradiance_in_plane); 
            diffuse_solar_irradiance_in_plane:$(unit.diffuse_solar_irradiance_in_plane); 
            spec_flow_rate:$(unit.spec_flow_rate);
            long_wave_irradiance:$(unit.long_wave_irradiance);
            collector_temperature_last_timestep:$(unit.last_average_temperature)"
            throw(InputError)
        end
    end
end

# for validation
function calc_thermal_power_fixed_flow_input_temp_old!(unit, sim_params, target_temperature)
    average_temperature = find_zero(
        t_avg -> spec_thermal_power_func_input_temp_old(t_avg, unit.input_temp, unit.last_average_temperature, unit.spec_flow_rate, unit, sim_params["time_step_seconds"]), 
        (unit.input_temp-500, unit.input_temp+500),
        Bisection())

    delta_T = (average_temperature - unit.input_temp) * 2
    # if delta_T > unit.delta_T_min
    if true # TODO for debugging
        p_spec_th = unit.spec_flow_rate * delta_T * unit.vol_heat_cap
        return (
            p_spec_th,
            average_temperature,
            unit.input_temp + delta_T
        )
    end
end

# for validation
"""
Function for calculating the thermal power output of a solarthermal collector.
Can be used to solve after t_avg and find average_temperature for specific used energy
"""
function spec_thermal_power_func_input_temp(t_avg, t_input, t_avg_last, spec_flow_rate, unit, time_step)     
    unit.eta_0_b * unit.K_b * unit.beam_solar_irradiance_in_plane + 
    unit.eta_0_b * unit.K_d * unit.diffuse_solar_irradiance_in_plane -
    unit.a_params[1] * (t_avg - unit.ambient_temperature) -
    unit.a_params[2] * (t_avg - unit.ambient_temperature)^2 -
    unit.a_params[3] * unit.reduced_wind_speed * (t_avg - unit.ambient_temperature) +
    unit.a_params[4] * (unit.long_wave_irradiance - unit.sigma * (unit.ambient_temperature + 273.15)^4) -
    unit.a_params[5] * ((t_avg-t_avg_last) / time_step) - 
    unit.a_params[6] * unit.reduced_wind_speed * (unit.beam_solar_irradiance_in_plane + unit.diffuse_solar_irradiance_in_plane) -
    unit.a_params[7] * unit.reduced_wind_speed * (unit.long_wave_irradiance - unit.sigma * (unit.ambient_temperature + 273.15)^4) -
    unit.a_params[8] * (t_avg - unit.ambient_temperature)^4 -
    spec_flow_rate * (t_avg - t_input) * 2 * unit.vol_heat_cap
end
# for validation
"""
Derivative of the thermal power output fuction to speed up solve function
"""
function derivate_spec_thermal_power_func_input_temp(t_avg, spec_flow_rate, unit, time_step)
    unit.a_params[1] * (-1) -
    unit.a_params[2] * 2 * (t_avg - unit.ambient_temperature) -
    unit.a_params[3] * unit.reduced_wind_speed -
    unit.a_params[5] * 1 / time_step -
    unit.a_params[8] * 4 * (t_avg - unit.ambient_temperature)^3 -
    spec_flow_rate * 2 * unit.vol_heat_cap
end

# for validation old benchmarks
"""
Function for calculating the thermal power output of a solarthermal collector.
Can be used to solve after t_avg and find average_temperature for specific used energy
"""
function spec_thermal_power_func_input_temp_old(t_avg, t_input, t_avg_last, spec_flow_rate, unit, time_step)   
    # eta_0_hem = 0.162
    # b_params = [10.04, 0.81, 0.152]
    eta_0_hem = 0.702
    b_params = [32.64, 3.59, 0.066]
    
    f = 
    eta_0_hem * (1-b_params[3] * (unit.reduced_wind_speed)) * ((unit.beam_solar_irradiance_in_plane + unit.diffuse_solar_irradiance_in_plane) + 0.85 * (unit.long_wave_irradiance - 5.670374419e-8 * (unit.ambient_temperature + 273.15)^4)) -
    (b_params[1] + b_params[2] * (unit.reduced_wind_speed)) * (t_avg - unit.ambient_temperature) -
    spec_flow_rate * (t_avg - t_input) * 2 * unit.vol_heat_cap
end

"""
Function for calculating the thermal power output of a solarthermal collector.
Can be used to solve after t_avg and find average_temperature for specific used energy
"""
function spec_thermal_power_func(t_avg, t_target, t_avg_last, spec_flow_rate, unit::SolarthermalCollectorVal, time_step)         
    unit.eta_0_b * unit.K_b * unit.beam_solar_irradiance_in_plane + 
    unit.eta_0_b * unit.K_d * unit.diffuse_solar_irradiance_in_plane -
    unit.a_params[1] * (t_avg - unit.ambient_temperature) -
    unit.a_params[2] * (t_avg - unit.ambient_temperature)^2 -
    unit.a_params[3] * unit.reduced_wind_speed * (t_avg - unit.ambient_temperature) +
    unit.a_params[4] * (unit.long_wave_irradiance - unit.sigma * (unit.ambient_temperature + 273.15)^4) -
    unit.a_params[5] * ((t_avg-t_avg_last) / time_step) -
    unit.a_params[6] * unit.reduced_wind_speed * (unit.beam_solar_irradiance_in_plane + unit.diffuse_solar_irradiance_in_plane) -
    unit.a_params[7] * unit.reduced_wind_speed * (unit.long_wave_irradiance - unit.sigma * (unit.ambient_temperature + 273.15)^4) -
    unit.a_params[8] * (t_avg - unit.ambient_temperature)^4 -
    spec_flow_rate * (t_target - t_avg) * 2 * unit.vol_heat_cap
end

"""
Derivative of the thermal power output fuction to speed up solve function
"""
function derivate_spec_thermal_power_func(t_avg, spec_flow_rate, unit::SolarthermalCollectorVal, time_step)
    unit.a_params[1] * (-1) -
    unit.a_params[2] * 2 * (t_avg - unit.ambient_temperature) -
    unit.a_params[3] * unit.reduced_wind_speed -
    unit.a_params[5] * 1 / time_step -
    unit.a_params[8] * 4 * (t_avg - unit.ambient_temperature)^3 +
    spec_flow_rate * 2 * unit.vol_heat_cap
end

function output_values(unit::SolarthermalCollectorVal)::Vector{String}

    return [string(unit.m_heat_out)*" OUT",
            "Temperature", 
            "Max_Energy",
            "Average_Temperature",
            "Ambient_Temperature",
            "Used_Energy",
            "direct_normal_irradiance",
            "beam_solar_irradiance_in_plane",
            "diffuse_solar_irradiance_in_plane",
            "delta_T",
            "spec_flow_rate",
            "spec_flow_rate_profile",
            "zenith_angle",
            "angle_of_incidence"
            ]
end

function output_value(unit::SolarthermalCollectorVal, key::OutputKey)::Float64
    if key.value_key == "OUT"
        return calculate_energy_flow(unit.output_interfaces[key.m_heat_out])
    elseif key.value_key == "Temperature"
        return unit.output_temperature
    elseif key.value_key == "Max_Energy"
        return unit.max_energy
    elseif key.value_key == "Average_Temperature"
        return unit.average_temperature
    elseif key.value_key == "Ambient_Temperature"
        return unit.ambient_temperature
    elseif key.value_key == "Used_Energy"
        return unit.used_energy
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
    elseif key.value_key == "spec_flow_rate_profile"
        return unit.spec_flow_rate_profile_value
    elseif key.value_key == "zenith_angle"
        return unit.zenith_angle
    elseif key.value_key == "angle_of_incidence"
        return unit.angle_of_incidence
    end
    throw(KeyError(key.value_key))
end

export SolarthermalCollectorVal
