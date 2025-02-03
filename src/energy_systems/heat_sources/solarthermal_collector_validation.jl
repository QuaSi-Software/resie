using Interpolations
using Roots

include("solar_irradiance.jl")

"""
Implementation of a solarthermal collector.
Works for flat plate collectors, vacuum tube collectors and PVT modules.
The calculation is based on EN ISO 9806:2017 for quasi-dynamic models.

Stagnation is ignored under the assumption that the system has either measures to prevent 
stagnation harming the collectors or the designed size is small enough in comparision to 
the global system for stagnation not to become a problem.

## ATTENTION: Geothermal heat collector is currently work in progress and not completed!!
"""

mutable struct SolarthermalCollectorVal <: Component
    uac::String
    controller::Controller
    sys_function::SystemFunction

    input_interfaces::InterfaceMap
    output_interfaces::InterfaceMap

    medium::Symbol

    collector_gross_area::Float64
    tilt_angle::Float32
    azimuth_angle::Float32

    eta_0_b::Float32
    K_b_array::Array{Union{Missing, Float32}, 2}
    K_b_itp
    K_b::Float32
    K_d::Float32
    a_params::Array{Float32, 1}
    sigma::Float64

    ambient_temperature_profile::Union{Profile,Nothing}
    global_solar_hor_irradiance_profile::Union{Profile,Nothing}
    diffuse_solar_hor_irradiance_profile::Union{Profile,Nothing}
    long_wave_irradiance_profile::Union{Profile,Nothing}
    wind_speed_profile::Union{Profile,Nothing}

    wind_speed_reduction::Float32

    ambient_temperature::Temperature
    beam_solar_irradiance_in_plane::Float32
    diffuse_solar_irradiance_in_plane::Float32
    reduced_wind_speed::Float32
    long_wave_irradiance::Float32
    direct_normal_irradiance::Float32

    delta_T::Union{Float32,Nothing}
    spec_flow_rate::Union{Float64,Nothing}
    delta_T_min::Union{Float32,Nothing}
    spec_flow_rate_min::Union{Float32,Nothing}

    spec_thermal_power::Float64
    max_energy::Union{Float64,Nothing}
    used_energy::Float64
    output_temperature::Temperature
    average_temperature::Temperature
    last_average_temperature::Temperature

    available_energies::Dict{String, Floathing}
    average_temperatures::Dict{String, Temperature}
    output_temperatures::Dict{String, Temperature}
    runtimes::Dict{String, Floathing}

    vol_heat_cap::Float64

    # for validation
    flow_rate_profile::Union{Profile,Nothing}
    input_temp_profile::Union{Profile,Nothing}
    input_temp::Float32

    function SolarthermalCollectorVal(uac::String, config::Dict{String,Any}, sim_params::Dict{String,Any})
        ambient_temperature_profile = 
            "ambient_temperature_profile_file_path" in keys(config) ?
            Profile(config["ambient_temperature_profile_file_path"], sim_params) :
            nothing
        global_solar_hor_irradiance_profile = 
            "global_solar_hor_irradiance_profile_file_path" in keys(config) ?
            Profile(config["global_solar_hor_irradiance_profile_file_path"], sim_params) :
            nothing
        diffuse_solar_hor_irradiance_profile = 
            "diffuse_solar_hor_irradiance_profile_file_path" in keys(config) ?
            Profile(config["diffuse_solar_hor_irradiance_profile_file_path"], sim_params) :
            nothing
        long_wave_irradiance_profile = 
            "long_wave_irradiance_profile_file_path" in keys(config) ?
            Profile(config["long_wave_irradiance_profile_file_path"], sim_params) :
            nothing
        wind_speed_profile = 
            "wind_speed_profile_file_path" in keys(config) ?
            Profile(config["wind_speed_profile_file_path"], sim_params) :
            nothing

        medium = Symbol(default(config, "medium", "m_h_w_ht1"))
        register_media([medium])

        return new(                                                
            uac, # uac
            controller_for_strategy( # controller
                config["strategy"]["name"], config["strategy"], sim_params
            ),
            sf_bounded_source, # sys_function
            InterfaceMap(), # input_interfaces
            InterfaceMap( # output_interfaces
                medium => nothing
            ),
            medium, # medium name of output interface
            
            config["collector_gross_area"], # gross area of collector
            config["tilt_angle"], # tilt angle
            config["azimuth_angle"], # azimuth angle or orientation angle
        
            config["eta_0_b"],
            vcat([10,20,30,40,50,60,70,80,90]', config["K_b_t_array"]', config["K_b_l_array"]'), 
                                # row1: angle 0° to 90°/vertical/west, 
                                # row2: transversal/tilt, 
                                # row3: longitudinal/azimuth/orientation
            (), # initialise Interpolation functions for K_b
                                # row1: angle 0° to 90°/vertical/west, 
                                # row2: transversal/tilt, 
                                # row3: longitudinal/azimuth/orientation
            1, # Calculate K_b from K_b_array in control()
            config["K_d"], # collector parameter for diffuse irradiance 
            config["a_params"], # collector sim_params a1 to a8 corresponding to EN ISO 9806:2017
            5.670374419*10^-8, # Stefan Bolzmann Constant
        
            ambient_temperature_profile,
            global_solar_hor_irradiance_profile,
            diffuse_solar_hor_irradiance_profile,
            long_wave_irradiance_profile,
            wind_speed_profile,

            default(config, "wind_speed_reduction", 0.5),
            
            0.0, # ambient temperature from profile
            0.0, # beam_solar_irradiance in collector plane, calculated from profile
            0.0, # diffuse_solar_irradiance_in_plane from profile
            0.0, # reduced wind_speed from profile multiplied with wind speed reduction
            0.0, # long_wave_irradiance from profile
            0.0, # direct_normal_irradiance calculated from sun position
            
            default(config, "delta_T", nothing), # delta_T between input and output temperature
            default(config, "spec_flow_rate", nothing), # nominal specific volume flow of the thermal collector
            default(config, "delta_T_min", 2), # minimal delta_T between input and output temperature for the collector to start producing energy; used together with spec_flow_rate
            default(config, "spec_flow_rate_min", 0.000002), # minimal specific volume flow of the thermal collector to start producing energy; used together with delta_T

            0.0, # specific thermal power of the solarthermal collector
            0.0, # maximum available energy
            0.0, # energy that was removed from the system in current time step
            0.0, # output temperature in current time step, calculated in control()
            0.0, # average temperature in current time step, calculated in control()
            0.0, # average temperature from last time step, calculated in control()

            Dict{String, Float64}(), # list of available_energies in for each component connected to the collector
            Dict{String, Temperature}(), # list of average_temperatures in for each component connected to the collector
            Dict{String, Temperature}(), # list of output_temperatures in for each component connected to the collector
            Dict{String, Float64}(), # list of runtimes at different temperature levels for each component connected to the collector
    
            default(config, "vol_heat_capacity", 4.2e6), # volumnetric heat capacity of the fluid in the collector in [J/(m³*K)]
            
            # for validation
            nothing,
            nothing,
            -1.5
            )
    end
end

function initialise!(unit::SolarthermalCollectorVal, sim_params::Dict{String,Any})
    set_storage_transfer!(
        unit.output_interfaces[unit.medium],
        default(
            unit.controller.parameter, "load_storages " * String(unit.medium), true
        )
    )
    unit.output_temperature = Profiles.value_at_time(
        unit.ambient_temperature_profile, sim_params["time"]
        )
    unit.average_temperature = Profiles.value_at_time(
        unit.ambient_temperature_profile, sim_params["time"]
        )
    unit.last_average_temperature = Profiles.value_at_time(
        unit.ambient_temperature_profile, sim_params["time"]
        )

    unit.K_b_itp = init_K_b(unit.K_b_array)
    # for validation
    unit.flow_rate_profile = Profile("./profiles/validation/RSTV/WMZ_flow_rate_m3_h.prf", sim_params)
    unit.input_temp_profile = Profile("./profiles/validation/RSTV/all_temperature_in.prf", sim_params)
end

function control(
    unit::SolarthermalCollectorVal,
    components::Grouping,
    sim_params::Dict{String,Any}
    )
    move_state(unit, components, sim_params)

    # for validation
    unit.spec_flow_rate = Profiles.value_at_time(
        unit.flow_rate_profile, sim_params["time"]
        ) / 3600 / unit.collector_gross_area
    # unit.spec_flow_rate = ifelse(unit.spec_flow_rate < 2.78e-5 * 0.8, 2.78e-5, unit.spec_flow_rate)
    # runtime_percent = 1
    runtime_percent = (Profiles.value_at_time(unit.flow_rate_profile, sim_params["time"]) / 3600 / unit.collector_gross_area
        ) / unit.spec_flow_rate

    unit.input_temp = Profiles.value_at_time(unit.input_temp_profile, sim_params["time"])

    # get values from profiles
    global_solar_hor_irradiance = Profiles.value_at_time(
        unit.global_solar_hor_irradiance_profile, sim_params["time"]
        )
        
    diffuse_solar_hor_irradiance = Profiles.value_at_time(
        unit.diffuse_solar_hor_irradiance_profile, sim_params["time"]
        )
    unit.diffuse_solar_irradiance_in_plane = diffuse_solar_hor_irradiance * ((1 + cosd(unit.tilt_angle)) / 2) 

    unit.ambient_temperature = Profiles.value_at_time(
        unit.ambient_temperature_profile, sim_params["time"]
        )

    unit.reduced_wind_speed = max(
        Profiles.value_at_time(unit.wind_speed_profile, sim_params["time"]) * unit.wind_speed_reduction - 3, 0)

    unit.long_wave_irradiance = Profiles.value_at_time(
        unit.long_wave_irradiance_profile, sim_params["time"]
        )

    start_date = DateTime(Dates.year(Dates.now()), 1, 1, 0, 30, 0) - Dates.Hour(1) # TODO centralise start_date
    time = start_date + Dates.Second(sim_params["time"])
    
    solar_zenith, solar_azimuth = sun_position(time, 9.18, 47.67, 1.0, unit.ambient_temperature)
    
    unit.beam_solar_irradiance_in_plane, unit.direct_normal_irradiance, angle_of_incidence, longitudinal_angle, transversal_angle = beam_irr_in_plane(
        unit.tilt_angle, unit.azimuth_angle, solar_zenith, solar_azimuth, 
        global_solar_hor_irradiance, diffuse_solar_hor_irradiance
        )
    
    # Calculate K_b based on a table with provided values for longitudinal and transversal angles.
    # Table with provided values is mirrored for negative angles.
    if abs(transversal_angle) < 90 && abs(longitudinal_angle) < 90
        unit.K_b = unit.K_b_itp[1](abs(transversal_angle)) * unit.K_b_itp[2](abs(longitudinal_angle))
    else
        unit.K_b = 0
    end
 
    unit.available_energies, has_calculated_all_maxima, unit.average_temperatures, unit.output_temperatures, unit.runtimes = calculate_energies(unit, sim_params)

    unit.max_energy = ifelse(all(isnothing.(values(unit.available_energies))), 0, _sum(collect(values(unit.available_energies))) * runtime_percent)
    unit.output_temperature = highest(collect(values(unit.output_temperatures)))
    
    set_max_energy!(unit.output_interfaces[unit.medium], collect(values(sort(unit.available_energies))) .* runtime_percent, collect(keys(sort(unit.available_energies))), has_calculated_all_maxima)
    set_temperature!(unit.output_interfaces[unit.medium], nothing, unit.output_temperature)
end

function process(unit::SolarthermalCollectorVal, sim_params::Dict{String,Any})

    exchanges = balance_on(unit.output_interfaces[unit.medium], unit.output_interfaces[unit.medium].target)
    energy_demands = [abs(e.balance + e.energy_potential) for e in exchanges]
    
    output_temperature_interface = lowest(unit.output_temperature, temp_min_highest(exchanges))

    if sum(energy_demands) > 0 && _sum(collect(values(unit.available_energies))) > 0
        for e in exchanges
            comp_uac = e.purpose_uac
            comp_energy_demand = abs(e.balance + e.energy_potential)
            used_energy_comp = min(comp_energy_demand, unit.available_energies[comp_uac])

            # recalculate the average_temperature if not all energy is used to keep a constant 
            # flow rate
            if unit.delta_T === nothing && unit.spec_flow_rate !== nothing && used_energy_comp < unit.available_energies[comp_uac]
                unit.average_temperatures[comp_uac] = e.temperature_min - wh_to_watts(used_energy_comp) / unit.collector_gross_area / (unit.spec_flow_rate * 2 * unit.vol_heat_cap)
            end

            unit.used_energy += used_energy_comp
            unit.runtimes[comp_uac] = comp_energy_demand / unit.available_energies[comp_uac] * unit.runtimes[comp_uac]
        end

        unit.runtimes[unit.uac] = 0
        add!(
            unit.output_interfaces[unit.medium],
            unit.used_energy,
            output_temperature_interface
        )

    else
        unit.used_energy = 0.0
        set_max_energy!(unit.output_interfaces[unit.medium], 0.0)
    end
    
    for uac in keys(unit.runtimes)
        if !(uac in [e.purpose_uac for e in exchanges] || uac == unit.uac) && !isnothing(unit.average_temperatures[uac])
            unit.runtimes[uac] = 0
        elseif isnothing(unit.average_temperatures[uac])
            unit.runtimes[uac] = nothing
        end
    end

    unit.output_temperature = _weighted_mean(collect(values(sort(unit.output_temperatures))), collect(values(sort(unit.runtimes))))
    unit.average_temperature = _weighted_mean(collect(values(sort(unit.average_temperatures))), collect(values(sort(unit.runtimes))))
    unit.last_average_temperature = unit.average_temperature
end

function calculate_energies(unit::SolarthermalCollectorVal, sim_params::Dict{String,Any})
    
    exchanges = balance_on(unit.output_interfaces[unit.medium], unit.output_interfaces[unit.medium].target)

    left_energy_factor = 1.0

    if sum([abs(e.balance + e.energy_potential) for e in exchanges]) == 0
        left_energy_factor = 0
    end

    available_energies = Dict{String, Floathing}()
    calculated_values = Dict()
    average_temperatures = Dict{String, Temperature}()
    output_temperatures = Dict{String, Temperature}()
    runtimes = Dict{String, Floathing}()

    uacs = [e.purpose_uac for e in exchanges]
    target_temperatures = temp_min_all(exchanges)

    if is_max_energy_nothing(unit.output_interfaces[unit.medium].max_energy)
        potential_energies = fill(Inf, length(exchanges))
    else
        potential_energies = [abs(e.balance + e.energy_potential) for e in exchanges]
    end

    component_idx = 1
   
    while component_idx <= length(exchanges) && left_energy_factor > 0
        
        # calculate specific thermal power depending on control mode; reuse already caluclated values if possible
        if target_temperatures[component_idx] in keys(calculated_values)
            p_spec_th, t_avg, t_target = calculated_values[target_temperatures[component_idx]]

        elseif unit.delta_T !== nothing && unit.spec_flow_rate === nothing
            p_spec_th, t_avg, t_target = calc_thermal_power_fixed_delta_T(unit, sim_params, target_temperatures[component_idx])

        elseif unit.delta_T === nothing && unit.spec_flow_rate !== nothing
            p_spec_th, t_avg, t_target = calc_thermal_power_fixed_flow_input_temp!(unit, sim_params, target_temperatures[component_idx])

        else
            @error "Error in config file: Exclusively delta_T OR spec_flow_rate must have a value"
            throw(InputError)
        end

        produced_energy = isnothing(p_spec_th) ? 0.0 : watt_to_wh(p_spec_th * unit.collector_gross_area) * left_energy_factor
        calculated_values[target_temperatures[component_idx]] = (p_spec_th, t_avg, t_target)

        # check if more energy can be produced for a certain component than is demanded
        if any(isinf, potential_energies)
            available_energies[uacs[component_idx]] = produced_energy
            runtimes[uacs[component_idx]] = 1

        elseif isnothing(p_spec_th)
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
        output_temperatures[uacs[component_idx]] =  t_target
        component_idx += 1
    end

    # if no energy at expected temperature levels is available then calculate 
    # average temperature with a used energy of 0
    if _sum(collect(values(available_energies))) == 0
        average_temperature = find_zero(
            (t_avg->spec_thermal_power_func(t_avg, t_avg, unit.last_average_temperature, 0, unit, sim_params["time_step_seconds"]), 
                t_avg->derivate_spec_thermal_power_func(t_avg, 0, unit, sim_params["time_step_seconds"])
            ), 
            unit.last_average_temperature, 
            Roots.Newton()
            )

        calc_max = false
        average_temperatures[unit.uac] = average_temperature
        map!(v->nothing, values(output_temperatures))
        output_temperatures[unit.uac] = average_temperature
        map!(v->nothing, values(runtimes))
        runtimes[unit.uac] = 1

    elseif any(isinf, potential_energies)
        average_temperature = find_zero(
            (t_avg->spec_thermal_power_func(t_avg, t_avg, unit.last_average_temperature, 0, unit, sim_params["time_step_seconds"]), 
                t_avg->derivate_spec_thermal_power_func(t_avg, 0, unit, sim_params["time_step_seconds"])
            ), 
            unit.last_average_temperature, 
            Roots.Newton()
            )

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

    return (
        available_energies,
        calc_max,
        average_temperatures,
        output_temperatures,
        runtimes
    )
end

"""
Calculate K_b based on a table with provided values for longitudinal and transversal angles 
and the angle of irradiance on the plane of the solarthermal collector.
Table with provided values is mirrored for negative angles.
"""
function init_K_b(K_b_array)
    K_b_array = hcat([0, 1, 1], K_b_array)
    K_b_array[:,10] = [90,0,0]
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

        itp_t = interpolate((K_b_filtered[1, :],), K_b_filtered[2, :], Gridded(Linear()))
        itp_l = interpolate((K_b_filtered[1, :],), K_b_filtered[3, :], Gridded(Linear()))

        for idx in missing_idx
            K_b_array[2, idx] = itp_t(angle_range[idx])
            K_b_array[3, idx] = itp_l(angle_range[idx])
        end
    end
    interp_transversal = scale(interpolate(K_b_array[2, :], BSpline(Cubic(Flat(OnGrid())))), angle_range)
    interp_longitudinal = scale(interpolate(K_b_array[2, :], BSpline(Cubic(Flat(OnGrid())))), angle_range)

    return (interp_transversal, interp_longitudinal)
end

"""
Calculate the specific thermal power of the solarthermal collector under the assumption that 
the collector has a fixed delta_T and the flow rate is variable 
"""
function calc_thermal_power_fixed_delta_T(unit::SolarthermalCollectorVal, sim_params::Dict{String,Any}, target_temperature)

    average_temperature = target_temperature - unit.delta_T / 2
    p_spec_th = spec_thermal_power_func(average_temperature, target_temperature, unit.last_average_temperature, 0, unit, sim_params["time_step_seconds"])
    spec_flow_rate = p_spec_th / (unit.delta_T * unit.vol_heat_cap)

    if spec_flow_rate >= unit.spec_flow_rate_min
        return (
            p_spec_th,
            average_temperature,
            target_temperature
        )
    else
        return (
            nothing,
            nothing,
            nothing
        )
    end
end

"""
Calculate the specific thermal power of the solarthermal collector under the assumption that 
the collector has a fixed flow rate and the delta_T is variable 
"""
function calc_thermal_power_fixed_flow_rate(unit::SolarthermalCollectorVal, sim_params::Dict{String,Any}, target_temperature)
    # go through the expected temperature levels and from high to low and check if the 
    # solarthermal collector can deliver energy at this level
    try
        average_temperature = find_zero(
            (t_avg->spec_thermal_power_func(t_avg, target_temperature, unit.last_average_temperature, unit.spec_flow_rate, unit, sim_params["time_step_seconds"]), 
            t_avg->derivate_spec_thermal_power_func(t_avg, unit.spec_flow_rate, unit, sim_params["time_step_seconds"])
            ), 
            target_temperature - unit.delta_T_min, 
            Roots.Newton()
            )
        delta_T = (target_temperature - average_temperature) * 2
        if delta_T >= unit.delta_T_min
            p_spec_th = unit.spec_flow_rate * (target_temperature - average_temperature) * 2 * unit.vol_heat_cap
            return (
                p_spec_th,
                average_temperature,
                target_temperature
            )
        else
            return (
                nothing,
                nothing,
                nothing
            )
        end
    catch
        return (
            nothing,
            nothing,
            nothing
        )
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
    eta_0_hem = 0.702
    b_params = [32.64, 3.59, 0.066]
    
    f = 
    eta_0_hem * (1-b_params[3]) * ((unit.beam_solar_irradiance_in_plane + unit.diffuse_solar_irradiance_in_plane) + 0.85 * (unit.long_wave_irradiance - 5.670374419e-8 * (unit.ambient_temperature + 273.15)^4)) -
    (b_params[1] + b_params[2] * unit.reduced_wind_speed) * (t_avg - unit.ambient_temperature) -
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

    return [string(unit.medium)*" OUT",
            "Temperature", 
            "Max_Energy",
            "Average_Temperature",
            "Ambient_Temperature",
            "Used_Energy",
            "direct_normal_irradiance",
            "beam_solar_irradiance_in_plane",
            "diffuse_solar_irradiance_in_plane",
            "delta_T",
            "spec_flow_rate"
            ]
end

function output_value(unit::SolarthermalCollectorVal, key::OutputKey)::Float64
    if key.value_key == "OUT"
        return calculate_energy_flow(unit.output_interfaces[key.medium])
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
        if unit.output_temperature == unit.average_temperature # TODO evtl. andere Flag setzen
            return 0
        else 
            return wh_to_watts(unit.used_energy) / unit.collector_gross_area / ((unit.output_temperature - unit.average_temperature) * 2 * unit.vol_heat_cap) 
        end
    elseif key.value_key == "spec_flow_rate_profile"
        return unit.spec_flow_rate
    end
    throw(KeyError(key.value_key))
end

export SolarthermalCollectorVal
