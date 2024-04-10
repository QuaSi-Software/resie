using Interpolations
using Impute
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

    target_temperatures::Array{Float64, 1}
    delta_T::Union{Float32,Nothing}
    spec_flow_rate::Union{Float64,Nothing}

    spec_thermal_power::Float64
    max_energy::Union{Float64,Nothing}
    used_energy::Float64
    output_temperature::Temperature
    average_temperature::Temperature
    last_average_temperature::Temperature

    c_p::Float64

    #TODO for debugging
    target_temperatures_profile::Union{Profile,Nothing}
    direct_normal_irradiance

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
            0.0, # long_wave_irradiance from profile
            0.0, # reduced wind_speed from profile multiplied with wind speed reduction
            
            [0.0], # Target temperature for the output of the thermal collector
            default(config, "delta_T", nothing), # delta_T between input and output temperature
            default(config, "spec_flow_rate", nothing), # nominal specific volume flow of the thermal collector

            0.0, # specific thermal power of the solarthermal collector
            0.0, # maximum available energy
            0.0, # energy that was removed from the system in current time step
            0.0, # output temperature in current time step, calculated in control()
            0.0, # average temperature in current time step, calculated in control()
                     # TODO: initilaize values with ambient_temperature
            0.0, # average temperature from last time step, calculated in control()

            default(config, "c_p", 4200), # specific thermal capacity of the fluid in the collector

            #TODO for debugging access to target_temperatures:
            Profile("profiles/validation/demand_temperature_st_validation.prf", sim_params),
            0.0
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
end

function control(
    unit::SolarthermalCollectorVal,
    components::Grouping,
    sim_params::Dict{String,Any}
    )
    move_state(unit, components, sim_params)
    
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

    unit.reduced_wind_speed = min(
        Profiles.value_at_time(unit.wind_speed_profile, sim_params["time"]) * unit.wind_speed_reduction, 0
        ) #TODO add -3 to wind speed after debugging

    unit.long_wave_irradiance = Profiles.value_at_time(
        unit.long_wave_irradiance_profile, sim_params["time"]
        )

    if sim_params["time"] isa DateTime 
        time = sim_params["time"]
    else
        start_date = DateTime(2005, 1, 1, 0, 30, 0) - Dates.Hour(1)
        # time = start_date + Dates.Second(sim_params["time"])
        time = start_date + Dates.Second(floor(sim_params["time"]/3600)*3600) #TODO for debugging
    end

    solar_zenith, solar_azimuth = sun_position(time, 9.18, 47.67, 1.0, unit.ambient_temperature)

    unit.beam_solar_irradiance_in_plane, unit.direct_normal_irradiance, angle_of_incidence, longitudinal_angle, transversal_angle = beam_irr_in_plane(
        unit.tilt_angle, unit.azimuth_angle, solar_zenith, solar_azimuth, 
        global_solar_hor_irradiance, diffuse_solar_hor_irradiance, global_solar_hor_irradiance
        )

    unit.K_b = calc_K_b!(unit.K_b_array, transversal_angle, longitudinal_angle)

    # TODO debugging access to target_temperatures
    unit.target_temperatures = [Profiles.value_at_time(unit.target_temperatures_profile, sim_params["time"])]

    # TODO debugging access to direct irradiance in plane
    unit.beam_solar_irradiance_in_plane = global_solar_hor_irradiance

    # calcualte specific thermal power depending on control mode
    if unit.delta_T !== nothing && unit.spec_flow_rate === nothing
        calc_thermal_power_fixed_delta_T!(unit, sim_params)
    elseif unit.delta_T === nothing && unit.spec_flow_rate !== nothing
        # calc_thermal_power_fixed_flow!(unit, sim_params)
        calc_thermal_power_fixed_flow_input_temp!(unit, sim_params) # TODO input_temp calculation for debugging
    else
        error("Error in config file: Either delta_T or spec_flow_rate must have a value") #TODO throw correct error
    end
    
    unit.max_energy = watt_to_wh(max(unit.spec_thermal_power * unit.collector_gross_area, 0))

    if unit.max_energy != 0
        set_max_energy!(unit.output_interfaces[unit.medium], unit.max_energy)
    end

    set_temperature!(unit.output_interfaces[unit.medium], nothing, unit.output_temperature)

    unit.last_average_temperature = unit.average_temperature
end

function process(unit::SolarthermalCollectorVal, sim_params::Dict{String,Any})
    outface = unit.output_interfaces[unit.medium]
    exchanges = balance_on(outface, outface.target)
    blnc = balance(exchanges)

    # # define needed output temperature depending on if a demand needs energy
    # demand_temperatures = []
    # storage_temperatures = []
    # for exchange in exchanges
    #     if exchange.balance < 0.0
    #         push!(demand_temperatures, exchange.temperature)
    #     elseif exchange.energy_potential < 0.0
    #         push!(storage_temperatures, exchange.temperature)
    #     end
    # end

    # if max(demand_temperatures) !== nothing
    #     unit.output_temperature = max(demand_temperatures)
    # elseif max(storage_temperatures) !== nothing
    #     unit.output_temperature = max(storage_temperatures)
    # end

    # loop through demands & storages that are connected to the solarthermal collector and
    # write the demands in a dictonary with the corresponding temperature as the key
    # additionaly write target_temperatures to a list for the next control() step
    unit.target_temperatures = []
    energy_demands = Dict()
    for exchange in exchanges
        energy_demand_exchange = exchange.balance + exchange.energy_potential
        if energy_demand_exchange < 0
            push!(unit.target_temperatures, exchange.temperature_min + 1)
            if haskey(energy_demands, exchange.temperature_min)
                energy_demands[floor(exchange.temperature_min; digits = 1)] += energy_demand_exchange
            else
                energy_demands[floor(exchange.temperature_min; digits = 1)] = energy_demand_exchange
            end
        end
    end
    unique!(unit.target_temperatures)
    sort!(unit.target_temperatures, rev=true)

    # check which temperature demands can be supplied by the solarthermal collector 
    # output_temperature and sum up all energy demands that are lower than that
    temps = []
    for temp in collect(keys(energy_demands))
        if temp <= unit.output_temperature
            push!(temps, round(temp; digits = 1))
        end
    end
    sort!(temps, rev=true)

    energy_demand = 0
    if length(temps) > 0
        for temp in temps
            energy_demand += energy_demands[temp]
        end
    end

    if energy_demand < 0.0
        unit.used_energy = min(abs(energy_demand), unit.max_energy)
        add!(
            outface,
            unit.used_energy,
            unit.output_temperature
        )

        # recalculate the average_temperature if not all energy is used to keep a constant 
        # flow rate
        # TODO for debugging diabled
        # if unit.delta_T === nothing && unit.spec_flow_rate !== nothing && unit.used_energy < unit.max_energy
        #     unit.average_temperature = unit.output_temperature - wh_to_watt(unit.used_energy) / unit.collector_gross_area / (unit.spec_flow_rate * 1000 * 2 * unit.c_p)
        # end
    else
        unit.used_energy = 0.0
        set_max_energy!(unit.output_interfaces[unit.medium], 0.0)
    end
end

"""
Calculate K_b based on a table with provided values for longitudinal and transversal angles 
and the angle of irradiance on the plane of the solarthermal collector.
Table with provided values is mirrored for negative angles.
"""
function calc_K_b!(K_b_array, transversal_angle, longitudinal_angle)
    # TODO prepare interpolations only once at initialise()

    # K_b_array_mirr = reverse(K_b_array, dims=2)
    # K_b_array_mirr[1, :] = K_b_array_mirr[1, :] .* -1
    # K_b_array = hcat(K_b_array_mirr, [0, 1, 1], K_b_array)

    transversal_angle = abs(transversal_angle)
    longitudinal_angle = abs(longitudinal_angle)

    K_b_array = hcat([0, 1, 1], K_b_array)
    angle_range = K_b_array[1, 1]:10:K_b_array[1, end]

    if any(ismissing, K_b_array) 
        K_b_array[:,19] = [90,0,0]
        # K_b_array[:,1] = [-90,0,0]
        K_b_array = Impute.interp(K_b_array; dims=1)
    end

    interp_transversal = cubic_spline_interpolation(angle_range, K_b_array[2,:])
    interp_longitudinal = cubic_spline_interpolation(angle_range, K_b_array[3,:])
    if transversal_angle < 90 && longitudinal_angle < 90
        K_b = interp_transversal(transversal_angle) * interp_longitudinal(longitudinal_angle)
    else
        K_b = 0
    end
    return K_b
end

"""
Calculate the specific thermal power of the solarthermal collector under the assumption that 
the collector has a fixed delta_T and the flow rate is variable 
"""
function calc_thermal_power_fixed_delta_T!(unit, sim_params)

    spec_flow_rate_min = 0.000002 #TODO: put into constructor & config
    energy_available = false

    # go through the needed output temperatures and from high to low and check if the 
    # solarthermal collector can deliver energy at this level
    for target_temperature in unit.target_temperatures
        t_avg = target_temperature - unit.delta_T / 2
        p_spec_th = spec_thermal_power_func(t_avg, target_temperature, unit.last_average_temperature, 0, unit, sim_params["time_step_seconds"])
        spec_flow_rate = p_spec_th / (unit.delta_T * unit.c_p * 1000) 
        if spec_flow_rate > spec_flow_rate_min
            energy_available = true
            unit.spec_thermal_power = p_spec_th
            unit.average_temperature = t_avg
            unit.output_temperature = target_temperature
            break
        end
    end

    # if no energy at expected temperature levels is available then calculate 
    # average temperature with a used energy of 0
    if !energy_available
        unit.spec_thermal_power = 0
        unit.average_temperature = find_zero(
            (t_avg->spec_thermal_power_func(t_avg, t_avg, unit.last_average_temperature, 0, unit, sim_params["time_step_seconds"]), 
                t_avg->derivate_spec_thermal_power_func(t_avg, 0, unit, sim_params["time_step_seconds"])
            ), 
            unit.last_average_temperature, 
            Roots.Newton()
            )

        unit.output_temperature = unit.average_temperature
    end
end

"""
Calculate the specific thermal power of the solarthermal collector under the assumption that 
the collector has a fixed flow rate and the delta_T is variable 
"""
function calc_thermal_power_fixed_flow!(unit, sim_params)
    delta_T_min = 0 #TODO: put into constructor & config
    energy_available = false

    # go through the expected temperature levels and from high to low and check if the 
    # solarthermal collector can deliver energy at this level
    for target_temperature in unit.target_temperatures
        try
            average_temperature = find_zero(
                (t_avg->spec_thermal_power_func(t_avg, target_temperature, unit.last_average_temperature, unit.spec_flow_rate, unit, sim_params["time_step_seconds"]), 
                t_avg->derivate_spec_thermal_power_func(t_avg, unit.spec_flow_rate, unit, sim_params["time_step_seconds"])
                ), 
                target_temperature - delta_T_min, 
                Roots.Newton()
                )
            delta_T = (target_temperature - average_temperature) * 2
            if delta_T > delta_T_min
                energy_available = true
                unit.spec_thermal_power = unit.spec_flow_rate * 1000 * (target_temperature - average_temperature) * 2 * unit.c_p
                unit.average_temperature = average_temperature
                unit.output_temperature = target_temperature
                break
            end
        catch
            continue
        end
    end
          
    # if no energy at expected temperature levels is available then calculate 
    # average temperature with a used energy of 0
    if !energy_available
        unit.spec_thermal_power = 0
        unit.average_temperature = find_zero(
            (t_avg->spec_thermal_power_func(t_avg, t_avg, unit.last_average_temperature, 0, unit, sim_params["time_step_seconds"]), 
             t_avg->derivate_spec_thermal_power_func(t_avg, 0, unit, sim_params["time_step_seconds"])
            ), 
            unit.last_average_temperature, 
            Roots.Newton()
            )
        unit.output_temperature = unit.average_temperature
    end

end

"""
Calculate the specific thermal power of the solarthermal collector under the assumption that 
the collector has a fixed flow rate and the delta_T is variable 
"""
function calc_thermal_power_fixed_flow_input_temp!(unit, sim_params)
    delta_T_min = 0 #TODO: put into constructor & config
    energy_available = false
    input_temp = 10

    # go through the expected temperature levels and from high to low and check if the 
    # solarthermal collector can deliver energy at this level

    try
        average_temperature = find_zero(
            (t_avg->spec_thermal_power_func_input_temp(t_avg, input_temp, unit.last_average_temperature, unit.spec_flow_rate, unit, sim_params["time_step_seconds"]), 
            t_avg->derivate_spec_thermal_power_func_input_temp(t_avg, unit.spec_flow_rate, unit, sim_params["time_step_seconds"])
            ), 
            input_temp + delta_T_min, 
            Roots.Newton()
            )
        delta_T = (average_temperature - input_temp) * 2
        # if delta_T > delta_T_min
        if true # TODO for debugging
            energy_available = true
            unit.spec_thermal_power = unit.spec_flow_rate * 1000 * delta_T * unit.c_p
            unit.average_temperature = average_temperature
            unit.output_temperature = input_temp + delta_T
        end
    catch 
    end

          
    # if no energy at expected temperature levels is available then calculate 
    # average temperature with a used energy of 0
    if !energy_available
        unit.spec_thermal_power = 0
        unit.average_temperature = find_zero(
            (t_avg->spec_thermal_power_func(t_avg, t_avg, unit.last_average_temperature, 0, unit, sim_params["time_step_seconds"]), 
             t_avg->derivate_spec_thermal_power_func(t_avg, 0, unit, sim_params["time_step_seconds"])
            ), 
            unit.last_average_temperature, 
            Roots.Newton()
            )
        unit.output_temperature = unit.average_temperature
    end

end

"""
Calculate the specific thermal power of the solarthermal collector under the assumption that 
the collector has a fixed flow rate and the delta_T is variable.
Simulates smaller timesteps to decrease influence of thermal capacity.
"""
function calc_thermal_power_fixed_flow_sub!(unit, sim_params)
    delta_T_min = 0 #TODO: put into constructor & config
    energy_available = false
    t_in_min = 0

    # go through the expected temperature levels and from high to low and check if the 
    # solarthermal collector can deliver energy at this level
    for target_temperature in unit.target_temperatures
        t_avg_last = unit.last_average_temperature
        sub_time_step = 900
        average_temperature = 0.0
        try
            for time_step in sub_time_step:sub_time_step:sim_params["time_step_seconds"]
                average_temperature = find_zero(
                    (t_avg->spec_thermal_power_func(t_avg, target_temperature, t_avg, unit.spec_flow_rate, unit, sub_time_step), 
                    t_avg->derivate_spec_thermal_power_func_static(t_avg, unit.spec_flow_rate, unit, sub_time_step)
                    ), 
                    target_temperature - delta_T_min, 
                    Roots.Newton()
                    )

                delta_T = (target_temperature - average_temperature) * 2
                if delta_T > delta_T_min && (target_temperature - delta_T) > t_in_min
                    t_avg_last = average_temperature
                else
                    average_temperature = find_zero(
                        (t_avg->spec_thermal_power_func(t_avg, t_avg, t_avg, 0, unit, sub_time_step), 
                         t_avg->derivate_spec_thermal_power_func_static(t_avg, 0, unit, sub_time_step)
                        ), 
                        t_avg_last, 
                        Roots.Newton()
                        )
                    t_avg_last = average_temperature
                end
            end
            delta_T = (target_temperature - average_temperature) * 2
            if delta_T > delta_T_min && (target_temperature - delta_T) > t_in_min
                energy_available = true
                unit.spec_thermal_power = unit.spec_flow_rate * 1000 * (target_temperature - average_temperature) * 2 * unit.c_p
                unit.average_temperature = average_temperature
                unit.output_temperature = target_temperature
                break
            end
        catch
            continue
        end
    end
          
    # if no energy at expected temperature levels is available then calculate 
    # average temperature with a used energy of 0
    if !energy_available
        unit.spec_thermal_power = 0
        unit.average_temperature = find_zero(
            (t_avg->spec_thermal_power_func(t_avg, t_avg, unit.last_average_temperature, 0, unit, sim_params["time_step_seconds"]), 
             t_avg->derivate_spec_thermal_power_func(t_avg, 0, unit, sim_params["time_step_seconds"])
            ), 
            unit.last_average_temperature, 
            Roots.Newton()
            )
        unit.output_temperature = unit.average_temperature
    end

end

"""
Function for calculating the thermal power output of a solarthermal collector.
Can be used to solve after t_avg and find average_temperature for specific used energy
"""
function spec_thermal_power_func(t_avg, t_target, t_avg_last, spec_flow_rate, unit, time_step)         
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
    spec_flow_rate * 1000 * (t_target - t_avg) * 2 * unit.c_p
end

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
    spec_flow_rate * 1000 * (t_avg - t_input) * 2 * unit.c_p
end

"""
Derivative of the thermal power output fuction to speed up solve function
"""
function derivate_spec_thermal_power_func(t_avg, spec_flow_rate, unit, time_step)
    unit.a_params[1] * (-1) -
    unit.a_params[2] * 2 * (t_avg - unit.ambient_temperature) -
    unit.a_params[3] * unit.reduced_wind_speed -
    unit.a_params[5] * 1 / time_step -
    unit.a_params[8] * 4 * (t_avg - unit.ambient_temperature)^3 +
    spec_flow_rate * 1000 * 2 * unit.c_p
end

"""
Derivative of the thermal power output fuction to speed up solve function
"""
function derivate_spec_thermal_power_func_input_temp(t_avg, spec_flow_rate, unit, time_step)
    unit.a_params[1] * (-1) -
    unit.a_params[2] * 2 * (t_avg - unit.ambient_temperature) -
    unit.a_params[3] * unit.reduced_wind_speed -
    unit.a_params[5] * 1 / time_step -
    unit.a_params[8] * 4 * (t_avg - unit.ambient_temperature)^3 -
    spec_flow_rate * 1000 * 2 * unit.c_p
end

"""
Derivative of the thermal power output fuction to speed up solve function
Derivative of the thermal power output fuction to speed up solve function
"""
function derivate_spec_thermal_power_func_static(t_avg, spec_flow_rate, unit, time_step)
    unit.a_params[1] * (-1) -
    unit.a_params[2] * 2 * (t_avg - unit.ambient_temperature) -
    unit.a_params[3] * unit.reduced_wind_speed -
    unit.a_params[8] * 4 * (t_avg - unit.ambient_temperature)^3 +
    spec_flow_rate * 1000 * 2 * unit.c_p
end

function output_values(unit::SolarthermalCollectorVal)::Vector{String}

    return [string(unit.medium)*" OUT",
            "Temperature", 
            "Max_Energy",
            "Average_Temperature",
            "Ambient_Temperature",
            "Used_Energy",
            "beam_solar_irradiance_in_plane",
            "diffuse_solar_irradiance_in_plane",
            "delta_T",
            "spec_flow_rate",
            "direct_normal_irradiance" # TODO for debugging
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
            return wh_to_watt(unit.used_energy) / unit.collector_gross_area / ((unit.output_temperature - unit.average_temperature) * 2 * unit.c_p * 1000) 
        end
    elseif key.value_key == "direct_normal_irradiance" # TODO for debugging
        return unit.direct_normal_irradiance
    end
    throw(KeyError(key.value_key))
end

export SolarthermalCollector
