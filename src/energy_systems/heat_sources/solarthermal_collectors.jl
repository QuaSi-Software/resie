using Interpolations
using Impute
using Roots

"""
Implementation of a solarthermal collector.
Works for flat plate collectors, vacuum tube collectors and PVT modules.
The calculation is based on EN ISO 9806:2017 for quasi-dynamic models.

Stagnation is ignored under the assumption that the system has either measures to prevent 
stagnation harming the collectors or the designed size is small enough in comparision to 
the total system for stagnation not to become a problem.

## ATTENTION: Geothermal heat collector is currently work in progress and not completed!!
"""

mutable struct SolarthermalCollector <: Component
    uac::String
    controller::Controller
    sys_function::SystemFunction

    input_interfaces::InterfaceMap
    output_interfaces::InterfaceMap

    medium::Symbol

    collector_gross_area::Float64
    tilt_angle::Float32
    azimut_angle::Float32

    eta_0_b::Float32
    K_b_array::Array{Union{Missing, Float32}, 2}
    K_b::Float32
    K_d::Float32
    a_params::Array{Float32, 1}
    sigma::Float64

    ambient_temperature_profile::Union{Profile,Nothing}
    beam_solar_irradiance_profile::Union{Profile,Nothing}
    diffuse_solar_irradiance_profile::Union{Profile,Nothing}
    long_wave_irradiance_profile::Union{Profile,Nothing}
    wind_speed_profile::Union{Profile,Nothing}

    ambient_temperature::Temperature
    beam_solar_irradiance::Float32
    diffuse_solar_irradiance::Float32
    reduced_wind_speed::Float32
    long_wave_irradiance::Float32

    target_temperatures::Array{Float64, 1}
    delta_T::Union{Float32,Nothing}
    spec_vol_flow::Union{Float64,Nothing}

    spec_thermal_power::Float64
    max_energy::Union{Float64,Nothing}
    used_energy::Float64
    output_temperature::Temperature
    average_temperature::Temperature
    last_average_temperature::Temperature

    function SolarthermalCollector(uac::String, config::Dict{String,Any}, sim_params::Dict{String,Any})
        ambient_temperature_profile = 
            "ambient_temperature_profile_file_path" in keys(config) ?
            Profile(config["ambient_temperature_profile_file_path"], sim_params) :
            nothing
        beam_solar_irradiance_profile = 
            "beam_solar_irradiance_profile_file_path" in keys(config) ?
            Profile(config["beam_solar_irradiance_profile_file_path"], sim_params) :
            nothing
        diffuse_solar_irradiance_profile = 
            "diffuse_solar_irradiance_profile_file_path" in keys(config) ?
            Profile(config["diffuse_solar_irradiance_profile_file_path"], sim_params) :
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
            config["tilt_angle"], # tilt or transversal angle
            config["azimut_angle"], # azimut, orientation or longitudinal angle 
        
            config["eta_0_b"],
            vcat([10,20,30,40,50,60,70,80,90]', config["K_b_t_array"]', config["K_b_l_array"]'), 
                                 # row1: angle 0° to 90°/vertical/west, 
                                 # row2: transversal/tilt, 
                                 # row3: longitudinal/azimut/orientation
            1, # Calculate K_b from K_b_array in control()
            config["K_d"], # collector parameter for diffuse irradiance 
            config["a_params"], # collector sim_params a1 to a8 EN ISO 9806:2017
            5.670374419*10^-8, # Stefan Bolzmann Constant
        
            ambient_temperature_profile,
            beam_solar_irradiance_profile,
            diffuse_solar_irradiance_profile,
            long_wave_irradiance_profile,
            wind_speed_profile,
            
            0.0, # ambient temperature from profile
            0.0, # beam_solar_irradiance from profile
            0.0, # diffuse_solar_irradiance from profile
            0.0, # long_wave_irradiance from profile
            0.0, # wind_speed from profile
            
            [0.0], # Target temperature for the output of the thermal collector
            default(config, "delta_T", nothing), # delta_T between input and output temperature
            default(config, "spec_vol_flow", nothing), # nominal specific volume flow of the thermal collector

            0.0, # specific thermal power of the solarthermal collector
            0.0, # maximum available energy
            0.0, # energy that was removed from the system in current time step
            0.0, # output temperature in current time step, calculated in control()
            0.0, # average temperature in current time step, calculated in control()
                     # TODO: initilaize values with ambient_temperature
            0.0, # average temperature from last time step, calculated in control()
        )
    end
end

function initialise!(unit::SolarthermalCollector, sim_params::Dict{String,Any})
    set_storage_transfer!(
        unit.output_interfaces[unit.medium],
        default(
            unit.controller.parameter, "load_storages " * String(unit.medium), true
        )
    )

    unit.K_b = find_K_b!(unit.K_b_array, unit.tilt_angle, unit.azimut_angle)

end

function control(
    unit::SolarthermalCollector,
    components::Grouping,
    sim_params::Dict{String,Any}
    )
    move_state(unit, components, sim_params)
    
    # get values from profiles
    unit.beam_solar_irradiance = Profiles.value_at_time(
        unit.beam_solar_irradiance_profile, sim_params["time"]
        )
    unit.diffuse_solar_irradiance = Profiles.value_at_time(
        unit.diffuse_solar_irradiance_profile, sim_params["time"]
        )
    unit.ambient_temperature = Profiles.value_at_time(
        unit.ambient_temperature_profile, sim_params["time"]
        )
    unit.reduced_wind_speed = min(
        Profiles.value_at_time(unit.wind_speed_profile, sim_params["time"]) - 3, 0
        )
    unit.long_wave_irradiance = Profiles.value_at_time(
        unit.long_wave_irradiance_profile, sim_params["time"]
        )

    if unit.delta_T !== nothing && unit.spec_vol_flow === nothing
        calc_thermal_power_fixed_delta_T!(unit, sim_params)
    elseif unit.delta_T === nothing && unit.spec_vol_flow !== nothing
        calc_thermal_power_fixed_flow!(unit, sim_params)
    else
        error("Error in config file: Either delta_T or spec_vol_flow must have a value") #TODO throw correct error
    end
    
    unit.max_energy = watt_to_wh(max(unit.spec_thermal_power * unit.collector_gross_area, 0))

    if unit.max_energy != 0
        set_max_energy!(unit.output_interfaces[unit.medium], unit.max_energy)
    end

    set_temperature!(unit.output_interfaces[unit.medium], nothing, unit.output_temperature)

    unit.last_average_temperature = unit.average_temperature
end

function process(unit::SolarthermalCollector, sim_params::Dict{String,Any})
    outface = unit.output_interfaces[unit.medium]
    exchanges = balance_on(outface, outface.target)
    blnc = balance(exchanges)

    # # define needed output temperature depending on if a demand needs energy
    # demand_temperatures = []
    # storage_temperatures = []
    # for exchange in exchanges
    #     if exchange.balance < 0.0
    #         push!(demand_temperatures, exchange.temperature)
    #     elseif exchange.storage_potential < 0.0
    #         push!(storage_temperatures, exchange.temperature)
    #     end
    # end

    # if max(demand_temperatures) !== nothing
    #     unit.output_temperature = max(demand_temperatures)
    # elseif max(storage_temperatures) !== nothing
    #     unit.output_temperature = max(storage_temperatures)
    # end

    # unit.last_average_temperature = unit.output_temperature - unit.delta_T / 2

    unit.target_temperatures = []
    energy_demands = Dict()
    for exchange in exchanges
        energy_demand_exchange = exchange.balance + exchange.energy_potential
        if energy_demand_exchange < 0
            push!(unit.target_temperatures, exchange.temperature_min + 1)
            if haskey(energy_demands, exchange.temperature_min)
                energy_demands[round(exchange.temperature_min; digits = 1)] += energy_demand_exchange
            else
                energy_demands[round(exchange.temperature_min; digits = 1)] = energy_demand_exchange
            end
        end
    end
    unique!(unit.target_temperatures)
    sort!(unit.target_temperatures, rev=true)

    temps = [-100.0]
    for temp in collect(keys(energy_demands))
        if temp <= unit.output_temperature
            push!(temps, round(temp; digits = 1))
        end
    end
    if haskey(energy_demands, maximum(temps))
        energy_demand = energy_demands[maximum(temps)]
    else
        energy_demand = 0
    end

    if energy_demand < 0.0
        unit.used_energy = min(abs(energy_demand), unit.max_energy)
        add!(
            outface,
            unit.used_energy,
            unit.output_temperature
        )
        if unit.delta_T === nothing && unit.spec_vol_flow !== nothing && unit.used_energy < unit.max_energy
            unit.average_temperature = unit.output_temperature - wh_to_watt(unit.used_energy) / unit.collector_gross_area / (unit.spec_vol_flow * 1000 * 2 * 4180)
        end
    else
        unit.used_energy = 0.0
    end
end

function find_K_b!(K_b_array, tilt_angle, azimut_angle)
    K_b_array = hcat([0, 1, 1], K_b_array)

    if any(ismissing, K_b_array) 
        K_b_array[:,10] = [90,0,0]
        K_b_array = Impute.interp(K_b_array; dims=1)
    end

    interp_tilt = linear_interpolation(K_b_array[1,:], K_b_array[2,:])
    interp_azimut = linear_interpolation(K_b_array[1,:], K_b_array[3,:])
    K_b = interp_tilt(tilt_angle) * interp_azimut(azimut_angle)

    return K_b
end

function calc_thermal_power_fixed_delta_T!(unit, sim_params)

    spec_vol_flow_min = 0.000002 #TODO: put into constructor & config
    energy_available = false
    for target_temperature in unit.target_temperatures
        t_avg = target_temperature - unit.delta_T / 2
        p_spec_th = spec_thermal_power_func(t_avg, target_temperature, 0, unit, sim_params)
        spec_vol_flow = p_spec_th / (unit.delta_T * 4180 * 1000) 
        if spec_vol_flow > spec_vol_flow_min
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
            (t_avg->spec_thermal_power_func(t_avg, t_avg, 0, unit, sim_params), 
             t_avg->derivate_spec_thermal_power_func(t_avg, 0, unit, sim_params)
            ), 
            unit.last_average_temperature, 
            Roots.Newton()
            )
        unit.output_temperature = unit.average_temperature
    end
end

function calc_thermal_power_fixed_flow!(unit, sim_params)
    delta_T_min = 2 #TODO: put into constructor & config
    energy_available = false
    for target_temperature in unit.target_temperatures
        average_temperature = find_zero(
            (t_avg->spec_thermal_power_func(t_avg, target_temperature, unit.spec_vol_flow, unit, sim_params), 
             t_avg->derivate_spec_thermal_power_func(t_avg, unit.spec_vol_flow, unit, sim_params)
            ), 
            target_temperature - delta_T_min, 
            Roots.Newton()
            )
        delta_T = (target_temperature - average_temperature) * 2
        if delta_T > delta_T_min
            energy_available = true
            unit.spec_thermal_power = unit.spec_vol_flow * 1000 * (target_temperature - average_temperature) * 2 * 4180
            unit.average_temperature = average_temperature
            unit.output_temperature = target_temperature
            break
        end
    end
          
    # if no energy at expected temperature levels is available then calculate 
    # average temperature with a used energy of 0
    if !energy_available
        unit.spec_thermal_power = 0
        unit.average_temperature = find_zero(
            (t_avg->spec_thermal_power_func(t_avg, t_avg, 0, unit, sim_params), 
             t_avg->derivate_spec_thermal_power_func(t_avg, 0, unit, sim_params)
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
function spec_thermal_power_func(t_avg, t_target, spec_vol_flow, unit, sim_params)         
    unit.eta_0_b * unit.K_b * unit.beam_solar_irradiance + 
    unit.eta_0_b * unit.K_d * unit.diffuse_solar_irradiance -
    unit.a_params[1] * (t_avg - unit.ambient_temperature) -
    unit.a_params[2] * (t_avg - unit.ambient_temperature)^2 -
    unit.a_params[3] * unit.reduced_wind_speed * (t_avg - unit.ambient_temperature) +
    unit.a_params[4] * (unit.long_wave_irradiance - unit.sigma * (unit.ambient_temperature + 273.15)^4) -
    unit.a_params[5] * ((t_avg-unit.last_average_temperature) / sim_params["time_step_seconds"]) -
    unit.a_params[6] * unit.reduced_wind_speed * (unit.beam_solar_irradiance + unit.diffuse_solar_irradiance) -
    unit.a_params[7] * unit.reduced_wind_speed * (unit.long_wave_irradiance - unit.sigma * (unit.ambient_temperature + 273.15)^4) -
    unit.a_params[8] * (t_avg - unit.ambient_temperature)^4 -
    spec_vol_flow * 1000 * (t_target - t_avg) * 2 * 4180
end

"""
Derivative of the thermal power output fuction to speed up solve function
"""
function derivate_spec_thermal_power_func(t_avg, spec_vol_flow, unit, sim_params)
    unit.a_params[1] * (-1) -
    unit.a_params[2] * 2 * (t_avg - unit.ambient_temperature) -
    unit.a_params[3] * unit.reduced_wind_speed -
    unit.a_params[5] * 1 / sim_params["time_step_seconds"] -
    unit.a_params[8] * 4 * (t_avg - unit.ambient_temperature)^3 +
    spec_vol_flow * 1000 * 2 * 4180
end

function output_values(unit::SolarthermalCollector)::Vector{String}

    return [string(unit.medium)*" OUT",
            "Temperature", 
            "Max_Energy",
            "Average_Temperature",
            "Ambient_Temperature",
            "Used_Energy",
            "beam_solar_irradiance",
            "delta_T",
            "spec_vol_flow"
            ]
end

function output_value(unit::SolarthermalCollector, key::OutputKey)::Float64
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
    elseif key.value_key == "beam_solar_irradiance"
        return unit.beam_solar_irradiance
    elseif key.value_key == "delta_T"
        return (unit.output_temperature - unit.average_temperature) * 2
    elseif key.value_key == "spec_vol_flow"
        if unit.output_temperature == unit.average_temperature
            return 0
        else
            return wh_to_watt(unit.used_energy) / unit.collector_gross_area / ((unit.output_temperature - unit.average_temperature) * 2 * 4180 * 1000) 
        end
    end
    throw(KeyError(key.value_key))
end

export SolarthermalCollector