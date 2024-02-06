using Interpolations
using Impute

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
    direct_solar_irradiance_profile::Union{Profile,Nothing}
    diffuse_solar_irradiance_profile::Union{Profile,Nothing}
    long_wave_irradiance_profile::Union{Profile,Nothing}
    wind_speed_profile::Union{Profile,Nothing}

    T_target::Float32
    delta_T::Float32

    max_energy::Float64
    output_temperature::Temperature
    input_temperature::Temperature
    average_temperature::Temperature
    last_average_temperature::Temperature
    ambient_temperature::Temperature

    function SolarthermalCollector(uac::String, config::Dict{String,Any})
        ambient_temperature_profile = 
            "ambient_temperature_profile_file_path" in keys(config) ?
            Profile(config["ambient_temperature_profile_file_path"]) :
            nothing
        direct_solar_irradiance_profile = 
            "direct_solar_irradiance_profile_file_path" in keys(config) ?
            Profile(config["direct_solar_irradiance_profile_file_path"]) :
            nothing
        diffuse_solar_irradiance_profile = 
            "diffuse_solar_irradiance_profile_file_path" in keys(config) ?
            Profile(config["diffuse_solar_irradiance_profile_file_path"]) :
            nothing
        long_wave_irradiance_profile = 
            "long_wave_irradiance_profile_file_path" in keys(config) ?
            Profile(config["long_wave_irradiance_profile_file_path"]) :
            nothing
        wind_speed_profile = 
            "wind_speed_profile_file_path" in keys(config) ?
            Profile(config["wind_speed_profile_file_path"]) :
            nothing

        medium = Symbol(default(config, "medium", "m_h_w_ht1"))
        register_media([medium])

        return new(                                                
            uac, # uac
            controller_for_strategy( # controller
                config["strategy"]["name"], config["strategy"]
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
            config["a_params"], # collector parameters a1 to a8 EN ISO 9806:2017
            5.670374419*10^-8, # Stefan Bolzmann Constant
        
            ambient_temperature_profile,
            direct_solar_irradiance_profile,
            diffuse_solar_irradiance_profile,
            long_wave_irradiance_profile,
            wind_speed_profile,
            
            config["T_target"], # Target temperature for the output of the thermal collector
            default(config, "delta_T", 8), # delta_T between input and output temperature

            0.0, # maximum available energy
            nothing, # output temperature in current time step, calculated in control()
            nothing, # input temperature in current time step, calculated in control()
            nothing, # average temperature in current time step, calculated in control()
                     # TODO: initilaize values with ambient_temperature
            0.0, # average temperature from last time step, calculated in control()
            0.0, # ambient temperature from profile
        )
    end
end

function output_values(unit::SolarthermalCollector)::Vector{String}

    return [string(unit.medium)*" OUT",
            "Temperature", 
            "Max_Energy",
            "Average_Temperature",
            "Ambient_Temperature"
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
    end
    throw(KeyError(key.value_key))
end

function control(
    unit::SolarthermalCollector,
    components::Grouping,
    parameters::Dict{String,Any}
)
    move_state(unit, components, parameters)

    unit.K_b = find_K_b!(unit.K_b_array, unit.tilt_angle, unit.azimut_angle)
    
    unit.output_temperature = unit.T_target
    unit.input_temperature = unit.output_temperature - unit.delta_T
    unit.average_temperature = (unit.input_temperature + unit.output_temperature) / 2

    direct_solar_irradiance = Profiles.value_at_time(
        unit.direct_solar_irradiance_profile, parameters["time"]
        )
    diffuse_solar_irradiance = Profiles.value_at_time(
        unit.diffuse_solar_irradiance_profile, parameters["time"]
        )
    unit.ambient_temperature = Profiles.value_at_time(
        unit.ambient_temperature_profile, parameters["time"]
        )
    reduced_wind_speed = Profiles.value_at_time(
        unit.wind_speed_profile, parameters["time"]
        ) - 3
    long_wave_irradiance = Profiles.value_at_time(
        unit.long_wave_irradiance_profile, parameters["time"]
        )

    spec_thermal_power = 
        unit.eta_0_b * unit.K_b * direct_solar_irradiance +
        unit.eta_0_b * unit.K_d * diffuse_solar_irradiance -
        unit.a_params[1] * (unit.average_temperature - unit.ambient_temperature) -
        unit.a_params[2] * (unit.average_temperature - unit.ambient_temperature)^2 -
        unit.a_params[3] * reduced_wind_speed * (unit.average_temperature - unit.ambient_temperature) -
        unit.a_params[4] * (long_wave_irradiance - unit.sigma * (unit.ambient_temperature + 273.15)^4) -
        unit.a_params[5] * ((unit.average_temperature-unit.last_average_temperature) / parameters["time_step_seconds"]) -
        unit.a_params[6] * reduced_wind_speed * (direct_solar_irradiance + diffuse_solar_irradiance) -
        unit.a_params[7] * reduced_wind_speed * (long_wave_irradiance - unit.sigma * (unit.ambient_temperature + 273.15)^4) -
        unit.a_params[8] * (unit.average_temperature - unit.ambient_temperature)^4
    
    spec_mass_flow = spec_thermal_power / (unit.delta_T * 4.18) #TODO Add C_p calculation
    unit.max_energy = watt_to_wh(spec_thermal_power * unit.collector_gross_area)
    set_max_energy!(unit.output_interfaces[unit.medium], unit.max_energy)

    unit.output_interfaces[unit.medium].temperature = unit.output_temperature

    unit.last_average_temperature = unit.average_temperature
end

#TODO: define function process 

function process(unit::SolarthermalCollector, parameters::Dict{String,Any})
    outface = unit.output_interfaces[unit.medium]
    exchange = balance_on(outface, outface.target)
    if exchange.balance < 0.0
        add!(
            outface,
            min(abs(exchange.balance), unit.max_energy),
            unit.output_temperature
        )
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

export SolarthermalCollector