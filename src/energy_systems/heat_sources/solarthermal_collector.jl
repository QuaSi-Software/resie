import Interpolations as ip
using Roots
using Resie.SolarIrradiance

"""
Implementation of a solarthermal collector.
Works for flat plate collectors, vacuum tube collectors and PVT modules.
The calculation is based on EN ISO 9806:2017 for quasi-dynamic models.
"""

#! format: off
const SOLARTHERMAL_COLLECTOR_PARAMETERS = Dict(
    "m_heat_out" => (
        default="m_h_w_ht1",
        description="Heat output medium",
        display_name="Medium heat_out",
        required=false,
        type=String,
        json_type="string",
        unit="-"
    ),
    "ambient_temperature_profile_file_path" => (
        default=nothing,
        description="Path to ambient temperature profile file",
        display_name="Ambient temp. profile",
        required=false,
        conditionals=[
            ("ambient_temperature_from_global_file", "mutex"),
            ("constant_ambient_temperature", "mutex")
        ],
        type=String,
        json_type="string",
        unit="-"
    ),
    "ambient_temperature_from_global_file" => (
        default=nothing,
        description="If given points to a key in the global weather data file with the ambient temperature profile",
        display_name="Global file amb. temp. key",
        required=false,
        conditionals=[
            ("ambient_temperature_profile_file_path", "mutex"),
            ("constant_ambient_temperature", "mutex")
        ],
        type=String,
        json_type="string",
        unit="-"
    ),
    "constant_ambient_temperature" => (
        default=nothing,
        description="Constant ambient temperature value",
        display_name="Constant ambient temp.",
        required=false,
        conditionals=[
            ("ambient_temperature_profile_file_path", "mutex"),
            ("ambient_temperature_from_global_file", "mutex")
        ],
        type=Temperature,
        json_type="number",
        unit="°C"
    ),
    "beam_solar_radiation_profile_file_path" => (
        default=nothing,
        description="Path to beam solar radiation profile file",
        display_name="Beam solar rad. profile",
        required=false,
        conditionals=[
            ("beam_solar_radiation_from_global_file", "mutex"),
            ("constant_beam_solar_radiation", "mutex")
        ],
        type=String,
        json_type="string",
        unit="-"
    ),
    "beam_solar_radiation_from_global_file" => (
        default=nothing,
        description="If given points to a key in the global weather data file with beam solar radiation",
        display_name="Global file beam rad. key",
        required=false,
        conditionals=[
            ("beam_solar_radiation_profile_file_path", "mutex"),
            ("constant_beam_solar_radiation", "mutex")
        ],
        type=String,
        json_type="string",
        unit="-"
    ),
    "constant_beam_solar_radiation" => (
        default=nothing,
        description="Constant beam solar radiation value",
        display_name="Constant beam solar rad.",
        required=false,
        conditionals=[
            ("beam_solar_radiation_profile_file_path", "mutex"),
            ("beam_solar_radiation_from_global_file", "mutex")
        ],
        type=Float64,
        json_type="number",
        unit="W/m^2"
    ),
    "diffuse_solar_radiation_profile_file_path" => (
        default=nothing,
        description="Path to diffuse solar radiation profile file",
        display_name="Diffuse solar rad. profile",
        required=false,
        conditionals=[
            ("diffuse_solar_radiation_from_global_file", "mutex"),
            ("constant_diffuse_solar_radiation", "mutex")
        ],
        type=String,
        json_type="string",
        unit="-"
    ),
    "diffuse_solar_radiation_from_global_file" => (
        default=nothing,
        description="If given points to a key in the global weather data file with diffuse solar radiation",
        display_name="Global file diffuse rad. key",
        required=false,
        conditionals=[
            ("diffuse_solar_radiation_profile_file_path", "mutex"),
            ("constant_diffuse_solar_radiation", "mutex")
        ],
        type=String,
        json_type="string",
        unit="-"
    ),
    "constant_diffuse_solar_radiation" => (
        default=nothing,
        description="Constant diffuse solar radiation value",
        display_name="Constant diffuse solar rad.",
        required=false,
        conditionals=[
            ("diffuse_solar_radiation_profile_file_path", "mutex"),
            ("diffuse_solar_radiation_from_global_file", "mutex")
        ],
        type=Float64,
        json_type="number",
        unit="W/m^2"
    ),
    "infrared_sky_radiation_profile_file_path" => (
        default=nothing,
        description="Path to long wave infrared sky radiation profile file",
        display_name="Infrared sky rad. profile",
        required=false,
        conditionals=[
            ("infrared_sky_radiation_from_global_file", "mutex"),
            ("constant_infrared_sky_radiation", "mutex")
        ],
        type=String,
        json_type="string",
        unit="-"
    ),
    "infrared_sky_radiation_from_global_file" => (
        default=nothing,
        description="If given points to a key in the global weather data file with infrared sky radiation",
        display_name="Global file infrared rad. key",
        required=false,
        conditionals=[
            ("infrared_sky_radiation_profile_file_path", "mutex"),
            ("constant_infrared_sky_radiation", "mutex")
        ],
        type=String,
        json_type="string",
        unit="-"
    ),
    "constant_infrared_sky_radiation" => (
        default=nothing,
        description="Constant long wave infrared sky radiation value",
        display_name="Constant infrared sky rad.",
        required=false,
        conditionals=[
            ("infrared_sky_radiation_profile_file_path", "mutex"),
            ("infrared_sky_radiation_from_global_file", "mutex")
        ],
        type=Float64,
        json_type="number",
        unit="W/m²"
    ),
    "wind_speed_profile_file_path" => (
        default=nothing,
        description="Path to wind speed profile file",
        display_name="Wind speed profile",
        required=false,
        conditionals=[
            ("wind_speed_from_global_file", "mutex"),
            ("constant_wind_speed", "mutex")
        ],
        type=String,
        json_type="string",
        unit="-"
    ),
    "wind_speed_from_global_file" => (
        default=nothing,
        description="If given points to a key in the global weather data file with wind speed",
        display_name="Global file wind speed key",
        required=false,
        conditionals=[
            ("wind_speed_profile_file_path", "mutex"),
            ("constant_wind_speed", "mutex")
        ],
        type=String,
        json_type="string",
        unit="-"
    ),
    "constant_wind_speed" => (
        default=nothing,
        description="Constant wind speed value",
        display_name="Constant wind speed",
        required=false,
        conditionals=[
            ("wind_speed_profile_file_path", "mutex"),
            ("wind_speed_from_global_file", "mutex")
        ],
        type=Float64,
        json_type="number",
        unit="m/s"
    ),
    "collector_gross_area" => (
        default=nothing,
        description="Gross area of the solarthermal collector",
        display_name="Collector gross area",
        required=true,
        type=Float64,
        json_type="number",
        unit="m^2"
    ),
    "tilt_angle" => (
        default=nothing,
        description="Tilt angle of the collector between 0° (horizontal) and 90° (vertical)",
        display_name="Tilt angle",
        required=true,
        type=Float64,
        json_type="number",
        unit="°"
    ),
    "azimuth_angle" => (
        default=nothing,
        description="Azimuth angle (orientation) of the collector between -180° and 180° " *
                    "with 0°=south, -90°=east, 90°=west",
        display_name="Azimuth angle",
        required=true,
        type=Float64,
        json_type="number",
        unit="°"
    ),
    "ground_reflectance" => (
        default=0.2,
        description="Reflectance (albedo) of the ground around the collector",
        display_name="Ground reflectance",
        required=false,
        type=Float64,
        json_type="number",
        unit="-"
    ),
    "eta_0_b" => (
        default=nothing,
        description="Zero-loss efficiency eta_0_b at (mean_temp - ambient_temp)=0 K based " *
                    "on beam solar irradiance. From Solar Keymark database.",
        display_name="Zero-loss eff. eta_0_b",
        required=true,
        type=Float64,
        json_type="number",
        unit="-"
    ),
    "K_b_t_array" => (
        default=nothing,
        description="Array with transversal incidence angle modifier values from 10° to " *
                    "90°. Use null for interpolation.",
        display_name="K_b transversal array",
        required=true,
        type=Vector{Float64},
        json_type="array",
        unit="-"
    ),
    "K_b_l_array" => (
        default=nothing,
        description="Array with longitudinal incidence angle modifier values from 10° to " *
                    "90°. Use null for interpolation.",
        display_name="K_b longitudinal array",
        required=true,
        type=Vector{Float64},
        json_type="array",
        unit="-"
    ),
    "K_d" => (
        default=nothing,
        description="Incidence angle modifier for diffuse irradiance. From Solar Keymark database.",
        display_name="Diffuse incidence angle modifier K_d",
        required=true,
        type=Float64,
        json_type="number",
        unit="-"
    ),
    "a_params" => (
        default=nothing,
        description="Parameters a1-a8 defining the collector according to EN ISO " *
                    "9806:2017. From Solar Keymark database.",
        display_name="ISO 9806:2017 a-parameters",
        required=true,
        type=Vector{Float64},
        json_type="array",
        unit="-"
    ),
    "vol_heat_capacity" => (
        default=4.2e6,
        description="Volumetric heat capacity of the heat transfer fluid in the collector",
        display_name="Volumetric heat capacity",
        required=false,
        type=Float64,
        json_type="number",
        unit="J/m^3*K"
    ),
    "wind_speed_reduction" => (
        default=1.0,
        description="Factor to adjust wind speed to account for different wind conditions vs. 10 m height reference",
        display_name="Wind speed reduction factor",
        required=false,
        type=Float64,
        json_type="number",
        unit="-"
    ),
    "delta_T" => (
        default=nothing,
        description="Constant temperature difference between collector input and output. " *
                    "Either delta_T or spec_flow_rate must be defined.",
        display_name="Temperature difference",
        required=false,
        conditionals=[("spec_flow_rate", "mutex")],
        type=Float64,
        json_type="number",
        unit="K"
    ),
    "spec_flow_rate" => (
        default=nothing,
        description="Constant specific flow rate per m² collector area. Either " *
                    "spec_flow_rate or delta_T must be defined.",
        display_name="Specific flow rate",
        required=false,
        conditionals=[("delta_T", "mutex")],
        type=Float64,
        json_type="number",
        unit="m^3/s*m^2"
    ),
    "spec_flow_rate_min" => (
        default=0.000002,
        description="Minimal spec_flow_rate to start producing energy; used together with delta_T",
        display_name="Min. specific flow rate",
        required=false,
        type=Float64,
        json_type="number",
        unit="m^3/s*m^2"
    ),
    "delta_T_min" => (
        default=2.0,
        description="Minimal delta_T to start producing energy; used together with spec_flow_rate",
        display_name="Min. temperature difference",
        required=false,
        type=Float64,
        json_type="number",
        unit="K"
    )
)
#! format: on

mutable struct SolarthermalCollector <: Component
    ## general
    uac::String
    controller::Controller
    sys_function::SystemFunction
    input_interfaces::InterfaceMap
    output_interfaces::InterfaceMap
    m_heat_out::Symbol

    ## collector installation
    collector_gross_area::Float64
    tilt_angle::Float64
    azimuth_angle::Float64
    ground_reflectance::Float64

    ## collector parameters
    eta_0_b::Float64
    # row1: angle 0° to 90°/vertical/west
    # row2: transversal/tilt
    # row3: longitudinal/azimuth/orientation
    K_b_array::Array{Union{Missing,Float64},2}
    # interpolation functions
    K_b_itp::Tuple{AbstractArray,AbstractArray}
    # calculated from K_b_array in control()
    K_b::Float64
    K_d::Float64
    a_params::Array{Float64,1}
    # equals Boltzmann-constant
    sigma::Float64
    vol_heat_cap::Float64

    ## weather profiles
    ambient_temperature_profile::Union{Profile,Nothing}
    beam_solar_hor_irradiance_profile::Union{Profile,Nothing}
    diffuse_solar_hor_irradiance_profile::Union{Profile,Nothing}
    long_wave_irradiance_profile::Union{Profile,Nothing}
    wind_speed_profile::Union{Profile,Nothing}

    ## weather parameters
    wind_speed_reduction::Float64
    ambient_temperature::Temperature
    # diffuse_solar_irradiance_in_plane in [W/m²]
    beam_solar_irradiance_in_plane::Float64
    # diffuse_solar_irradiance_in_plane [W/m²]
    diffuse_solar_irradiance_in_plane::Float64
    reduced_wind_speed::Floathing
    beam_solar_hor_irradiance::Floathing
    diffuse_solar_hor_irradiance::Floathing
    long_wave_irradiance::Floathing
    # direct_normal_irradiance calculated from sun position [W/m²]
    direct_normal_irradiance::Float64

    ## operation parameters
    delta_T::Union{Float64,Nothing}
    spec_flow_rate::Floathing
    delta_T_min::Union{Float64,Nothing}
    spec_flow_rate_min::Floathing

    ## results
    # energy that was removed from the system in current time step
    used_energy::Array{Floathing}
    # output temperature in current time step, calculated in control()
    output_temperature::Temperature
    # average temperature in current time step, calculated in control()
    average_temperature::Temperature
    # average temperature from last time step, calculated in control()
    last_average_temperature::Temperature
    # actual flow rate resulting from simulation and used energy
    spec_flow_rate_actual::Floathing
    # actual delta_T resulting from simulation and written to output
    delta_T_actual::Float64

    ## internal parameters for temperature layers
    # list of available_energies in for each component connected to the collector
    available_energies::Array{Floathing}
    # list of average_temperatures in for each component connected to the collector
    average_temperatures::Array{Floathing}
    # list of output_temperatures in for each component connected to the collector
    output_temperatures::Array{Floathing}
    # list of runtimes at different temperature levels for each component connected to the collector
    runtimes::Array{Floathing}
    # list of all calculated output temperature and output energy pairs
    temperature_energy_pairs::Dict{Temperature,Tuple}
    # uacs for which energy was calculated
    calc_uacs::Array{Stringing}
    mean_ambient_temperature::Float64

    function SolarthermalCollector(uac::String, config::Dict{String,Any}, sim_params::Dict{String,Any})
        new(SSOT_parameter_constructor(SolarthermalCollector, uac, config, sim_params)...)
    end
end

function component_parameters(x::Type{SolarthermalCollector})::Dict{String,NamedTuple}
    return deepcopy(SOLARTHERMAL_COLLECTOR_PARAMETERS)
end

function extract_parameter(x::Type{SolarthermalCollector}, config::Dict{String,Any}, param_name::String,
                           param_def::NamedTuple, sim_params::Dict{String,Any}, uac::String)
    if param_name in ("ambient_temperature_from_global_file", "beam_solar_radiation_from_global_file",
                      "infrared_sky_radiation_from_global_file", "diffuse_solar_radiation_from_global_file",
                      "wind_speed_from_global_file")
        return load_profile_from_global_weather_file(config, param_name, sim_params, uac)
    elseif param_name in ("ambient_temperature_profile_file_path", "beam_solar_radiation_profile_file_path",
                          "infrared_sky_radiation_profile_file_path", "diffuse_solar_radiation_profile_file_path",
                          "wind_speed_profile_file_path")
        return load_optional_profile(config, param_name, sim_params)
    elseif param_name in ("constant_ambient_temperature", "constant_beam_solar_radiation",
                          "constant_infrared_sky_radiation", "constant_diffuse_solar_radiation", "constant_wind_speed")
        if occursin("temperature", param_name)
            return convert(Temperature, default(config, param_name, nothing))
        else
            return convert(Floathing, default(config, param_name, nothing))
        end
    end

    return extract_parameter(Component, config, param_name, param_def, sim_params, uac)
end

function validate_config(x::Type{SolarthermalCollector}, config::Dict{String,Any}, extracted::Dict{String,Any},
                         uac::String, sim_params::Dict{String,Any})
    validate_config(Component, extracted, uac, sim_params, component_parameters(SolarthermalCollector))
end

function init_from_params(x::Type{SolarthermalCollector}, uac::String, params::Dict{String,Any},
                          raw_params::Dict{String,Any}, sim_params::Dict{String,Any})::Tuple
    m_heat_out = Symbol(params["m_heat_out"])
    register_media([m_heat_out])

    const_wind_speed = params["constant_wind_speed"]
    if const_wind_speed === nothing
        const_wind_speed = 0.0
    else
        const_wind_speed = max(0.0, const_wind_speed * params["wind_speed_reduction"] - 3.0)
    end

    return (uac,
            Controller(params["control_parameters"]),
            sf_flexible_source,
            InterfaceMap(),
            InterfaceMap(m_heat_out => nothing),
            m_heat_out,
            params["collector_gross_area"],
            params["tilt_angle"],
            params["azimuth_angle"],
            params["ground_reflectance"],
            params["eta_0_b"],
            # K_b_array
            vcat([10, 20, 30, 40, 50, 60, 70, 80, 90]', params["K_b_t_array"]', params["K_b_l_array"]'),
            ([], []), # K_b_itp
            1,        # K_b
            params["K_d"],
            params["a_params"],
            5.670374419 * 10^-8, # sigma
            params["vol_heat_capacity"],
            some_or_none(params["ambient_temperature_profile_file_path"],
                         params["ambient_temperature_from_global_file"]),
            some_or_none(params["beam_solar_radiation_profile_file_path"],
                         params["beam_solar_radiation_from_global_file"]),
            some_or_none(params["diffuse_solar_radiation_profile_file_path"],
                         params["diffuse_solar_radiation_from_global_file"]),
            some_or_none(params["infrared_sky_radiation_profile_file_path"],
                         params["infrared_sky_radiation_from_global_file"]),
            some_or_none(params["wind_speed_profile_file_path"],
                         params["wind_speed_from_global_file"]),
            params["wind_speed_reduction"],
            params["constant_ambient_temperature"],
            0.0, # beam_solar_irradiance_in_plane
            0.0, # diffuse_solar_irradiance_in_plane
            const_wind_speed,
            params["constant_beam_solar_radiation"],
            params["constant_diffuse_solar_radiation"],
            params["constant_infrared_sky_radiation"],
            0.0, # direct_normal_irradiance
            params["delta_T"],
            params["spec_flow_rate"],
            params["delta_T_min"],
            params["spec_flow_rate_min"],
            [],   # used_energy
            0.0,  # output_temperature
            0.0,  # average_temperature
            0.0,  # last_average_temperature
            0.0,  # spec_flow_rate_actual
            0.0,  # delta_T_actual
            [],   # available_energies
            [],   # average_temperatures
            [],   # output_temperatures
            [],   # runtimes
            Dict{Temperature,Tuple}(), # temperature_energy_pairs
            [],   # calc_uacs
            10.0) # mean_ambient_temperature
end

function initialise!(unit::SolarthermalCollector, sim_params::Dict{String,Any})
    set_storage_transfer!(unit.output_interfaces[unit.m_heat_out],
                          load_storages(unit.controller, unit.m_heat_out))

    if sim_params["longitude"] === nothing || sim_params["latitude"] === nothing
        @error "Longitude and latitude must be provided through a weather file or in the 
        simulation parameters to calculate the sun position."
        throw(InputError())
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
        unit.reduced_wind_speed = max(Profiles.value_at_time(unit.wind_speed_profile, sim_params) *
                                      unit.wind_speed_reduction - 3, 0)
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

    unit.available_energies,
    has_calculated_all_maxima,
    unit.average_temperatures,
    unit.output_temperatures,
    unit.calc_uacs = calculate_energies(unit, components, sim_params)

    set_max_energy!(unit.output_interfaces[unit.m_heat_out],
                    unit.available_energies,
                    [nothing for _ in 1:(length(unit.available_energies))],
                    unit.output_temperatures,
                    unit.calc_uacs,
                    has_calculated_all_maxima,
                    false)
end

function process(unit::SolarthermalCollector, sim_params::Dict{String,Any})
    exchanges = balance_on(unit.output_interfaces[unit.m_heat_out], unit.output_interfaces[unit.m_heat_out].target)
    energy_demands = [abs(e.balance + e.energy_potential) for e in exchanges]

    unit.used_energy = zeros(length(unit.calc_uacs))
    unit.runtimes = zeros(length(unit.calc_uacs))
    purpose_uac = Array{Stringing}(nothing, length(unit.calc_uacs))

    if sum(energy_demands) > 0 && _sum(unit.available_energies) > 0
        temperature_min = Array{Temperature}(nothing, length(unit.calc_uacs))
        temperature_max = Array{Temperature}(nothing, length(unit.calc_uacs))

        for e in exchanges
            uac_idx = findfirst(==(e.purpose_uac), unit.calc_uacs)
            comp_energy_demand = abs(e.balance + e.energy_potential)
            used_energy_comp = min(comp_energy_demand, unit.available_energies[uac_idx])

            unit.used_energy[uac_idx] += used_energy_comp
            temperature_min[uac_idx] = unit.output_temperatures[uac_idx]
            temperature_max[uac_idx] = unit.output_temperatures[uac_idx]
            purpose_uac[uac_idx] = e.purpose_uac
            if unit.available_energies[uac_idx] > 0
                unit.runtimes[uac_idx] = used_energy_comp / unit.available_energies[uac_idx]
            end
        end

        mask = abs.(unit.used_energy) .> sim_params["epsilon"]
        if any(mask)
            add!(unit.output_interfaces[unit.m_heat_out],
                 unit.used_energy[mask],
                 temperature_min[mask],
                 temperature_max[mask],
                 purpose_uac[mask])
        else
            set_max_energy!(unit.output_interfaces[unit.m_heat_out], 0.0)
        end
    else
        set_max_energy!(unit.output_interfaces[unit.m_heat_out], 0.0)
    end

    # correct runtimes to make sure the amount of non nothing runtimes is the same as non 
    # nothing temperatures and the sum of runtimes is normalized to cover the whole time step
    runtimes_sum = _sum(unit.runtimes)
    for idx in eachindex(unit.runtimes)
        if isnothing(unit.average_temperatures[idx])
            unit.runtimes[idx] = nothing
            unit.used_energy[idx] = nothing
        elseif runtimes_sum > 0
            unit.runtimes[idx] = unit.runtimes[idx] / runtimes_sum
        end
    end

    if isnothing(_sum(unit.runtimes)) || _sum(unit.runtimes) == 0
        stagnation_temp = get_stagnation_temperature(unit, sim_params)
        unit.output_temperature = stagnation_temp
        unit.average_temperature = stagnation_temp
        unit.spec_flow_rate_actual = 0
        unit.delta_T_actual = 0
    else
        # recalculate the average_temperature to keep a constant flow rate if necessary
        if unit.spec_flow_rate !== nothing
            for uac_idx in eachindex(unit.runtimes)
                if !isnothing(unit.runtimes[uac_idx]) && unit.runtimes[uac_idx] != 0
                    unit.average_temperatures[uac_idx] = unit.output_temperatures[uac_idx] -
                                                         (sim_params["wh_to_watts"](unit.used_energy[uac_idx]) /
                                                          unit.runtimes[uac_idx]) /
                                                         unit.collector_gross_area /
                                                         (unit.spec_flow_rate * 2 * unit.vol_heat_cap)
                end
            end
        end

        unit.output_temperature = _weighted_mean(unit.output_temperatures, unit.runtimes)
        unit.average_temperature = _weighted_mean(unit.average_temperatures, unit.runtimes)

        unit.spec_flow_rate_actual = sim_params["wh_to_watts"](_sum(unit.used_energy)) /
                                     unit.collector_gross_area /
                                     ((unit.output_temperature - unit.average_temperature) * 2 * unit.vol_heat_cap)
        unit.delta_T_actual = (unit.output_temperature - unit.average_temperature) * 2

        unit.output_temperature = _weighted_mean(unit.output_temperatures, unit.used_energy)
    end

    unit.last_average_temperature = unit.average_temperature
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
- `available_energies::Array{Floathing}`: The available energy to each connected component with uac as keys.
- `calc_max::Bool`: A Bool indicating if max energy was calculated or not.
- `average_temperatures::Array{Temperature}`: The average collector temperatures for each connected component with uac as keys.
- `output_temperatures::Array{Temperature}`: The output collector temperatures for each connected component with uac as keys.
- `runtimes::Array{Floathing}`: The partial time as part of 1 each connected component gets delivered energy with uac as keys.
"""
function calculate_energies(unit::SolarthermalCollector, components, sim_params::Dict{String,Any})
    exchanges = balance_on(unit.output_interfaces[unit.m_heat_out], unit.output_interfaces[unit.m_heat_out].target)

    average_temperatures = Array{Floathing}(nothing, length(exchanges))
    output_temperatures = Array{Floathing}(nothing, length(exchanges))
    max_energies = zeros(Float64, length(exchanges))
    target_uac = Array{Stringing}(nothing, length(exchanges))
    for (ex_idx, exchange) in enumerate(exchanges)
        # no check for possibly written max energy --> is done by set_max_energy!() for 
        # 1-to-1 connections and for a connection to a bus, the interface is only used by 
        # the current component
        success, exchange_temperature, max_energy = determine_temperature_and_energy(unit.controller,
                                                                                     components,
                                                                                     unit.uac,
                                                                                     exchange.purpose_uac,
                                                                                     sim_params)
        if !success
            # no control module is provided between source and target
            # set temperature to exchange.temperature_min --> can also be nothing, but in 
            # case it is given we will use it if no control module is used
            exchange_temperature, max_energy = check_temperature_and_get_max_energy(unit,
                                                                                    sim_params,
                                                                                    exchange.temperature_min,
                                                                                    false)
        end

        if max_energy == 0
            t_avg = nothing
            exchange_temperature = nothing
        else
            _, t_avg, exchange_temperature = unit.temperature_energy_pairs[exchange_temperature]
        end

        target_uac[ex_idx] = exchange.purpose_uac
        output_temperatures[ex_idx] = exchange_temperature
        max_energies[ex_idx] = max_energy
        average_temperatures[ex_idx] = t_avg
    end

    if length(exchanges) > 1
        has_calculated_all_maxima = true
    else
        has_calculated_all_maxima = false
    end

    return (max_energies,
            has_calculated_all_maxima,
            average_temperatures,
            output_temperatures,
            target_uac)
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
                                                                        sim_params["time_step_seconds"], 0)),
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
                                                                        sim_params["time_step_seconds"], 0)),
                             unit.last_average_temperature,
                             Roots.Newton())
        current_max_temperature = avg_temp + unit.delta_T_min / 2
        energy_at_max_temp = unit.spec_flow_rate * unit.delta_T_min * unit.vol_heat_cap

    else
        @error "Error in config file: Exclusively delta_T OR spec_flow_rate must have a value"
        throw(InputError())
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
                                                                           sim_params["time_step_seconds"], 0)),
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
        throw(InputError())
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

    if temperature_output === nothing || temperature_output > source_max_out_temperature
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
                            unit::SolarthermalCollector, time_step:Number, 
                            independent_target_temp::Int64=1)

Derivative of the thermal power output function to speed up solve function for 
spec_thermal_power_func

# Arguments
- `t_avg::Temperature`: Average temperature of the collector. Value to solve after.
- `spec_flow_rate::Float64`: Specific flow rate through the collector in m³/(s*m²)
- `unit::SolarthermalCollector`: The solarthermal collector unit for which the calculation 
                                 is performed.
- `time_step::Number`: time_step of the simulation in seconds
- `independent_target_temp::Int64`: If the target temperature is dependent on t_avg set this 
                                    value to 0 otherwise its 1. Corrects the derivate.
 
# Returns
- 'Float64': Derivative of difference between generated energy, losses and power output
"""
function derivate_spec_thermal_power_func(t_avg, spec_flow_rate, unit::SolarthermalCollector,
                                          time_step, independent_target_temp::Int64=1)
    unit.a_params[1] * (-1) -
    unit.a_params[2] * 2 * (t_avg - unit.ambient_temperature) -
    unit.a_params[3] * unit.reduced_wind_speed -
    unit.a_params[5] * 1 / time_step -
    unit.a_params[8] * 4 * (t_avg - unit.ambient_temperature)^3 +
    spec_flow_rate * 2 * unit.vol_heat_cap * independent_target_temp
end

function output_values(unit::SolarthermalCollector)::Vector{String}
    return [string(unit.m_heat_out) * ":OUT",
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
        return unit.delta_T_actual
    elseif key.value_key == "spec_flow_rate"
        return unit.spec_flow_rate_actual
    end
    throw(KeyError(key.value_key))
end

export SolarthermalCollector
