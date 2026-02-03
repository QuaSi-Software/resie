using GLMakie
using Roots
using Plots: plot, savefig
using Plots: Plots
using SparseArrays
using LinearAlgebra

#! format: off
const GEOTHERMAL_HEAT_COLLECTOR_PARAMETERS = Dict(
    "m_heat_in" => (
        default="m_h_w_ht1",
        description="Heat input medium (for regeneration/loading)",
        display_name="Medium heat_in",
        required=false,
        type=String,
        json_type="string",
        unit="-"
    ),
    "m_heat_out" => (
        default="m_h_w_lt1",
        description="Heat output medium (for extraction/unloading)",
        display_name="Medium heat_out",
        required=false,
        type=String,
        json_type="string",
        unit="-"
    ),
    "model_type" => (
        default="simplified",
        description="Operation model: currently 'simplified', with constant fluid-to-soil " *
                    "resistance, and 'detailed' with calculated fluid-to-soil resistance " *
                    "in every time step, are available.",
        display_name="Model type",
        required=false,
        type=String,
        json_type="string",
        options=["simplified", "detailed"],
        unit="-"
    ),
    "accuracy_mode" => (
        default="normal",
        description="Simulation accuracy mode",
        display_name="Accuracy mode",
        required=false,
        type=String,
        json_type="string",
        options=["very_rough", "rough", "normal", "high", "very_high"],
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
        description="If given points to a key in the global weather data file with the " *
                    "ambient temperature profile to be used",
        display_name="Global file amb. temp. key",
        required=false,
        conditionals=[
            ("consider_losses", "is_true"),
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
        type=Float64,
        json_type="number",
        unit="°C"
    ),
    "global_solar_radiation_profile_file_path" => (
        default=nothing,
        description="Path to global solar radiation profile file",
        display_name="Global solar radiation profile",
        required=false,
        conditionals=[
            ("global_solar_radiation_from_global_file", "mutex"),
            ("constant_global_solar_radiation", "mutex")
        ],
        type=String,
        json_type="string",
        unit="-"
    ),
    "global_solar_radiation_from_global_file" => (
        default=nothing,
        description="If given points to a key in the global weather data file with the " *
                    "global solar radiation profile to be used",
        display_name="Global file glob. rad. key",
        required=false,
        conditionals=[
            ("global_solar_radiation_profile_file_path", "mutex"),
            ("constant_global_solar_radiation", "mutex")
        ],
        type=String,
        json_type="string",
        unit="-"
    ),
    "constant_global_solar_radiation" => (
        default=nothing,
        description="Constant global solar radiation value",
        display_name="Constant global solar radiation",
        required=false,
        conditionals=[
            ("global_solar_radiation_profile_file_path", "mutex"),
            ("global_solar_radiation_from_global_file", "mutex")
        ],
        type=Float64,
        json_type="number",
        unit="Wh/m^2"
    ),
    "infrared_sky_radiation_profile_file_path" => (
        default=nothing,
        description="Path to infrared sky radiation profile file",
        display_name="Infrared radiation profile",
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
        description="If given points to a key in the global weather data file with the " *
                    "infrared sky radiation profile to be used",
        display_name="Global file infr. sky rad. key",
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
        description="Constant infrared sky radiation value",
        display_name="Constant infrared radiation",
        required=false,
        conditionals=[
            ("infrared_sky_radiation_profile_file_path", "mutex"),
            ("infrared_sky_radiation_from_global_file", "mutex")
        ],
        type=Float64,
        json_type="number",
        unit="Wh/m^2"
    ),
    "pipe_radius_outer" => (
        default=0.016,
        description="Outer radius of the collector pipe",
        display_name="Pipe outer radius",
        required=false,
        validations=[
            ("self", "value_gt_num", 0.0)
        ],
        type=Float64,
        json_type="number",
        unit="m"
    ),
    "pipe_thickness" => (
        default=0.003,
        description="Thickness of the collector pipe",
        display_name="Pipe thickness",
        required=false,
        validations=[
            ("self", "value_gt_num", 0.0)
        ],
        type=Float64,
        json_type="number",
        unit="m"
    ),
    "pipe_heat_conductivity" => (
        default=0.4,
        description="Heat conductivity of pipe material",
        display_name="Pipe heat conductivity",
        required=false,
        validations=[
            ("self", "value_gt_num", 0.0)
        ],
        type=Float64,
        json_type="number",
        unit="W/m*K"
    ),
    "pipe_laying_depth" => (
        default=1.5,
        description="Depth of pipe system below ground surface",
        display_name="Pipe laying depth",
        required=false,
        validations=[
            ("self", "value_gt_num", 0.0)
        ],
        type=Float64,
        json_type="number",
        unit="m"
    ),
    "pipe_length" => (
        default=100,
        description="Length of one collector pipe",
        display_name="Pipe length",
        required=false,
        validations=[
            ("self", "value_gt_num", 0.0)
        ],
        type=Float64,
        json_type="number",
        unit="m"
    ),
    "number_of_pipes" => (
        default=1,
        description="Number of parallel pipes, each with a length of 'pipe_length'",
        display_name="Number of pipes",
        required=false,
        validations=[
            ("self", "value_gte_num", 1.0)
        ],
        type=Float64,
        json_type="number",
        unit="-"
    ),
    "pipe_spacing" => (
        default=0.5,
        description="Distance between pipes",
        display_name="Pipe spacing",
        required=false,
        validations=[
            ("self", "value_gt_num", 0.0)
        ],
        type=Float64,
        json_type="number",
        unit="m"
    ),
    "pipe_soil_thermal_resistance" => (
        default=0.1,
        description="Thermal resistance between pipe and soil",
        display_name="Pipe-soil thermal resistance",
        required=false,
        validations=[
            ("self", "value_gt_num", 0.0)
        ],
        type=Float64,
        json_type="number",
        unit="m*K/W"
    ),
    "considered_soil_depth" => (
        default=10.0,
        description="Depth of soil considered in simulation",
        display_name="Considered soil depth",
        required=false,
        validations=[
            ("self", "value_gt_rel", "pipe_laying_depth")
        ],
        type=Float64,
        json_type="number",
        unit="m"
    ),
    "soil_density" => (
        default=2000,
        description="Density of unfrozen soil",
        display_name="Soil density",
        required=false,
        validations=[
            ("self", "value_gt_num", 0.0)
        ],
        type=Float64,
        json_type="number",
        unit="kg/m^3"
    ),
    "soil_specific_heat_capacity" => (
        default=1000,
        description="Specific heat capacity of unfrozen soil",
        display_name="Soil specific heat capacity",
        required=false,
        validations=[
            ("self", "value_gt_num", 0.0)
        ],
        type=Float64,
        json_type="number",
        unit="J/kg*K"
    ),
    "soil_specific_heat_capacity_frozen" => (
        default=900,
        description="Specific heat capacity of frozen soil",
        display_name="Soil specific heat capacity (frozen)",
        required=false,
        validations=[
            ("self", "value_gt_num", 0.0)
        ],
        type=Float64,
        json_type="number",
        unit="J/kg*K"
    ),
    "soil_heat_conductivity" => (
        default=1.5,
        description="Heat conductivity of unfrozen soil (lambda)",
        display_name="Soil heat conductivity",
        required=false,
        validations=[
            ("self", "value_gt_num", 0.0)
        ],
        type=Float64,
        json_type="number",
        unit="W/m*K"
    ),
    "soil_heat_conductivity_frozen" => (
        default=2.0,
        description="Heat conductivity of frozen soil (lambda)",
        display_name="Soil heat conductivity (frozen)",
        required=false,
        validations=[
            ("self", "value_gt_num", 0.0)
        ],
        type=Float64,
        json_type="number",
        unit="W/m*K"
    ),
    "soil_specific_enthalpy_of_fusion" => (
        default=90000,
        description="Specific enthalpy of fusion of soil",
        display_name="Soil enthalpy of fusion",
        required=false,
        validations=[
            ("self", "value_gt_num", 0.0)
        ],
        type=Float64,
        json_type="number",
        unit="J/kg"
    ),
    "phase_change_upper_boundary_temperature" => (
        default=-0.25,
        description="Upper boundary temperature for phase change",
        display_name="Phase change upper boundary",
        required=false,
        type=Float64,
        json_type="number",
        unit="°C"
    ),
    "phase_change_lower_boundary_temperature" => (
        default=-1,
        description="Lower boundary temperature for phase change",
        display_name="Phase change lower boundary",
        required=false,
        type=Float64,
        json_type="number",
        unit="°C"
    ),
    "surface_convective_heat_transfer_coefficient" => (
        default=14.7,
        description="Convective heat transfer coefficient at surface",
        display_name="Surface convective heat transfer",
        required=false,
        validations=[
            ("self", "value_gt_num", 0.0)
        ],
        type=Float64,
        json_type="number",
        unit="W/m^2*K"
    ),
    "surface_reflection_factor" => (
        default=0.25,
        description="Surface reflection factor (albedo)",
        display_name="Surface reflection factor",
        required=false,
        validations=[
            ("self", "value_gte_num", 0.0),
            ("self", "value_lte_num", 1.0)
        ],
        type=Float64,
        json_type="number",
        unit="-"
    ),
    "surface_emissivity" => (
        default=0.9,
        description="Surface emissivity on ground surface",
        display_name="Surface emissivity",
        required=false,
        validations=[
            ("self", "value_gte_num", 0.0),
            ("self", "value_lte_num", 1.0)
        ],
        type=Float64,
        json_type="number",
        unit="-"
    ),
    "fluid_specific_heat_capacity" => (
        default=3800,
        description="Specific heat capacity of heat transfer fluid",
        display_name="Fluid specific heat capacity",
        required=false,
        validations=[
            ("self", "value_gt_num", 0.0)
        ],
        type=Float64,
        json_type="number",
        unit="J/kg*K"
    ),
    "fluid_density" => (
        default=1045,
        description="Density of heat transfer fluid; default for 30 % glycol at 0 °C",
        display_name="Fluid density",
        required=false,
        validations=[
            ("self", "value_gt_num", 0.0)
        ],
        type=Float64,
        json_type="number",
        unit="kg/m^3"
    ),
    "fluid_prandtl_number" => (
        default=30,
        description="Prandtl number of heat transfer fluid; default for 30 % glycol at 0 °C",
        display_name="Fluid Prandtl number",
        required=false,
        validations=[
            ("self", "value_gt_num", 0.0)
        ],
        type=Float64,
        json_type="number",
        unit="-"
    ),
    "fluid_kinematic_viscosity" => (
        default=3.9e-6,
        description="Kinematic viscosity of heat transfer fluid; default for 30 % glycol at 0 °C",
        display_name="Fluid kinematic viscosity",
        required=false,
        validations=[
            ("self", "value_gt_num", 0.0)
        ],
        type=Float64,
        json_type="number",
        unit="m^2/s"
    ),
    "fluid_heat_conductivity" => (
        default=0.5,
        description="Heat conductivity of heat transfer fluid; default for 30 % glycol at 0 °C",
        display_name="Fluid heat conductivity",
        required=false,
        validations=[
            ("self", "value_gt_num", 0.0)
        ],
        type=Float64,
        json_type="number",
        unit="W/m*K"
    ),
    "use_dynamic_fluid_properties" => (
        default=false,
        description="Use temperature-dependent fluid properties; false for constant, true " *
                    "for temperature-dependent fluid properties accoring to TRNSYS Type 710",
        display_name="Use dynamic fluid properties",
        required=false,
        type=Bool,
        json_type="boolean",
        unit="-"
    ),
    "nusselt_approach" => (
        default="Stephan",
        description="Approach for Nusselt number calculation",
        display_name="Nusselt approach",
        required=false,
        type=String,
        json_type="string",
        options=["Stephan", "Ramming"],
        unit="-"
    ),
    "start_temperature_fluid_and_pipe" => (
        default=15.5,
        description="Starting temperature of fluid and soil near pipe during initialisation",
        display_name="Start temperature",
        required=false,
        type=Float64,
        json_type="number",
        unit="°C"
    ),
    "undisturbed_ground_temperature" => (
        default=9.0,
        description="Undisturbed ground temperature at lower simulation boundary",
        display_name="Undisturbed ground temperature",
        required=false,
        type=Float64,
        json_type="number",
        unit="°C"
    ),
    "max_output_power" => (
        default=20,
        description="Maximum output power per unit area. Depends on ground and climate " *
                    "localization. See VDI 4640-2.",
        display_name="Maximum output power",
        required=false,
        validations=[
            ("self", "value_gt_num", 0.0)
        ],
        type=Float64,
        json_type="number",
        unit="W/m^2"
    ),
    "max_input_power" => (
        default=20,
        description="Maximum input power per unit area. Depends on ground and climate " *
                    "localization. See VDI 4640-2",
        display_name="Maximum input power",
        required=false,
        validations=[
            ("self", "value_gt_num", 0.0)
        ],
        type=Float64,
        json_type="number",
        unit="W/m^2"
    ),
    "unloading_temperature_spread" => (
        default=3.0,
        description="Temperature spread between forward and return during unloading",
        display_name="Unloading temperature spread",
        required=false,
        validations=[
            ("self", "value_gte_num", 0.0)
        ],
        type=Float64,
        json_type="number",
        unit="K"
    ),
    "loading_temperature_spread" => (
        default=3.0,
        description="Temperature spread between forward and return during loading",
        display_name="Loading temperature spread",
        required=false,
        validations=[
            ("self", "value_gte_num", 0.0)
        ],
        type=Float64,
        json_type="number",
        unit="K"
    ),
    "fluid_min_output_temperature" => (
        default=nothing,
        description="Minimum output temperature of fluid during unloading",
        display_name="Min output temperature",
        required=false,
        type=Temperature,
        json_type="number",
        unit="°C"
    ),
    "fluid_max_input_temperature" => (
        default=nothing,
        description="Maximum input temperature of fluid during loading",
        display_name="Max input temperature",
        required=false,
        type=Temperature,
        json_type="number",
        unit="°C"
    ),
    "regeneration" => (
        default=true,
        description="Enable regeneration (loading) of the collector",
        display_name="Regeneration enabled",
        required=false,
        type=Bool,
        json_type="boolean",
        unit="-"
    ),
    "max_picard_iter" => (
        default=3,
        description="Maximum number of iterations for implicit solving (update of k(T))",
        display_name="Max Picard iterations",
        required=false,
        validations=[
            ("self", "value_gte_num", 1.0)
        ],
        type=Int,
        json_type="integer",
        unit="-"
    ),
    "picard_tol" => (
        default=1e-3,
        description="Absolute tolerance for Picard iteration",
        display_name="Picard tolerance",
        required=false,
        validations=[
            ("self", "value_gt_num", 0.0)
        ],
        type=Float64,
        json_type="number",
        unit="K"
    ),
)
#! format: on

"""
Implementation of a geothermal heat collector.
This implementation acts as storage as is can produce and load energy.

Possible improvements in the future:
- Addition of convective heat transport at the earth's surface as a function of wind speed
  and ambient temperature instead of a constant factor
- Adaptation of the model to be able to simulate a single pipe (requires adaptation of the
  boundary conditions)
"""
mutable struct GeothermalHeatCollector <: Component
    uac::String
    controller::Controller
    sys_function::SystemFunction
    input_interfaces::InterfaceMap
    output_interfaces::InterfaceMap
    m_heat_in::Symbol
    m_heat_out::Symbol

    ambient_temperature_profile::Union{Profile,Nothing}
    constant_ambient_temperature::Temperature
    global_radiation_profile::Union{Profile,Nothing}
    constant_global_radiation::Union{Float64,Nothing}
    infrared_sky_radiation_profile::Union{Profile,Nothing}
    constant_infrared_sky_radiation::Union{Float64,Nothing}
    unloading_temperature_spread::Temperature
    fluid_min_output_temperature::Temperature
    fluid_max_input_temperature::Temperature
    loading_temperature_spread::Temperature

    # depends on ground and climate localization. [VDI 4640-2.]
    max_output_power::Union{Nothing,Float64}
    # depends on ground and climate localization. [VDI 4640-2.]
    max_input_power::Union{Nothing,Float64}
    regeneration::Bool
    max_output_energy::Float64
    current_max_output_energy::Float64
    max_input_energy::Float64
    current_max_input_energy::Float64

    current_output_temperature::Temperature
    current_input_temperature::Temperature
    ambient_temperature::Temperature

    soil_specific_heat_capacity::Float64
    soil_specific_heat_capacity_frozen::Float64
    soil_density::Float64
    soil_heat_conductivity::Float64
    soil_heat_conductivity_frozen::Float64
    soil_density_vector::Vector
    soil_specific_enthalpy_of_fusion::Float64

    phase_change_upper_boundary_temperature::Temperature
    phase_change_lower_boundary_temperature::Temperature

    surface_convective_heat_transfer_coefficient::Float64
    surface_reflection_factor::Float64

    # holds the temperature of the last timestep
    t1::Array{Float64}
    # holds the temperature of the current timestep
    t2::Array{Float64}
    # holds the specific heat capacity for each node
    cp::Array{Float64}
    # precalculated density * volume of the soil around each node
    soil_weight::Array{Float64}
    # identifies the fluid node index in y direction
    fluid_node_y_idx::Int
    # identifies nodes surrounding the fluid-node: pipe-surrounding: true; else: false
    is_pipe_surrounding::Array{Bool}

    pipe_radius_outer::Float64
    pipe_thickness::Float64
    # pipe inner diameter
    pipe_d_i::Float64
    # pipe outer diameter
    pipe_d_o::Float64
    pipe_laying_depth::Float64
    pipe_length::Float64
    number_of_pipes::Float64
    pipe_spacing::Float64
    considered_soil_depth::Float64
    accuracy_mode::String
    model_type::String
    pipe_soil_thermal_resistance::Floathing

    # solar global radiation on horizontal surface, to be read from weather-profile
    global_radiation_power::Float64
    # boltzmann_constant [W/(m^2 K^4)]: Stefan-Boltzmann-Constant
    boltzmann_constant::Float64
    surface_emissivity::Float64

    # horizontal dimension parallel to ground surface and orthogonal to pipe
    dx::Vector{Float64}
    # vertical dimension orthogonal to ground surface and orthogonal to pipe
    dy::Vector{Float64}
    # this is dx between the nodes (while dx is the x-width assigned to each node)
    dx_mesh::Vector{Float64}
    # this is dy between the nodes (while dy is the y-width assigned to each node)
    dy_mesh::Vector{Float64}
    # horizontal dimension parallel to ground surface and parallel to pipe (equals the length of one pipe, constant)
    dz::Float64

    # is set to fluid_start_temperature at the beginning
    fluid_temperature::Temperature

    # total heat flux in or out of collector
    collector_total_heat_energy_in_out::Float64

    pipe_heat_conductivity::Float64
    fluid_specific_heat_capacity::Float64
    fluid_prandtl_number::Float64
    fluid_density::Float64
    fluid_kinematic_viscosity::Float64
    fluid_heat_conductivity::Float64
    use_dynamic_fluid_properties::Bool
    nusselt_approach::String

    # to be calculated by a function
    fluid_reynolds_number::Float64
    # convective heat transfer coefficient between fluid and pipe
    alpha_fluid_pipe::Float64
    average_temperature_adjacent_to_pipe::Float64
    undisturbed_ground_temperature::Float64

    # holds temperature field of nodes for output plot
    temp_field_output::Array{Float64}
    # precalculated parameter for freezing function
    sigma_lat::Float64
    # precalculated parameter for freezing function
    t_lat::Float64
    # precalculated parameter for freezing function
    delta_t_lat::Float64
    # precalculated parameter
    volume_adjacent_to_pipe::Float64

    # indicating if the process step has already been performed in the current time step
    process_done::Bool
    # indicating if the load step has already been performed in the current time step
    load_done::Bool

    max_picard_iter::Int
    picard_tol::Float64

    function GeothermalHeatCollector(uac::String, config::Dict{String,Any}, sim_params::Dict{String,Any})
        new(SSOT_parameter_constructor(GeothermalHeatCollector, uac, config, sim_params)...)
    end
end

function component_parameters(x::Type{GeothermalHeatCollector})::Dict{String,NamedTuple}
    return deepcopy(GEOTHERMAL_HEAT_COLLECTOR_PARAMETERS) # return a copy to prevent external modification
end

function extract_parameter(x::Type{GeothermalHeatCollector}, config::Dict{String,Any}, param_name::String,
                           param_def::NamedTuple, sim_params::Dict{String,Any}, uac::String)
    if param_name in ("ambient_temperature_from_global_file", "global_solar_radiation_from_global_file",
                      "infrared_sky_radiation_from_global_file")
        return load_profile_from_global_weather_file(config, param_name, sim_params, uac)
    elseif param_name in ("ambient_temperature_profile_file_path", "global_solar_radiation_profile_file_path",
                          "infrared_sky_radiation_profile_file_path")
        return load_optional_profile(config, param_name, sim_params)
    elseif param_name in ("constant_ambient_temperature", "constant_global_solar_radiation",
                          "constant_infrared_sky_radiation")
        if occursin("temperature", param_name)
            return convert(Temperature, default(config, param_name, nothing))
        else
            return convert(Floathing, default(config, param_name, nothing))
        end
    end

    return extract_parameter(Component, config, param_name, param_def, sim_params, uac)
end

function validate_config(x::Type{GeothermalHeatCollector}, config::Dict{String,Any}, extracted::Dict{String,Any},
                         uac::String, sim_params::Dict{String,Any})
    validate_config(Component, extracted, uac, sim_params, component_parameters(GeothermalHeatCollector))
end

function init_from_params(x::Type{GeothermalHeatCollector}, uac::String, params::Dict{String,Any},
                          raw_params::Dict{String,Any}, sim_params::Dict{String,Any})::Tuple
    m_heat_in = Symbol(params["m_heat_in"])
    m_heat_out = Symbol(params["m_heat_out"])
    register_media([m_heat_in, m_heat_out])

    return (uac,
            Controller(params["control_parameters"]),
            sf_storage,
            InterfaceMap(m_heat_in => nothing),
            InterfaceMap(m_heat_out => nothing),
            m_heat_in,
            m_heat_out,
            some_or_none(params["ambient_temperature_profile_file_path"],
                         params["ambient_temperature_from_global_file"]),
            params["constant_ambient_temperature"],
            some_or_none(params["global_solar_radiation_profile_file_path"],
                         params["global_solar_radiation_from_global_file"]),
            params["constant_global_solar_radiation"],
            some_or_none(params["infrared_sky_radiation_profile_file_path"],
                         params["infrared_sky_radiation_from_global_file"]),
            params["constant_infrared_sky_radiation"],
            params["unloading_temperature_spread"],
            params["fluid_min_output_temperature"],
            params["fluid_max_input_temperature"],
            params["loading_temperature_spread"],
            params["max_output_power"],
            params["max_input_power"],
            params["regeneration"],
            0.0, # max_output_energy
            0.0, # current_max_output_energy
            0.0, # max_input_energy
            0.0, # current_max_input_energy
            0.0, # current_output_temperature
            0.0, # current_input_temperature
            0.0, # ambient_temperature
            params["soil_specific_heat_capacity"],
            params["soil_specific_heat_capacity_frozen"],
            params["soil_density"],
            params["soil_heat_conductivity"],
            params["soil_heat_conductivity_frozen"],
            Array{Float64}(undef, 0), # soil_density_vector
            params["soil_specific_enthalpy_of_fusion"],
            params["phase_change_upper_boundary_temperature"],
            params["phase_change_lower_boundary_temperature"],
            params["surface_convective_heat_transfer_coefficient"],
            params["surface_reflection_factor"],
            Array{Float64}(undef, 0, 0), # t1
            Array{Float64}(undef, 0, 0), # t2
            Array{Float64}(undef, 0, 0), # cp
            Array{Float64}(undef, 0, 0), # soil_weight
            0,                           # fluid_node_y_idx
            Array{Bool}(undef, 0, 0),    # is_pipe_surrounding
            params["pipe_radius_outer"],
            params["pipe_thickness"],
            0.0,                         # pipe_d_i
            0.0,                         # pipe_d_o
            params["pipe_laying_depth"],
            params["pipe_length"],
            params["number_of_pipes"],
            params["pipe_spacing"],
            params["considered_soil_depth"],
            params["accuracy_mode"],
            params["model_type"],
            params["pipe_soil_thermal_resistance"],
            0.0,                      # global_radiation_power
            5.6697e-8,                # boltzmann_constant
            params["surface_emissivity"],
            Array{Float64}(undef, 0), # dx
            Array{Float64}(undef, 0), # dy
            Array{Float64}(undef, 0), # dx_mesh
            Array{Float64}(undef, 0), # dy_mesh
            0.0,                      # dz
            0.0,                      # fluid_temperature 
            0.0,                      # collector_total_heat_energy_in_out
            params["pipe_heat_conductivity"],
            params["fluid_specific_heat_capacity"],
            params["fluid_prandtl_number"],
            params["fluid_density"],
            params["fluid_kinematic_viscosity"],
            params["fluid_heat_conductivity"],
            params["use_dynamic_fluid_properties"],
            params["nusselt_approach"],
            0.0,                            # fluid_reynolds_number
            0.0,                            # alpha_fluid_pipe
            params["start_temperature_fluid_and_pipe"],
            params["undisturbed_ground_temperature"],
            Array{Float64}(undef, 0, 0, 0), # temp_field_output
            0.0,                            # sigma_lat
            0.0,                            # t_lat
            0.0,                            # delta_t_lat
            0.0,                            # volume_adjacent_to_pipe
            false,                          # process_done,
            false,                          # load_done
            params["max_picard_iter"],
            params["picard_tol"])
end
function initialise!(unit::GeothermalHeatCollector, sim_params::Dict{String,Any})
    if unit.regeneration
        set_storage_transfer!(unit.input_interfaces[unit.m_heat_in],
                              unload_storages(unit.controller, unit.m_heat_in))
    end
    set_storage_transfer!(unit.output_interfaces[unit.m_heat_out],
                          load_storages(unit.controller, unit.m_heat_out))

    # calculate diameters of pipe
    unit.pipe_d_i = 2 * unit.pipe_radius_outer - (2 * unit.pipe_thickness)
    unit.pipe_d_o = 2 * unit.pipe_radius_outer

    # calculate max_energy
    A_collector = unit.pipe_length * unit.pipe_spacing * unit.number_of_pipes
    unit.max_output_energy = sim_params["watt_to_wh"](unit.max_output_power * A_collector)
    if unit.regeneration
        unit.max_input_energy = sim_params["watt_to_wh"](unit.max_input_power * A_collector)
    else
        unit.max_input_energy = 0.0
    end

    # calculate coefficients for freezing function 
    unit.sigma_lat = (unit.phase_change_upper_boundary_temperature - unit.phase_change_lower_boundary_temperature) / 5
    unit.t_lat = (unit.phase_change_upper_boundary_temperature + unit.phase_change_lower_boundary_temperature) / 2
    unit.delta_t_lat = unit.phase_change_upper_boundary_temperature - unit.phase_change_lower_boundary_temperature

    # set fluid start temperature
    unit.fluid_temperature = unit.average_temperature_adjacent_to_pipe

    # calculate simulation mesh
    if unit.accuracy_mode == "very_rough"
        min_mesh_width = unit.pipe_d_o
        max_mesh_width = unit.pipe_d_o * 256
        expansion_factor = 2.0
    elseif unit.accuracy_mode == "rough"
        min_mesh_width = unit.pipe_d_o / 2
        max_mesh_width = unit.pipe_d_o * 128
        expansion_factor = 2.0
    elseif unit.accuracy_mode == "normal"
        min_mesh_width = unit.pipe_d_o / 4
        max_mesh_width = unit.pipe_d_o * 64
        expansion_factor = 2.0
    elseif unit.accuracy_mode == "high"
        min_mesh_width = unit.pipe_d_o / 8
        max_mesh_width = unit.pipe_d_o * 32
        expansion_factor = 2.0
    elseif unit.accuracy_mode == "very_high"
        min_mesh_width = unit.pipe_d_o / 16
        max_mesh_width = unit.pipe_d_o * 16
        expansion_factor = 2.0
    else
        @error "In geothermal collector $(unit.uac), the accuracy_mode has to be one of: very_rough, rough, " *
               "normal, high or very_high"
        throw(InputError())
    end

    # dy_mesh holds the delta between the nodes, while dy is the y-width assigned to each node
    unit.dy, y_pipe_node_num = create_mesh_y(min_mesh_width,
                                             max_mesh_width,
                                             expansion_factor,
                                             unit.pipe_laying_depth,
                                             unit.pipe_radius_outer,
                                             unit.considered_soil_depth)

    for y in 1:(length(unit.dy) - 1)
        push!(unit.dy_mesh, (unit.dy[y] + unit.dy[y + 1]) / 2)
    end

    # dx_mesh holds the delta between the nodes, while dx is the x-width assigned to each node 
    unit.dx = create_mesh_x(min_mesh_width,
                            max_mesh_width,
                            expansion_factor,
                            unit.pipe_radius_outer,
                            unit.pipe_spacing)

    for x in 1:(length(unit.dx) - 1)
        push!(unit.dx_mesh, (unit.dx[x] + unit.dx[x + 1]) / 2)
    end

    unit.dz = copy(unit.pipe_length)

    n_nodes_x = length(unit.dx)
    n_nodes_y = length(unit.dy)

    # localize fluid and adjacent nodes as they will be calculated separately later
    unit.fluid_node_y_idx = y_pipe_node_num

    # set pipe-surrounding nodes assuming that the pipe has one central
    # and 8 surrounding nodes, of which 5 are considered in the axisymmetric grid used here.
    unit.is_pipe_surrounding = fill(false, n_nodes_y, n_nodes_x)  # pipe-surrounding nodes are true, all others are false
    unit.is_pipe_surrounding[y_pipe_node_num - 1, 1] = true
    unit.is_pipe_surrounding[y_pipe_node_num - 1, 2] = true
    unit.is_pipe_surrounding[y_pipe_node_num + 0, 2] = true
    unit.is_pipe_surrounding[y_pipe_node_num + 1, 1] = true
    unit.is_pipe_surrounding[y_pipe_node_num + 1, 2] = true

    # set starting temperature distribution as follows:
    # top:         average ambient temperature, linearly interpolated to pipe surroundings
    # pipe:        6 nodes with the same starting temperature for pipe and pipe surroundings
    # lower bound: undisturbed ground temperature
    unit.t1 = fill(0.0, n_nodes_y, n_nodes_x)
    n_nodes_surface_to_pipe_surrounding = y_pipe_node_num - 3
    n_nodes_pipe_surrounding_to_ground = n_nodes_y - n_nodes_surface_to_pipe_surrounding - 5

    # set temperatures from surface to pipe surroundings with linear interpolation
    if unit.constant_ambient_temperature !== nothing
        average_ambient_temperature = unit.constant_ambient_temperature
    else
        average_ambient_temperature = sum(values(unit.ambient_temperature_profile.data)) /
                                      length(values(unit.ambient_temperature_profile.data))
    end
    step = (unit.average_temperature_adjacent_to_pipe - average_ambient_temperature) /
           sum(unit.dy_mesh[1:n_nodes_surface_to_pipe_surrounding])
    for i in 1:n_nodes_surface_to_pipe_surrounding
        unit.t1[i, :] .= average_ambient_temperature + sum(unit.dy_mesh[1:(i - 1)]) * step
    end

    # set temperatures for pipe and pipe surroundings
    unit.t1[(n_nodes_surface_to_pipe_surrounding + 1):(n_nodes_surface_to_pipe_surrounding + 5), :] .= unit.average_temperature_adjacent_to_pipe

    # set temperatures from pipe surroundings to ground with linear interpolation
    step = (unit.undisturbed_ground_temperature - unit.average_temperature_adjacent_to_pipe) /
           sum(unit.dy_mesh[(n_nodes_surface_to_pipe_surrounding + 5):(n_nodes_y - 1)])
    for i in 1:n_nodes_pipe_surrounding_to_ground
        dy_idx = i + n_nodes_surface_to_pipe_surrounding + 5
        unit.t1[dy_idx, :] = unit.t1[dy_idx - 1, :] .+ unit.dy_mesh[dy_idx - 1] * step
    end

    unit.t2 = copy(unit.t1)

    # calculate volume around the pipe
    unit.volume_adjacent_to_pipe = ((unit.dx[1] + unit.dx[2]) *                                             # area adjacent to pipe in x-direction
                                    (unit.dy[unit.fluid_node_y_idx - 1] + unit.dy[unit.fluid_node_y_idx] +  # area adjacent to pipe in y-direction
                                     unit.dy[unit.fluid_node_y_idx + 1]) -
                                    (0.5 * pi * unit.pipe_radius_outer^2)) * unit.dz                        # 1/2 cross section of pipe

    # set soil density. currently only homogenous soil is considered, but more layers are possible
    unit.soil_density_vector = fill(unit.soil_density, n_nodes_y)

    # calculate soil weight of the soil around each node
    unit.soil_weight = fill(0.0, n_nodes_y, n_nodes_x)
    for h in 1:n_nodes_y
        for i in 1:n_nodes_x
            unit.soil_weight[h, i] = unit.soil_density_vector[h] * unit.dz * unit.dx[i] * unit.dy[h]
        end
    end

    # specific heat capacity for each node needed, because of apparent heat capacity method.
    unit.cp = fill(unit.soil_specific_heat_capacity, n_nodes_y, n_nodes_x)

    # vector to hold the results of the temperatures for each node in each simulation time step
    unit.temp_field_output = zeros(Float64, sim_params["number_of_time_steps_output"], n_nodes_y, n_nodes_x)
end

""" 
    create_mesh_y(min_mesh_width, max_mesh_width, expansion_factor, pipe_laying_depth,
                  pipe_radius_outer, total_depth_simulation_domain)

Creates a non-uniform mesh for a geothermal heat collector in y-direction (orthogonal to
ground surface and orthogonal to pipe).

| -------- (ground surface)
| O        (pipe cross section)
V ........ (lower bound)

The arrow shows the y-direction.

The mesh starts at the ground surface with an initial width of `min_mesh_width`. The width
expands by a factor of `expansion_factor` for each node until reaching halfway to the
pipe-laying depth. At that point, the mesh width begins to decrease. The pipe itself is
represented by one node, surrounded by one node each in the positive and negative
y-directions. Below the pipe area, the mesh width starts again at `min_mesh_width`,
increases up to halfway to the lower bound, and then decreases symmetrically.

The function determines the width of each volume element around a node point, starting with
the surface layer. The node itself is located in the center of the respective width.

Note: The node spacing at the turning point between increasing and decreasing mesh width may
deviate from the value calculated using the expansion factor.
"""
function create_mesh_y(min_mesh_width::Float64,
                       max_mesh_width::Float64,
                       expansion_factor::Float64,
                       pipe_laying_depth::Float64,
                       pipe_radius_outer::Float64,
                       total_depth_simulation_domain::Float64)
    function calculate_increasing_decreasing_distances(midpoint::Float64,
                                                       min_mesh_width::Float64,
                                                       max_mesh_width::Float64,
                                                       expansion_factor::Float64)
        distances = []
        current_width = min_mesh_width

        while distances == [] || sum(distances) + current_width < midpoint
            push!(distances, current_width)
            current_width = min(max_mesh_width, current_width * expansion_factor)
        end
        distance_to_midpoint = midpoint - sum(distances)
        if distance_to_midpoint > distances[end]
            push!(distances, distance_to_midpoint)
            push!(distances, distance_to_midpoint)
            mirrored_values = reverse(distances[1:(end - 2)])
        else
            push!(distances, 2 * distance_to_midpoint)
            mirrored_values = reverse(distances[1:(end - 1)])
        end
        append!(distances, mirrored_values)

        return distances
    end

    # define segments of the computing grid in y-direction 
    sy1 = 0                                         # surface
    sy2 = pipe_laying_depth - pipe_radius_outer * 3 / 2     # node above fluid node
    sy3 = pipe_laying_depth + pipe_radius_outer * 3 / 2     # node below fluid node
    sy4 = total_depth_simulation_domain             # lower simulation boundary

    # segment 1: surface to pipe 
    midpoint = (sy2 - sy1) / 2
    dy_1 = calculate_increasing_decreasing_distances(midpoint, min_mesh_width, max_mesh_width, expansion_factor)

    # detect node number of pipe 
    y_pipe_node_num = length(dy_1) + 2

    # segment 2: pipe
    dy_2 = [pipe_radius_outer, pipe_radius_outer, pipe_radius_outer]

    # segment 3: pipe to lower boundary 
    midpoint = (sy4 - sy3) / 2
    dy_3 = calculate_increasing_decreasing_distances(midpoint, min_mesh_width, max_mesh_width, expansion_factor)

    return [dy_1..., dy_2..., dy_3...], y_pipe_node_num
end

""" 
    create_mesh_x(min_mesh_width, max_mesh_width, expansion_factor, pipe_radius_outer, pipe_spacing)

Creates a non-uniform mesh for a geothermal collector in x-direction (parallel to ground
surface and orthogonal to pipe).

-----------------------  (ground surface)
O ------->               (pipe cross section) 
.......................  (lower bound)

The arrow shows the x-direction.

Starts with half the pipe diameter to represent the pipe nodes and increases the mesh width
by the expansion factor for each node until half the pipe spacing is reached.

The function determines the width of each volume element around a node point, starting with
the vertical mirror axis in which the pipe lies. The node itself is located in the center of
the respective width, except of the first node which is located right on the mirror axis.

Note that the last dx can be smaller than calculated to meet the given boundary by
`pipe_spacing`.
"""
function create_mesh_x(min_mesh_width::Float64,
                       max_mesh_width::Float64,
                       expansion_factor::Float64,
                       pipe_radius_outer::Float64,
                       pipe_spacing::Float64)
    # maximum width of grid in x direction
    x_bound = pipe_spacing / 2     # [m]

    # set first two dx to half the pipe radius
    dx = [pipe_radius_outer / 2, pipe_radius_outer]       # [m]

    # set next dx to the minimum width
    append!(dx, min_mesh_width)

    # add dx while x_bound is not exeeded
    while sum(dx) < x_bound
        append!(dx, min(dx[end] * expansion_factor, max_mesh_width))
    end

    # limit to x_bound
    if sum(dx) > x_bound
        dx[end] -= sum(dx) - x_bound
    end

    return dx
end

function control(unit::GeothermalHeatCollector,
                 components::Grouping,
                 sim_params::Dict{String,Any})
    update(unit.controller)

    # reset energy summarizer
    unit.collector_total_heat_energy_in_out = 0.0

    # get ambient temperature and global radiation from profile for current time step if needed
    if unit.constant_ambient_temperature !== nothing
        unit.ambient_temperature = unit.constant_ambient_temperature
    else
        unit.ambient_temperature = Profiles.value_at_time(unit.ambient_temperature_profile, sim_params)
    end
    if unit.constant_global_radiation !== nothing
        unit.global_radiation_power = unit.constant_global_radiation  # W/m^2
    else
        unit.global_radiation_power = Profiles.power_at_time(unit.global_radiation_profile, sim_params) # W/m^2
    end

    unit.current_output_temperature = unit.fluid_temperature + unit.unloading_temperature_spread / 2
    unit.current_output_temperature = highest(unit.fluid_min_output_temperature, unit.current_output_temperature)

    # limit max_energy if current_output_temperature undercuts fluid_min_output_temperature
    if unit.fluid_min_output_temperature === nothing
        unit.current_max_output_energy = unit.max_output_energy
    else
        unit.current_max_output_energy = unit.current_output_temperature <= unit.fluid_min_output_temperature ?
                                         0.0 : unit.max_output_energy
    end
    set_max_energy!(unit.output_interfaces[unit.m_heat_out], unit.current_max_output_energy, nothing,
                    unit.current_output_temperature)

    # get input temperature for energy input (regeneration) and set temperature and max_energy to input interface
    if unit.regeneration
        unit.current_input_temperature = unit.fluid_temperature - unit.loading_temperature_spread / 2
        unit.current_input_temperature = lowest(unit.fluid_max_input_temperature, unit.current_input_temperature)

        # limit max_energy if current_input_temperature exceeds fluid_max_input_temperature
        if unit.fluid_max_input_temperature === nothing
            unit.current_max_input_energy = unit.max_input_energy
        else
            unit.current_max_input_energy = unit.current_input_temperature >= unit.fluid_max_input_temperature ?
                                            0.0 : unit.max_input_energy
        end
        set_max_energy!(unit.input_interfaces[unit.m_heat_in], unit.current_max_input_energy,
                        unit.current_input_temperature, nothing)
    end
end

function calculate_new_temperature_field!(unit::GeothermalHeatCollector, q_in_out::Float64, sim_params)
    # Calculate temperature field using implicit Euler
    # --- helpers ---
    function idx(nx::Int, h::Int, i::Int)
        return (h - 1) * nx + i
    end

    function liquid_fraction(unit::GeothermalHeatCollector, T::Float64)
        T0 = unit.phase_change_lower_boundary_temperature
        T1 = unit.phase_change_upper_boundary_temperature
        if T <= T0
            return 0.0
        elseif T >= T1
            return 1.0
        else
            s = (T - T0) / (T1 - T0)   # in (0,1)
            return s * s * (3.0 - 2.0 * s)   # smoothstep
        end
    end

    function dliquid_dT(unit::GeothermalHeatCollector, T::Float64)
        T0 = unit.phase_change_lower_boundary_temperature
        T1 = unit.phase_change_upper_boundary_temperature
        dT = T1 - T0
        if T <= T0 || T >= T1
            return 0.0
        else
            s = (T - T0) / dT
            return (6.0 * s * (1.0 - s)) / dT    # d/dT smoothstep
        end
    end

    function c_sensible(unit::GeothermalHeatCollector, T::Float64)
        # linear blend of sensible cp across the interval
        T0 = unit.phase_change_lower_boundary_temperature
        T1 = unit.phase_change_upper_boundary_temperature
        if T <= T0
            return unit.soil_specific_heat_capacity_frozen
        elseif T >= T1
            return unit.soil_specific_heat_capacity
        else
            w = (T - T0) / (T1 - T0)
            return (1.0 - w) * unit.soil_specific_heat_capacity_frozen + w * unit.soil_specific_heat_capacity
        end
    end

    # evaluate conductivity at given temperature
    function k_from_T(unit::GeothermalHeatCollector, T::Float64)
        T0 = unit.phase_change_lower_boundary_temperature
        T1 = unit.phase_change_upper_boundary_temperature
        if T >= T1
            return unit.soil_heat_conductivity
        elseif T <= T0
            return unit.soil_heat_conductivity_frozen
        else
            f = (T1 - T) / (T1 - T0)
            return unit.soil_heat_conductivity * (1.0 - f) + unit.soil_heat_conductivity_frozen * f
        end
    end

    function k_between(unit::GeothermalHeatCollector, Tmat::AbstractMatrix, h1::Int, i1::Int, h2::Int, i2::Int)
        return (k_from_T(unit, Tmat[h1, i1]) + k_from_T(unit, Tmat[h2, i2])) / 2.0
    end

    # write output
    if sim_params["current_date"] >= sim_params["start_date_output"]
        unit.temp_field_output[Int(sim_params["time_since_output"] / sim_params["time_step_seconds"]) + 1, :, :] = copy(unit.t1)
    end

    # geometry
    nx = length(unit.dx)
    ny = length(unit.dy)
    N = nx * ny

    # ambient & radiation inputs
    ambient_T = unit.constant_ambient_temperature === nothing ?
                Profiles.value_at_time(unit.ambient_temperature_profile, sim_params) :
                unit.constant_ambient_temperature

    global_rad = unit.constant_global_radiation === nothing ?
                 Profiles.power_at_time(unit.global_radiation_profile, sim_params) :
                 unit.constant_global_radiation

    infrared_sky_radiation = unit.constant_infrared_sky_radiation === nothing ?
                             Profiles.power_at_time(unit.infrared_sky_radiation_profile, sim_params) :
                             unit.constant_infrared_sky_radiation
    sky_temperature = (infrared_sky_radiation / unit.boltzmann_constant)^0.25  # [K]

    # pipe source strength
    specific_heat_flux_pipe = sim_params["wh_to_watts"](q_in_out) / (unit.pipe_length * unit.number_of_pipes) # [W/m]
    q_in_out_surrounding = specific_heat_flux_pipe * unit.pipe_length                                      # [W per pipe]
    hpipe = unit.fluid_node_y_idx

    # initial state
    T_old = copy(unit.t1)

    # Picard iterations (update k(T))
    T_guess = copy(T_old)

    for _ in 1:(unit.max_picard_iter)
        # Assemble A * T^{n+1} = b
        I = Int[]
        J = Int[]
        V = Float64[]
        b = zeros(Float64, N)

        for h in 1:ny
            for i in 1:nx
                p = idx(nx, h, i)

                # Bottom boundary (Dirichlet: undisturbed ground temperature)
                if h == ny
                    push!(I, p)
                    push!(J, p)
                    push!(V, 1.0)
                    b[p] = unit.undisturbed_ground_temperature
                    continue
                end

                # Capacity term (apparent cp)
                c_eff = c_sensible(unit, T_guess[h, i]) +
                        unit.soil_specific_enthalpy_of_fusion * dliquid_dT(unit, T_guess[h, i])
                C = c_eff * unit.soil_weight[h, i]  # [J/K] for the control volume
                aP = C / sim_params["time_step_seconds"]
                rhs = aP * T_old[h, i]

                # Conduction (finite-volume), k at current Picard iterate
                # East
                if i < nx
                    kE = k_between(unit, T_guess, h, i, h, i + 1)
                    aE = unit.dz * unit.dy[h] * kE / unit.dx_mesh[i]
                    aP += aE
                    push!(I, p)
                    push!(J, idx(nx, h, i + 1))
                    push!(V, -aE)
                end
                # West (omit link into the fluid node)
                if i > 1
                    if !(i == 2 && h == hpipe)
                        kW = k_between(unit, T_guess, h, i, h, i - 1)
                        aW = unit.dz * unit.dy[h] * kW / unit.dx_mesh[i - 1]
                        aP += aW
                        push!(I, p)
                        push!(J, idx(nx, h, i - 1))
                        push!(V, -aW)
                    end
                end
                # South (+y)
                if h < ny
                    kS = k_between(unit, T_guess, h, i, h + 1, i)
                    aS = unit.dz * unit.dx[i] * kS / unit.dy_mesh[h]
                    aP += aS
                    push!(I, p)
                    push!(J, idx(nx, h + 1, i))
                    push!(V, -aS)
                end
                # North (-y) (interior only; surface handled below)
                if h > 1
                    kN = k_between(unit, T_guess, h, i, h - 1, i)
                    aN = unit.dz * unit.dx[i] * kN / unit.dy_mesh[h - 1]
                    aP += aN
                    push!(I, p)
                    push!(J, idx(nx, h - 1, i))
                    push!(V, -aN)
                end

                # Surface cell: add convection implicitly; radiation & solar explicitly
                if h == 1
                    A_surf = unit.dz * unit.dx[i]
                    hconv = unit.surface_convective_heat_transfer_coefficient
                    aP += hconv * A_surf
                    rhs += hconv * A_surf * ambient_T

                    # explicit sources
                    rhs += A_surf * (1.0 - unit.surface_reflection_factor) * global_rad
                    rhs += A_surf * unit.surface_emissivity * unit.boltzmann_constant *
                           (sky_temperature^4 - (T_guess[h, i] + 273.15)^4)
                end

                # finalize central coeff & RHS
                push!(I, p)
                push!(J, p)
                push!(V, aP)
                b[p] += rhs
            end
        end

        # Pipe source distribution (5 surrounding nodes)
        b[idx(nx, hpipe - 1, 1)] += (1.0 / 16.0) * q_in_out_surrounding
        b[idx(nx, hpipe + 1, 1)] += (1.0 / 16.0) * q_in_out_surrounding
        b[idx(nx, hpipe, 2)] += (1.0 / 8.0) * q_in_out_surrounding
        b[idx(nx, hpipe - 1, 2)] += (1.0 / 8.0) * q_in_out_surrounding
        b[idx(nx, hpipe + 1, 2)] += (1.0 / 8.0) * q_in_out_surrounding

        # Solve sparse linear system
        A = sparse(I, J, V, N, N)
        T_vec = A \ b

        # Map vector back to (ny, nx)
        T_new = similar(unit.t1)
        for hh in 1:ny
            for ii in 1:nx
                T_new[hh, ii] = T_vec[idx(nx, hh, ii)]
            end
        end

        # Convergence on successive iterates
        max_dT = 0.0
        for hh in 1:ny
            for ii in 1:nx
                max_dT = max(max_dT, abs(T_guess[hh, ii] - T_new[hh, ii]))
            end
        end

        if max_dT < unit.picard_tol
            break
        else
            # light damping helps if many cells cross the phase band in one go
            T_guess .= 0.5 .* T_guess .+ 0.5 .* T_new
        end
    end

    # Commit field
    unit.t2 .= T_guess

    # Pipe-adjacent average
    avg_adj = (unit.t2[hpipe - 1, 1] +
               unit.t2[hpipe + 1, 1] +
               2.0 * unit.t2[hpipe, 2] +
               2.0 * unit.t2[hpipe + 1, 2] +
               2.0 * unit.t2[hpipe - 1, 2]) / 8.0
    unit.average_temperature_adjacent_to_pipe = avg_adj

    # Fluid temperature using pipe/soil resistance (compute with final local k)
    if unit.model_type == "simplified"
        pipe_R_len = unit.pipe_soil_thermal_resistance
    else
        unit.alpha_fluid_pipe, unit.fluid_reynolds_number = calculate_alpha_pipe(unit, q_in_out,
                                                                                 sim_params["wh_to_watts"])
        k_loc = k_from_T(unit, unit.t2[hpipe, 2])
        k_eff = 1.0 / (unit.pipe_d_o / (unit.alpha_fluid_pipe * unit.pipe_d_i) +
                       (log(unit.pipe_d_o / unit.pipe_d_i) * unit.pipe_d_o) / (2.0 * unit.pipe_heat_conductivity) +
                       unit.dx_mesh[2] / (2.0 * k_loc))
        pipe_R_len = 1.0 / (k_eff * pi * unit.pipe_d_o)
    end

    unit.fluid_temperature = unit.average_temperature_adjacent_to_pipe +
                             pipe_R_len * specific_heat_flux_pipe
    unit.t2[hpipe, 1] = unit.fluid_temperature
    unit.t2[hpipe - 1, 1] = avg_adj
    unit.t2[hpipe + 1, 1] = avg_adj
    unit.t2[hpipe, 2] = avg_adj
    unit.t2[hpipe - 1, 2] = avg_adj
    unit.t2[hpipe + 1, 2] = avg_adj

    # advance
    unit.t1 .= unit.t2
end

function plot_optional_figures_begin(unit::GeothermalHeatCollector,
                                     output_path::String,
                                     output_formats::Vector{String},
                                     sim_params::Dict{String,Any})
    # plot mesh of geothermal collector
    x_positions = cumsum([0; unit.dx])
    y_positions = cumsum([0; unit.dy]) .* (-1)

    plt = plot(;
               title="Mesh of the collector at \"$(unit.accuracy_mode)\"\n" *
                     "accuracy (dx_min = $(round(minimum(unit.dx)*100;digits=1)) mm) ",
               xlabel="horizontal dimension [m]",
               ylabel="vertical dimension [m]",
               legend=false,
               linewidth=6,
               gridlinewidth=1,
               size=(900, 1200),
               titlefontsize=28,
               guidefontsize=24,
               tickfontsize=24,
               legendfontsize=24,
               margin=15Plots.mm,
               aspect_ratio=:equal)

    # Draw vertical lines
    for x in x_positions
        Plots.plot!(plt, [x, x], [y_positions[1], y_positions[end]]; color=:black, linewidth=1)
    end

    # Draw horizontal lines
    for y in y_positions
        Plots.plot!(plt, [x_positions[1], x_positions[end]], [y, y]; color=:black, linewidth=1)
    end

    fig_name = "collector_simulation_mesh_$(unit.uac)"
    for output_format in output_formats
        savefig(output_path * "/" * fig_name * "." * output_format)
    end

    return true
end

function plot_optional_figures_end(unit::GeothermalHeatCollector, sim_params::Dict{String,Any}, output_path::String)
    # plot temperature field as 3D mesh with time-slider
    @info "Plotting time-shiftable temperature distribution of geothermal collector $(unit.uac). " *
          "Close figure to continue..."

    f = Figure()
    ax = Axis3(f[1, 1])

    ax.zlabel = "Temperature [°C]"
    ax.xlabel = "Vertical expansion (depth) [m]"
    ax.ylabel = "Horizontal expansion [m]"
    min_temp = minimum(unit.temp_field_output)
    min_temp = min_temp < 0.0 ? 1.1 * min_temp : 0.9 * min_temp
    max_temp = maximum(unit.temp_field_output)
    max_temp = max_temp < 0.0 ? 0.9 * max_temp : 1.1 * max_temp
    Makie.zlims!(ax, min_temp, max_temp)

    x_abs = [0; cumsum(unit.dx_mesh)]               # Absolute x coordinates
    y_abs = [unit.dy[1] / 2; cumsum(unit.dy_mesh)]  # Absolute y coordinates

    ## activate for equal axis ratio
    # xlims!(ax, 0, max(x_abs[end], y_abs[end]))
    # ylims!(ax, 0, max(x_abs[end], y_abs[end]))

    time = Observable(1)
    surfdata = @lift(unit.temp_field_output[$time, :, :])
    GLMakie.surface!(ax, y_abs, x_abs, surfdata)
    GLMakie.scatter!(ax, y_abs, x_abs, surfdata)
    slg = SliderGrid(f[2, 1], (; range=1:1:sim_params["number_of_time_steps_output"], label="Time"))

    on(slg.sliders[1].value) do v
        time[] = v
    end
    wait(display(f))

    return false
end

# function to calculate heat transfer coefficient alpha.
function calculate_alpha_pipe(unit::GeothermalHeatCollector, q_in_out::Float64, wh2w::Function)

    # calculate mass flow in pipe
    collector_power_in_out_per_pipe = wh2w(abs(q_in_out)) / unit.number_of_pipes  # W/pipe
    temperature_spread = q_in_out > 0 ? unit.loading_temperature_spread : unit.unloading_temperature_spread
    collector_mass_flow_per_pipe = collector_power_in_out_per_pipe /
                                   (unit.fluid_specific_heat_capacity * temperature_spread)  # kg/s

    if unit.use_dynamic_fluid_properties
        # calculate reynolds-number based on dynamic viscosity using dynamic temperature-dependend fluid properties, 
        # adapted from TRNSYS Type 710, for 30 Vol-% ethylene glycol mix:
        fluid_dynamic_viscosity = 0.0000017158 * unit.fluid_temperature^2 -
                                  0.0001579079 * unit.fluid_temperature + 0.0048830621
        unit.fluid_heat_conductivity = 0.0010214286 * unit.fluid_temperature + 0.447
        unit.fluid_prandtl_number = fluid_dynamic_viscosity * unit.fluid_specific_heat_capacity /
                                    unit.fluid_heat_conductivity
        fluid_reynolds_number = (4 * collector_mass_flow_per_pipe) / (pi * unit.pipe_d_i * fluid_dynamic_viscosity)
    else
        # calculate reynolds-number, based on kinematic viscosity with constant fluid properties.
        fluid_reynolds_number = (4 * collector_mass_flow_per_pipe) /
                                (unit.fluid_density * unit.fluid_kinematic_viscosity * unit.pipe_d_i * pi)
    end

    if fluid_reynolds_number <= 2300  # laminar
        nusselt = calculate_Nu_laminar(unit, fluid_reynolds_number)
    elseif fluid_reynolds_number > 2300 && fluid_reynolds_number <= 1e4 # transitional
        # Gielinski 1995
        factor = (fluid_reynolds_number - 2300) / (1e4 - 2300)
        nusselt = (1 - factor) * calculate_Nu_laminar(unit, 2300.0) +
                  factor * calculate_Nu_turbulent(unit, 10_000.0)
    else  # turbulent
        nusselt = calculate_Nu_turbulent(unit, fluid_reynolds_number)
    end

    alpha = nusselt * unit.fluid_heat_conductivity / unit.pipe_d_i

    return alpha, fluid_reynolds_number
end

function calculate_Nu_laminar(unit::GeothermalHeatCollector, fluid_reynolds_number::Float64)
    if unit.nusselt_approach == "Ramming"
        # Approach used in Ramming 2007 from Elsner, Norbert; Fischer, Siegfried; Huhn, Jörg; "Grundlagen der
        # Technischen Thermodynamik",  Band 2 Wärmeübertragung, Akademie Verlag, Berlin 1993. 
        k_a = 1.1 - 1 / (3.4 + 0.0667 * unit.fluid_prandtl_number)
        k_n = 0.35 + 1 / (7.825 + 2.6 * sqrt(unit.fluid_prandtl_number))

        # calculate Nu-Number
        nusselt_laminar = ((k_a / (1 - k_n) *
                            (unit.fluid_prandtl_number * unit.pipe_d_i * fluid_reynolds_number / unit.pipe_length)^k_n)^3 +
                           4.364^3)^(1 / 3)
    elseif unit.nusselt_approach == "Stephan"
        # Stephan
        pr_water = 13.44                # Pr Number Water 0 °C as reference
        nusselt_laminar = 3.66 +
                          (0.0677 *
                           (fluid_reynolds_number * unit.fluid_prandtl_number * unit.pipe_d_i / unit.pipe_length)^1.33) /
                          (1 +
                           0.1 * unit.fluid_prandtl_number *
                           (fluid_reynolds_number * unit.pipe_d_i / unit.pipe_length)^0.83) *
                          (unit.fluid_prandtl_number / pr_water)^0.11
    else
        @error "In geothermal collector $(unit.uac), the nusselt_approach has to be one of: Ramming, Stephan."
        throw(InputError())
    end
    return nusselt_laminar
end

function calculate_Nu_turbulent(unit::GeothermalHeatCollector, fluid_reynolds_number::Float64)
    # Approach used from Gnielinski in: V. Gnielinski: Ein neues Berechnungsverfahren für die Wärmeübertragung 
    # im Übergangsbereich zwischen laminarer und turbulenter Rohrströmung. Forsch im Ing Wes 61:240-248, 1995. 
    zeta = (1.8 * log(fluid_reynolds_number) - 1.5)^-2
    nusselt_turbulent = (zeta / 8 * fluid_reynolds_number * unit.fluid_prandtl_number) /
                        (1 + 12.7 * sqrt(zeta / 8) * (unit.fluid_prandtl_number^(2 / 3) - 1))
    return nusselt_turbulent
end

# process function that provides energy from the geothermal heat collector and calculates new temperatures 
# according to actual delivered or received energy
function process(unit::GeothermalHeatCollector, sim_params::Dict{String,Any})
    # get actual required energy from output interface
    outface = unit.output_interfaces[unit.m_heat_out]
    exchanges = balance_on(outface, outface.target)
    energy_demanded = balance(exchanges) + energy_potential(exchanges)
    energy_available = unit.current_max_output_energy  # is positive

    # shortcut if there is no energy demanded
    if energy_demanded >= -sim_params["epsilon"]
        handle_component_update!(unit, "process", sim_params)
        set_max_energy!(unit.output_interfaces[unit.m_heat_out], 0.0)
        return
    end

    for exchange in exchanges
        demanded_on_interface = exchange.balance + exchange.energy_potential

        if demanded_on_interface >= -sim_params["epsilon"]
            continue
        end

        if (exchange.temperature_min !== nothing &&
            exchange.temperature_min > unit.current_output_temperature)
            # we can only supply energy at a temperature at or below the collector's current
            # output temperature
            continue
        end

        used_heat = abs(demanded_on_interface)

        if energy_available > used_heat
            energy_available -= used_heat
            add!(outface, used_heat, nothing, unit.current_output_temperature)
        else
            add!(outface, energy_available, nothing, unit.current_output_temperature)
            energy_available = 0.0
        end
    end

    # write output heat flux into vector
    energy_delivered = -(unit.current_max_output_energy - energy_available)
    unit.collector_total_heat_energy_in_out = energy_delivered
    handle_component_update!(unit, "process", sim_params)
end

function handle_component_update!(unit::GeothermalHeatCollector, step::String, sim_params::Dict{String,Any})
    if step == "process"
        unit.process_done = true
    elseif step == "load"
        unit.load_done = true
    end
    if unit.process_done && unit.load_done
        # calculate new temperatures of field to account for possible ambient effects
        calculate_new_temperature_field!(unit, unit.collector_total_heat_energy_in_out, sim_params)
        # reset 
        unit.process_done = false
        unit.load_done = false
    end
end

function load(unit::GeothermalHeatCollector, sim_params::Dict{String,Any})
    if !unit.regeneration
        handle_component_update!(unit, "load", sim_params)
        return
    end

    inface = unit.input_interfaces[unit.m_heat_in]
    exchanges = balance_on(inface, inface.source)
    energy_available = balance(exchanges) + energy_potential(exchanges)
    energy_demand = unit.current_max_input_energy  # is positive

    if energy_available <= sim_params["epsilon"]
        # shortcut if there is no energy to be used
        handle_component_update!(unit, "load", sim_params)
        set_max_energy!(unit.input_interfaces[unit.m_heat_in], 0.0)
        return
    end

    for exchange in exchanges
        exchange_energy_available = exchange.balance + exchange.energy_potential

        if exchange_energy_available < sim_params["epsilon"]
            continue
        end

        if exchange.temperature_max !== nothing &&
           exchange.temperature_max < unit.current_input_temperature
            # we can only take energy if it's at a higher/equal temperature than the
            # collector's current input temperature
            continue
        end

        if energy_demand > exchange_energy_available
            energy_demand -= exchange_energy_available
            sub!(inface, exchange_energy_available, unit.current_input_temperature, nothing)
        else
            sub!(inface, energy_demand, unit.current_input_temperature, nothing)
            energy_demand = 0.0
        end
    end

    energy_taken = unit.current_max_input_energy - energy_demand
    unit.collector_total_heat_energy_in_out += energy_taken
    handle_component_update!(unit, "load", sim_params)
end

function output_values(unit::GeothermalHeatCollector)::Vector{String}
    output_vals = []
    if unit.regeneration
        push!(output_vals, string(unit.m_heat_in) * ":IN")
    end
    if unit.model_type == "detailed"
        push!(output_vals, "fluid_reynolds_number")
        push!(output_vals, "alpha_fluid_pipe")
    end
    append!(output_vals,
            [string(unit.m_heat_out) * ":OUT",
             "fluid_temperature",
             "ambient_temperature",
             "global_radiation_power"])
    # push!(output_vals, "TEMPERATURE_xNodeNum_yNodeNum")

    return output_vals
end

function output_value(unit::GeothermalHeatCollector, key::OutputKey)::Float64
    if key.value_key == "IN"
        return calculate_energy_flow(unit.input_interfaces[key.medium])
    elseif key.value_key == "OUT"
        return calculate_energy_flow(unit.output_interfaces[key.medium])
    elseif startswith(key.value_key, "Temperature_")
        splitted = split(key.value_key, "_")
        x_idx = parse(Int, splitted[2])
        y_idx = parse(Int, splitted[3])
        if !(1 <= y_idx <= length(unit.dy) + 1) || !(1 <= x_idx <= length(unit.dx) + 1)
            throw(ArgumentError("The indexes ($(x_idx), $(y_idx)) of the requested temperature-output of the " *
                                "geothermal collector $(unit.uac) exeed the number of nodes of the mesh. The maximum " *
                                "is ($(length(unit.dx)), $(length(unit.dy)))."))
        else
            return unit.t2[x_idx, y_idx]
        end
    elseif key.value_key == "fluid_temperature"
        return unit.fluid_temperature
    elseif key.value_key == "fluid_reynolds_number"
        return unit.fluid_reynolds_number
    elseif key.value_key == "ambient_temperature"
        return unit.ambient_temperature
    elseif key.value_key == "global_radiation_power"
        return unit.global_radiation_power
    elseif key.value_key == "alpha_fluid_pipe"
        return unit.alpha_fluid_pipe
    end
    throw(KeyError(key.value_key))
end

export GeothermalHeatCollector
