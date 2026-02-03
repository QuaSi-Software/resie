#! format: off
const FIXED_SINK_PARAMETERS = Dict(
    "medium" => (
        description="Medium of the sink (e.g. electricity, heat, gas, etc.)",
        display_name="Medium",
        required=true,
        type=String,
        json_type="string",
        unit="-"
    ),
    "temperature_profile_file_path" => (
        default=nothing,
        description="Path to a temperature profile file",
        display_name="Temperature profile file",
        required=false,
        conditionals=[
            ("temperature_from_global_file", "mutex"),
            ("constant_temperature", "mutex")
        ],
        type=String,
        json_type="string",
        unit="-"
    ),
    "temperature_from_global_file" => (
        default=nothing,
        description="If given points to a key in the global weather data file with the " *
                    "temperature profile to be used",
        display_name="Global file temp. key",
        required=false,
        conditionals=[
            ("temperature_profile_file_path", "mutex"),
            ("constant_temperature", "mutex")
        ],
        type=String,
        json_type="string",
        unit="-"
    ),
    "constant_temperature" => (
        default=nothing,
        description="Constant temperature value",
        display_name="Constant temperature",
        required=false,
        conditionals=[
            ("temperature_profile_file_path", "mutex"),
            ("temperature_from_global_file", "mutex")
        ],
        type=Float64,
        json_type="number",
        unit="°C"
    ),
    "energy_profile_file_path" => (
        default=nothing,
        description="Path to a profile file with energy values",
        display_name="Energy profile file",
        required=false,
        conditionals=[
            ("constant_demand", "mutex"),
        ],
        type=String,
        json_type="string",
        unit="-"
    ),
    "constant_demand" => (
        default=nothing,
        description="Constant demand (power, not work)",
        display_name="Constant demand",
        required=false,
        conditionals=[
            ("energy_profile_file_path", "mutex"),
        ],
        type=Float64,
        json_type="number",
        unit="W"
    ),
    "scale" => (
        default=1.0,
        description="Scaling factor for the energy profile",
        display_name="Energy scale",
        required=false,
        conditionals=[("energy_profile_file_path", "is_not_nothing")],
        type=Float64,
        json_type="number",
        unit="W"
    ),
)
#! format: on

"""
Implementation of a component modeling a generic fixed sink of a chosen medium.

This is particularly useful for testing, but can also be used to model any fixed
component or other equipment unit that consumes energy of a given medium. This amount of
energy ought to be provided by other components in the energy system.
Note that "fixed" in this context means that the amount of energy the unit consumes is
fixed within a timestep, but can vary over multiple timesteps. No calculation other than
scaling of profile values is performed in each timestep.
"""
mutable struct FixedSink <: Component
    uac::String
    controller::Controller
    sys_function::SystemFunction
    medium::Symbol

    input_interfaces::InterfaceMap
    output_interfaces::InterfaceMap

    energy_profile::Union{Profile,Nothing}
    temperature_profile::Union{Profile,Nothing}
    scaling_factor::Float64

    demand::Float64
    temperature::Temperature

    constant_demand::Union{Nothing,Float64}
    constant_temperature::Temperature

    function FixedSink(uac::String, config::Dict{String,Any}, sim_params::Dict{String,Any})
        return new(SSOT_parameter_constructor(FixedSink, uac, config, sim_params)...)
    end
end

"""
A component that models the demand consumers in a building require.

As the simulation does not encompass demand calculations, this is usually taken from other
tools that calculate the demand before an energy system simulation is performed. These
profiles usually are normalized to some degree, therefore Demand instances require a scaling
factor to turn the relative values to absolute values of required energy.

This is an alias to the generic implementation of a fixed sink.
"""
const Demand = FixedSink

function component_parameters(x::Type{FixedSink})::Dict{String,NamedTuple}
    return deepcopy(FIXED_SINK_PARAMETERS) # return a copy to prevent external modification
end

function extract_parameter(x::Type{FixedSink}, config::Dict{String,Any}, param_name::String, param_def::NamedTuple,
                           sim_params::Dict{String,Any}, uac::String)
    if param_name == "temperature_from_global_file"
        return load_profile_from_global_weather_file(config, param_name, sim_params, uac)
    elseif param_name == "temperature_profile_file_path"
        return load_optional_profile(config, param_name, sim_params)
    elseif param_name == "constant_temperature"
        return convert(Temperature, default(config, param_name, nothing))
    end

    return extract_parameter(Component, config, param_name, param_def, sim_params, uac)
end

function validate_config(x::Type{FixedSink}, config::Dict{String,Any}, extracted::Dict{String,Any}, uac::String,
                         sim_params::Dict{String,Any})
    validate_config(Component, extracted, uac, sim_params, component_parameters(FixedSink))
end

function init_from_params(x::Type{FixedSink}, uac::String, params::Dict{String,Any},
                          raw_params::Dict{String,Any}, sim_params::Dict{String,Any})::Tuple
    medium = Symbol(params["medium"])
    register_media([medium])

    energy_profile = params["energy_profile_file_path"] !== nothing ?
                     Profile(params["energy_profile_file_path"], sim_params) :
                     nothing

    # return tuple in the order expected by new()
    return (uac,                                     # uac
            Controller(params["control_parameters"]),
            sf_fixed_sink,                           # sys_function
            medium,                                  # medium
            InterfaceMap(medium => nothing),         # input_interfaces
            InterfaceMap(medium => nothing),         # output_interfaces
            energy_profile,                          # energy_profile
            some_or_none(params["temperature_profile_file_path"], params["temperature_from_global_file"]),
            params["scale"],                         # scaling_factor
            0.0,                                     # demand
            nothing,                                 # temperature
            params["constant_demand"],               # constant_demand (power, not work!)
            params["constant_temperature"])          # constant_temperature
end

function initialise!(unit::FixedSink, sim_params::Dict{String,Any})
    set_storage_transfer!(unit.input_interfaces[unit.medium],
                          unload_storages(unit.controller, unit.medium))
end

function control(unit::FixedSink,
                 components::Grouping,
                 sim_params::Dict{String,Any})
    update(unit.controller)

    if unit.constant_demand !== nothing
        unit.demand = sim_params["watt_to_wh"](unit.constant_demand)
    elseif unit.energy_profile !== nothing
        unit.demand = unit.scaling_factor * Profiles.work_at_time(unit.energy_profile, sim_params)
    else
        unit.demand = 0.0
    end

    if unit.constant_temperature !== nothing
        unit.temperature = unit.constant_temperature
    elseif unit.temperature_profile !== nothing
        unit.temperature = Profiles.value_at_time(unit.temperature_profile, sim_params)
    end
    set_max_energy!(unit.input_interfaces[unit.medium], unit.demand, unit.temperature, nothing)
end

function process(unit::FixedSink, sim_params::Dict{String,Any})
    inface = unit.input_interfaces[unit.medium]
    sub!(inface, unit.demand, unit.temperature, nothing)
end

function output_values(unit::FixedSink)::Vector{String}
    if unit.temperature_profile === nothing && unit.constant_temperature === nothing
        return [string(unit.medium) * ":IN",
                "Demand"]
    else
        return [string(unit.medium) * ":IN",
                "Demand",
                "Temperature"]
    end
end

function output_value(unit::FixedSink, key::OutputKey)::Float64
    if key.value_key == "IN"
        return calculate_energy_flow(unit.input_interfaces[key.medium])
    elseif key.value_key == "Demand"
        return unit.demand
    elseif key.value_key == "Temperature"
        return unit.temperature
    end
    throw(KeyError(key.value_key))
end

export FixedSink, Demand
