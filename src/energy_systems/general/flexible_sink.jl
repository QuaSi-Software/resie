#! format: off
const FLEXIBLE_SINK_COMPONENT_PARAMETERS = Dict(
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
                    "temperature profile to be used. Use `temp_ambient_air` as key.",
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
    "max_power_profile_file_path" => (
        default=nothing,
        description="Path to a profile file with maximum power values",
        display_name="Max. power profile file",
        required=false,
        conditionals=[("constant_power", "mutex")],
        validations=[
            ("at_least_one", "max_power_profile_file_path", "constant_power")
        ],
        type=String,
        json_type="string",
        unit="-"
    ),
    "constant_power" => (
        default=nothing,
        description="Constant maximum power to take in",
        display_name="Constant max. power",
        required=false,
        conditionals=[("max_power_profile_file_path", "mutex")],
        validations=[
            ("self", "value_gte_num_or_nothing", 0.0),
            ("at_least_one", "max_power_profile_file_path", "constant_power")
        ],
        type=Float64,
        json_type="number",
        unit="W"
    ),
    "scale" => (
        default=1.0,
        description="Scaling factor for the max. power profile",
        display_name="Max. power scale",
        required=false,
        conditionals=[("max_power_profile_file_path", "is_not_nothing")],
        type=Float64,
        json_type="number",
        unit="W"
    ),
)

const FLEXIBLE_SINK_ECONOMIC_PARAMETERS = get_economic_standard_params("connection", 
    Dict{String,Any}(
        "energy_price_profile_file_path" => nothing,
        "energy_price_profile_scale" => 1.0,
        "constant_energy_price" => nothing,
        "energy_price_change_rate_per_year" =>  0.02,
        "base_cost_per_year" => 0.0,
        "base_cost_change_rate_per_year" => 0.0,

        "lifetime_years" => 20,
        "capex_specific" => "const:0.0",
        "capex_price_change_rate_per_year" => 0.0,
        "maintenance_inspection_rate_per_year" => 0.0,
        "maintenance_inspection_price_change_rate_per_year" =>  0.0,
        "repair_rate_per_year" => 0.0,
        "repair_price_change_rate_per_year" =>  0.0,
        "operational_labour_hours_per_year" =>  0.0,
        "subsidy_rate_of_capex" => nothing,
        "subsidy_max" => nothing
    ),
    Dict{String,Any}(            
        "capex_specific" => "€/(constant_power or scale)"
    ),
)

const FLEXIBLE_SINK_EMISSIONS_PARAMETERS = get_emissions_standard_params("connection_sink", 
    Dict{String,Any}(
        "energy_emission_credits_profile_file_path" => nothing,
        "energy_emission_credits_profile_scale" => 1.0,
        "constant_energy_emission_credits" => nothing,
        "energy_emission_credits_change_rate_per_year" =>  0.0,
    
        "lifetime_years" => 20,
        "embodied_emissions_specific" => "const:0.0",
        "embodied_emissions_change_rate_per_year" => 0.0
    ),
    Dict{String,Any}(
        "embodied_emissions_specific" => "g CO2/(constant_power or scale)"
    )
)
#! format: on

"""
Implementation of a component modeling a generic flexible sink of a chosen medium.

This is particularly useful for testing, but can also be used to model any flexible
component or other equipment unit that consumes energy of a given medium. The component
might still have a maximum power intake in a single time step, but can consume any fraction
of this.
"""
mutable struct FlexibleSink <: Component
    uac::String
    controller::Controller
    sys_function::SystemFunction
    medium::Symbol

    input_interfaces::InterfaceMap
    output_interfaces::InterfaceMap

    economic_parameters::Dict{String,Any}
    emissions_parameters::Dict{String,Any}

    max_power_profile::Union{Profile,Nothing}
    temperature_profile::Union{Profile,Nothing}
    scaling_factor::Float64

    max_energy::Float64
    temperature::Temperature
    constant_power::Union{Nothing,Float64}
    constant_temperature::Temperature

    function FlexibleSink(uac::String, config::Dict{String,Any}, sim_params::Dict{String,Any})
        return new(SSOT_parameter_constructor(FlexibleSink, uac, config, sim_params)...)
    end
end

function component_parameters(x::Type{FlexibleSink})::Dict{String,Any}
    return deepcopy(FLEXIBLE_SINK_COMPONENT_PARAMETERS)
end

function economic_parameters(x::Type{FlexibleSink})::Dict{String,Any}
    return deepcopy(FLEXIBLE_SINK_ECONOMIC_PARAMETERS)
end

function emissions_parameters(x::Type{FlexibleSink})::Dict{String,Any}
    return deepcopy(FLEXIBLE_SINK_EMISSIONS_PARAMETERS)
end

function extract_parameter(x::Type{FlexibleSink}, config::Dict{String,Any}, param_name::String, param_def::NamedTuple,
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

function validate_config(x::Type{FlexibleSink}, config::Dict{String,Any}, extracted::Dict{String,Any}, uac::String,
                         sim_params::Dict{String,Any}, param_type::String)
    if param_type == "economy"
        parameter = economic_parameters(FlexibleSink)
        uac = uac * " - economic_parameters"
    elseif param_type == "emissions"
        parameter = emissions_parameters(FlexibleSink)
        uac = uac * " - emissions_parameters"
    elseif param_type == "component"
        parameter = component_parameters(FlexibleSink)
    end
    validate_config(Component, extracted, uac, sim_params, parameter)
end

function init_from_params(x::Type{FlexibleSink}, uac::String, params::Dict{String,Any},
                          raw_params::Dict{String,Any}, sim_params::Dict{String,Any})::Tuple
    medium = Symbol(params["medium"])

    max_power_profile = params["max_power_profile_file_path"] !== nothing ?
                        Profile(params["max_power_profile_file_path"], sim_params) :
                        nothing

    # return tuple in the order expected by new()
    return (uac,                                     # uac
            Controller(params["control_parameters"]),
            sf_flexible_sink,                        # sys_function
            medium,                                  # medium
            InterfaceMap(medium => nothing),         # input_interfaces
            InterfaceMap(medium => nothing),         # output_interfaces
            params["economic_parameters"],
            params["emissions_parameters"],
            max_power_profile,                       # max_power_profile
            some_or_none(params["temperature_profile_file_path"], params["temperature_from_global_file"]),
            params["scale"],                         # scaling_factor
            0.0,                                     # max_energy
            nothing,                                 # temperature
            params["constant_power"],                # constant_power
            params["constant_temperature"])          # constant_temperature
end

function initialise!(unit::FlexibleSink, sim_params::Dict{String,Any})
    set_storage_transfer!(unit.input_interfaces[unit.medium],
                          unload_storages(unit.controller, unit.medium))
end

function control(unit::FlexibleSink,
                 components::Grouping,
                 sim_params::Dict{String,Any})
    update(unit.controller)

    if unit.constant_power !== nothing
        unit.max_energy = sim_params["watt_to_wh"](unit.constant_power)
    elseif unit.max_power_profile !== nothing
        unit.max_energy = unit.scaling_factor * Profiles.work_at_time(unit.max_power_profile, sim_params)
    else
        unit.max_energy = 0.0
    end

    if unit.constant_temperature !== nothing
        unit.temperature = unit.constant_temperature
    elseif unit.temperature_profile !== nothing
        unit.temperature = Profiles.value_at_time(unit.temperature_profile, sim_params)
    end
    set_max_energy!(unit.input_interfaces[unit.medium], unit.max_energy, unit.temperature, nothing)
end

function process(unit::FlexibleSink, sim_params::Dict{String,Any})
    inface = unit.input_interfaces[unit.medium]
    exchanges = balance_on(inface, inface.source)

    # if we get the exchanges from a bus, the temperature check has already been performed
    if inface.source.sys_function == EnergySystems.sf_bus
        energy_supply = [e.balance + e.energy_potential for e in exchanges]
        temperature_min = [e.temperature_min for e in exchanges]
        temperature_max = [e.temperature_max for e in exchanges]
    else # check temperature
        energy_supply,
        temperature_min,
        temperature_max = check_temperatures_sink(exchanges, unit.temperature, unit.max_energy)
    end
    if sum(energy_supply; init=0.0) > 0.0
        sub!(inface, energy_supply, temperature_min, temperature_max)
    end
end

function get_reference_for_capex_and_embodied_emissions(unit::FlexibleSink)
    if unit.constant_power !== nothing
        return unit.max_energy
    else
        return unit.scaling_factor
    end
end

function output_values(unit::FlexibleSink)::Vector{String}
    if unit.temperature_profile === nothing && unit.constant_temperature === nothing
        return [string(unit.medium) * ":IN",
                "Max_Energy"]
    else
        return [string(unit.medium) * ":IN",
                "Max_Energy",
                "Temperature"]
    end
end

function output_value(unit::FlexibleSink, key::OutputKey)::Float64
    if key.value_key == "IN"
        return calculate_energy_flow(unit.input_interfaces[key.medium])
    elseif key.value_key == "Max_Energy"
        return unit.max_energy
    elseif key.value_key == "Temperature"
        return unit.temperature
    end
    throw(KeyError(key.value_key))
end

export FlexibleSink
