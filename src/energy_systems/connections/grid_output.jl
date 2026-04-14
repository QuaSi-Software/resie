#! format: off
const GRID_OUTPUT_PARAMETERS = Dict(
    "medium" => (
        description="Medium of the grid output (e.g. electrical or thermal medium)",
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
        description="Constant temperature value for the grid output",
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
)

const GRID_OUTPUT_ECONOMIC_PARAMETERS = get_economic_standard_params("connection",
    Dict{String,Any}(
        "energy_price_profile_file_path" => nothing,
        "energy_price_profile_scale" => 1.0,
        "constant_energy_price" => nothing,
        "energy_price_change_rate_per_year" =>  0.02,
        "base_cost_per_year" => 0.0,
        "base_cost_change_rate_per_year" => 0.0
    ),
    Dict{String,Any}(),
)

const GRID_OUTPUT_EMISSIONS_PARAMETERS = get_emissions_standard_params("connection",
    Dict{String,Any}(
        "energy_emissions_profile_file_path" => nothing,
        "energy_emissions_profile_scale" => 1.0,
        "constant_energy_emissions" => nothing,
        "energy_emissions_change_rate_per_year" =>  0.0
    ),
    Dict{String,Any}(),
)
#! format: on

"""
Implementation of a component modelling the output connection to a public grid of a certain medium.

Public grids are considered to have an unlimited amount of energy they can provide. They
exist to model real connections to a public grid that can take in any remaining surplus of
energy. A grid output can also be configured to take in energy (if a thermal medium is
modelled) at a given temperature.
"""
mutable struct GridOutput <: Component
    uac::String
    controller::Controller
    sys_function::SystemFunction
    medium::Symbol

    input_interfaces::InterfaceMap
    output_interfaces::InterfaceMap

    economic_parameters::Dict{String,Any}
    emissions_parameters::Dict{String,Any}

    temperature_profile::Union{Profile,Nothing}
    constant_temperature::Temperature
    temperature::Temperature

    input_sum::Float64

    function GridOutput(uac::String, config::Dict{String,Any}, sim_params::Dict{String,Any})
        return new(SSOT_parameter_constructor(GridOutput, uac, config, sim_params)...)
    end
end

function component_parameters(x::Type{GridOutput})::Dict{String,Any}
    return deepcopy(GRID_OUTPUT_PARAMETERS)
end

function economic_parameters(x::Type{GridOutput})::Dict{String,Any}
    return deepcopy(GRID_OUTPUT_ECONOMIC_PARAMETERS)
end

function emissions_parameters(x::Type{GridOutput})::Dict{String,Any}
    return deepcopy(GRID_OUTPUT_EMISSIONS_PARAMETERS)
end

function extract_parameter(x::Type{GridOutput}, config::Dict{String,Any}, param_name::String,
                           param_def::NamedTuple, sim_params::Dict{String,Any}, uac::String)
    if param_name == "temperature_from_global_file"
        return load_profile_from_global_weather_file(config, param_name, sim_params, uac)
    elseif param_name == "temperature_profile_file_path"
        return load_optional_profile(config, param_name, sim_params)
    elseif param_name == "constant_temperature"
        return convert(Temperature, default(config, param_name, nothing))
    end

    return extract_parameter(Component, config, param_name, param_def, sim_params, uac)
end

function validate_config(x::Type{GridOutput}, config::Dict{String,Any}, extracted::Dict{String,Any},
                         uac::String, sim_params::Dict{String,Any}, param_type::String)
    if param_type == "economy"
        parameter = economic_parameters(GridOutput)
        uac = uac * " - economic_parameters"
    elseif param_type == "emissions"
        parameter = emissions_parameters(GridOutput)
        uac = uac * " - emissions_parameters"
    elseif param_type == "component"
        parameter = component_parameters(GridOutput)
    end
    validate_config(Component, extracted, uac, sim_params, parameter)
end

function init_from_params(x::Type{GridOutput}, uac::String, params::Dict{String,Any},
                          raw_params::Dict{String,Any}, sim_params::Dict{String,Any})::Tuple
    medium = Symbol(params["medium"])

    # return tuple in the order expected by new()
    return (uac,                                     # uac
            Controller(params["control_parameters"]),#
            sf_flexible_sink,                        # sys_function
            medium,                                  # medium
            InterfaceMap(medium => nothing),         # input_interfaces
            InterfaceMap(medium => nothing),         # output_interfaces
            params["economic_parameters"],
            params["emissions_parameters"],
            some_or_none(params["temperature_profile_file_path"], params["temperature_from_global_file"]),
            params["constant_temperature"],          # constant_temperature
            nothing,                                 # temperature
            0.0)                                     # input_sum
end

function initialise!(unit::GridOutput, sim_params::Dict{String,Any})
    set_storage_transfer!(unit.input_interfaces[unit.medium], unload_storages(unit.controller, unit.medium))
end

function control(unit::GridOutput, components::Grouping, sim_params::Dict{String,Any})
    update(unit.controller)

    if unit.constant_temperature !== nothing
        unit.temperature = unit.constant_temperature
    elseif unit.temperature_profile !== nothing
        unit.temperature = Profiles.value_at_time(unit.temperature_profile, sim_params)
    end

    set_max_energy!(unit.input_interfaces[unit.medium], Inf, unit.temperature, nothing)
end

function process(unit::GridOutput, sim_params::Dict{String,Any})
    inface = unit.input_interfaces[unit.medium]
    exchanges = balance_on(inface, inface.source)

    # if we get multiple exchanges from balance_on, a bus is involved, which means the
    # temperature check has already been performed. we only need to check the case for
    # a single input which can happen for direct 1-to-1 connections or if the bus has
    # filtered inputs down to a single entry, which works the same as the 1-to-1 case
    if length(exchanges) > 1
        energy_supply = balance(exchanges) + energy_potential(exchanges)
    else
        e = first(exchanges)
        if (unit.temperature === nothing ||
            (e.temperature_min === nothing || e.temperature_min <= unit.temperature) &&
            (e.temperature_max === nothing || e.temperature_max >= unit.temperature))
            # end of condition
            energy_supply = e.balance + e.energy_potential
        else
            energy_supply = 0.0
        end
    end

    if energy_supply > 0.0
        unit.input_sum += energy_supply
        sub!(inface, abs(energy_supply), unit.temperature, nothing)
    end
end

function output_values(unit::GridOutput)::Vector{String}
    output_vals = [string(unit.medium) * ":IN", "Input_sum"]
    if unit.temperature !== nothing
        push!(output_vals, "Temperature")
    end
    return output_vals
end

function output_value(unit::GridOutput, key::OutputKey)::Float64
    if key.value_key == "IN"
        return calculate_energy_flow(unit.input_interfaces[key.medium])
    elseif key.value_key == "Input_sum"
        return unit.input_sum
    elseif key.value_key == "Temperature"
        return unit.temperature
    end
    throw(KeyError(key.value_key))
end

export GridOutput
