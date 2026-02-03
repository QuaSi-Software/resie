#! format: off
const GRID_CONNECTION_PARAMETERS = Dict(
    "medium" => (
        description="Medium of the grid connection (e.g. electrical or thermal medium)",
        display_name="Medium",
        required=true,
        type=String,
        json_type="string",
        unit="-"
    ),
    "is_source" => (
        default=false,
        description="If true, the connection provides energy (source); if false, it " *
                    "receives energy (sink). Can be omitted if aliases GridInput and " *
                    "GridOutput are used.",
        display_name="Is source?",
        required=true,
        type=Bool,
        json_type="boolean",
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
        description="Constant temperature value for the grid connection",
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
#! format: on

"""
Implementation of a component modeling the connection to a public grid of a certain medium.

Public grids are considered to have an unlimited amount of energy they can provide, but
might be limited in the power they can provide, although this behaviour is not yet
implemented. They exist to model real connections to a public grid that can provide any
remaining demand of energy or take in any excess of energy. To make it possible to model a
one-way connection they are split into two instances for providing or receiving energy and
must be handled as such in the input for constructing a project. A grid can also require or 
provide a given temperature.
"""
mutable struct GridConnection{IsSource} <: Component
    uac::String
    controller::Controller
    sys_function::SystemFunction
    medium::Symbol

    input_interfaces::InterfaceMap
    output_interfaces::InterfaceMap

    temperature_profile::Union{Profile,Nothing}
    constant_temperature::Temperature
    temperature::Temperature

    output_sum::Float64
    input_sum::Float64
end

# Outer constructor that uses the type parameter `IsSource` to determine whether this
# connection is a source or a sink.
function GridConnection{IsSource}(uac::String, config::Dict{String,Any}, sim_params::Dict{String,Any}) where {IsSource}
    constructor_errored = false

    # extract all parameters using the parameter dictionary as the source of truth
    extracted_params = Dict{String,Any}()
    for (param_name, param_def) in GRID_CONNECTION_PARAMETERS
        try
            extracted_params[param_name] = extract_parameter(GridConnection{IsSource}, config, param_name,
                                                             param_def, sim_params, uac)
        catch e
            @error "$(sprint(showerror, e))"
            constructor_errored = true
        end
    end

    try
        # extract control_parameters, which is essentially a subconfig
        extracted_params["control_parameters"] = extract_control_parameters(Component, config)

        # validate configuration, e.g. for interdependencies and allowed values
        validate_config(GridConnection{IsSource}, config, extracted_params, uac, sim_params)
    catch e
        @error "$(sprint(showerror, e))"
        constructor_errored = true
    end

    # we delayed throwing the errors upwards so that many errors are caught at once
    if constructor_errored
        throw(InputError("Can't construct component $uac because of errors parsing " *
                         "the parameter config. Check the error log for more " *
                         "information on each error."))
    end

    # initialize and construct the object
    init_values = init_from_params(GridConnection{IsSource}, uac, extracted_params, config, sim_params)
    return GridConnection{IsSource}(init_values...)
end

# Backwards-compatible constructor: if no type parameter is provided, use the
# `is_source` entry in `config` (if present) to pick the correct parametric type.
function GridConnection(uac::String, config::Dict{String,Any}, sim_params::Dict{String,Any})
    is_src = Bool(default(config, "is_source", false))
    if is_src
        return GridConnection{true}(uac, config, sim_params)
    else
        return GridConnection{false}(uac, config, sim_params)
    end
end

# type-aliases for more convenient construction
const GridInput = GridConnection{true}
const GridOutput = GridConnection{false}

function component_parameters(x::Type{GridConnection{IsSource}})::Dict{String,NamedTuple} where {IsSource}
    return deepcopy(GRID_CONNECTION_PARAMETERS) # return a copy to prevent external modification
end

function component_parameters(x::Type{GridInput})::Dict{String,NamedTuple}
    params = deepcopy(GRID_CONNECTION_PARAMETERS) # return a copy to prevent external modification
    delete!(params, "is_source") # remove is_source as this is true for GridInput by definition
    return params
end

function component_parameters(x::Type{GridOutput})::Dict{String,NamedTuple}
    params = deepcopy(GRID_CONNECTION_PARAMETERS) # return a copy to prevent external modification
    delete!(params, "is_source") # remove is_source as this is false for GridOutput by definition
    return params
end

function extract_parameter(x::Type{GridConnection{IsSource}}, config::Dict{String,Any}, param_name::String,
                           param_def::NamedTuple, sim_params::Dict{String,Any}, uac::String) where {IsSource}
    if param_name == "is_source"
        return IsSource
    end

    if param_name == "temperature_from_global_file"
        return load_profile_from_global_weather_file(config, param_name, sim_params, uac)
    elseif param_name == "temperature_profile_file_path"
        return load_optional_profile(config, param_name, sim_params)
    elseif param_name == "constant_temperature"
        return convert(Temperature, default(config, param_name, nothing))
    end

    return extract_parameter(Component, config, param_name, param_def, sim_params, uac)
end

function validate_config(x::Type{GridConnection{IsSource}}, config::Dict{String,Any}, extracted::Dict{String,Any},
                         uac::String, sim_params::Dict{String,Any}) where {IsSource}
    validate_config(Component, extracted, uac, sim_params, component_parameters(GridConnection{IsSource}))
end

function init_from_params(x::Type{GridConnection{IsSource}}, uac::String, params::Dict{String,Any},
                          raw_params::Dict{String,Any}, sim_params::Dict{String,Any})::Tuple where {IsSource}
    medium = Symbol(params["medium"])
    register_media([medium])

    sys_function = params["is_source"] ? sf_flexible_source : sf_flexible_sink

    # return tuple in the order expected by new()
    return (uac,                                     # uac
            Controller(params["control_parameters"]),
            sys_function,                            # sys_function
            medium,                                  # medium
            InterfaceMap(medium => nothing),         # input_interfaces
            InterfaceMap(medium => nothing),         # output_interfaces
            some_or_none(params["temperature_profile_file_path"], params["temperature_from_global_file"]),
            params["constant_temperature"],          # constant_temperature
            nothing,                                 # temperature
            0.0,                                     # output_sum
            0.0)                                     # input_sum
end

function initialise!(unit::GridConnection{IsSource}, sim_params::Dict{String,Any}) where {IsSource}
    if unit.sys_function === sf_flexible_source
        set_storage_transfer!(unit.output_interfaces[unit.medium],
                              load_storages(unit.controller, unit.medium))
    else
        set_storage_transfer!(unit.input_interfaces[unit.medium],
                              unload_storages(unit.controller, unit.medium))
    end
end

function control(unit::GridConnection{IsSource},
                 components::Grouping,
                 sim_params::Dict{String,Any}) where {IsSource}
    update(unit.controller)

    if unit.constant_temperature !== nothing
        unit.temperature = unit.constant_temperature
    elseif unit.temperature_profile !== nothing
        unit.temperature = Profiles.value_at_time(unit.temperature_profile, sim_params)
    end

    if unit.sys_function === sf_flexible_source
        set_max_energy!(unit.output_interfaces[unit.medium], Inf, nothing, unit.temperature)
    else
        set_max_energy!(unit.input_interfaces[unit.medium], Inf, unit.temperature, nothing)
    end
end

function process(unit::GridConnection{IsSource}, sim_params::Dict{String,Any}) where {IsSource}
    if unit.sys_function === sf_flexible_source
        outface = unit.output_interfaces[unit.medium]
        exchanges = balance_on(outface, outface.target)

        # if we get multiple exchanges from balance_on, a bus is involved, which means the
        # temperature check has already been performed. we only need to check the case for
        # a single output which can happen for direct 1-to-1 connections or if the bus has
        # filtered outputs down to a single entry, which works the same as the 1-to-1 case
        if length(exchanges) > 1
            energy_demand = balance(exchanges) + energy_potential(exchanges)
            temp_out = temp_min_highest(exchanges)
        else
            e = first(exchanges)
            if (unit.temperature === nothing ||
                (e.temperature_min === nothing || e.temperature_min <= unit.temperature) &&
                (e.temperature_max === nothing || e.temperature_max >= unit.temperature))
                # end of condition
                energy_demand = e.balance + e.energy_potential
                temp_out = lowest(e.temperature_min, unit.temperature)
            else
                energy_demand = 0.0
            end
        end
        if energy_demand < 0.0
            unit.output_sum += energy_demand
            add!(outface, abs(energy_demand), nothing, temp_out)
        end
    else
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
end

function output_values(unit::GridConnection{IsSource})::Vector{String} where {IsSource}
    output_vals = []
    if unit.sys_function == sf_flexible_source
        append!(output_vals, [string(unit.medium) * ":OUT", "Output_sum"])
    elseif unit.sys_function == sf_flexible_sink
        append!(output_vals, [string(unit.medium) * ":IN", "Input_sum"])
    end
    if unit.temperature !== nothing
        push!(output_vals, "Temperature")
    end
    return output_vals
end

function output_value(unit::GridConnection{IsSource}, key::OutputKey)::Float64 where {IsSource}
    if key.value_key == "IN"
        return calculate_energy_flow(unit.input_interfaces[key.medium])
    elseif key.value_key == "OUT"
        return calculate_energy_flow(unit.output_interfaces[key.medium])
    elseif key.value_key == "Output_sum"
        return unit.output_sum
    elseif key.value_key == "Input_sum"
        return unit.input_sum
    elseif key.value_key == "Temperature"
        return unit.temperature
    end
    throw(KeyError(key.value_key))
end

export GridConnection, GridInput, GridOutput
