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
    medium = Symbol(config["medium"])
    register_media([medium])

    constant_temperature,
    temperature_profile = get_parameter_profile_from_config(config,
                                                            sim_params,
                                                            "temperature",
                                                            "temperature_profile_file_path",
                                                            "temperature_from_global_file",
                                                            "constant_temperature",
                                                            uac)

    sys_function = IsSource ? sf_flexible_source : sf_flexible_sink

    return GridConnection{IsSource}(uac,                        # uac
                                    Controller(default(config, "control_parameters", nothing)),
                                    sys_function,               # sys_function
                                    medium,                     # medium
                                    InterfaceMap(medium => nothing), # input_interfaces
                                    InterfaceMap(medium => nothing), # output_interfaces
                                    temperature_profile,        # temperature_profile
                                    constant_temperature,       # constant_temperature
                                    nothing,                    # temperature
                                    0.0,                        # output_sum
                                    0.0)                        # input_sum
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

const GridInput = GridConnection{true}

const GridOutput = GridConnection{false}

export GridConnection, GridInput, GridOutput
