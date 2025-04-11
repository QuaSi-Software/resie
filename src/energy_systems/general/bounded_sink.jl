"""
Implementation of a component modeling a generic bounded sink of a chosen medium.

This is particularly useful for testing, but can also be used to model any bounded
component or other equipment unit that consumes energy of a given medium. The component
might still have a maximum power intake in a single time step, but can consume any fraction
of this.
"""
mutable struct BoundedSink <: Component
    uac::String
    controller::Controller
    sys_function::SystemFunction
    medium::Symbol

    input_interfaces::InterfaceMap
    output_interfaces::InterfaceMap

    max_power_profile::Union{Profile,Nothing}
    temperature_profile::Union{Profile,Nothing}
    scaling_factor::Float64

    max_energy::Float64
    temperature::Temperature
    constant_power::Union{Nothing,Float64}
    constant_temperature::Temperature

    function BoundedSink(uac::String, config::Dict{String,Any}, sim_params::Dict{String,Any})
        max_power_profile = "max_power_profile_file_path" in keys(config) ?
                            Profile(config["max_power_profile_file_path"], sim_params) :
                            nothing

        constant_temperature,
        temperature_profile = get_parameter_profile_from_config(config,
                                                                sim_params,
                                                                "temperature",
                                                                "temperature_profile_file_path",
                                                                "temperature_from_global_file",
                                                                "constant_temperature",
                                                                uac)
        medium = Symbol(config["medium"])
        register_media([medium])

        return new(uac, # uac
                   Controller(default(config, "control_parameters", nothing)),
                   sf_bounded_sink,                 # sys_function
                   medium,                          # medium
                   InterfaceMap(medium => nothing), # input_interfaces
                   InterfaceMap(medium => nothing), # output_interfaces
                   max_power_profile,               # max_power_profile
                   temperature_profile,             # temperature_profile
                   default(config, "scale", 1.0),   # scaling_factor
                   0.0,                             # max_energy
                   nothing,                         # temperature
                   default(config, "constant_power", nothing),    # constant_power
                   constant_temperature)            # constant_temperature
    end
end

function initialise!(unit::BoundedSink, sim_params::Dict{String,Any})
    set_storage_transfer!(unit.input_interfaces[unit.medium],
                          unload_storages(unit.controller, unit.medium))
end

function control(unit::BoundedSink,
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

function process(unit::BoundedSink, sim_params::Dict{String,Any})
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
        sub!(inface, min(energy_supply, unit.max_energy), unit.temperature)
    end
end

function output_values(unit::BoundedSink)::Vector{String}
    if unit.temperature_profile === nothing && unit.constant_temperature === nothing
        return [string(unit.medium) * " IN",
                "Max_Energy"]
    else
        return [string(unit.medium) * " IN",
                "Max_Energy",
                "Temperature"]
    end
end

function output_value(unit::BoundedSink, key::OutputKey)::Float64
    if key.value_key == "IN"
        return calculate_energy_flow(unit.input_interfaces[key.medium])
    elseif key.value_key == "Max_Energy"
        return unit.max_energy
    elseif key.value_key == "Temperature"
        return unit.temperature
    end
    throw(KeyError(key.value_key))
end

export BoundedSink
