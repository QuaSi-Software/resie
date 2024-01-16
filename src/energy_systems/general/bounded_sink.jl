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

        temperature_profile = get_temperature_profile_from_config(config, sim_params, uac)

        medium = Symbol(config["medium"])
        register_media([medium])

        return new(
            uac, # uac
            controller_for_strategy( # controller
                config["strategy"]["name"], config["strategy"], sim_params
            ),
            sf_bounded_sink, # sys_function
            medium, # medium
            InterfaceMap( # input_interfaces
                medium => nothing
            ),
            InterfaceMap( # output_interfaces
                medium => nothing
            ),
            max_power_profile, # max_power_profile
            temperature_profile, #temperature_profile
            default(config, "scale", 1.0), # scaling_factor
            0.0, # max_energy
            nothing, # temperature
            default(config, "constant_power", nothing), # constant_power
            default(config, "constant_temperature", nothing), # constant_temperature
        )
    end
end

function control(
    unit::BoundedSink,
    components::Grouping,
    sim_params::Dict{String,Any}
)
    move_state(unit, components, sim_params)

    if unit.constant_power !== nothing
        unit.max_energy = watt_to_wh(unit.constant_power)
    elseif unit.max_power_profile !== nothing
        unit.max_energy = unit.scaling_factor * Profiles.work_at_time(
            unit.max_power_profile, sim_params["time"]
        )
    else
        unit.max_energy = 0.0
    end
    set_max_energy!(unit.input_interfaces[unit.medium], unit.max_energy)

    if unit.constant_temperature !== nothing
        unit.temperature = unit.constant_temperature
    elseif unit.temperature_profile !== nothing
        unit.temperature = Profiles.value_at_time(
            unit.temperature_profile, sim_params["time"]
        )
    end
    set_temperature!(unit.input_interfaces[unit.medium], highest(
        unit.temperature,
        unit.input_interfaces[unit.medium].temperature
    ))
end

function process(unit::BoundedSink, sim_params::Dict{String,Any})
    inface = unit.input_interfaces[unit.medium]
    exchanges = balance_on(inface, inface.source)
    blnc = balance(exchanges)
    if blnc > 0.0
        sub!(
            inface,
            min(abs(blnc), unit.max_energy)
        )
    end
end

function output_values(unit::BoundedSink)::Vector{String}
    if unit.temperature_profile === nothing && unit.constant_temperature === nothing
        return [string(unit.medium)*" IN",
                "Max_Energy"]
    else
        return [string(unit.medium)*" IN",
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