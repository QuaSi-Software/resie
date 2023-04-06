"""
Implementation of an energy system modeling a generic bounded sink of a chosen medium.

This is particularly useful for testing, but can also be used to model any bounded
energy system or other equipment unit that consumes energy of a given medium. The system
might still have a maximum power intake in a single time step, but can consume any fraction
of this.
"""
mutable struct BoundedSink <: ControlledSystem
    uac::String
    controller::Controller
    sys_function::SystemFunction
    medium::Symbol

    input_interfaces::InterfaceMap
    output_interfaces::InterfaceMap

    max_power_profile::Profile
    temperature_profile::Union{Profile,Nothing}
    scaling_factor::Float64

    max_energy::Float64
    temperature::Temperature
    static_temperature::Temperature

    function BoundedSink(uac::String, config::Dict{String,Any})
        max_power_profile = Profile(config["max_power_profile_file_path"])
        temperature_profile = "temperature_profile_file_path" in keys(config) ?
                               Profile(config["temperature_profile_file_path"]) :
                               nothing
        medium = Symbol(config["medium"])
        register_media([medium])

        return new(
            uac, # uac
            controller_for_strategy( # controller
                config["strategy"]["name"], config["strategy"]
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
            config["scale"], # scaling_factor
            0.0, # max_energy
            nothing, # temperature
            default(config, "static_temperature", nothing) # static_temperature
        )
    end
end

function output_values(unit::BoundedSink)::Vector{String}
    return ["IN", "Max_Energy"]
end

function output_value(unit::BoundedSink, key::OutputKey)::Float64
    if key.value_key == "IN"
        return calculate_energy_flow(unit.input_interfaces[key.medium])
    elseif key.value_key == "Max_Energy"
        return unit.max_energy
    end
    throw(KeyError(key.value_key))
end

function control(
    unit::BoundedSink,
    systems::Grouping,
    parameters::Dict{String,Any}
)
    move_state(unit, systems, parameters)
    unit.max_energy = unit.scaling_factor * Profiles.work_at_time(unit.max_power_profile, parameters["time"])
    set_max_energy!(unit.input_interfaces[unit.medium], unit.max_energy)

    if unit.static_temperature !== nothing
        unit.temperature = unit.static_temperature
    elseif unit.temperature_profile !== nothing
        unit.temperature = Profiles.value_at_time(unit.temperature_profile, parameters["time"])
    end
    unit.input_interfaces[unit.medium].temperature = highest_temperature(unit.temperature, unit.input_interfaces[unit.medium].temperature)

end

function produce(unit::BoundedSink, parameters::Dict{String,Any}, watt_to_wh::Function)
    inface = unit.input_interfaces[unit.medium]
    exchange = balance_on(inface, inface.source)
    if exchange.balance > 0.0
        sub!(
            inface,
            min(abs(exchange.balance), unit.max_energy)
        )
    end
end

export BoundedSink