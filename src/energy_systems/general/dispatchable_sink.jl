"""
Implementation of an energy system modeling a generic dispatchable sink of a chosen medium.

This is particularly useful for testing, but can also be used to model any dispatchable
energy system or other equipment unit that consumes energy of a given medium. The system
might still have a maximum power intake in a single time step, but can consume any fraction
of this.
"""
mutable struct DispatchableSink <: ControlledSystem
    uac::String
    controller::Controller
    sys_function::SystemFunction
    medium::Symbol

    input_interfaces::InterfaceMap
    output_interfaces::InterfaceMap

    max_power_profile::Profile
    scaling_factor::Float64

    max_power::Float64

    function DispatchableSink(uac::String, config::Dict{String,Any})
        max_power_profile = Profile(config["max_power_profile_file_path"])
        medium = getproperty(EnergySystems, Symbol(config["medium"]))

        return new(
            uac, # uac
            controller_for_strategy( # controller
                config["strategy"]["name"], config["strategy"]
            ),
            sf_dispatchable_sink, # sys_function
            medium, # medium
            InterfaceMap( # input_interfaces
                medium => nothing
            ),
            InterfaceMap( # output_interfaces
                medium => nothing
            ),
            max_power_profile, # max_power_profile
            config["scale"], # scaling_factor
            0.0, # max_power
        )
    end
end

function output_values(unit::DispatchableSink)::Vector{String}
    return ["IN", "Max_Power"]
end

function output_value(unit::DispatchableSink, key::OutputKey)::Float64
    if key.value_key == "IN"
        return unit.input_interfaces[key.medium].sum_abs_change * 0.5
    elseif key.value_key == "Max_Power"
        return unit.max_power
    end
    throw(KeyError(key.value_key))
end

function control(
    unit::DispatchableSink,
    systems::Grouping,
    parameters::Dict{String,Any}
)
    move_state(unit, systems, parameters)
    unit.max_power = unit.scaling_factor * Profiles.power_at_time(unit.max_power_profile, parameters["time"])
    set_max_power!(unit.input_interfaces[unit.medium], unit.max_power)
end

function produce(unit::DispatchableSink, parameters::Dict{String,Any}, watt_to_wh::Function)
    inface = unit.input_interfaces[unit.medium]
    balance, _, _ = balance_on(inface, inface.source)
    if balance > 0.0
        sub!(
            inface,
            min(abs(balance), watt_to_wh(unit.max_power))
        )
    end
end

export DispatchableSink