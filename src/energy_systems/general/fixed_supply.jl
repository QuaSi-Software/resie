"""
Implementation of an energy system modeling an abstract fixed supply of some medium.

This is particularly useful for testing, but can also be used to model any energy system
or other equipment unit that produces energy in a medium, all of which has to be consumed
as the energy system cannot be dispatched like a grid connection can.
Note that "fixed" in this context means that the amount of energy the unit produces is
fixed within a timestep, but can vary over multiple timesteps. No calculation other than
scaling of profile values is performed in each timestep.
"""
mutable struct FixedSupply <: ControlledSystem
    uac::String
    controller::Controller
    sys_function::SystemFunction
    medium::MediumCategory

    input_interfaces::InterfaceMap
    output_interfaces::InterfaceMap

    energy_profile::Profile
    temperature_profile::Union{Profile,Nothing}
    scaling_factor::Float64

    supply::Float64
    temperature::Temperature

    function FixedSupply(uac::String, config::Dict{String,Any})
        energy_profile = Profile(config["energy_profile_file_path"])
        temperature_profile = "temperature_profile_file_path" in keys(config) ?
                              Profile(config["temperature_profile_file_path"]) :
                              nothing
        medium = getproperty(EnergySystems, Symbol(config["medium"]))

        return new(
            uac, # uac
            controller_for_strategy( # controller
                config["strategy"]["name"], config["strategy"]
            ),
            sf_fixed_source, # sys_function
            medium, # medium
            InterfaceMap( # input_interfaces
                medium => nothing
            ),
            InterfaceMap( # output_interfaces
                medium => nothing
            ),
            energy_profile, # energy_profile
            temperature_profile, #temperature_profile
            config["scale"], # scaling_factor
            0.0, # supply
            nothing # temperature
        )
    end
end

function output_values(unit::FixedSupply)::Vector{String}
    return ["OUT", "Supply", "Temperature"]
end

function output_value(unit::FixedSupply, key::OutputKey)::Float64
    if key.value_key == "OUT"
        return unit.output_interfaces[key.medium].sum_abs_change * 0.5
    elseif key.value_key == "Supply"
        return unit.supply
    elseif key.value_key == "Temperature"
        return unit.temperature
    end
    throw(KeyError(key.value_key))
end

function control(
    unit::FixedSupply,
    systems::Grouping,
    parameters::Dict{String,Any}
)
    move_state(unit, systems, parameters)
    unit.supply = unit.scaling_factor * Profiles.work_at_time(unit.energy_profile, parameters["time"])
    if unit.temperature_profile !== nothing
        unit.temperature = Profiles.value_at_time(unit.temperature_profile, parameters["time"])
        unit.output_interfaces[unit.medium].temperature = unit.temperature
    end
end

function produce(unit::FixedSupply, parameters::Dict{String,Any}, watt_to_wh::Function)
    outface = unit.output_interfaces[unit.medium]
    add!(outface, unit.supply, unit.temperature)
end

export FixedSupply