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

    function FixedSink(uac::String, config::Dict{String,Any}, parameters::Dict{String,Any})
        energy_profile = "energy_profile_file_path" in keys(config) ?
                         Profile(config["energy_profile_file_path"], parameters) :
                         nothing
        temperature_profile = "temperature_profile_file_path" in keys(config) ?
                              Profile(config["temperature_profile_file_path"], parameters) :
                              nothing
        medium = Symbol(config["medium"])
        register_media([medium])

        return new(
            uac, # uac
            controller_for_strategy( # controller
                config["strategy"]["name"], config["strategy"], parameters
            ),
            sf_fixed_sink, # sys_function
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
            0.0, # demand
            nothing, # temperature
            default(config, "constant_demand", nothing), # constant_demand (power, not work!)
            default(config, "constant_temperature", nothing), # constant_temperature
        )
    end
end

function control(
    unit::FixedSink,
    components::Grouping,
    parameters::Dict{String,Any}
)
    move_state(unit, components, parameters)

    if unit.constant_demand !== nothing
        unit.demand = watt_to_wh(unit.constant_demand)
    elseif unit.energy_profile !== nothing
        unit.demand = unit.scaling_factor * Profiles.work_at_time(
            unit.energy_profile, parameters["time"]
        )
    else
        unit.demand = 0.0
    end
    set_max_energy!(unit.input_interfaces[unit.medium], unit.demand)

    if unit.constant_temperature !== nothing
        unit.temperature = unit.constant_temperature
    elseif unit.temperature_profile !== nothing
        unit.temperature = Profiles.value_at_time(
            unit.temperature_profile, parameters["time"]
        )
    end
    unit.input_interfaces[unit.medium].temperature = highest_temperature(
        unit.temperature,
        unit.input_interfaces[unit.medium].temperature
    )
end

function process(unit::FixedSink, parameters::Dict{String,Any})
    inface = unit.input_interfaces[unit.medium]
    sub!(inface, unit.demand, unit.temperature)
end

function output_values(unit::FixedSink)::Vector{String}
    if unit.temperature_profile === nothing && unit.constant_temperature === nothing
        return [string(unit.medium)*" IN",
                "Demand"]
    else
        return [string(unit.medium)*" IN",
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

"""
A component that models the demand consumers in a building require.

As the simulation does not encompass demand calculations, this is usually taken from other
tools that calculate the demand before an energy system simulation is performed. These
profiles usually are normalized to some degree, therefore Demand instances require a scaling
factor to turn the relative values to absolute values of required energy.

This is an alias to the generic implementation of a fixed sink.
"""
const Demand = FixedSink

export FixedSink, Demand