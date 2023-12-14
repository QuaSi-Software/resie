"""
Implementation of a component modeling an abstract fixed supply of some medium.

This is particularly useful for testing, but can also be used to model any component
or other equipment unit that processes energy in a medium, all of which has to be consumed
as the component cannot be dispatched like a grid connection can.
Note that "fixed" in this context means that the amount of energy the unit processes is
fixed within a timestep, but can vary over multiple timesteps. No calculation other than
scaling of profile values is performed in each timestep.
"""
mutable struct FixedSupply <: Component
    uac::String
    controller::Controller
    sys_function::SystemFunction
    medium::Symbol

    input_interfaces::InterfaceMap
    output_interfaces::InterfaceMap

    energy_profile::Union{Profile,Nothing}
    temperature_profile::Union{Profile,Nothing}
    scaling_factor::Float64

    supply::Float64
    temperature::Temperature

    constant_supply::Union{Nothing,Float64}
    constant_temperature::Temperature

    function FixedSupply(uac::String, config::Dict{String,Any}, sim_params::Dict{String,Any})
        energy_profile = "energy_profile_file_path" in keys(config) ?
                         Profile(config["energy_profile_file_path"], sim_params) :
                         nothing

        temperature_profile = get_temperature_profile_from_config(config, sim_params, uac)

        medium = Symbol(config["medium"])
        register_media([medium])

        return new(
            uac, # uac
            controller_for_strategy( # controller
                config["strategy"]["name"], config["strategy"], sim_params
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
            nothing, # temperature
            default(config, "constant_supply", nothing), # constant_supply (power, not work!)
            default(config, "constant_temperature", nothing), # constant_temperature
        )
    end
end

function control(
    unit::FixedSupply,
    components::Grouping,
    sim_params::Dict{String,Any}
)
    move_state(unit, components, sim_params)

    if unit.constant_supply !== nothing
        unit.supply = watt_to_wh(unit.constant_supply)
    elseif unit.energy_profile !== nothing
        unit.supply = unit.scaling_factor * Profiles.work_at_time(unit.energy_profile, sim_params["time"])
    else
        unit.supply = 0.0
    end
    set_max_energy!(unit.output_interfaces[unit.medium], unit.supply)

    if unit.constant_temperature !== nothing
        unit.temperature = unit.constant_temperature
    elseif unit.temperature_profile !== nothing
        unit.temperature = Profiles.value_at_time(unit.temperature_profile, sim_params["time"])
    end
    unit.output_interfaces[unit.medium].temperature = highest_temperature(unit.temperature, unit.output_interfaces[unit.medium].temperature)
end

function process(unit::FixedSupply, sim_params::Dict{String,Any})
    outface = unit.output_interfaces[unit.medium]
    add!(outface, unit.supply, unit.temperature)
end

function output_values(unit::FixedSupply)::Vector{String}
    if unit.temperature_profile === nothing && unit.constant_temperature === nothing
        return [string(unit.medium)*" OUT",
                "Supply"]
    else
        return [string(unit.medium)*" OUT",
                "Supply",
                "Temperature"]
    end
end

function output_value(unit::FixedSupply, key::OutputKey)::Float64
    if key.value_key == "OUT"
        return calculate_energy_flow(unit.output_interfaces[key.medium])
    elseif key.value_key == "Supply"
        return unit.supply
    elseif key.value_key == "Temperature"
        return unit.temperature
    end
    throw(KeyError(key.value_key))
end

export FixedSupply