"""
Implementation of a component modeling an abstract bounded supply of some medium.

This is particularly useful for testing, but can also be used to model any bounded
component or other equipment unit that processes energy in a given medium. The component
might still have a maximum power draw in a single time step, but can provide any fraction
of this to connected components.
"""
mutable struct BoundedSupply <: Component
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

    function BoundedSupply(uac::String, config::Dict{String,Any}, sim_params::Dict{String,Any})
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
            sf_bounded_source, # sys_function
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

function initialise!(unit::BoundedSupply, sim_params::Dict{String,Any})
    set_storage_transfer!(
        unit.output_interfaces[unit.medium],
        default(
            unit.controller.parameter, "load_storages " * String(unit.medium), true
        )
    )
end

function control(
    unit::BoundedSupply,
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
    set_max_energy!(unit.output_interfaces[unit.medium], unit.max_energy)

    if unit.constant_temperature !== nothing
        unit.temperature = unit.constant_temperature
    elseif unit.temperature_profile !== nothing
        unit.temperature = Profiles.value_at_time(
            unit.temperature_profile, sim_params["time"]
        )
    end
    set_temperature!(
        unit.output_interfaces[unit.medium],
        nothing,
        unit.temperature
    )
end

function process(unit::BoundedSupply, sim_params::Dict{String,Any})
    outface = unit.output_interfaces[unit.medium]
    exchanges = balance_on(outface, outface.target)
    blnc = balance(exchanges) + energy_potential(exchanges)

    if (
        unit.controller.parameter["name"] == "extended_storage_control"
        && unit.controller.parameter["load_any_storage"]
    )
        energy_demand = blnc + storage_potential(exchanges)
    else
        energy_demand = blnc
    end

    if energy_demand < 0.0
        add!(outface, min(abs(energy_demand), unit.max_energy), unit.temperature)
    end
end

function output_values(unit::BoundedSupply)::Vector{String}
    if unit.temperature_profile === nothing && unit.constant_temperature === nothing
        return [string(unit.medium)*" OUT",
                "Max_Energy"]
    else
        return [string(unit.medium)*" OUT",
                "Max_Energy",
                "Temperature"]
    end
end

function output_value(unit::BoundedSupply, key::OutputKey)::Float64
    if key.value_key == "OUT"
        return calculate_energy_flow(unit.output_interfaces[key.medium])
    elseif key.value_key == "Max_Energy"
        return unit.max_energy
    elseif key.value_key == "Temperature"
        return unit.temperature
    end
    throw(KeyError(key.value_key))
end

export BoundedSupply