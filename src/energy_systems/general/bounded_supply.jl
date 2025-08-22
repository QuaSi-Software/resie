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
                   sf_bounded_source,   # sys_function
                   medium,              # medium
                   InterfaceMap(medium => nothing), # input_interfaces
                   InterfaceMap(medium => nothing), # output_interfaces
                   max_power_profile,               # max_power_profile
                   temperature_profile,             # temperature_profile
                   default(config, "scale", 1.0),   # scaling_factor
                   0.0,                             # max_energy
                   nothing,                         # temperature
                   default(config, "constant_power", nothing),       # constant_power
                   constant_temperature)            # constant_temperature
    end
end

function initialise!(unit::BoundedSupply, sim_params::Dict{String,Any})
    set_storage_transfer!(unit.output_interfaces[unit.medium],
                          load_storages(unit.controller, unit.medium))
end

function control(unit::BoundedSupply,
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
    set_max_energy!(unit.output_interfaces[unit.medium], unit.max_energy, nothing, unit.temperature)
end

function process(unit::BoundedSupply, sim_params::Dict{String,Any})
    outface = unit.output_interfaces[unit.medium]
    exchanges = balance_on(outface, outface.target)

    # if we get the exchanges from a bus, the temperature check has already been performed
    if outface.target.sys_function == EnergySystems.sf_bus
        energy_demand = [e.balance + e.energy_potential for e in exchanges]
        temperature_min = [nothing for _ in exchanges]
        temperature_max = [unit.temperature for _ in exchanges]
    else # check temperature
        energy_demand,
        temperature_min,
        temperature_max, _ = check_temperatures_source(exchanges, unit.temperature, unit.max_energy)
    end

    if sum(energy_demand; init=0.0) < 0.0
        add!(outface, abs.(energy_demand), temperature_min, temperature_max)
    end
end

function output_values(unit::BoundedSupply)::Vector{String}
    if unit.temperature_profile === nothing && unit.constant_temperature === nothing
        return [string(unit.medium) * " OUT",
                "Max_Energy"]
    else
        return [string(unit.medium) * " OUT",
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
