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

    treat_profile_as_volume_flow_in_qm_per_hour::Bool
    rho_medium::Floathing
    cp_medium::Floathing

    function FixedSupply(uac::String, config::Dict{String,Any}, sim_params::Dict{String,Any})
        energy_profile = "energy_profile_file_path" in keys(config) ?
                         Profile(config["energy_profile_file_path"], sim_params) :
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
                   sf_fixed_source,                 # sys_function
                   medium,                          # medium
                   InterfaceMap(medium => nothing), # input_interfaces
                   InterfaceMap(medium => nothing), # output_interfaces
                   energy_profile,                  # energy_profile
                   temperature_profile,             # temperature_profile
                   default(config, "scale", 1.0),   # scaling_factor
                   0.0,                             # supply
                   nothing,                         # temperature
                   default(config, "constant_supply", nothing),      # constant_supply (power, not work!)
                   constant_temperature,            # constant_temperature
                   default(config, "treat_profile_as_volume_flow_in_qm_per_hour", false),
                   default(config, "rho_medium", nothing),  # [kg/m^3]
                   default(config, "cp_medium", nothing))   # [J/kgK]
    end
end

function initialise!(unit::FixedSupply, sim_params::Dict{String,Any})
    set_storage_transfer!(unit.output_interfaces[unit.medium],
                          load_storages(unit.controller, unit.medium))

    if unit.treat_profile_as_volume_flow_in_qm_per_hour
        if unit.rho_medium === nothing
            @error "In fixed supply $(unit.uac), the profile should be treated as volume flow. Please provide the medium density rho_medium."
        end
        if unit.cp_medium === nothing
            @error "In fixed supply $(unit.uac), the profile should be treated as volume flow. Please provide the medium thermal capacity cp_medium."
        end
    end
end

function control(unit::FixedSupply,
                 components::Grouping,
                 sim_params::Dict{String,Any})
    update(unit.controller)

    if unit.constant_temperature !== nothing
        unit.temperature = unit.constant_temperature
    elseif unit.temperature_profile !== nothing
        unit.temperature = Profiles.value_at_time(unit.temperature_profile, sim_params)
    end
    if unit.constant_supply !== nothing
        unit.supply = sim_params["watt_to_wh"](unit.constant_supply)
    elseif unit.energy_profile !== nothing
        if unit.treat_profile_as_volume_flow_in_qm_per_hour &&
           isa(unit.output_interfaces[unit.medium].target, LimitCoolingInputTemperatureTarget)
            t_low = unit.output_interfaces[unit.medium].target.current_energy_input_return_temperature # [°C]
            t_high = unit.temperature # [°C]
            volume_flow = unit.scaling_factor * Profiles.power_at_time(unit.energy_profile, sim_params) # [m^3/h]
            power = unit.rho_medium * volume_flow * unit.cp_medium / 3600 * (t_high - t_low) # [W]
            unit.supply = sim_params["watt_to_wh"](power)
        else
            unit.supply = unit.scaling_factor * Profiles.work_at_time(unit.energy_profile, sim_params)
        end
    else
        unit.supply = 0.0
    end

    set_max_energy!(unit.output_interfaces[unit.medium], unit.supply, nothing, unit.temperature)
end

function process(unit::FixedSupply, sim_params::Dict{String,Any})
    outface = unit.output_interfaces[unit.medium]
    add!(outface, unit.supply, nothing, unit.temperature)
end

function output_values(unit::FixedSupply)::Vector{String}
    if unit.temperature_profile === nothing && unit.constant_temperature === nothing
        return [string(unit.medium) * ":OUT",
                "Supply"]
    else
        return [string(unit.medium) * ":OUT",
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
