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

    treat_profile_as_volume_flow_in_qm_per_hour::Bool
    rho_medium::Floathing
    cp_medium::Floathing

    function FixedSink(uac::String, config::Dict{String,Any}, sim_params::Dict{String,Any})
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
                   sf_fixed_sink,   # sys_function
                   medium,          # medium
                   InterfaceMap(medium => nothing),  # input_interfaces
                   InterfaceMap(medium => nothing),  # output_interfaces
                   energy_profile,                   # energy_profile
                   temperature_profile,              # temperature_profile
                   default(config, "scale", 1.0),    # scaling_factor
                   0.0,     # demand
                   nothing, # temperature
                   default(config, "constant_demand", nothing),       # constant_demand (power, not work!)
                   constant_temperature,   # constant_temperature
                   default(config, "treat_profile_as_volume_flow_in_qm_per_hour", false),
                   default(config, "rho_medium", nothing),  # [kg/m^3]
                   default(config, "cp_medium", nothing))   # [J/kgK]
    end
end

function initialise!(unit::FixedSink, sim_params::Dict{String,Any})
    set_storage_transfer!(unit.input_interfaces[unit.medium],
                          unload_storages(unit.controller, unit.medium))
    if unit.treat_profile_as_volume_flow_in_qm_per_hour
        if unit.rho_medium === nothing
            @error "In fixed sink $(unit.uac), the profile should be treated as volume flow. Please provide the medium density rho_medium."
        end
        if unit.cp_medium === nothing
            @error "In fixed sink $(unit.uac), the profile should be treated as volume flow. Please provide the medium thermal capacity cp_medium."
        end
    end
end

function control(unit::FixedSink,
                 components::Grouping,
                 sim_params::Dict{String,Any})
    update(unit.controller)

    if unit.constant_temperature !== nothing
        unit.temperature = unit.constant_temperature
    elseif unit.temperature_profile !== nothing
        unit.temperature = Profiles.value_at_time(unit.temperature_profile, sim_params)
    end

    if unit.constant_demand !== nothing
        unit.demand = sim_params["watt_to_wh"](unit.constant_demand)
    elseif unit.energy_profile !== nothing
        if unit.treat_profile_as_volume_flow_in_qm_per_hour &&
           isa(unit.input_interfaces[unit.medium].source, LimitCoolingInputTemperatureTarget)
            t_low = unit.temperature # [°C]
            t_high = unit.input_interfaces[unit.medium].source.current_max_output_temperature # [°C]
            volume_flow = unit.scaling_factor * Profiles.power_at_time(unit.energy_profile, sim_params) # [m^3/h]
            power = unit.rho_medium * volume_flow * unit.cp_medium / 3600 * (t_high - t_low) # [W]
            unit.demand = sim_params["watt_to_wh"](power)
        else
            unit.demand = unit.scaling_factor * Profiles.work_at_time(unit.energy_profile, sim_params)
        end
    else
        unit.demand = 0.0
    end

    set_max_energy!(unit.input_interfaces[unit.medium], unit.demand, unit.temperature, nothing)
end

function process(unit::FixedSink, sim_params::Dict{String,Any})
    inface = unit.input_interfaces[unit.medium]
    sub!(inface, unit.demand, unit.temperature, nothing)
end

function output_values(unit::FixedSink)::Vector{String}
    if unit.temperature_profile === nothing && unit.constant_temperature === nothing
        return [string(unit.medium) * ":IN",
                "Demand"]
    else
        return [string(unit.medium) * ":IN",
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
