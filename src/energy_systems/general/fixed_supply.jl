#! format: off
const FIXED_SUPPLY_PARAMETERS = Dict(
    "medium" => (
        description="Medium of the supply (e.g. electricity, heat, gas, etc.)",
        display_name="Medium",
        required=true,
        type=String,
        json_type="string",
        unit="-"
    ),
    "temperature_profile_file_path" => (
        default=nothing,
        description="Path to a temperature profile file",
        display_name="Temperature profile file",
        required=false,
        conditionals=[
            ("temperature_from_global_file", "mutex"),
            ("constant_temperature", "mutex")
        ],
        type=String,
        json_type="string",
        unit="-"
    ),
    "temperature_from_global_file" => (
        default=false,
        description="If true, take the temperature profile from the global weather data file",
        display_name="Temperature from global file",
        required=false,
        conditionals=[
            ("temperature_profile_file_path", "mutex"),
            ("constant_temperature", "mutex")
        ],
        type=Bool,
        json_type="boolean",
        unit="-"
    ),
    "constant_temperature" => (
        default=nothing,
        description="Constant temperature value",
        display_name="Constant temperature",
        required=false,
        conditionals=[
            ("temperature_profile_file_path", "mutex"),
            ("temperature_from_global_file", "mutex")
        ],
        type=Float64,
        json_type="number",
        unit="°C"
    ),
    "energy_profile_file_path" => (
        default=nothing,
        description="Path to a profile file with energy values",
        display_name="Energy profile file",
        required=false,
        conditionals=[
            ("constant_supply", "mutex"),
        ],
        type=String,
        json_type="string",
        unit="-"
    ),
    "constant_supply" => (
        default=nothing,
        description="Constant supply (power, not work)",
        display_name="Constant supply",
        required=false,
        conditionals=[
            ("energy_profile_file_path", "mutex"),
        ],
        type=Float64,
        json_type="number",
        unit="W"
    ),
    "scale" => (
        default=1.0,
        description="Scaling factor for the energy profile",
        display_name="Energy scale",
        required=false,
        conditionals=[("energy_profile_file_path", "is_not_nothing")],
        type=Float64,
        json_type="number",
        unit="-"
    ),
)
#! format: on

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
        return new(SSOT_parameter_constructor(FixedSupply, uac, config, sim_params)...)
    end
end

function component_parameters(x::Type{FixedSupply})::Dict{String,NamedTuple}
    return deepcopy(FIXED_SUPPLY_PARAMETERS) # return a copy to prevent external modification
end

function extract_parameter(x::Type{FixedSupply}, config::Dict{String,Any}, param_name::String, param_def::NamedTuple,
                           sim_params::Dict{String,Any}, uac::String)
    if param_name == "constant_temperature" || param_name == "temperature_profile_file_path"
        constant_temperature,
        temperature_profile = get_parameter_profile_from_config(config,
                                                                sim_params,
                                                                "temperature",
                                                                "temperature_profile_file_path",
                                                                "temperature_from_global_file",
                                                                "constant_temperature",
                                                                uac)
        return param_name == "constant_temperature" ? constant_temperature : temperature_profile
    end

    return extract_parameter(Component, config, param_name, param_def, sim_params, uac)
end

function validate_config(x::Type{FixedSupply}, config::Dict{String,Any}, extracted::Dict{String,Any}, uac::String,
                         sim_params::Dict{String,Any})
    validate_config(Component, extracted, uac, sim_params, component_parameters(FixedSupply))
end

function init_from_params(x::Type{FixedSupply}, uac::String, params::Dict{String,Any},
                          raw_params::Dict{String,Any}, sim_params::Dict{String,Any})::Tuple
    medium = Symbol(params["medium"])
    register_media([medium])

    energy_profile = params["energy_profile_file_path"] !== nothing ?
                     Profile(params["energy_profile_file_path"], sim_params) :
                     nothing

    # return tuple in the order expected by new()
    return (uac,                                     # uac
            Controller(params["control_parameters"]),
            sf_fixed_source,                         # sys_function
            medium,                                  # medium
            InterfaceMap(medium => nothing),         # input_interfaces
            InterfaceMap(medium => nothing),         # output_interfaces
            energy_profile,                          # energy_profile
            params["temperature_profile_file_path"], # temperature_profile, might be from global weather data
            params["scale"],                         # scaling_factor
            0.0,                                     # supply
            nothing,                                 # temperature
            params["constant_supply"],               # constant_supply (power, not work!)
            params["constant_temperature"])          # constant_temperature
end

function initialise!(unit::FixedSupply, sim_params::Dict{String,Any})
    set_storage_transfer!(unit.output_interfaces[unit.medium],
                          load_storages(unit.controller, unit.medium))
end

function control(unit::FixedSupply,
                 components::Grouping,
                 sim_params::Dict{String,Any})
    update(unit.controller)

    if unit.constant_supply !== nothing
        unit.supply = sim_params["watt_to_wh"](unit.constant_supply)
    elseif unit.energy_profile !== nothing
        unit.supply = unit.scaling_factor * Profiles.work_at_time(unit.energy_profile, sim_params)
    else
        unit.supply = 0.0
    end

    if unit.constant_temperature !== nothing
        unit.temperature = unit.constant_temperature
    elseif unit.temperature_profile !== nothing
        unit.temperature = Profiles.value_at_time(unit.temperature_profile, sim_params)
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
