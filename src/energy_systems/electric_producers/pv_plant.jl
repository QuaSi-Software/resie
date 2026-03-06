#! format: off
const PV_PLANT_PARAMETERS = Dict(
    "m_el_out" => (
        default="m_e_ac_230v",
        description="Medium of the output electricity",
        display_name="Medium el_out",
        required=false,
        type=String,
        json_type="string",
        unit="-"
    ),
    "energy_profile_file_path" => (
        description="Path to a profile file with the energy values",
        display_name="Energy profile file",
        required=true,
        type=String,
        json_type="string",
        unit="-"
    ),
    "scale" => (
        default=1.0,
        description="Scaling factor for the energy profile",
        display_name="Energy scale",
        required=false,
        type=Float64,
        json_type="number",
        unit="-"
    ),
)
#! format: on

"""
Implementation of a photovoltaic (PV) power plant.

No calculation of power is happening here, this is mostly just a wrapper around a yield
profile that must be calculated before a simulation and be imported.
"""
mutable struct PVPlant <: Component
    uac::String
    controller::Controller
    sys_function::SystemFunction

    input_interfaces::InterfaceMap
    output_interfaces::InterfaceMap

    m_el_out::Symbol

    energy_profile::Profile
    scaling_factor::Float64

    supply::Float64

    function PVPlant(uac::String, config::Dict{String,Any}, sim_params::Dict{String,Any})
        return new(SSOT_parameter_constructor(PVPlant, uac, config, sim_params)...)
    end
end

function component_parameters(x::Type{PVPlant})::Dict{String,NamedTuple}
    return deepcopy(PV_PLANT_PARAMETERS) # return a copy to prevent external modification
end

function extract_parameter(x::Type{PVPlant}, config::Dict{String,Any}, param_name::String, param_def::NamedTuple,
                           sim_params::Dict{String,Any}, uac::String)
    return extract_parameter(Component, config, param_name, param_def, sim_params, uac)
end

function validate_config(x::Type{PVPlant}, config::Dict{String,Any}, extracted::Dict{String,Any}, uac::String,
                         sim_params::Dict{String,Any})
    validate_config(Component, extracted, uac, sim_params, component_parameters(PVPlant))
end

function init_from_params(x::Type{PVPlant}, uac::String, params::Dict{String,Any},
                          raw_params::Dict{String,Any}, sim_params::Dict{String,Any})::Tuple
    # turn media names into Symbol and register them
    m_el_out = Symbol(params["m_el_out"])
    register_media([m_el_out])

    # load energy profile from path
    energy_profile = Profile(params["energy_profile_file_path"], sim_params)

    return (uac,
            Controller(params["control_parameters"]),
            sf_fixed_source,
            InterfaceMap(),
            InterfaceMap(m_el_out => nothing),
            m_el_out,
            energy_profile,
            params["scale"],
            0.0)  # supply
end

function initialise!(unit::PVPlant, sim_params::Dict{String,Any})
    set_storage_transfer!(unit.output_interfaces[unit.m_el_out],
                          load_storages(unit.controller, unit.m_el_out))
end

function control(unit::PVPlant,
                 components::Grouping,
                 sim_params::Dict{String,Any})
    update(unit.controller)
    unit.supply = unit.scaling_factor * Profiles.work_at_time(unit.energy_profile, sim_params)
    set_max_energy!(unit.output_interfaces[unit.m_el_out], unit.supply)
end

function process(unit::PVPlant, sim_params::Dict{String,Any})
    outface = unit.output_interfaces[unit.m_el_out]
    add!(outface, unit.supply)
end

function output_values(unit::PVPlant)::Vector{String}
    return [string(unit.m_el_out) * ":OUT",
            "Supply"]
end

function output_value(unit::PVPlant, key::OutputKey)::Float64
    if key.value_key == "OUT"
        return calculate_energy_flow(unit.output_interfaces[key.medium])
    elseif key.value_key == "Supply"
        return unit.supply
    end
    throw(KeyError(key.value_key))
end

export PVPlant
