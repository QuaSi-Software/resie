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
        m_el_out = Symbol(default(config, "m_el_out", "m_e_ac_230v"))
        register_media([m_el_out])

        # load energy profile from path
        energy_profile = Profile(config["energy_profile_file_path"], sim_params)

        return new(uac, # uac
                   Controller(default(config, "control_parameters", nothing)),
                   sf_fixed_source,                   # sys_function
                   InterfaceMap(),                    # input_interfaces
                   InterfaceMap(m_el_out => nothing), # output_interfaces
                   m_el_out,
                   energy_profile,  # energy_profile
                   config["scale"], # scaling_factor
                   0.0)             # supply
    end
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
    return [string(unit.m_el_out) * " OUT",
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
