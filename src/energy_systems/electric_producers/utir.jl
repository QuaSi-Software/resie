"""
Unified implementation of electric transformers, inverters and rectifiers (UTIR).

This model can be used as a transformer-type component to connect electric components
requiring or providing both AC and DC current at different voltage levels.

Implements traits: PLRDEComponent
"""
mutable struct UTIR <: Component
    uac::String
    controller::Controller
    sys_function::SystemFunction

    input_interfaces::InterfaceMap
    output_interfaces::InterfaceMap

    m_el_in::Symbol
    m_el_out::Symbol

    power::Float64
    linear_interface::Symbol
    min_power_fraction::Float64
    # efficiency functions by input/output
    efficiencies::Dict{Symbol,Function}
    # list of names of input and output interfaces, used internally only
    interface_list::Tuple{Symbol,Symbol}
    # lookup tables for conversion of energy values to PLR
    energy_to_plr::Dict{Symbol,Vector{Tuple{Float64,Float64}}}
    discretization_step::Float64

    losses::Float64

    function UTIR(uac::String, config::Dict{String,Any}, sim_params::Dict{String,Any})
        m_el_in = Symbol(config["m_el_in"])
        m_el_out = Symbol(config["m_el_out"])
        register_media([m_el_in, m_el_out])
        interface_list = (Symbol("el_in"), Symbol("el_out"))

        linear_interface = Symbol(replace(default(config, "linear_interface", "el_in"), "m_" => ""))
        if !(linear_interface in interface_list)
            @error "Given unknown interface name $linear_interface designated as linear " *
                   "for component $uac"
        end

        efficiencies = Dict{Symbol,Function}(
            Symbol("el_in") => parse_efficiency_function(default(config, "efficiency_el_in", "const:1.0")),
            Symbol("el_out") => parse_efficiency_function(default(config, "efficiency_el_out", "const:1.0")),
        )

        return new(uac,
                   Controller(default(config, "control_parameters", nothing)),
                   sf_transformer,
                   InterfaceMap(m_el_in => nothing),
                   InterfaceMap(m_el_out => nothing),
                   m_el_in,
                   m_el_out,
                   config["power"] / efficiencies[Symbol("el_out")](1.0),
                   linear_interface,
                   default(config, "min_power_fraction", 0.1),
                   efficiencies,
                   interface_list,
                   Dict{Symbol,Vector{Tuple{Float64,Float64}}}(),        # energy_to_plr
                   1.0 / default(config, "nr_discretization_steps", 30), # discretization_step
                   0.0)  # losses
    end
end

function initialise!(unit::UTIR, sim_params::Dict{String,Any})
    set_storage_transfer!(unit.input_interfaces[unit.m_el_in],
                          unload_storages(unit.controller, unit.m_el_in))
    set_storage_transfer!(unit.output_interfaces[unit.m_el_out],
                          load_storages(unit.controller, unit.m_el_out))

    unit.energy_to_plr = create_plr_lookup_tables(unit, sim_params)
end

function control(unit::UTIR,
                 components::Grouping,
                 sim_params::Dict{String,Any})
    update(unit.controller)
    set_max_energy!(unit.output_interfaces[unit.m_el_out], nothing)
end

"""
Set maximum energies that can be taken in and put out by the unit
"""
function set_max_energies!(unit::UTIR, el_in::Float64, el_out::Float64)
    set_max_energy!(unit.input_interfaces[unit.m_el_in], el_in)
    set_max_energy!(unit.output_interfaces[unit.m_el_out], el_out)
end

function calculate_energies(unit::UTIR, sim_params::Dict{String,Any})::Tuple{Bool,Vector{Floathing}}
    # get maximum PLR from control modules
    max_plr = upper_plr_limit(unit.controller, sim_params)
    if max_plr <= 0.0
        return (false, [])
    end

    return calculate_energies_for_plrde(unit, sim_params, unit.min_power_fraction, max_plr)
end

function potential(unit::UTIR, sim_params::Dict{String,Any})
    success, energies = calculate_energies(unit, sim_params)

    if !success || sum(energies[1]; init=0.0) < sim_params["epsilon"]
        set_max_energies!(unit, 0.0, 0.0)
    else
        set_max_energies!(unit, energies[1], energies[2])
    end
end

function process(unit::UTIR, sim_params::Dict{String,Any})
    success, energies = calculate_energies(unit, sim_params)

    if !success || sum(energies[1]; init=0.0) < sim_params["epsilon"]
        unit.losses = 0.0
        set_max_energies!(unit, 0.0, 0.0)
        return
    end

    sub!(unit.input_interfaces[unit.m_el_in], energies[1])
    add!(unit.output_interfaces[unit.m_el_out], energies[2])

    unit.losses = check_epsilon(energies[1] - energies[2], sim_params)
end

function component_has_minimum_part_load(unit::UTIR)
    return unit.min_power_fraction > 0.0
end

function output_values(unit::UTIR)::Vector{String}
    return [string(unit.m_el_in) * ":IN",
            string(unit.m_el_out) * ":OUT",
            "LossesGains"]
end

function output_value(unit::UTIR, key::OutputKey)::Float64
    if key.value_key == "IN"
        return calculate_energy_flow(unit.input_interfaces[key.medium])
    elseif key.value_key == "OUT"
        return calculate_energy_flow(unit.output_interfaces[key.medium])
    elseif key.value_key == "LossesGains"
        return -unit.losses
    end
    throw(KeyError(key.value_key))
end

export UTIR
