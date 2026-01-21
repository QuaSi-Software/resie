
#! format: off
const UTIR_PARAMETERS = Dict(
    "m_el_in" => (
        description="Electricity input medium",
        display_name="Medium el_in",
        required=true,
        type=String,
        json_type="string",
        unit="-"
    ),
    "m_el_out" => (
        description="Electricity output medium",
        display_name="Medium el_out",
        required=true,
        type=String,
        json_type="string",
        unit="-"
    ),
    "linear_interface" => (
        default="el_in",
        description="Which interface is considered linear relative to the part-load-ratio",
        display_name="Linear interface",
        options=["el_in", "el_out"],
        required=false,
        type=String,
        json_type="string",
        unit="-"
    ),
    "efficiency_el_in" => (
        default="const:1.0",
        description="Efficiency function for the electricity input",
        display_name="Efficiency el_in",
        required=false,
        type=String,
        json_type="string",
        function_type="1dim",
        unit="-"
    ),
    "efficiency_el_out" => (
        default="const:1.0",
        description="Efficiency function for the electricity output",
        display_name="Efficiency el_out",
        required=false,
        type=String,
        json_type="string",
        function_type="1dim",
        unit="-"
    ),
    "power" => (
        default=nothing,
        description="Max. power rating",
        display_name="Max. power",
        required=true,
        validations=[
            ("self", "value_gt_num", 0.0)
        ],
        type=Float64,
        json_type="number",
        unit="W"
    ),
    "min_power_fraction" => (
        default=0.0,
        description="Minimum part-load ratio to operate",
        display_name="Min. power fraction",
        required=false,
        validations=[
            ("self", "value_gte_num", 0.0),
            ("self", "value_lte_num", 1.0)
        ],
        type=Float64,
        json_type="number",
        unit="-"
    ),
    "nr_discretization_steps" => (
        default=30,
        description="Number of intervals for interpolated efficiency functions",
        display_name="Nr. discretization steps",
        required=false,
        validations=[
            ("self", "value_gte_num", 1.0)
        ],
        type=UInt,
        json_type="number",
        unit="-"
    ),
)
#! format: on

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
        return new(SSOT_parameter_constructor(UTIR, uac, config, sim_params)...)
    end
end

function component_parameters(x::Type{UTIR})::Dict{String,NamedTuple}
    return deepcopy(UTIR_PARAMETERS) # Return a copy to prevent external modification
end

function extract_parameter(x::Type{UTIR}, config::Dict{String,Any}, param_name::String,
                           param_def::NamedTuple, sim_params::Dict{String,Any}, uac::String)
    return extract_parameter(Component, config, param_name, param_def, sim_params, uac)
end

function validate_config(x::Type{UTIR}, config::Dict{String,Any}, extracted::Dict{String,Any},
                         uac::String, sim_params::Dict{String,Any})
    validate_config(Component, extracted, uac, sim_params, component_parameters(UTIR))
end

function init_from_params(x::Type{UTIR}, uac::String, params::Dict{String,Any},
                          raw_params::Dict{String,Any}, sim_params::Dict{String,Any})::Tuple
    # turn media names into Symbol and register them
    m_el_in = Symbol(params["m_el_in"])
    m_el_out = Symbol(params["m_el_out"])
    register_media([m_el_in, m_el_out])
    interface_list = (Symbol("el_in"), Symbol("el_out"))
    linear_interface = Symbol(params["linear_interface"])

    efficiencies = Dict{Symbol,Function}(
        Symbol("el_in") => params["efficiency_el_in"],
        Symbol("el_out") => params["efficiency_el_out"],
    )

    # return tuple in the order expected by new()
    return (uac,
            Controller(params["control_parameters"]),
            sf_transformer,
            InterfaceMap(m_el_in => nothing),
            InterfaceMap(m_el_out => nothing),
            m_el_in,
            m_el_out,
            params["power"] / efficiencies[Symbol("el_out")](1.0),
            linear_interface,
            params["min_power_fraction"],
            efficiencies,
            interface_list,
            Dict{Symbol,Vector{Tuple{Float64,Float64}}}(), # energy_to_plr
            1.0 / params["nr_discretization_steps"], # discretization_step
            0.0) # losses
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
