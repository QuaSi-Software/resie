#! format: off
const FUEL_BOILER_PARAMETERS = Dict(
    "m_fuel_in" => (
        description="Fuel input medium",
        display_name="Medium fuel_in",
        required=true,
        type=String,
        json_type="string",
        unit="-"
    ),
    "m_heat_out" => (
        default="m_h_w_ht1",
        description="Heat output medium",
        display_name="Medium heat_out",
        required=false,
        type=String,
        json_type="string",
        unit="-"
    ),
    "linear_interface" => (
        default="heat_out",
        description="Which interface is considered linear relative to the part-load-ratio",
        display_name="Linear interface",
        options=["fuel_in", "heat_out"],
        required=false,
        type=String,
        json_type="string",
        unit="-"
    ),
    "efficiency_fuel_in" => (
        default="const:1.1",
        description="Efficiency function for the fuel input",
        display_name="Efficiency fuel_in",
        required=false,
        type=String,
        json_type="string",
        function_type="1dim",
        unit="-"
    ),
    "efficiency_heat_out" => (
        default="const:1.0",
        description="Efficiency function for the heat output",
        display_name="Efficiency heat_out",
        required=false,
        type=String,
        json_type="string",
        function_type="1dim",
        unit="-"
    ),
    "power_th" => (
        default=nothing,
        description="Design thermal power",
        display_name="Thermal power",
        required=true,
        validations=[
            ("self", "value_gt_num", 0.0)
        ],
        type=Float64,
        json_type="number",
        unit="W"
    ),
    "min_power_fraction" => (
        default=0.1,
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
    "output_temperature" => (
        default=nothing,
        description="Fixed output temperature, or nothing for auto-detection",
        display_name="Output temperature",
        required=false,
        type=Floathing,
        json_type="number",
        unit="°C"
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
Implementation of a boiler producing heat (as hot water) from the chemical energy in a fuel.

While not technically correct, this implementation can also be used to model an electric
boiler, as the input medium can be freely chosen and therefore can be chosen to be
electricity instead of a chemical fuel.

Concerning the efficiency of fuel boilers special care should be taken in regards to which
fuel value the efficiencies correlate. There is a lot of difference between a condensing
gas boiler supplying a 45 °C heating demand and a wood boiler supplying a 80 °C DHW demand.
You can find more information on this topic in the online documentation.

Implements traits: PLRDEComponent
"""
mutable struct FuelBoiler <: Component
    uac::String
    controller::Controller
    sys_function::SystemFunction

    input_interfaces::InterfaceMap
    output_interfaces::InterfaceMap

    m_fuel_in::Symbol
    m_heat_out::Symbol

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

    output_temperature::Temperature
    losses::Float64

    function FuelBoiler(uac::String, config::Dict{String,Any}, sim_params::Dict{String,Any})
        return new(SSOT_parameter_constructor(FuelBoiler, uac, config, sim_params)...)
    end
end

function component_parameters(x::Type{FuelBoiler})::Dict{String,NamedTuple}
    return deepcopy(FUEL_BOILER_PARAMETERS) # Return a copy to prevent external modification
end

function extract_parameter(x::Type{FuelBoiler}, config::Dict{String,Any}, param_name::String,
                           param_def::NamedTuple, sim_params::Dict{String,Any}, uac::String)
    return extract_parameter(Component, config, param_name, param_def, sim_params, uac)
end

function validate_config(x::Type{FuelBoiler}, config::Dict{String,Any}, extracted::Dict{String,Any},
                         uac::String, sim_params::Dict{String,Any})
    validate_config(Component, extracted, uac, sim_params, component_parameters(FuelBoiler))
end

function init_from_params(x::Type{FuelBoiler}, uac::String, params::Dict{String,Any},
                          raw_params::Dict{String,Any}, sim_params::Dict{String,Any})::Tuple
    # turn media names into Symbol and register them
    m_fuel_in = Symbol(params["m_fuel_in"])
    m_heat_out = Symbol(params["m_heat_out"])
    register_media([m_fuel_in, m_heat_out])
    interface_list = (Symbol("fuel_in"), Symbol("heat_out"))
    linear_interface = Symbol(params["linear_interface"])

    efficiencies = Dict{Symbol,Function}(
        Symbol("fuel_in") => params["efficiency_fuel_in"],
        Symbol("heat_out") => params["efficiency_heat_out"],
    )

    # return tuple in the order expected by new()
    return (uac,
            Controller(params["control_parameters"]),
            sf_transformer,
            InterfaceMap(m_fuel_in => nothing),
            InterfaceMap(m_heat_out => nothing),
            m_fuel_in,
            m_heat_out,
            params["power_th"] / efficiencies[Symbol("heat_out")](1.0),
            linear_interface,
            params["min_power_fraction"],
            efficiencies,
            interface_list,
            Dict{Symbol,Vector{Tuple{Float64,Float64}}}(), # energy_to_plr
            1.0 / params["nr_discretization_steps"], # discretization_step
            params["output_temperature"],
            0.0) # losses
end

function initialise!(unit::FuelBoiler, sim_params::Dict{String,Any})
    set_storage_transfer!(unit.input_interfaces[unit.m_fuel_in],
                          unload_storages(unit.controller, unit.m_fuel_in))
    set_storage_transfer!(unit.output_interfaces[unit.m_heat_out],
                          load_storages(unit.controller, unit.m_heat_out))

    unit.energy_to_plr = create_plr_lookup_tables(unit, sim_params)
end

function control(unit::FuelBoiler,
                 components::Grouping,
                 sim_params::Dict{String,Any})
    update(unit.controller)
    set_max_energy!(unit.output_interfaces[unit.m_heat_out], nothing, nothing, unit.output_temperature)
end

"""
Set maximum energies that can be taken in and put out by the unit
"""
function set_max_energies!(unit::FuelBoiler, fuel_in::Float64, heat_out::Float64)
    set_max_energy!(unit.input_interfaces[unit.m_fuel_in], fuel_in)
    set_max_energy!(unit.output_interfaces[unit.m_heat_out], heat_out, nothing, unit.output_temperature)
end

function calculate_energies(unit::FuelBoiler,
                            sim_params::Dict{String,Any})::Tuple{Bool,Vector{Floathing}}
    # get maximum PLR from control modules
    max_plr = upper_plr_limit(unit.controller, sim_params)
    if max_plr <= 0.0
        return (false, [])
    end

    return calculate_energies_for_plrde(unit, sim_params, unit.min_power_fraction, max_plr)
end

function potential(unit::FuelBoiler, sim_params::Dict{String,Any})
    success, energies = calculate_energies(unit, sim_params)

    if !success || sum(energies[1]; init=0.0) < sim_params["epsilon"]
        set_max_energies!(unit, 0.0, 0.0)
    else
        set_max_energies!(unit, energies[1], energies[2])
    end
end

function process(unit::FuelBoiler, sim_params::Dict{String,Any})
    success, energies = calculate_energies(unit, sim_params)

    if !success || sum(energies[1]; init=0.0) < sim_params["epsilon"]
        unit.losses = 0.0
        set_max_energies!(unit, 0.0, 0.0)
        return
    end

    sub!(unit.input_interfaces[unit.m_fuel_in], energies[1])
    add!(unit.output_interfaces[unit.m_heat_out], energies[2], nothing, unit.output_temperature)

    unit.losses = check_epsilon(energies[1] - energies[2], sim_params)
end

function component_has_minimum_part_load(unit::FuelBoiler)
    return unit.min_power_fraction > 0.0
end

function output_values(unit::FuelBoiler)::Vector{String}
    return [string(unit.m_fuel_in) * ":IN",
            string(unit.m_heat_out) * ":OUT",
            "LossesGains"]
end

function output_value(unit::FuelBoiler, key::OutputKey)::Float64
    if key.value_key == "IN"
        return calculate_energy_flow(unit.input_interfaces[key.medium])
    elseif key.value_key == "OUT"
        return calculate_energy_flow(unit.output_interfaces[key.medium])
    elseif key.value_key == "LossesGains"
        return -unit.losses
    end
    throw(KeyError(key.value_key))
end

export FuelBoiler, plr_from_energy
