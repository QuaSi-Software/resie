"""
Implementation of a combined-heat-and-power plant (CHPP) component.

The implementation is flexible enough to cover any type of CHPP that takes in chemical fuel
and produces electricity and heat. The archetypical CHPP chosen uses natural gas in a
reprocating internal combustion engine and utilises heat extraction via jacketing and
exhaust cooling. Heat extraction via engine oil cooling is possible, but typically not done
with plants of less than around 500 kW electrical output power.

The heat is typically extracted at temperatures higher than those common in heating and
cooling systems and thus has been chosen as not having an upper limit to the temperatures.

Implements traits: PLRDEComponent
"""
mutable struct CHPP <: Component
    uac::String
    controller::Controller
    sys_function::SystemFunction

    input_interfaces::InterfaceMap
    output_interfaces::InterfaceMap

    m_fuel_in::Symbol
    m_heat_out::Symbol
    m_el_out::Symbol

    power::Float64
    linear_interface::Symbol
    min_power_fraction::Float64
    # efficiency functions by input/output
    efficiencies::Dict{Symbol,Function}
    # list of names of input and output interfaces, used internally only
    interface_list::Tuple{Symbol,Symbol,Symbol}
    # lookup tables for conversion of energy values to PLR
    energy_to_plr::Dict{Symbol,Vector{Tuple{Float64,Float64}}}
    discretization_step::Float64

    min_run_time::UInt
    output_temperature::Temperature

    losses::Float64

    function CHPP(uac::String, config::Dict{String,Any}, sim_params::Dict{String,Any})
        m_fuel_in = Symbol(default(config, "m_fuel_in", "m_c_g_natgas"))
        m_heat_out = Symbol(default(config, "m_heat_out", "m_h_w_ht1"))
        m_el_out = Symbol(default(config, "m_el_out", "m_e_ac_230v"))
        register_media([m_fuel_in, m_heat_out, m_el_out])
        interface_list = (Symbol("fuel_in"), Symbol("el_out"), Symbol("heat_out"))

        linear_interface = Symbol(
            replace(
                default(config, "linear_interface", "fuel_in"),
                "m_" => ""
            )
        )
        if !(linear_interface in interface_list)
            @error "Given unknown interface name $linear_interface designated as linear " *
                "for component $uac"
        end

        efficiencies = Dict{Symbol,Function}(
            Symbol("fuel_in") => parse_efficiency_function(default(config,
                "efficiency_fuel_in", "const:1.0"
            )),
            Symbol("el_out") => parse_efficiency_function(default(config,
                "efficiency_el_out",
                "pwlin:0.01,0.17,0.25,0.31,0.35,0.37,0.38,0.38,0.38"
            )),
            Symbol("heat_out") => parse_efficiency_function(default(config,
                "efficiency_heat_out",
                "pwlin:0.8,0.69,0.63,0.58,0.55,0.52,0.5,0.49,0.49"
            )),
        )

        return new(
            uac, # uac
            controller_for_strategy( # controller
                config["strategy"]["name"], config["strategy"], sim_params
            ),
            sf_transformer, # sys_function
            InterfaceMap( # input_interfaces
                m_fuel_in => nothing
            ),
            InterfaceMap( # output_interfaces
                m_heat_out => nothing,
                m_el_out => nothing
            ),
            m_fuel_in,
            m_heat_out,
            m_el_out,
            config["power_el"] / efficiencies[Symbol("el_out")](1.0),
            linear_interface,
            default(config, "min_power_fraction", 0.2),
            efficiencies,
            interface_list,
            Dict{Symbol,Vector{Tuple{Float64,Float64}}}(), # energy_to_plr
            1.0 / default(config, "nr_discretization_steps", 8), # discretization_step
            default(config, "min_run_time", 1800),
            default(config, "output_temperature", nothing),
            0.0, # losses
        )
    end
end

function initialise!(unit::CHPP, sim_params::Dict{String,Any})
    set_storage_transfer!(
        unit.input_interfaces[unit.m_fuel_in],
        default(
            unit.controller.parameter, "unload_storages " * String(unit.m_fuel_in), true
        )
    )
    set_storage_transfer!(
        unit.output_interfaces[unit.m_heat_out],
        default(
            unit.controller.parameter, "load_storages " * String(unit.m_heat_out), true
        )
    )
    set_storage_transfer!(
        unit.output_interfaces[unit.m_el_out],
        default(
            unit.controller.parameter, "load_storages " * String(unit.m_el_out), true
        )
    )

    unit.energy_to_plr = create_plr_lookup_tables(unit, sim_params)
end

function control(
    unit::CHPP,
    components::Grouping,
    sim_params::Dict{String,Any}
)
    move_state(unit, components, sim_params)
    set_temperature!(
        unit.output_interfaces[unit.m_heat_out],
        nothing,
        unit.output_temperature,
    )
end

function set_max_energies!(unit::CHPP, fuel_in::Float64, el_out::Float64, heat_out::Float64)
    set_max_energy!(unit.input_interfaces[unit.m_fuel_in], fuel_in)
    set_max_energy!(unit.output_interfaces[unit.m_el_out], el_out)
    set_max_energy!(unit.output_interfaces[unit.m_heat_out], heat_out)
end

function calculate_energies(
    unit::CHPP,
    sim_params::Dict{String,Any},
)::Tuple{Bool, Vector{Floathing}}
    # check operational state for strategy storage_driven
    if (
        unit.controller.strategy == "storage_driven"
        && unit.controller.state_machine.state != 2
    )
        return (false, [])
    end

    # get max PLR of external profile, if any
    max_plr = (
        unit.controller.parameter["operation_profile_path"] === nothing
        ? 1.0
        : value_at_time(unit.controller.parameter["operation_profile"], sim_params["time"])
    )
    if max_plr <= 0.0
        return (false, [])
    end

    return calculate_energies_for_plrde(
        unit, sim_params, unit.min_power_fraction, max_plr
    )
end

function potential(
    unit::CHPP,
    sim_params::Dict{String,Any}
)
    success, energies = calculate_energies(unit, sim_params)

    if !success
        set_max_energies!(unit, 0.0, 0.0, 0.0)
    else
        set_max_energies!(unit, energies[1], energies[2], energies[3])
    end
end

function process(unit::CHPP, sim_params::Dict{String,Any})
    success, energies = calculate_energies(unit, sim_params)

    if !success
        set_max_energies!(unit, 0.0, 0.0, 0.0)
        return
    end

    sub!(unit.input_interfaces[unit.m_fuel_in], energies[1])
    add!(unit.output_interfaces[unit.m_el_out], energies[2])
    add!(unit.output_interfaces[unit.m_heat_out], energies[3])

    unit.losses = energies[1] - energies[2] - energies[3]
end

function output_values(unit::CHPP)::Vector{String}
    return [string(unit.m_fuel_in)*" IN",
            string(unit.m_el_out)*" OUT",
            string(unit.m_heat_out)*" OUT",
            "Losses"]
end

function output_value(unit::CHPP, key::OutputKey)::Float64
    if key.value_key == "IN"
        return calculate_energy_flow(unit.input_interfaces[key.medium])
    elseif key.value_key == "OUT"
        return calculate_energy_flow(unit.output_interfaces[key.medium])
    elseif key.value_key == "Losses"
        return unit.losses
    end
    throw(KeyError(key.value_key))
end

export CHPP