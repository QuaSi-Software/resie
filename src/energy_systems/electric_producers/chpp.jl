"""
Implementation of a combined-heat-power-plant (CHPP) component.

For the moment this remains a simple implementation that converts natural gas into
electricity and heat (as medium m_h_w_ht1) at a defined ratio of 1:0.4:0.6. Has a minimum
run time of 1800s taken into consideration in its control behaviour and a minimum power
fraction of 20%. The power_gas is considered the maximum amount of both heat and electricity
that the CHPP can produce.

The only currently implemented operation strategy involves checking the load of a linked
buffer tank and en-/disabling the CHPP when a threshold is reached, in addition to an
overfill shutoff condition.
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
    design_power_medium::Symbol
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

        design_power_medium = Symbol(
            replace(
                default(config, "design_power_medium", "el_out"),
                "m_" => ""
            )
        )
        if !(design_power_medium in interface_list)
            @error "Given unknown design power medium $design_power_medium for $uac"
        end

        efficiencies = Dict{Symbol,Function}(
            Symbol("fuel_in") => parse_efficiency_function(default(config,
                "efficiency_fuel_in",
                "pwlin:100,5.89,4,3.23,2.86,2.70,2.63,2.63,2.63"
            )),
            Symbol("el_out") => parse_efficiency_function(default(config,
                "efficiency_el_out", "const:1.0"
            )),
            Symbol("heat_out") => parse_efficiency_function(default(config,
                "efficiency_heat_out",
                "pwlin:80,4.06,2.52,1.87,1.57,1.41,1.32,1.29,1.29"
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
            design_power_medium,
            default(config, "min_power_fraction", 0.2),
            efficiencies,
            interface_list,
            Dict{Symbol,Vector{Tuple{Float64,Float64}}}(), # energy_to_plr
            1.0 / default(config, "nr_discretization_steps", 30), # discretization_step
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