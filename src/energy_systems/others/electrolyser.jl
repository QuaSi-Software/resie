"""
Implementation of an electrolyser, turning electricity and water into H2, O2 and heat.

For the moment this remains a simple implementation that converts electricity into
the gases and heat (as medium m_h_w_ht1) at a defined ratio (default 1:0.6:0.4). Has a
minimum run time taken into consideration in its control behaviour and a minimum power
fraction in its processing. The power_el is considered the maximum amount of electricity
that the electrolyser can consume.
"""
mutable struct Electrolyser <: Component
    uac::String
    controller::Controller
    sys_function::SystemFunction

    input_interfaces::InterfaceMap
    output_interfaces::InterfaceMap

    m_el_in::Symbol
    m_heat_ht_out::Symbol
    m_heat_lt_out::Symbol
    m_h2_out::Symbol
    m_o2_out::Symbol

    power::Float64
    linear_interface::Symbol
    min_power_fraction::Float64
    # efficiency functions by input/output
    efficiencies::Dict{Symbol,Function}
    # list of names of input and output interfaces, used internally only
    interface_list::Tuple{Symbol,Symbol,Symbol,Symbol,Symbol}
    # lookup tables for conversion of energy values to PLR
    energy_to_plr::Dict{Symbol,Vector{Tuple{Float64,Float64}}}
    discretization_step::Float64

    min_run_time::UInt

    heat_lt_is_usable::Bool
    output_temperature_ht::Temperature
    output_temperature_lt::Temperature

    losses::Float64
    losses_heat::Float64
    losses_hydrogen::Float64

    function Electrolyser(uac::String, config::Dict{String,Any}, sim_params::Dict{String,Any})
        heat_lt_is_usable = default(config, "heat_lt_is_usable", false)

        m_el_in = Symbol(default(config, "m_el_in", "m_e_ac_230v"))
        m_heat_ht_out = Symbol(default(config, "m_heat_ht_out", "m_h_w_ht1"))
        m_heat_lt_out = Symbol(default(config, "m_heat_lt_out", "m_h_w_lt1"))
        m_h2_out = Symbol(default(config, "m_h2_out", "m_c_g_h2"))
        m_o2_out = Symbol(default(config, "m_o2_out", "m_c_g_o2"))
        register_media([m_el_in, m_heat_ht_out, m_heat_lt_out, m_h2_out, m_o2_out])
        interface_list = (
            Symbol("el_in"), Symbol("heat_ht_out"), Symbol("heat_lt_out"),
            Symbol("h2_out"), Symbol("o2_out")
        )

        linear_interface = Symbol(
            replace(
                default(config, "linear_interface", "el_in"),
                "m_" => ""
            )
        )
        if !(linear_interface in interface_list)
            @error "Given unknown interface name $linear_interface designated as linear " *
                "for component $uac"
        end

        efficiencies = Dict{Symbol,Function}(
            Symbol("el_in") => parse_efficiency_function(default(config,
                "efficiency_el_in", "const:1.0"
            )),
            Symbol("heat_ht_out") => parse_efficiency_function(default(config,
                "efficiency_heat_ht_out", "const:0.15"
            )),
            Symbol("heat_lt_out") => parse_efficiency_function(default(config,
                "efficiency_heat_lt_out", "const:0.07"
            )),
            Symbol("h2_out") => parse_efficiency_function(default(config,
                "efficiency_h2_out", "const:0.6"
            )),
            Symbol("o2_out") => parse_efficiency_function(default(config,
                "efficiency_o2_out", "const:0.6"
            )),
        )

        output_interfaces = InterfaceMap(
            m_heat_ht_out => nothing,
            m_h2_out => nothing,
            m_o2_out => nothing,
        )
        if heat_lt_is_usable
            output_interfaces[m_heat_lt_out] = nothing
        end

        return new(
            uac, # uac
            controller_for_strategy( # controller
                config["strategy"]["name"], config["strategy"], sim_params
            ),
            sf_transformer, # sys_function
            InterfaceMap( # input_interfaces
                m_el_in => nothing
            ),
            output_interfaces,
            m_el_in,
            m_heat_ht_out,
            m_heat_lt_out,
            m_h2_out,
            m_o2_out,
            config["power_el"] / efficiencies[Symbol("el_in")](1.0),
            linear_interface,
            default(config, "min_power_fraction", 0.2),
            efficiencies,
            interface_list,
            Dict{Symbol,Vector{Tuple{Float64,Float64}}}(), # energy_to_plr
            1.0 / default(config, "nr_discretization_steps", 8), # discretization_step
            default(config, "min_run_time", 3600),
            heat_lt_is_usable,
            default(config, "output_temperature_ht", 55.0),
            default(config, "output_temperature_lt", 25.0),
            0.0, # losses
            0.0, # losses_heat
            0.0  # losses_hydrogen
        )
    end
end

function initialise!(unit::Electrolyser, sim_params::Dict{String,Any})
    set_storage_transfer!(
        unit.input_interfaces[unit.m_el_in],
        default(
            unit.controller.parameter, "unload_storages " * String(unit.m_el_in), true
        )
    )

    set_storage_transfer!(
        unit.output_interfaces[unit.m_heat_ht_out],
        default(
            unit.controller.parameter, "load_storages " * String(unit.m_heat_ht_out), true
        )
    )

    if unit.heat_lt_is_usable
        set_storage_transfer!(
            unit.output_interfaces[unit.m_heat_lt_out],
            default(
                unit.controller.parameter,
                "load_storages " * String(unit.m_heat_lt_out),
                true
            )
        )
    else
        unit.controller.parameter["consider_m_heat_lt_out"] = false
    end

    set_storage_transfer!(
        unit.output_interfaces[unit.m_h2_out],
        default(
            unit.controller.parameter, "load_storages " * String(unit.m_h2_out), true
        )
    )

    set_storage_transfer!(
        unit.output_interfaces[unit.m_o2_out],
        default(
            unit.controller.parameter, "load_storages " * String(unit.m_o2_out), true
        )
    )

    unit.energy_to_plr = create_plr_lookup_tables(unit, sim_params)
end

function control(
    unit::Electrolyser,
    components::Grouping,
    sim_params::Dict{String,Any}
)
    move_state(unit, components, sim_params)
    set_temperature!(
        unit.output_interfaces[unit.m_heat_ht_out],
        nothing,
        unit.output_temperature_ht
    )
    if unit.heat_lt_is_usable
        set_temperature!(
            unit.output_interfaces[unit.m_heat_lt_out],
            nothing,
            unit.output_temperature_lt
        )
    end
end

function set_max_energies!(
    unit::Electrolyser, el_in::Float64, heat_ht_out::Float64, heat_lt_out::Float64,
    h2_out::Float64, o2_out::Float64
)
    set_max_energy!(unit.input_interfaces[unit.m_el_in], el_in)
    set_max_energy!(unit.output_interfaces[unit.m_heat_ht_out], heat_ht_out)
    if unit.heat_lt_is_usable
        set_max_energy!(unit.output_interfaces[unit.m_heat_lt_out], heat_lt_out)
    end
    set_max_energy!(unit.output_interfaces[unit.m_h2_out], h2_out)
    set_max_energy!(unit.output_interfaces[unit.m_o2_out], o2_out)
end

function calculate_energies(
    unit::Electrolyser,
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
    unit::Electrolyser,
    sim_params::Dict{String,Any}
)
    success, energies = calculate_energies(unit, sim_params)

    if !success
        set_max_energies!(unit, 0.0, 0.0, 0.0, 0.0, 0.0)
    else
        set_max_energies!(
            unit, energies[1], energies[2], energies[3], energies[4], energies[5]
        )
    end
end

function process(unit::Electrolyser, sim_params::Dict{String,Any})
    success, energies = calculate_energies(unit, sim_params)

    if !success
        set_max_energies!(unit, 0.0, 0.0, 0.0, 0.0, 0.0)
        return
    end

    sub!(unit.input_interfaces[unit.m_el_in], energies[1])
    add!(unit.output_interfaces[unit.m_heat_ht_out], energies[2])
    if unit.heat_lt_is_usable
        add!(unit.output_interfaces[unit.m_heat_lt_out], energies[3])
    end
    add!(unit.output_interfaces[unit.m_h2_out], energies[4])
    add!(unit.output_interfaces[unit.m_o2_out], energies[5])

    unit.losses_heat = energies[1] - energies[2] - energies[4] +
        (unit.heat_lt_is_usable ? -1 : 0) * energies[3]
    unit.losses_hydrogen = 0.0
    unit.losses = unit.losses_heat + unit.losses_hydrogen
end

# has its own reset function as here more losses are present that need to be reset in every timestep
function reset(unit::Electrolyser)
    for inface in values(unit.input_interfaces)
        if inface !== nothing
            reset!(inface)
        end
    end
    for outface in values(unit.output_interfaces)
        if outface !== nothing
            reset!(outface)
        end
    end

    # reset losses
    unit.losses = 0.0
    unit.losses_hydrogen = 0.0
    unit.losses_heat = 0.0
end

function output_values(unit::Electrolyser)::Vector{String}
    channels = [
        string(unit.m_el_in)*" IN",
        string(unit.m_h2_out)*" OUT",
        string(unit.m_o2_out)*" OUT",
        string(unit.m_heat_ht_out)*" OUT",
        "Losses",
        "Losses_heat",
        "Losses_hydrogen"
    ]

    if unit.heat_lt_is_usable
        append!(channels, [string(unit.m_heat_lt_out)*" OUT"])
        return channels
    else
        return channels
    end
end

function output_value(unit::Electrolyser, key::OutputKey)::Float64
    if key.value_key == "IN"
        return calculate_energy_flow(unit.input_interfaces[key.medium])
    elseif key.value_key == "OUT"
        return calculate_energy_flow(unit.output_interfaces[key.medium])
    elseif key.value_key == "Losses_heat"
        return unit.losses_heat
    elseif key.value_key == "Losses_hydrogen"
        return unit.losses_hydrogen
    elseif key.value_key == "Losses"
        return unit.losses_hydrogen + unit.losses_heat
    end
    throw(KeyError(key.value_key))
end

export Electrolyser