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
    efficiency_fuel_in::Function
    efficiency_el_out::Function
    efficiency_heat_out::Function
    # lookup tables for conversion of energy values to PLR
    fuel_in_to_plr::Vector{Tuple{Float64,Float64}}
    el_out_to_plr::Vector{Tuple{Float64,Float64}}
    heat_out_to_plr::Vector{Tuple{Float64,Float64}}
    discretization_step::Float64

    min_run_time::UInt
    output_temperature::Temperature

    losses::Float64

    function CHPP(uac::String, config::Dict{String,Any}, sim_params::Dict{String,Any})
        m_fuel_in = Symbol(default(config, "m_fuel_in", "m_c_g_natgas"))
        m_heat_out = Symbol(default(config, "m_heat_out", "m_h_w_ht1"))
        m_el_out = Symbol(default(config, "m_el_out", "m_e_ac_230v"))
        register_media([m_fuel_in, m_heat_out, m_el_out])

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
            config["power"],
            Symbol(default(config, "design_power_medium", m_el_out)),
            default(config, "min_power_fraction", 0.2),
            parse_efficiency_function(default(config,
                "efficiency_fuel_in",
                "pwlin:100,5.89,4,3.23,2.86,2.70,2.63,2.63,2.63"
            )),
            parse_efficiency_function(default(config,"efficiency_el_out", "const:1.0")),
            parse_efficiency_function(default(config,
                "efficiency_heat_out",
                "pwlin:80,4.06,2.52,1.87,1.57,1.41,1.32,1.29,1.29"
            )),
            [], # fuel_in_to_plr
            [], # el_out_to_plr
            [], # heat_out_to_plr
            1.0 / default(config, "nr_discretization_steps", 30),
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

    # fill energy_to_plr lookup tables
    media_def = (
        (unit.efficiency_fuel_in, unit.fuel_in_to_plr),
        (unit.efficiency_el_out, unit.el_out_to_plr),
        (unit.efficiency_heat_out, unit.heat_out_to_plr)
    )

    for (eff_func, lookup_table) in media_def
        for plr in collect(0.0:unit.discretization_step:1.0)
            push!(lookup_table, (watt_to_wh(unit.power) * plr * eff_func(plr), plr))
        end

        # check if inverse function (as lookup table) is monotonically increasing
        last_energy = 0.0
        for (energy, plr) in lookup_table
            if energy > sim_params["epsilon"] && energy <= last_energy
                @warn "PLR-from-energy function of component $(unit.uac) at PLR $plr " *
                    "is not monotonic"
            end
            last_energy = energy
        end
    end
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

function check_fuel_in(
    unit::CHPP,
    sim_params::Dict{String,Any}
)
    if !unit.controller.parameter["consider_m_fuel_in"]
        return Inf
    end

    if (
        unit.input_interfaces[unit.m_fuel_in].source.sys_function == sf_transformer
        && unit.input_interfaces[unit.m_fuel_in].max_energy === nothing
    )
        return Inf
    else
        exchanges = balance_on(
            unit.input_interfaces[unit.m_fuel_in],
            unit.input_interfaces[unit.m_fuel_in].source
        )
        potential_energy_fuel = balance(exchanges) + energy_potential(exchanges)
        if potential_energy_fuel <= sim_params["epsilon"]
            return nothing
        end
        return potential_energy_fuel
    end
end

function check_el_out(
    unit::CHPP,
    sim_params::Dict{String,Any}
)
    if !unit.controller.parameter["consider_m_el_out"]
        return -Inf
    end

    exchanges = balance_on(
        unit.output_interfaces[unit.m_el_out],
        unit.output_interfaces[unit.m_el_out].target
    )
    potential_energy_el = balance(exchanges) + energy_potential(exchanges)
    if potential_energy_el >= -sim_params["epsilon"]
        return nothing
    end

    return potential_energy_el
end

function check_heat_out(
    unit::CHPP,
    sim_params::Dict{String,Any}
)
    if !unit.controller.parameter["consider_m_heat_out"]
        return -Inf
    end

    exchanges = balance_on(
        unit.output_interfaces[unit.m_heat_out],
        unit.output_interfaces[unit.m_heat_out].target
    )

    # if we get multiple exchanges from balance_on, a bus is involved, which means the
    # temperature check has already been performed. we only need to check the case for
    # a single input which can happen for direct 1-to-1 connections or if the bus has
    # filtered inputs down to a single entry, which works the same as the 1-to-1 case
    if length(exchanges) > 1
        potential_energy_heat_out = balance(exchanges) + energy_potential(exchanges)
    else
        e = first(exchanges)
        if (
            unit.output_temperature === nothing
            || (e.temperature_min === nothing || e.temperature_min <= unit.output_temperature)
            && (e.temperature_max === nothing || e.temperature_max >= unit.output_temperature)
        )
            potential_energy_heat_out = e.balance + e.energy_potential
        else
            potential_energy_heat_out = 0.0
        end
    end

    if potential_energy_heat_out >= -sim_params["epsilon"]
        return nothing
    end
    return potential_energy_heat_out
end

"""
    plr_from_energy(unit, energy_value)

Calculate part load ratio as inverse function from the given energy value.

The efficiency for the CHPP is defined as relative to one input or output being linear in
respect to the PLR, while the other energy values must be calculated as the inverse of the
efficiency.

As the efficiency function can be a variety of functions and is not necessarily (easily)
invertable, the inverse is calculated numerically at initialisation as a piece-wise linear
function from a customizable number of support values from an even distribution of a PLR
from 0.0 to 1.0. When calling plr_from_energy this approximated function is evaluated for
the given energy value and a linear interpolation between the two surrounding support values
is performed to calculate the corresponding PLR.

# Arguments
- `unit::CHPP`: The CHPP
- `medium::Symbol`: The medium to which the energy value corresponds. Should be one of the
    media defined in the input of the CHPP.
- `energy_value::Float64`: The energy value
# Returns
- `Float64`: The part load ratio (from 0.0 to 1.0) as inverse from the energy value
"""
function plr_from_energy(unit::CHPP, medium::Symbol, energy_value::Float64)::Float64
    # shortcut if the given medium is the design power medium as it is linear to PLR
    if medium === unit.design_power_medium
        return energy_value / watt_to_wh(unit.power)
    end

    lookup_table = unit.el_out_to_plr
    if medium == unit.m_fuel_in
        lookup_table = unit.fuel_in_to_plr
    elseif medium == unit.m_heat_out
        lookup_table = unit.heat_out_to_plr
    end
    energy_at_max = last(lookup_table)[1]

    if energy_value <= 0.0
        return 0.0
    elseif energy_value >= energy_at_max
        return 1.0
    end

    nr_iter = 0
    candidate_idx = floor(Int64, length(lookup_table) * energy_value / energy_at_max)

    while (
        nr_iter < length(lookup_table)
        && candidate_idx < length(lookup_table)
        && candidate_idx >= 1
    )
        (energy_lb, plr_lb) = lookup_table[candidate_idx]
        (energy_ub, plr_ub) = lookup_table[candidate_idx+1]
        if energy_lb <= energy_value && energy_value < energy_ub
            return plr_lb + (plr_ub - plr_lb) *
                (energy_value - energy_lb) / (energy_ub - energy_lb)
        elseif energy_value < energy_lb
            candidate_idx -= 1
        elseif energy_value >= energy_ub
            candidate_idx += 1
        end

        nr_iter += 1
    end

    @warn "The energy_value of medium $(medium) in component $(unit.uac) is not within " *
        "the range of the lookup table."
    return 0.0
end

function calculate_energies(
    unit::CHPP,
    sim_params::Dict{String,Any},
)::Tuple{Bool, Floathing, Floathing, Floathing}
    # get max PLR of external profile, if any
    max_plr = (
        unit.controller.parameter["operation_profile_path"] === nothing
        ? 1.0
        : value_at_time(unit.controller.parameter["operation_profile"], sim_params["time"])
    )
    if max_plr <= 0.0
        return (false, nothing, nothing, nothing)
    end

    available_fuel_in = check_fuel_in(unit, sim_params)
    available_el_out = check_el_out(unit, sim_params)
    available_heat_out = check_heat_out(unit, sim_params)

    # shortcut if we're limited by zero input/output
    if (
        available_fuel_in === nothing
        || available_el_out === nothing
        || available_heat_out === nothing
    )
        return (false, nothing, nothing, nothing)
    end

    # in the following we want to work with positive values as it is easier
    available_el_out = abs(available_el_out)
    available_heat_out = abs(available_heat_out)

    # limit input/output to design power. for the medium designated as the design power
    # medium the efficiency at PLR of 1.0 is 1.0 while the others take efficiency into
    # account
    energy_at_max = watt_to_wh(unit.power)
    available_fuel_in = min(
        available_fuel_in,
        energy_at_max * unit.efficiency_fuel_in(1.0)
    )
    available_el_out = min(
        available_el_out,
        energy_at_max * unit.efficiency_el_out(1.0)
    )
    available_heat_out = min(
        available_heat_out,
        energy_at_max * unit.efficiency_heat_out(1.0)
    )

    plr_from_fuel_in = plr_from_energy(unit, unit.m_fuel_in, available_fuel_in)
    plr_from_el_out = plr_from_energy(unit, unit.m_el_out, available_el_out)
    plr_from_heat_out = plr_from_energy(unit, unit.m_heat_out, available_heat_out)
    used_plr = min(plr_from_fuel_in, plr_from_el_out, plr_from_heat_out, max_plr, 1.0)

    # check minimum power fraction limit
    if used_plr < unit.min_power_fraction
        return (false, nothing, nothing, nothing)
    end

    return (
        true,
        used_plr * energy_at_max * unit.efficiency_fuel_in(used_plr),
        used_plr * energy_at_max * unit.efficiency_el_out(used_plr),
        used_plr * energy_at_max * unit.efficiency_heat_out(used_plr),
    )
end

function potential(
    unit::CHPP,
    sim_params::Dict{String,Any}
)
    energies = calculate_energies(unit, sim_params)

    if !energies[1]
        set_max_energies!(unit, 0.0, 0.0, 0.0)
    else
        set_max_energies!(unit, energies[2], energies[3], energies[4])
    end
end

function process(unit::CHPP, sim_params::Dict{String,Any})
    energies = calculate_energies(unit, sim_params)

    if !energies[1]
        set_max_energies!(unit, 0.0, 0.0, 0.0)
        return
    end

    sub!(unit.input_interfaces[unit.m_fuel_in], energies[2])
    add!(unit.output_interfaces[unit.m_el_out], energies[3])
    add!(unit.output_interfaces[unit.m_heat_out], energies[4])

    unit.losses = energies[2] - energies[3] - energies[4]
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