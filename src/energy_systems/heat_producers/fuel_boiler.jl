"""
Implementation of a boiler producing heat (as hot water) from the chemical energy in a fuel.

While not technically correct, this implementation can also be used to model an electric
boiler, as the input medium can be freely chosen and therefore can be chosen to be
electricity instead of a chemical fuel.
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
    design_power_medium::Symbol
    min_power_fraction::Float64
    efficiency_fuel_in::Function
    efficiency_heat_out::Function
    # lookup tables for conversion of energy values to PLR
    fuel_in_to_plr::Vector{Tuple{Float64,Float64}}
    heat_out_to_plr::Vector{Tuple{Float64,Float64}}
    discretization_step::Float64

    min_run_time::UInt
    output_temperature::Temperature
    losses::Float64

    function FuelBoiler(uac::String, config::Dict{String,Any}, sim_params::Dict{String,Any})
        m_fuel_in = Symbol(config["m_fuel_in"])
        m_heat_out = Symbol(default(config, "m_heat_out", "m_h_w_ht1"))
        register_media([m_fuel_in, m_heat_out])

        return new(
            uac,
            controller_for_strategy(
                config["strategy"]["name"], config["strategy"], sim_params
            ),
            sf_transformer,
            InterfaceMap(
                m_fuel_in => nothing
            ),
            InterfaceMap(
                m_heat_out => nothing
            ),
            m_fuel_in,
            m_heat_out,
            config["power"],
            Symbol(default(config, "design_power_medium", m_heat_out)),
            default(config, "min_power_fraction", 0.1),
            parse_efficiency_function(default(config, "efficiency_fuel_in", "const:1.11")),
            parse_efficiency_function(default(config, "efficiency_heat_out", "const:1.0")),
            [], # fuel_in_to_plr
            [], # heat_out_to_plr
            1.0 / default(config, "nr_discretization_steps", 30),
            default(config, "min_run_time", 0),
            default(config, "output_temperature", nothing),
            0.0, # losses
        )
    end
end

function initialise!(unit::FuelBoiler, sim_params::Dict{String,Any})
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

    # fill energy_to_plr lookup tables
    media_def = (
        (unit.efficiency_fuel_in, unit.fuel_in_to_plr),
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
    unit::FuelBoiler,
    components::Grouping,
    sim_params::Dict{String,Any}
)
    move_state(unit, components, sim_params)
    set_temperature!(
        unit.output_interfaces[unit.m_heat_out],
        nothing,
        unit.output_temperature
    )
end

"""
Set maximum energies that can be taken in and put out by the unit
"""
function set_max_energies!(unit::FuelBoiler, fuel_in::Float64, heat_out::Float64)
    set_max_energy!(unit.input_interfaces[unit.m_fuel_in], fuel_in)
    set_max_energy!(unit.output_interfaces[unit.m_heat_out], heat_out)
end

function check_fuel_in(
    unit::FuelBoiler,
    sim_params::Dict{String,Any}
)::Floathing
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

function check_heat_out(
    unit::FuelBoiler,
    sim_params::Dict{String,Any}
)::Floathing
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

The efficiency for the FuelBoiler is defined as relative to one input or output being linear
in respect to the PLR, while the other energy values must be calculated as the inverse of
the efficiency.

As the efficiency function can be a variety of functions and is not necessarily (easily)
invertable, the inverse is calculated numerically at initialisation as a piece-wise linear
function from a customizable number of support values from an even distribution of a PLR
from 0.0 to 1.0. When calling plr_from_energy this approximated function is evaluated for
the given energy value and a linear interpolation between the two surrounding support values
is performed to calculate the corresponding PLR.

# Arguments
- `unit::FuelBoiler`: The fuel boiler
- `medium::Symbol`: The medium to which the energy value corresponds. Should be one of the
    media defined in the input of the fuel boiler.
- `energy_value::Float64`: The energy value
# Returns
- `Float64`: The part load ratio (from 0.0 to 1.0) as inverse from the energy value
"""
function plr_from_energy(unit::FuelBoiler, medium::Symbol, energy_value::Float64)::Float64
    # shortcut if the given medium is the design power medium as it is linear to PLR
    if medium === unit.design_power_medium
        return energy_value / watt_to_wh(unit.power)
    end

    lookup_table = unit.heat_out_to_plr
    if medium == unit.m_fuel_in
        lookup_table = unit.fuel_in_to_plr
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
    unit::FuelBoiler,
    sim_params::Dict{String,Any}
)::Tuple{Bool, Floathing, Floathing}
    # get max PLR of external profile, if any
    max_plr = (
        unit.controller.parameter["operation_profile_path"] === nothing
        ? 1.0
        : value_at_time(unit.controller.parameter["operation_profile"], sim_params["time"])
    )
    if max_plr <= 0.0
        return (false, nothing, nothing)
    end

    available_fuel_in = check_fuel_in(unit, sim_params)
    available_heat_out = check_heat_out(unit, sim_params)

    # shortcut if we're limited by zero input/output
    if available_fuel_in === nothing || available_heat_out === nothing
        return (false, nothing, nothing)
    end

    # in the following we want to work with positive values as it is easier
    available_heat_out = abs(available_heat_out)

    # limit input/output to design power. for the medium designated as the design power
    # medium the efficiency at PLR of 1.0 is 1.0 while the others take efficiency into
    # account
    energy_at_max = watt_to_wh(unit.power)
    available_fuel_in = min(
        available_fuel_in,
        energy_at_max * unit.efficiency_fuel_in(1.0)
    )
    available_heat_out = min(
        available_heat_out,
        energy_at_max * unit.efficiency_heat_out(1.0)
    )

    plr_from_fuel_in = plr_from_energy(unit, unit.m_fuel_in, available_fuel_in)
    plr_from_heat_out = plr_from_energy(unit, unit.m_heat_out, available_heat_out)
    used_plr = min(plr_from_fuel_in, plr_from_heat_out, max_plr, 1.0)

    # check minimum power fraction limit
    if used_plr < unit.min_power_fraction
        return (false, nothing, nothing)
    end

    return (
        true,
        used_plr * energy_at_max * unit.efficiency_fuel_in(used_plr),
        used_plr * energy_at_max * unit.efficiency_heat_out(used_plr),
    )
end

function potential(
    unit::FuelBoiler,
    sim_params::Dict{String,Any}
)
    energies = calculate_energies(unit, sim_params)

    if !energies[1]
        set_max_energies!(unit, 0.0, 0.0)
    else
        set_max_energies!(unit, energies[2], energies[3])
    end
end

function process(unit::FuelBoiler, sim_params::Dict{String,Any})
    energies = calculate_energies(unit, sim_params)

    if !energies[1]
        set_max_energies!(unit, 0.0, 0.0)
        return
    end

    sub!(unit.input_interfaces[unit.m_fuel_in], energies[2])
    add!(unit.output_interfaces[unit.m_heat_out], energies[3])

    unit.losses = energies[2] - energies[3]
end

function output_values(unit::FuelBoiler)::Vector{String}
    return [string(unit.m_fuel_in)*" IN", 
            string(unit.m_heat_out)*" OUT",
            "Losses"]
end

function output_value(unit::FuelBoiler, key::OutputKey)::Float64
    if key.value_key == "IN"
        return calculate_energy_flow(unit.input_interfaces[key.medium])
    elseif key.value_key == "OUT"
        return calculate_energy_flow(unit.output_interfaces[key.medium])
    elseif key.value_key == "Losses"
        return unit.losses
    end
    throw(KeyError(key.value_key))
end

export FuelBoiler, plr_from_energy