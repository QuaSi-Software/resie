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

    power_th::Float64
    min_power_fraction::Float64
    efficiency::Function
    # lookup table for conversion of part load ratio to expended energy
    plr_to_expended_energy::Vector{Tuple{Float64,Float64}}

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
            config["power_th"],
            default(config, "min_power_fraction", 0.1),
            parse_efficiency_function(default(config, "efficiency", "const:0.9")),
            [], # plr_to_expended_energy
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

    # fill plr_to_expended_energy lookup table
    start_value = 0.0
    end_value = 1.0
    step_size = 0.1
    for plr in collect(start_value:step_size:end_value)
        # append tuple (expended energy, part load ratio value) to the lookup table
        push!(unit.plr_to_expended_energy, (
            watt_to_wh(unit.power_th) * plr / unit.efficiency(plr), plr
        ))
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

function plr_from_expended_energy(
    unit::FuelBoiler,
    intake_fuel::Float64
)::Float64
    energy_at_max = last(unit.plr_to_expended_energy)[1]

    if intake_fuel <= 0.0
        return 0.0
    elseif intake_fuel >= energy_at_max
        return 1.0
    end

    nr_iter = 0
    candidate_idx = floor(
        Int64,
        length(unit.plr_to_expended_energy) * intake_fuel / energy_at_max
    )

    while (
        nr_iter < length(unit.plr_to_expended_energy)
        && candidate_idx < length(unit.plr_to_expended_energy)
        && candidate_idx >= 1
    )
        (energy_lb, plr_lb) = unit.plr_to_expended_energy[candidate_idx]
        (energy_ub, plr_ub) = unit.plr_to_expended_energy[candidate_idx+1]
        if energy_lb <= intake_fuel && intake_fuel < energy_ub
            return plr_lb + (plr_ub - plr_lb) *
                (intake_fuel - energy_lb) / (energy_ub - energy_lb)
        elseif intake_fuel < energy_lb
            candidate_idx -= 1
        elseif intake_fuel >= energy_ub
            candidate_idx += 1
        end

        nr_iter += 1
    end

    @warn "The intake_fuel value in component $(unit.uac) is not within the range of the "
        "lookup table."
    return 0.0
end

function plr_from_produced_energy(unit::FuelBoiler, produced_heat::Float64)::Float64
    return produced_heat / watt_to_wh(unit.power_th)
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

    # limit input/output to design power
    available_fuel_in = min(
        available_fuel_in,
        watt_to_wh(unit.power_th) / unit.efficiency(1.0)
    )
    available_heat_out = min(available_heat_out, watt_to_wh(unit.power_th))

    plr_from_demand = plr_from_produced_energy(unit, available_heat_out)
    plr_from_supply = plr_from_expended_energy(unit, available_fuel_in)
    used_plr = min(plr_from_demand, plr_from_supply, 1.0)

    # check minimum power fraction limit
    if used_plr < unit.min_power_fraction
        return (false, nothing, nothing)
    end

    return (
        true,
        used_plr * watt_to_wh(unit.power_th) / unit.efficiency(used_plr),
        used_plr * watt_to_wh(unit.power_th),
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

export FuelBoiler, plr_from_expended_energy