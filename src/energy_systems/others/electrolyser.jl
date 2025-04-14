"""
Implementation of an electrolyser, turning electricity and pure water into H2, O2 and heat.

At time of writing only electrolysers splitting pure water are supported as they are the
most relevant technology at time of writing. The produced heat has a high temperature output
(depending on technology 50-65 °C) and an optional low temperature output (25-35 °C), which
is more difficult to utilise in reality as this is waste heat from cooling power electronics
and the heat the stacks lose to the equipment housing. If the low temperature heat is not
used, the output is counted towards the heat losses.

The electrolyser consists of a customizable number of subunits (which can just be one),
each of which has its own power equipment. The efficiencies are considered to apply for each
unit. As the minimum PLR for the units is typically higher than the total minimum PLR of the
electrolyser there is a dispatch mechanism for the units taking this into account.

Implements traits: PLRDEComponent
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
    power_total::Float64
    nr_units::Integer
    dispatch_strategy::String
    optimal_unit_plr::Float64

    linear_interface::Symbol
    min_power_fraction::Float64
    min_power_fraction_total::Float64
    # efficiency functions by input/output
    efficiencies::Dict{Symbol,Function}
    # list of names of input and output interfaces, used internally only
    interface_list::Tuple{Symbol,Symbol,Symbol,Symbol,Symbol}
    # lookup tables for conversion of energy values to PLR
    energy_to_plr::Dict{Symbol,Vector{Tuple{Float64,Float64}}}
    discretization_step::Float64

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
        interface_list = (Symbol("el_in"),
                          Symbol("heat_ht_out"),
                          Symbol("heat_lt_out"),
                          Symbol("h2_out"),
                          Symbol("o2_out"))

        linear_interface = Symbol(replace(default(config, "linear_interface", "el_in"), "m_" => ""))
        if !(linear_interface in interface_list)
            @error "Given unknown interface name $linear_interface designated as linear " *
                   "for component $uac"
        end

        efficiencies = Dict{Symbol,Function}(
            Symbol("el_in") => parse_efficiency_function(default(config,
                                                                 "efficiency_el_in", "const:1.0")),
            Symbol("heat_ht_out") => parse_efficiency_function(default(config,
                                                                       "efficiency_heat_ht_out", "const:0.15")),
            Symbol("heat_lt_out") => parse_efficiency_function(default(config,
                                                                       "efficiency_heat_lt_out", "const:0.07")),
            Symbol("h2_out") => parse_efficiency_function(default(config,
                                                                  "efficiency_h2_out", "const:0.57")),
            Symbol("h2_out_lossless") => parse_efficiency_function(default(config,
                                                                           "efficiency_h2_out_lossless", "const:0.6")),
            Symbol("o2_out") => parse_efficiency_function(default(config,
                                                                  "efficiency_o2_out", "const:0.6")),
        )

        output_interfaces = InterfaceMap(m_heat_ht_out => nothing,
                                         m_h2_out => nothing,
                                         m_o2_out => nothing)
        if heat_lt_is_usable
            output_interfaces[m_heat_lt_out] = nothing
        else
            # if low temperature heat output is not used, make sure it is not limiting
            efficiencies[Symbol("heat_lt_out")] => parse_efficiency_function("const:1.0")
        end

        nr_units = default(config, "nr_switchable_units", 1)
        power_total = config["power_el"] / efficiencies[Symbol("el_in")](1.0)
        power = power_total / nr_units

        dispatch_strategy = default(config, "dispatch_strategy", "equal_with_mpf")
        if !(dispatch_strategy in ("all_equal", "try_optimal", "equal_with_mpf"))
            @error "Unknown dispatch strategy $dispatch_strategy for electrolyser $uac"
        end

        return new(uac, # uac
                   Controller(default(config, "control_parameters", nothing)),
                   sf_transformer,                    # sys_function
                   InterfaceMap(m_el_in => nothing),  # input_interfaces
                   output_interfaces,
                   m_el_in,
                   m_heat_ht_out,
                   m_heat_lt_out,
                   m_h2_out,
                   m_o2_out,
                   power,
                   power_total,
                   nr_units,
                   dispatch_strategy,
                   default(config, "optimal_unit_plr", 0.65),
                   linear_interface,
                   default(config, "min_power_fraction", 0.4),
                   default(config, "min_power_fraction_total", 0.2),
                   efficiencies,
                   interface_list,
                   Dict{Symbol,Vector{Tuple{Float64,Float64}}}(),   # energy_to_plr
                   1.0 / default(config, "nr_discretization_steps", 1), # discretization_step
                   heat_lt_is_usable,
                   default(config, "output_temperature_ht", 55.0),
                   default(config, "output_temperature_lt", 25.0),
                   0.0, # losses
                   0.0, # losses_heat
                   0.0) # losses_hydrogen
    end
end

function initialise!(unit::Electrolyser, sim_params::Dict{String,Any})
    set_storage_transfer!(unit.input_interfaces[unit.m_el_in],
                          unload_storages(unit.controller, unit.m_el_in))

    set_storage_transfer!(unit.output_interfaces[unit.m_heat_ht_out],
                          load_storages(unit.controller, unit.m_heat_ht_out))

    if unit.heat_lt_is_usable
        set_storage_transfer!(unit.output_interfaces[unit.m_heat_lt_out],
                              load_storages(unit.controller, unit.m_heat_lt_out))
    else
        unit.controller.parameters["consider_m_heat_lt_out"] = false
    end

    set_storage_transfer!(unit.output_interfaces[unit.m_h2_out],
                          load_storages(unit.controller, unit.m_h2_out))

    set_storage_transfer!(unit.output_interfaces[unit.m_o2_out],
                          load_storages(unit.controller, unit.m_o2_out))

    unit.energy_to_plr = create_plr_lookup_tables(unit, sim_params)
end

function control(unit::Electrolyser,
                 components::Grouping,
                 sim_params::Dict{String,Any})
    update(unit.controller)

    set_max_energy!(unit.output_interfaces[unit.m_heat_ht_out], nothing, nothing, unit.output_temperature_ht)
    if unit.heat_lt_is_usable
        set_max_energy!(unit.output_interfaces[unit.m_heat_lt_out], nothing, nothing, unit.output_temperature_lt)
    end
end

function set_max_energies!(unit::Electrolyser,
                           el_in::Float64,
                           heat_ht_out::Float64,
                           heat_lt_out::Float64,
                           h2_out::Float64,
                           o2_out::Float64)
    set_max_energy!(unit.input_interfaces[unit.m_el_in], el_in)
    set_max_energy!(unit.output_interfaces[unit.m_heat_ht_out], heat_ht_out, nothing, unit.output_temperature_ht)
    if unit.heat_lt_is_usable
        set_max_energy!(unit.output_interfaces[unit.m_heat_lt_out], heat_lt_out, nothing, unit.output_temperature_lt)
    end
    set_max_energy!(unit.output_interfaces[unit.m_h2_out], h2_out)
    set_max_energy!(unit.output_interfaces[unit.m_o2_out], o2_out)
end

"""
    dispatch_units(ely::Electrolyser, plr::Float64, limit_name::Symbol, limit_value)

Calculate the number of active units and the PLR for each in order to meet the given limit.

# Arguments
- `ely::Electrolyser`: The electrolyser
- `plr::Float64`: The total PLR over the whole electrolyser assembly
- `limit_name::Symbol`: The name of the interface that is limiting. Should be on of the
    values in the `interface_list` field.
- `limit_value::Float64`: The limiting value to meet. Can be larger than the total power of
    the electroylser, in which case all units are utilised to their full extent.
# Returns
- `Integer`: Number of active units
- `Float64`: PLR of each active unit
"""
function dispatch_units(ely::Electrolyser,
                        plr::Float64,
                        limit_name::Symbol,
                        limit_value::Float64,
                        w2wh::Function)::Tuple{Integer,Float64}
    if limit_value == Inf
        return ely.nr_units, 1.0
    end

    if ely.dispatch_strategy == "try_optimal"
        optimal_val_per_unit = ely.optimal_unit_plr * w2wh(ely.power) *
                               ely.efficiencies[limit_name](ely.optimal_unit_plr)
        nr_units = max(1, min(ceil(limit_value / optimal_val_per_unit - 0.5), ely.nr_units))
        plr_per_unit = plr_from_energy(ely, limit_name, limit_value / nr_units, w2wh)

    elseif ely.dispatch_strategy == "equal_with_mpf"
        if plr >= ely.min_power_fraction
            nr_units = ely.nr_units
            plr_per_unit = plr
        else
            min_val_per_unit = ely.min_power_fraction * w2wh(ely.power) *
                               ely.efficiencies[limit_name](ely.min_power_fraction)
            nr_units = max(1, min(floor(limit_value / min_val_per_unit), ely.nr_units))
            plr_per_unit = plr_from_energy(ely, limit_name, limit_value / nr_units, w2wh)
        end

    elseif ely.dispatch_strategy == "all_equal"
        nr_units = ely.nr_units
        plr_per_unit = plr
    end

    return nr_units, plr_per_unit
end

function calculate_energies(unit::Electrolyser, sim_params::Dict{String,Any})::Tuple{Bool,Vector{Floathing}}
    # get maximum PLR from control modules
    max_plr = upper_plr_limit(unit.controller, sim_params)
    if max_plr <= 0.0
        return (false, [])
    end

    # calculate limiting interfaces and the total PLR that meets the limit
    limiting_plr = 1.0
    limiting_energy = Inf
    limiting_interface = Symbol("h2_out")
    plr_from_nrg = []

    for name in unit.interface_list
        availability = getproperty(EnergySystems, Symbol("check_" * String(name)))
        energy = availability(unit, sim_params)

        # shortcut if we're limited by zero input/output
        if energy === nothing
            return (false, [])
        end

        # in the following we want to work with positive values as it is easier
        energy = abs(energy)

        # limit to total design power
        energy = min(sim_params["watt_to_wh"](unit.power_total) * unit.efficiencies[name](1.0), energy)

        # we can get the total PLR by assuming all units are activated equally, even if
        # dispatch happens differently later
        plr = plr_from_energy(unit, name, energy / unit.nr_units, sim_params["watt_to_wh"])
        push!(plr_from_nrg, plr)

        # keep track which was the limiting interface and how much energy is on that
        # interface. if all interfaces are infinite, we're limited by the design power or
        # some external condition, in which case hydrogen will be the limiting interface
        if plr < limiting_plr
            limiting_plr = plr
            limiting_energy = energy
            limiting_interface = name
        end
    end

    # the operation point of the electrolyser is the minimum of the PLR from all inputs or
    # outputs plus additional constraints and full load
    used_plr = min(minimum(x -> x, plr_from_nrg), max_plr, 1.0)

    # check total minimum PLR before dispatching units, which might have their own minimum
    # PLR which is typically different from the total
    if used_plr < unit.min_power_fraction_total
        return (false, [])
    end

    # we now have the PLR of the entire assembly and the limiting energy (which might be
    # the design power) and need to decide how to dispatch the units to meet the target
    # limiting energy
    nr_active, plr_of_unit = dispatch_units(unit,
                                            used_plr,
                                            limiting_interface,
                                            limiting_energy,
                                            sim_params["watt_to_wh"])

    # now the total energies can be calculated from the number and PLR of utilised units
    energies = []
    for name in unit.interface_list
        push!(energies, nr_active * energy_from_plr(unit, name, plr_of_unit, sim_params["watt_to_wh"]))
    end
    return (true, energies)
end

function potential(unit::Electrolyser, sim_params::Dict{String,Any})
    success, energies = calculate_energies(unit, sim_params)

    if !success || sum(energies[1]; init=0.0) < sim_params["epsilon"]
        set_max_energies!(unit, 0.0, 0.0, 0.0, 0.0, 0.0)
    else
        set_max_energies!(unit, energies[1], energies[2], energies[3], energies[4], energies[5])
    end
end

function process(unit::Electrolyser, sim_params::Dict{String,Any})
    success, energies = calculate_energies(unit, sim_params)

    if !success
        set_max_energies!(unit, 0.0, 0.0, 0.0, 0.0, 0.0)
        return
    end

    plr = energies[1] / sim_params["watt_to_wh"](unit.power_total * unit.efficiencies[Symbol("el_in")](1.0))
    h2_out_lossless = energies[1] * unit.efficiencies[Symbol("h2_out_lossless")](plr)
    unit.losses_hydrogen = h2_out_lossless - energies[4]
    unit.losses_heat = energies[1] - energies[2] +
                       (unit.heat_lt_is_usable ? -1 : 0) * energies[3] -
                       h2_out_lossless
    unit.losses_hydrogen = check_epsilon(unit.losses_hydrogen, sim_params)
    unit.losses_heat = check_epsilon(unit.losses_heat, sim_params)

    unit.losses = unit.losses_heat + unit.losses_hydrogen

    sub!(unit.input_interfaces[unit.m_el_in], energies[1])
    add!(unit.output_interfaces[unit.m_heat_ht_out], energies[2], nothing, unit.output_temperature_ht)
    if unit.heat_lt_is_usable
        add!(unit.output_interfaces[unit.m_heat_lt_out], energies[3], nothing, unit.output_temperature_lt)
    end
    add!(unit.output_interfaces[unit.m_h2_out], energies[4])
    add!(unit.output_interfaces[unit.m_o2_out], energies[5])
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
    channels = [string(unit.m_el_in) * " IN",
                string(unit.m_h2_out) * " OUT",
                string(unit.m_o2_out) * " OUT",
                string(unit.m_heat_ht_out) * " OUT",
                "Losses",
                "Losses_heat",
                "Losses_hydrogen"]

    if unit.heat_lt_is_usable
        append!(channels, [string(unit.m_heat_lt_out) * " OUT"])
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
        return unit.losses
    end
    throw(KeyError(key.value_key))
end

export Electrolyser
