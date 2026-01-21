#! format: off
const ELECTROLYSER_PARAMETERS = Dict(
    "m_el_in" => (
        default="m_e_ac_230v",
        description="Electricity input medium",
        display_name="Medium el_in",
        required=false,
        type=String,
        json_type="string",
        unit="-"
    ),
    "m_h2_out" => (
        default="m_c_g_h2",
        description="Hydrogen output medium",
        display_name="Medium h2_out",
        required=false,
        type=String,
        json_type="string",
        unit="-"
    ),
    "m_o2_out" => (
        default="m_c_g_o2",
        description="Oxygen output medium",
        display_name="Medium o2_out",
        required=false,
        type=String,
        json_type="string",
        unit="-"
    ),
    "m_heat_ht_out" => (
        default="m_h_w_ht1",
        description="High-temperature heat output medium",
        display_name="Medium heat_ht_out",
        required=false,
        type=String,
        json_type="string",
        unit="-"
    ),
    "m_heat_lt_out" => (
        default="m_h_w_lt1",
        description="Low-temperature heat output medium",
        display_name="Medium heat_lt_out",
        required=false,
        conditionals=[("heat_lt_is_usable", "is_true")],
        type=String,
        json_type="string",
        unit="-"
    ),
    "heat_lt_is_usable" => (
        default=false,
        description="Toggle if the low temperature heat output is usable",
        display_name="LT-heat is usable?",
        required=false,
        type=Bool,
        json_type="boolean",
        unit="-"
    ),
    "linear_interface" => (
        default="el_in",
        description="Which interface is considered linear relative to the part-load-ratio",
        display_name="Linear interface",
        options=["el_in", "h2_out", "o2_out", "heat_ht_out", "heat_lt_out"],
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
    "efficiency_h2_out" => (
        default="const:0.57",
        description="Efficiency function for the hydrogen output including cycling and drying losses",
        display_name="Efficiency h2_out",
        required=false,
        type=String,
        json_type="string",
        function_type="1dim",
        unit="-"
    ),
    "efficiency_h2_out_lossless" => (
        default="const:0.6",
        description="Efficiency function for the hydrogen output without cycling and drying losses",
        display_name="Efficiency h2_out_lossless",
        required=false,
        type=String,
        json_type="string",
        function_type="1dim",
        unit="-"
    ),
    "efficiency_o2_out" => (
        default="const:0.6",
        description="Efficiency function for the oxygen output, conversion to energy via " *
                    "stochiometric hydrogen equivalent",
        display_name="Efficiency o2_out",
        required=false,
        type=String,
        json_type="string",
        function_type="1dim",
        unit="-"
    ),
    "efficiency_heat_ht_out" => (
        default="const:0.15",
        description="Efficiency function for the high-temperature heat output",
        display_name="Efficiency heat_ht_out",
        required=false,
        type=String,
        json_type="string",
        function_type="1dim",
        unit="-"
    ),
    "efficiency_heat_lt_out" => (
        default="const:0.07",
        description="Efficiency function for the low-temperature heat output",
        display_name="Efficiency heat_lt_out",
        required=false,
        conditionals=[("heat_lt_is_usable", "is_true")],
        type=String,
        json_type="string",
        function_type="1dim",
        unit="-"
    ),
    "power_el" => (
        description="Electric design power",
        display_name="Electric power",
        required=true,
        validations=[
            ("self", "value_gt_num", 0.0)
        ],
        type=Float64,
        json_type="number",
        unit="W"
    ),
    "nr_switchable_units" => (
        default=1,
        description="Number of switchable units that make up the assembly",
        display_name="Nr. units",
        required=false,
        validations=[
            ("self", "value_gte_num", 1.0)
        ],
        type=UInt,
        json_type="number",
        unit="-"
    ),
    "dispatch_strategy" => (
        default="equal_with_mpf",
        description="Strategy to use for dispatching the units of the assembly to achieve desired part-load ratio.",
        display_name="Dispatch strategy",
        required=false,
        options=("all_equal", "try_optimal", "equal_with_mpf"),
        type=String,
        json_type="string",
        unit="-"
    ),
    "optimal_unit_plr" => (
        default=0.65,
        description="Part-load ratio of a single unit when operation is optimal",
        display_name="Optimal unit PLR",
        required=false,
        conditionals=[("dispatch_strategy", "has_value", "try_optimal")],
        validations=[
            ("self", "value_gt_num", 0.0),
            ("self", "value_lte_num", 1.0)
        ],
        type=Float64,
        json_type="number",
        unit="-"
    ),
    "min_power_fraction" => (
        default=0.4,
        description="Minimum part-load ratio to operate of a single unit",
        display_name="Min. power fraction unit",
        required=false,
        validations=[
            ("self", "value_gte_num", 0.0),
            ("self", "value_lte_num", 1.0)
        ],
        type=Float64,
        json_type="number",
        unit="-"
    ),
    "min_power_fraction_total" => (
        default=0.2,
        description="Minimum part-load ratio to operate of the whole assembly",
        display_name="Min. power fraction total",
        required=false,
        validations=[
            ("self", "value_gte_num", 0.0),
            ("self", "value_lte_num", 1.0)
        ],
        type=Float64,
        json_type="number",
        unit="-"
    ),
    "output_temperature_ht" => (
        default=55.0,
        description="Fixed output temperature for high-temperature, or nothing for auto-detection",
        display_name="Output high-temperature",
        required=false,
        type=Floathing,
        json_type="number",
        unit="°C"
    ),
    "output_temperature_lt" => (
        default=25.0,
        description="Fixed output temperature for low-temperature, or nothing for auto-detection",
        display_name="Output low-temperature",
        required=false,
        conditionals=[("heat_lt_is_usable", "is_true")],
        type=Floathing,
        json_type="number",
        unit="°C"
    ),
    "nr_discretization_steps" => (
        default=1,
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
        return new(SSOT_parameter_constructor(Electrolyser, uac, config, sim_params)...)
    end
end

function component_parameters(x::Type{Electrolyser})::Dict{String,NamedTuple}
    return deepcopy(ELECTROLYSER_PARAMETERS) # Return a copy to prevent external modification
end

function extract_parameter(x::Type{Electrolyser}, config::Dict{String,Any}, param_name::String,
                           param_def::NamedTuple, sim_params::Dict{String,Any}, uac::String)
    return extract_parameter(Component, config, param_name, param_def, sim_params, uac)
end

function validate_config(x::Type{Electrolyser}, config::Dict{String,Any}, extracted::Dict{String,Any},
                         uac::String, sim_params::Dict{String,Any})
    validate_config(Component, extracted, uac, sim_params, component_parameters(Electrolyser))
end

function init_from_params(x::Type{Electrolyser}, uac::String, params::Dict{String,Any},
                          raw_params::Dict{String,Any}, sim_params::Dict{String,Any})::Tuple
    # turn media names into Symbol and register them
    heat_lt_is_usable = params["heat_lt_is_usable"]
    m_el_in = Symbol(params["m_el_in"])
    m_heat_ht_out = Symbol(params["m_heat_ht_out"])
    m_heat_lt_out = Symbol(params["m_heat_lt_out"])
    m_h2_out = Symbol(params["m_h2_out"])
    m_o2_out = Symbol(params["m_o2_out"])
    register_media([m_el_in, m_heat_ht_out, m_heat_lt_out, m_h2_out, m_o2_out])
    interface_list = (Symbol("el_in"),
                      Symbol("heat_ht_out"),
                      Symbol("heat_lt_out"),
                      Symbol("h2_out"),
                      Symbol("o2_out"))

    efficiencies = Dict{Symbol,Function}(
        Symbol("el_in") => params["efficiency_el_in"],
        Symbol("heat_ht_out") => params["efficiency_heat_ht_out"],
        Symbol("heat_lt_out") => params["efficiency_heat_lt_out"],
        Symbol("h2_out") => params["efficiency_h2_out"],
        Symbol("h2_out_lossless") => params["efficiency_h2_out_lossless"],
        Symbol("o2_out") => params["efficiency_o2_out"],
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

    nr_units = params["nr_switchable_units"]
    power_total = params["power_el"] / efficiencies[Symbol("el_in")](1.0)
    power = power_total / nr_units

    # return tuple in the order expected by new()
    return (uac,
            Controller(params["control_parameters"]),
            sf_transformer,
            InterfaceMap(m_el_in => nothing),
            output_interfaces,
            m_el_in,
            m_heat_ht_out,
            m_heat_lt_out,
            m_h2_out,
            m_o2_out,
            power,
            power_total,
            nr_units,
            params["dispatch_strategy"],
            params["optimal_unit_plr"],
            Symbol(params["linear_interface"]),
            params["min_power_fraction"],
            params["min_power_fraction_total"],
            efficiencies,
            interface_list,
            Dict{Symbol,Vector{Tuple{Float64,Float64}}}(), # energy_to_plr
            1.0 / params["nr_discretization_steps"], # discretization_step
            heat_lt_is_usable,
            params["output_temperature_ht"],
            params["output_temperature_lt"],
            0.0, # losses
            0.0, # losses_heat
            0.0) # losses_hydrogen
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
                           heat_ht_out::Union{Floathing,Vector{<:Floathing}},
                           heat_lt_out::Float64,
                           h2_out::Float64,
                           o2_out::Float64,
                           purpose_uac_ht::Union{Stringing,Vector{<:Stringing}}=nothing)
    set_max_energy!(unit.input_interfaces[unit.m_el_in], el_in)
    set_max_energy!(unit.output_interfaces[unit.m_heat_ht_out], heat_ht_out, nothing, unit.output_temperature_ht,
                    purpose_uac_ht)
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
    the electrolyser, in which case all units are utilized to their full extent.
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

function calculate_energies(unit::Electrolyser,
                            sim_params::Dict{String,Any})::Tuple{Bool,
                                                                 Vector{Union{Floathing,
                                                                              Vector{<:Floathing},
                                                                              Vector{<:Stringing}}}}
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

    # save the energies and purpose_uac for the heat ht output interfaces
    energies_heat_ht_out = Float64[]
    purpose_uac_ht_out = Stringing[]

    for name in unit.interface_list
        availability = getproperty(EnergySystems, Symbol("check_" * String(name)))
        if name == Symbol("heat_ht_out")
            energies_heat_ht_out, purpose_uac_ht_out = availability(unit, sim_params)
            energy = sum(energies_heat_ht_out; init=0.0)
        else
            energy = availability(unit, sim_params)
        end

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

    # now the total energies can be calculated from the number and PLR of utilized units
    energies = Union{Floathing,Vector{Floathing},Vector{Stringing}}[]
    purpose_uac = Stringing[]
    for name in unit.interface_list
        if name == Symbol("heat_ht_out")
            # split energy for the ht heat output interface into the different targets
            used_energy = nr_active * energy_from_plr(unit, name, plr_of_unit, sim_params["watt_to_wh"])

            energy_per_uac = Floathing[]
            while used_energy > 0.0 && !isempty(energies_heat_ht_out)
                current_energy = min(used_energy, abs(popfirst!(energies_heat_ht_out)))
                push!(energy_per_uac, current_energy)
                push!(purpose_uac, popfirst!(purpose_uac_ht_out))
                used_energy -= current_energy
            end
            if isempty(energy_per_uac)
                energy_per_uac = Floathing[0.0]
            end
            if isempty(purpose_uac)
                purpose_uac = Stringing[nothing]
            end
            push!(energies, energy_per_uac)
        else
            push!(energies, nr_active * energy_from_plr(unit, name, plr_of_unit, sim_params["watt_to_wh"]))
        end
    end
    # add purpose_uac of ht heat output at the end
    push!(energies, purpose_uac)

    return (true, energies)
end

function potential(unit::Electrolyser, sim_params::Dict{String,Any})
    success, energies = calculate_energies(unit, sim_params)

    if !success || sum(energies[1]; init=0.0) < sim_params["epsilon"]
        set_max_energies!(unit, 0.0, 0.0, 0.0, 0.0, 0.0)
    else
        set_max_energies!(unit, energies[1], energies[2], energies[3], energies[4], energies[5], energies[6])
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
    unit.losses_heat = energies[1] - sum(energies[2]; init=0.0) +
                       (unit.heat_lt_is_usable ? -1 : 0) * energies[3] -
                       h2_out_lossless
    unit.losses_hydrogen = check_epsilon(unit.losses_hydrogen, sim_params)
    unit.losses_heat = check_epsilon(unit.losses_heat, sim_params)

    unit.losses = unit.losses_heat + unit.losses_hydrogen

    sub!(unit.input_interfaces[unit.m_el_in], energies[1])
    add!(unit.output_interfaces[unit.m_heat_ht_out],
         energies[2],
         fill(nothing, length(energies[2])),
         fill(unit.output_temperature_ht, length(energies[2])),
         energies[6])
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

function component_has_minimum_part_load(unit::Electrolyser)
    return (unit.dispatch_strategy == "equal_with_mpf" && unit.min_power_fraction > 0.0) ||
           unit.min_power_fraction_total > 0.0
end

function output_values(unit::Electrolyser)::Vector{String}
    channels = [string(unit.m_el_in) * ":IN",
                string(unit.m_h2_out) * ":OUT",
                string(unit.m_o2_out) * ":OUT",
                string(unit.m_heat_ht_out) * ":OUT",
                "LossesGains",
                "Losses_heat",
                "Losses_hydrogen"]

    if unit.heat_lt_is_usable
        append!(channels, [string(unit.m_heat_lt_out) * ":OUT"])
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
        return -unit.losses_heat
    elseif key.value_key == "Losses_hydrogen"
        return -unit.losses_hydrogen
    elseif key.value_key == "LossesGains"
        return -unit.losses
    end
    throw(KeyError(key.value_key))
end

export Electrolyser
