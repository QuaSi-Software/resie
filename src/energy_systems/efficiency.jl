"""
Trait-like type of components that implement the load ratio dependent efficiency (PLRDE)
functionality.
"""
const PLRDEComponent = Union{CHPP, Electrolyser, FuelBoiler}

"""
    parse_efficiency_function(eff_def)

Parse the given definition of an efficiency function and return it as a callable function.

The function should return an efficiency factor e (on [0,1]) given an input of the part load
ratio (PLR) value (on [0,1]).

The definition looks like this: `<function_model>:<list_of_numbers>` with `function_model`
being a string (see below) and `list_of_numbers` being a comma-seperated list of numbers
with a period as decimal seperator and no thousands-seperator. The meaning of the numbers
depends on the function model.

Three different function models are implemented:
    * const: Takes one number and uses it as a constant efficiency factor.
    * poly: Takes a list of numbers as uses them as the coefficients of a polynomial with
        order n-1 where n is the length of coefficients. The list starts with coefficients
        of the highest order. E.g. `poly:0.5,2.0,0.1` means e(x)=0.5x²+2x+0.1
    * pwlin: A piece-wise linear interpolation. Takes a list of numbers and uses them as
        support values for an even distribution of linear sections on the interval [0,1].
        The PLR-values (on the x axis) are implicit with a step size equal to the inverse of
        the length of support values minus 1. The first and last support values are used as
        the values for a PLR of 0.0 and 1.0 respectively. E.g. `pwlin:0.6,0.8,0.9` means
        two sections of step size 0.5 with a value of e(0.0)==0.6, e(0.5)==0.8, e(1.0)=0.9
        and linear interpolation inbetween.

# Arguments
- `eff_def::String`: The efficiency definition as described above
# Returns
- `Function`: A callable function which returns an efficiency value when given a part load
    ratio value (from 0.0 to 1.0) as argument
"""
function parse_efficiency_function(eff_def::String)::Function
    splitted = split(eff_def, ":")

    if length(splitted) > 1
        method = lowercase(splitted[1])
        data = splitted[2]

        if method == "const"
            c = parse(Float64, data)
            return plr -> c

        elseif method == "poly"
            params = map(x -> parse(Float64,x), split(data, ","))
            return function(plr)
                return sum(p * plr^(length(params)-i) for (i, p) in enumerate(params))
            end

        elseif method == "pwlin"
            params = map(x -> parse(Float64,x), split(data, ","))
            step = 1.0 / (length(params)-1)
            return function(plr)
                bracket_nr = floor(Int64, plr / step) + 1
                lower_bound = params[bracket_nr]
                upper_bound = params[min(bracket_nr+1,length(params))]
                return lower_bound + (upper_bound-lower_bound)*(plr%step) / step
            end
        end
    end

    @warn "Cannot parse efficiency function from: $eff_def"
    return plr -> plr
end

"""
    create_plr_lookup_tables(unit)

Approximates the inverse of the PLRDE from energy values as lookup table.

This is required to calculate the part load ratio of a component from an energy value of
consumed or produced energy. Because the inverse is not always analytically known, it is
numerically approximated at fixed intervals of the PLR. The lookup table can then be
searched for the bounds around a known energy value and linear interpolation between these
bounds results in the approximation of the inverse function.

A warning is logged if the calculated inverse, as a piece-wise linear function, is not
monotonically increasing. This can happen if the efficiency function is not linear and has
strong gradients. Simulation can continue, but may result in inconsistent behaviour as
multiple solutions exist.

# Arguments
-`unit::PLRDEComponent`: The component for which to calculate the inverse
-`sim_params::Dict{String,Any}`: Simulation parameters
# Returns
-`Dict{Symbol,Vector{Tuple{Float64,Float64}}}`: The created lookup tables as dictionary with
    keys taken from the field `interface_list` of the component.
"""
function create_plr_lookup_tables(unit::PLRDEComponent, sim_params::Dict{String,Any})
    tables = Dict{Symbol,Vector{Tuple{Float64,Float64}}}()

    for name in unit.interface_list
        lookup_table = []
        for plr in collect(0.0:unit.discretization_step:1.0)
            push!(lookup_table, (
                watt_to_wh(unit.power) * plr * unit.efficiencies[name](plr),
                plr
            ))
        end
        tables[name] = lookup_table

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

    return tables
end


"""
    plr_from_energy(unit, energy_value)

Calculate the part load ratio as inverse function from the given energy value.

The PLRDE is defined such that one input or output is linear in respect to the PLR, while
the other inputs and outputs need to be calculated as the inverse of their efficiency
function, where that efficiency is relative to the linear input/output.

As the efficiency function can be a variety of functions and is not necessarily (easily)
invertable, the inverse is calculated numerically at initialisation as a piece-wise linear
function from a customizable number of support values from an even distribution of a PLR
from 0.0 to 1.0. When calling plr_from_energy this approximated function is evaluated for
the given energy value and a linear interpolation between the two surrounding support values
is performed to calculate the corresponding PLR.

# Arguments
- `unit::PLRDEComponent`: The component
- `interface::Symbol`: The interface to which the energy value corresponds. This should be
    one of the names defined in field `interface_list` of the component.
- `energy_value::Float64`: The energy value
# Returns
- `Float64`: The part load ratio (from 0.0 to 1.0) as inverse from the energy value
"""
function plr_from_energy(
    unit::PLRDEComponent,
    interface::Symbol,
    energy_value::Float64
)::Float64
    # shortcut if the given interface is designated as linear in respect to PLR
    if interface === unit.linear_interface
        return energy_value / watt_to_wh(unit.power)
    end

    lookup_table = unit.energy_to_plr[interface]
    energy_at_max = last(lookup_table)[1]

    if energy_value <= 0.0
        return 0.0
    elseif energy_value >= energy_at_max
        return 1.0
    end

    nr_iter = 0
    candidate_idx = max(
        1,
        floor(Int64, length(lookup_table) * energy_value / energy_at_max)
    )

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

    @warn "The energy_value on interface $(interface) in component $(unit.uac) is not " *
        "within the range of the lookup table."
    return 0.0
end

"""
    energy_from_plr(unit, interface, plr)

Calculates the energy value for the given interface from the given PLR.

# Arguments
- `unit::PLRDEComponent`: The component
- `interface::Symbol`: The interface to which the energy value corresponds. This should be
    one of the names defined in field `interface_list` of the component.
- `plr::Float64`: The part load ratio, should be between 0.0 and 1.0
# Returns
- `Float64`: The energy value
"""
function energy_from_plr(
    unit::PLRDEComponent,
    interface::Symbol,
    plr::Float64
)::Float64
    # the intuitive solution would be to multiply the power with the part load ratio and the
    # the efficiency at that PLR, but for unknown reasons this introduces an error in the
    # linear interpolation that appears to be quadratic in the distance to the support
    # values
    # @TODO: Investigate this further for improved performance
    # return plr * watt_to_wh(unit.power) * unit.efficiencies[interface](plr)

    # shortcut if the given interface is designated as linear in respect to PLR
    if interface === unit.linear_interface
        return plr * watt_to_wh(unit.power)
    end

    lookup_table = unit.energy_to_plr[interface]
    if plr >= 1.0
        return last(lookup_table)[1]
    end

    nr_iter = 0
    candidate_idx = max(1, floor(Int64, length(lookup_table) * plr))

    while (
        nr_iter < length(lookup_table)
        && candidate_idx < length(lookup_table)
        && candidate_idx >= 1
    )
        (energy_lb, plr_lb) = lookup_table[candidate_idx]
        (energy_ub, plr_ub) = lookup_table[candidate_idx+1]
        if plr_lb <= plr && plr < plr_ub
            return energy_lb + (energy_ub - energy_lb) * (plr - plr_lb) / (plr_ub - plr_lb)
        elseif plr < plr_lb
            candidate_idx -= 1
        elseif plr >= plr_ub
            candidate_idx += 1
        end

        nr_iter += 1
    end

    @warn "The PLR on interface $(interface) in component $(unit.uac) is not " *
        "within the range of the lookup table."
    return 0.0
end

"""
    calculate_energies_for_plrde(unit, sim_params, min_plr, max_plr)

Calculates the energy values for inputs and outputs at the determined operation point.

The operation point takes the available energy on the inputs and outputs into consideration
and chooses the part load ratio such that it fulfills the constraints through available
energies. It also considers additional given constraints and the design power of the
component.

# Arguments
- `unit::PLRDEComponent`: The component
- `sim_params::Dict{String,Any}`: Simulation parameters
- `min_plr::Float64`: Minimum PLR constraint
- `max_plr::Float64`: Maximum PLR contsraint
# Returns
- `Bool`: If the calculation was a success or if a constraint was not fulfilled
- `Vector{Floathing}`: The energy values in the same order as the field `interface_list` of
    the component
"""
function calculate_energies_for_plrde(
    unit::PLRDEComponent,
    sim_params::Dict{String,Any},
    min_plr::Float64,
    max_plr::Float64
)::Tuple{Bool, Vector{Floathing}}
    # calculate available energy and PLR as inverse from that for each input/output
    available_energy = []
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

        # limit to design power
        energy = min(
            watt_to_wh(unit.power) * unit.efficiencies[name](1.0),
            energy
        )

        push!(available_energy, energy)
        push!(plr_from_nrg, plr_from_energy(unit, name, energy))
    end

    # the operation point of the component is the minimum of the PLR from all inputs/outputs
    # plus additional constraints and the design power (if all available energy is infinite)
    used_plr = min(minimum(x->x, plr_from_nrg), max_plr, 1.0)

    if used_plr < min_plr
        return (false, [])
    end

    energies = []
    for name in unit.interface_list
        push!(energies, energy_from_plr(unit, name, used_plr))
    end
    return (true, energies)
end

"""
    check_fuel_in(unit, sim_params)

Checks the available energy on the input fuel interface.

# Arguments
- `unit::Union{CHPP,FuelBoiler}`: The component
- `sim_params::Dict{String,Any}`: Simulation parameters
# Returns
- `Floathing`: The available energy on the interface. If the value is nothing, that means
    no energy is available on this interface. The value can be `Inf`, which is a special
    floating point value signifying an infinite value
"""
function check_fuel_in(
    unit::Union{CHPP,FuelBoiler},
    sim_params::Dict{String,Any}
)
    if !unit.controller.base_module.parameters["consider_m_fuel_in"]
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

"""
    check_el_in(unit, sim_params)

Checks the available energy on the input electricity interface.

# Arguments
- `unit::Electrolyser`: The component
- `sim_params::Dict{String,Any}`: Simulation parameters
# Returns
- `Floathing`: The available energy on the interface. If the value is nothing, that means
    no energy is available on this interface. The value can be `Inf`, which is a special
    floating point value signifying an infinite value
"""
function check_el_in(
    unit::Electrolyser,
    sim_params::Dict{String,Any}
)
    if !unit.controller.base_module.parameters["consider_m_el_in"]
        return Inf
    end

    if (
        unit.input_interfaces[unit.m_el_in].source.sys_function == sf_transformer
        && unit.input_interfaces[unit.m_el_in].max_energy === nothing
    )
        return Inf
    else
        exchanges = balance_on(
            unit.input_interfaces[unit.m_el_in],
            unit.input_interfaces[unit.m_el_in].source
        )
        potential_energy_el = balance(exchanges) + energy_potential(exchanges)
        if potential_energy_el <= sim_params["epsilon"]
            return nothing
        end
        return potential_energy_el
    end
end

"""
    check_el_out(unit, sim_params)

Checks the available energy on the electricity output interface.

# Arguments
- `unit::CHPP`: The component
- `sim_params::Dict{String,Any}`: Simulation parameters
# Returns
- `Floathing`: The available energy on the interface. If the value is nothing, that means
    no energy is available on this interface. The value can be `-Inf`, which is a special
    floating point value signifying an infinite value
"""
function check_el_out(
    unit::CHPP,
    sim_params::Dict{String,Any}
)
    if !unit.controller.base_module.parameters["consider_m_el_out"]
        return -Inf
    end

    if (
        unit.output_interfaces[unit.m_el_out].target.sys_function == sf_transformer
        && unit.output_interfaces[unit.m_el_out].max_energy === nothing
    )
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

"""
    check_h2_out(unit, sim_params)

Checks the available energy on the hydrogen output interface.

# Arguments
- `unit::Electrolyser`: The component
- `sim_params::Dict{String,Any}`: Simulation parameters
# Returns
- `Floathing`: The available energy on the interface. If the value is nothing, that means
    no energy is available on this interface. The value can be `-Inf`, which is a special
    floating point value signifying an infinite value
"""
function check_h2_out(
    unit::Electrolyser,
    sim_params::Dict{String,Any}
)
    if !unit.controller.base_module.parameters["consider_m_h2_out"]
        return -Inf
    end

    if (
        unit.output_interfaces[unit.m_h2_out].target.sys_function == sf_transformer
        && unit.output_interfaces[unit.m_h2_out].max_energy === nothing
    )
        return -Inf
    end

    exchanges = balance_on(
        unit.output_interfaces[unit.m_h2_out],
        unit.output_interfaces[unit.m_h2_out].target
    )
    potential_energy_h2 = balance(exchanges) + energy_potential(exchanges)
    if potential_energy_h2 >= -sim_params["epsilon"]
        return nothing
    end

    return potential_energy_h2
end

"""
    check_o2_out(unit, sim_params)

Checks the available energy on the oxygen output interface.

# Arguments
- `unit::Electrolyser`: The component
- `sim_params::Dict{String,Any}`: Simulation parameters
# Returns
- `Floathing`: The available energy on the interface. If the value is nothing, that means
    no energy is available on this interface. The value can be `-Inf`, which is a special
    floating point value signifying an infinite value
"""
function check_o2_out(
    unit::Electrolyser,
    sim_params::Dict{String,Any}
)
    if !unit.controller.base_module.parameters["consider_m_o2_out"]
        return -Inf
    end

    if (
        unit.output_interfaces[unit.m_o2_out].target.sys_function == sf_transformer
        && unit.output_interfaces[unit.m_o2_out].max_energy === nothing
    )
        return -Inf
    end

    exchanges = balance_on(
        unit.output_interfaces[unit.m_o2_out],
        unit.output_interfaces[unit.m_o2_out].target
    )
    potential_energy_o2 = balance(exchanges) + energy_potential(exchanges)
    if potential_energy_o2 >= -sim_params["epsilon"]
        return nothing
    end

    return potential_energy_o2
end

"""
    check_heat_out_impl(unit, sim_params)

Checks the available energy on a heat output interface.

Please note that this method is intended to be used by various internal aliases to heat
outputs of different temperature regimes. The methods actually called are named after the
interfaces in the corresponding components. For example a CHPP would call check_heat_out,
which uses check_heat_out_impl internally.

This also checks the temperatures of the exchanges that are returned from a balance_on call,
which can be more than one exchange in the case of a bus, and checks if the output
temperature of the component falls into the minimum and maximum temperature range of the
exchange, if any is given at all.

# Arguments
- `unit::Union{CHPP,Electrolyser,FuelBoiler}`: The component
- `sim_params::Dict{String,Any}`: Simulation parameters
# Returns
- `Floathing`: The available energy on the interface. If the value is nothing, that means
    no energy is available on this interface. The value can be `-Inf`, which is a special
    floating point value signifying an infinite value
"""
function check_heat_out_impl(
    unit::Union{CHPP,Electrolyser,FuelBoiler},
    interface_name::String,
    medium::Symbol,
    output_temperature::Floathing,
    sim_params::Dict{String,Any}
)
    if !unit.controller.base_module.parameters["consider_m_"*interface_name]
        return -Inf
    end

    if (
        unit.output_interfaces[medium].target.sys_function == sf_transformer
        && unit.output_interfaces[medium].max_energy === nothing
    )
        return -Inf
    end

    exchanges = balance_on(
        unit.output_interfaces[medium],
        unit.output_interfaces[medium].target
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
            output_temperature === nothing
            || (e.temperature_min === nothing || e.temperature_min <= output_temperature)
            && (e.temperature_max === nothing || e.temperature_max >= output_temperature)
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
Alias to check_heat_out_impl for a heat output interface called heat_out.
"""
function check_heat_out(
    unit::Union{CHPP,FuelBoiler},
    sim_params::Dict{String,Any}
)
    return check_heat_out_impl(
        unit,
        "heat_out",
        unit.m_heat_out,
        unit.output_temperature,
        sim_params
    )
end

"""
Alias to check_heat_out_impl for a heat output interface called heat_ht_out.
"""
function check_heat_ht_out(
    unit::Electrolyser,
    sim_params::Dict{String,Any}
)
    return check_heat_out_impl(
        unit,
        "heat_ht_out",
        unit.m_heat_ht_out,
        unit.output_temperature_ht,
        sim_params
    )
end

"""
Alias to check_heat_out_impl for a heat output interface called heat_lt_out.
"""
function check_heat_lt_out(
    unit::Electrolyser,
    sim_params::Dict{String,Any}
)
    return check_heat_out_impl(
        unit,
        "heat_lt_out",
        unit.m_heat_lt_out,
        unit.output_temperature_lt,
        sim_params
    )
end