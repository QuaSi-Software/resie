"""
Trait-like type of components that implement the load ratio dependent efficiency (PLRDE)
functionality.
"""
const PLRDEComponent = Union{CHPP, FuelBoiler}

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
- `medium::Symbol`: The medium to which the energy value corresponds. This should be one of
    the names defined in field `interface_list` of the component.
- `energy_value::Float64`: The energy value
# Returns
- `Float64`: The part load ratio (from 0.0 to 1.0) as inverse from the energy value
"""
function plr_from_energy(
    unit::PLRDEComponent,
    medium::Symbol,
    energy_value::Float64
)::Float64
    # shortcut if the given medium is the design power medium as it is linear to PLR
    if medium === unit.design_power_medium
        return energy_value / watt_to_wh(unit.power)
    end

    lookup_table = unit.energy_to_plr[medium]
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
