"""
Trait-like type of components that implement the load ratio dependent efficiency (PLRDE)
functionality.
"""
const PLRDEComponent = Union{CHPP,Electrolyser,FuelBoiler,UTIR}

"""
    parse_efficiency_function(eff_def)

Parse the given definition of an efficiency function and return it as a callable function.

The function should return an efficiency factor e (on [0,1]) given an input of the part load
ratio (PLR) value (on [0,1]).

The definition looks like this: `function_prototype:list_of_numbers` with
`function_prototype` being a string (see below) and `list_of_numbers` being a comma-
seperated list of numbers with a period as decimal seperator and no thousands-seperator. The
meaning of the numbers depends on the function prototype.

Three different function prototypes are implemented:
    * const: Takes one number and uses it as a constant efficiency factor.
    * poly: Takes a list of numbers and uses them as the coefficients of a polynomial with
        order n-1 where n is the length of coefficients. The list starts with coefficients
        of the highest order. E.g. `poly:0.5,2.0,0.1` means e(x)=0.5x²+2x+0.1
    * pwlin: A piece-wise linear interpolation. Takes a list of numbers and uses them as
        support values for an even distribution of linear sections on the interval [0,1].
        The PLR-values (on the x axis) are implicit with a step size equal to the inverse of
        the length of support values minus 1. The first and last support values are used as
        the values for a PLR of 0.0 and 1.0 respectively. E.g. `pwlin:0.6,0.8,0.9` means
        two sections of step size 0.5 with a value of e(0.0)==0.6, e(0.5)==0.8, e(1.0)=0.9
        and linear interpolation inbetween.
    * offset_lin: Takes one number and uses as the slope of a linear function with an offset
        of its complement (in regards to 1.0). E.g. `offset_lin:0.5` means
        e(x)=1.0-0.5*(1.0-x)
    * logarithmic: Takes two numbers and uses them as the coefficients of a
        quasi-logarithmic function. E.g. `logarithmic:0.5,0.3` means
        e(x)=0.5x/(0.3x+(1-0.3))
    * inv_poly: Takes a list of numbers and uses them as the coefficients of a polynomial
        with order n-1 where n is the length of coefficients. The list starts with
        coefficients of the highest order. The inverse of the polynomial multiplied with the
        PLR is used as the efficiency function. E.g. `inv_poly:0.5,2.0,0.1`
        means e(x)=x/(0.5x²+2x+0.1)
    * exp: Takes three numbers and uses them as the coefficients of an exponential function.
        E.g. `exp:0.1,0.2,0.3` means e(x)=0.1+0.2*exp(0.3x)
    * unified_plf: Takes four numbers and uses them as the coefficients of a composite
        function of a logarithmic and linear function as described in the documentation on
        the unified formulation for PLR-dependent efficiencies of heat pumps. The first two
        numbers are the optimal PLR and the PLF at that PLR. The third number is a scaling
        factor for the logarithmic part and the fourth number is the PLF at PLR=1.0.

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
            params = map(x -> parse(Float64, x), split(data, ","))
            return function (plr)
                return sum(p * plr^(length(params) - i) for (i, p) in enumerate(params))
            end

        elseif method == "pwlin"
            params = map(x -> parse(Float64, x), split(data, ","))
            step = 1.0 / (length(params) - 1)
            return function (plr)
                bracket_nr = floor(Int64, plr / step) + 1
                lower_bound = params[bracket_nr]
                upper_bound = params[min(bracket_nr + 1, length(params))]
                return lower_bound + (upper_bound - lower_bound) * (plr % step) / step
            end

        elseif method == "offset_lin"
            c = parse(Float64, data)
            return plr -> 1.0 - c * (1.0 - plr)

        elseif method == "logarithmic"
            params = map(x -> parse(Float64, x), split(data, ","))
            return plr -> params[1] * plr / (params[2] * plr + (1 - params[2]))

        elseif method == "inv_poly"
            params = map(x -> parse(Float64, x), split(data, ","))
            return plr -> plr / sum(p * plr^(length(params) - i) for (i, p) in enumerate(params))

        elseif method == "exp"
            params = map(x -> parse(Float64, x), split(data, ","))
            return plr -> params[1] + params[2] * exp(params[3] * plr)

        elseif method == "unified_plf"
            params = map(x -> parse(Float64, x), split(data, ","))
            return function (plr)
                if plr < params[1]
                    a = params[2] * (params[3] * (params[1] - 1) + 1) / params[1]
                    return a * plr / (params[3] * plr + 1 - params[3])
                else
                    return (params[4] - params[2]) * (plr - params[1]) / (1 - params[1]) + params[2]
                end
            end
        end
    end

    @error "Cannot parse efficiency function from: $eff_def"
    return plr -> plr
end

"""
    parse_2dim_function(eff_def)

Parse the given definition of an general 2D function and return it as a callable function.

The definition looks like this: `function_prototype:list_of_numbers` with
`function_prototype` being a string (see below) and `list_of_numbers` being a comma-
seperated list of numbers with a period as decimal seperator and no thousands-seperator. The
meaning of the numbers depends on the function prototype.

Two different function prototypes are implemented:
    * const: Takes one number and uses it as a constant value.
    * poly-2: Takes a list of numbers and uses them as the coefficients of a 2D polynomial
        with order three. The coefficients are the numbered constants in the following
        formula: f(x,y) = c_1 + c_2*x + c_3*y + c_4*x² + c_5*x*y + c_6*y²
            + c_7*x³ + c_8*x²*y + c_9*x*y² + c_10*y³

# Arguments
- `eff_def::String`: The function definition as described above
# Returns
- `Function`: A callable function which returns a scalar value given values (usually
    temperatures) as input.
"""
function parse_2dim_function(eff_def::String)::Function
    splitted = split(eff_def, ":")

    if length(splitted) > 1
        method = lowercase(splitted[1])
        data = splitted[2]

        if method == "const"
            c = parse(Float64, data)
            return (x, y) -> c

        elseif method == "poly-2"
            params = parse.(Float64, split(data, ","))
            return function (x, y)
                return params[1] +
                       params[2] * x +
                       params[3] * y +
                       params[4] * x * x +
                       params[5] * x * y +
                       params[6] * y * y +
                       params[7] * x * x * x +
                       params[8] * x * x * y +
                       params[9] * x * y * y +
                       params[10] * y * y * y
            end
        end
    end

    @error "Cannot parse 2-dimensional function from: $eff_def"
    return (x, y) -> 0.0
end

"""
Performs bilinear interpolation between four support values.

The interpolation is on the surface of two triangles with a common diagonal between the two
points (x1,y1) and (x3,y3). The four values are, in order, on the points (x1,y1), (x3,y1),
(x1,y3) and (x3,y3). The point (x2,y2) is the point at which interpolation happens.

# Arguments
- `x1::Float64`: x-coordinate of points of v1 and v3
- `x2::Float64`: x-coordinate to interpolate
- `x3::Float64`: x-coordinate of points of v2 and v4
- `y1::Float64`: y-coordinate of points of v1 and v2
- `y2::Float64`: y-coordinate to interpolate
- `y3::Float64`: y-coordinate of points of v3 and v4
- `v1::Float64`: value at (x1,y1)
- `v2::Float64`: value at (x3,y1)
- `v3::Float64`: value at (x1,y3)
- `v4::Float64`: value at (x3,y3)
# Returns
- `Float64`: Interpolated value at (x2,y2)
"""
function bilinear_interpolate(x1, x2, x3, y1, y2, y3, v1, v2, v3, v4)
    lin_x = (x2 - x1) / (x3 - x1)
    lin_y = (y2 - y1) / (y3 - y1)
    # value in plane spanned by v1, v3 and v4
    u = v1 + lin_x * (v4 - v3)
    v = v3 + lin_x * (v4 - v3)
    val_134 = u + lin_y * (v - u)
    # value in plane spanned by v1, v2 and v4
    u = v1 + lin_x * (v2 - v1)
    v = v4 - (1.0 - lin_x) * (v2 - v1)
    val_124 = u + lin_y * (v - u)
    # if the crease is convex, we take the min, otherwise the max
    if v2 < v1 + (v4 - v3)
        return min(val_134, val_124)
    else
        return max(val_134, val_124)
    end
end

"""
    parse_cop_function(eff_def)

Parse the given definition of a COP function and return it as a callable function.

The definition looks like this: `function_prototype:list_of_numbers` with
`function_prototype` being a string (see below) and `list_of_numbers` being a comma-
seperated list of numbers with a period as decimal seperator and no thousands-seperator and
semicolon as row seperator. The meaning of the numbers depends on the function prototype.

Four different function prototypes are implemented:
    * const: Takes one number and uses it as a constant COP.
    * carnot: Takes one number as uses it as the scaling factor for the Carnot-COP.
    * poly-2: Uses the same functionality for 2D polynomials of order 3 as described in
        function `parse_2dim_function`.
    * field: Takes a 2D array of values and performs bilinear interpolation between the
        surrounding four support values. An example with additional line breaks and spaces
        added for clarity:
        "field:
         0, 0,10,20,30;
         0,15, 9, 6, 4;
        10,15,15,10, 7;
        20,15,15,15,11"
        The first row are the grid points along the T_sink_out dimension, with the first
        value being ignored. The points cover a range from 0 °C to 30 °C with a spacing of
        10 K. The first column are the grid points along the T_source_in dimension, with the
        first value being ignored. The points cover a range from 0 °C to 20 °C with a
        spacing of 10 K. The support values are the COP (at PLR=1), for example a value of
        10 for T_source_in=10 and T_sink_out=20.

# Arguments
- `eff_def::String`: The function definition as described above
# Returns
- `Function`: A callable function which, when called with the input and output temperatures
    as arguments, returns the COP value at a PLR of 1.0.
"""
function parse_cop_function(eff_def::String)::Function
    splitted = split(eff_def, ":")

    if length(splitted) > 1
        method_cop = lowercase(splitted[1])
        data_cop = splitted[2]

        if method_cop == "const"
            c = parse(Float64, data_cop)
            return function (src, snk)
                return c
            end

        elseif method_cop == "carnot"
            c = parse(Float64, data_cop)
            return function (src, snk)
                if src === nothing || snk === nothing
                    return nothing
                end
                return c * (273.15 + snk) / (snk - src)
            end

        elseif method_cop == "poly-2"
            poly = parse_2dim_function(method_cop * ":" * data_cop)
            return poly

        elseif method_cop == "field"
            rows = split(data_cop, ';')
            cells = [parse.(Float64, split(row, ",")) for row in rows]
            dim_src = length(cells)
            dim_snk = length(cells[1])
            values = Array{Float64,2}(undef, dim_src, dim_snk)
            for (idx_src, row) in pairs(cells)
                for (idx_snk, val) in pairs(row)
                    values[idx_src, idx_snk] = val
                end
            end

            # first dim of values is source, second is sink. source is vertical (one row is
            # at the same source temp), sink is horizontal (one column is at the same sink
            # temp). first row is sink temperatures, first column is source temperatures
            return function (src, snk)
                src_idx = nothing
                for idx in 2:(dim_src - 1)
                    if src >= values[idx, 1] && src < values[idx + 1, 1]
                        src_idx = idx
                    end
                end
                snk_idx = nothing
                for idx in 2:(dim_snk - 1)
                    if snk >= values[1, idx] && snk < values[1, idx + 1]
                        snk_idx = idx
                    end
                end
                if snk_idx === nothing || src_idx === nothing
                    @error "Given temperatures $src and $snk outside of COP field."
                    throw(BoundsError(values, (src_idx, snk_idx)))
                end
                return bilinear_interpolate(values[1, snk_idx],
                                            snk,
                                            values[1, snk_idx + 1],
                                            values[src_idx, 1],
                                            src,
                                            values[src_idx + 1, 1],
                                            values[src_idx, snk_idx],
                                            values[src_idx, snk_idx + 1],
                                            values[src_idx + 1, snk_idx],
                                            values[src_idx + 1, snk_idx + 1])
            end
        end
    end

    @error "Cannot parse COP function from: $eff_def"
    return (x, y) -> 0.0
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
        for plr in collect(0.0:(unit.discretization_step):1.0)
            push!(lookup_table, (sim_params["watt_to_wh"](unit.power) * plr *
                                 unit.efficiencies[name](plr),
                                 plr))
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
- `w2wh::Function`: Function to convert power to work
# Returns
- `Float64`: The part load ratio (from 0.0 to 1.0) as inverse from the energy value
"""
function plr_from_energy(unit::PLRDEComponent,
                         interface::Symbol,
                         energy_value::Float64,
                         w2wh::Function)::Float64
    # shortcut if the given interface is designated as linear in respect to PLR
    if interface === unit.linear_interface
        return energy_value / w2wh(unit.power)
    end

    lookup_table = unit.energy_to_plr[interface]
    energy_at_max = last(lookup_table)[1]

    if energy_value <= 0.0
        return 0.0
    elseif energy_value >= energy_at_max
        return 1.0
    end

    nr_iter = 0
    candidate_idx = min(length(lookup_table) - 1,
                        max(1, floor(Int64, length(lookup_table) * energy_value / energy_at_max)))

    while (nr_iter < length(lookup_table)
           && candidate_idx < length(lookup_table)
           && candidate_idx >= 1)
        (energy_lb, plr_lb) = lookup_table[candidate_idx]
        (energy_ub, plr_ub) = lookup_table[candidate_idx + 1]
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
- `w2wh::Function`: Function to convert power to work
# Returns
- `Float64`: The energy value
"""
function energy_from_plr(unit::PLRDEComponent,
                         interface::Symbol,
                         plr::Float64,
                         w2wh::Function)::Float64
    # the intuitive solution would be to multiply the power with the part load ratio and the
    # the efficiency at that PLR, but for unknown reasons this introduces an error in the
    # linear interpolation that appears to be quadratic in the distance to the support
    # values
    # @TODO: Investigate this further for improved performance
    # return plr * w2wh(unit.power) * unit.efficiencies[interface](plr)

    # shortcut if the given interface is designated as linear in respect to PLR
    if interface === unit.linear_interface
        return plr * w2wh(unit.power)
    end

    lookup_table = unit.energy_to_plr[interface]
    if plr >= 1.0
        return last(lookup_table)[1]
    end

    nr_iter = 0
    candidate_idx = max(1, floor(Int64, length(lookup_table) * plr))

    while (nr_iter < length(lookup_table)
           && candidate_idx < length(lookup_table)
           && candidate_idx >= 1)
        (energy_lb, plr_lb) = lookup_table[candidate_idx]
        (energy_ub, plr_ub) = lookup_table[candidate_idx + 1]
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
function calculate_energies_for_plrde(unit::PLRDEComponent,
                                      sim_params::Dict{String,Any},
                                      min_plr::Float64,
                                      max_plr::Float64)::Tuple{Bool,Vector{Floathing}}
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
        energy = min(sim_params["watt_to_wh"](unit.power) * unit.efficiencies[name](1.0),
                     energy)

        push!(available_energy, energy)
        push!(plr_from_nrg, plr_from_energy(unit, name, energy, sim_params["watt_to_wh"]))
    end

    # the operation point of the component is the minimum of the PLR from all inputs/outputs
    # plus additional constraints and the design power (if all available energy is infinite)
    used_plr = min(minimum(x -> x, plr_from_nrg), max_plr, 1.0)

    if used_plr < min_plr
        return (false, [])
    end

    energies = []
    for name in unit.interface_list
        push!(energies, energy_from_plr(unit, name, used_plr, sim_params["watt_to_wh"]))
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
function check_fuel_in(unit::Union{CHPP,FuelBoiler}, sim_params::Dict{String,Any})
    if !unit.controller.parameters["consider_m_fuel_in"]
        return Inf
    end

    if (unit.input_interfaces[unit.m_fuel_in].source.sys_function == sf_transformer
        &&
        is_max_energy_nothing(unit.input_interfaces[unit.m_fuel_in].max_energy))
        # end of condition
        return Inf
    else
        exchanges = balance_on(unit.input_interfaces[unit.m_fuel_in],
                               unit.input_interfaces[unit.m_fuel_in].source)
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
- `unit::Union{Electrolyser,HeatPump,UTIR}`: The component
- `sim_params::Dict{String,Any}`: Simulation parameters
# Returns
- `Floathing`: The available energy on the interface. If the value is nothing, that means
    no energy is available on this interface. The value can be `Inf`, which is a special
    floating point value signifying an infinite value
"""
function check_el_in(unit::Union{Electrolyser,HeatPump,UTIR},
                     sim_params::Dict{String,Any})
    if !unit.controller.parameters["consider_m_el_in"]
        return Inf
    end

    if (unit.input_interfaces[unit.m_el_in].source.sys_function == sf_transformer
        &&
        is_max_energy_nothing(unit.input_interfaces[unit.m_el_in].max_energy))
        # end of condition
        return Inf
    else
        exchanges = balance_on(unit.input_interfaces[unit.m_el_in],
                               unit.input_interfaces[unit.m_el_in].source)
        potential_energy_el = balance(exchanges) + energy_potential(exchanges)
        if potential_energy_el <= sim_params["epsilon"]
            return nothing
        end
        return potential_energy_el
    end
end

"""
    check_el_in_layered(unit, sim_params)

Checks the available energy on the input electricity interface and returns energies and The
uacs of the sources as vectors.

# Arguments
- `unit::Union{Electrolyser,HeatPump,UTIR}`: The component
- `sim_params::Dict{String,Any}`: Simulation parameters
# Returns
- `Vector{Floathing}`: The available energy on the interface. If the value is nothing, that means
    no energy is available on this interface. The value can be `Inf`, which is a special
    floating point value signifying an infinite value
- `Vector{Stringing}`: The UACs of the sources on the interface.
"""
function check_el_in_layered(unit::Union{Electrolyser,HeatPump,UTIR},
                             sim_params::Dict{String,Any})
    if (unit.input_interfaces[unit.m_el_in].source.sys_function == sf_transformer
        &&
        is_max_energy_nothing(unit.input_interfaces[unit.m_el_in].max_energy))
        # direct connection to transformer that has not had its potential
        return ([Inf],
                [unit.input_interfaces[unit.m_el_in].source.uac])
    else
        exchanges = balance_on(unit.input_interfaces[unit.m_el_in],
                               unit.input_interfaces[unit.m_el_in].source)
        if unit.controller.parameters["consider_m_el_in"]
            return ([e.balance + e.energy_potential for e in exchanges],
                    [e.purpose_uac for e in exchanges])
        else # ignore the energy of the heat input, but purpose_uac is still required
            return ([Inf for _ in exchanges],
                    [e.purpose_uac for e in exchanges])
        end
    end
end

"""
    check_heat_in_layered(unit, sim_params)

Checks the available energy on the input heat interface.

# Arguments
- `unit::HeatPump`: The component
- `sim_params::Dict{String,Any}`: Simulation parameters
# Returns
- `Vector{Floathing}`: The available energy on the interface as one layer per source. The
    value can be `Inf`, which is a special floating point value signifying an infinite value
- `Vector{Temperature}`: The minimum temperatures on the interface as one layer per source.
- `Vector{Temperature}`: The maximum temperatures on the interface as one layer per source.
- `Vector{Stringing}`: The UACs of the sources on the interface.
"""
function check_heat_in_layered(unit::HeatPump, sim_params::Dict{String,Any})
    if (unit.input_interfaces[unit.m_heat_in].source.sys_function == sf_transformer
        &&
        is_max_energy_nothing(unit.input_interfaces[unit.m_heat_in].max_energy))
        # direct connection to transformer that has not had its potential
        return ([Inf],
                [unit.input_interfaces[unit.m_heat_in].max_energy.temperature_min[1]],
                [unit.input_interfaces[unit.m_heat_in].max_energy.temperature_max[1]],
                [unit.input_interfaces[unit.m_heat_in].source.uac])
    else
        exchanges = balance_on(unit.input_interfaces[unit.m_heat_in],
                               unit.input_interfaces[unit.m_heat_in].source)

        if unit.controller.parameters["consider_m_heat_in"]
            return ([e.balance + e.energy_potential for e in exchanges],
                    temp_min_all(exchanges),
                    temp_max_all(exchanges),
                    [e.purpose_uac for e in exchanges])
        else # ignore the energy of the heat input, but temperature is still required
            return ([Inf for _ in exchanges],
                    temp_min_all(exchanges),
                    temp_max_all(exchanges),
                    [e.purpose_uac for e in exchanges])
        end
    end
end

"""
    check_el_out(unit, sim_params)

Checks the available energy on the electricity output interface.

# Arguments
- `unit::Union{CHPP,UTIR}`: The component
- `sim_params::Dict{String,Any}`: Simulation parameters
# Returns
- `Floathing`: The available energy on the interface. If the value is nothing, that means
    no energy is available on this interface. The value can be `-Inf`, which is a special
    floating point value signifying an infinite value
"""
function check_el_out(unit::Union{CHPP,UTIR},
                      sim_params::Dict{String,Any})
    if !unit.controller.parameters["consider_m_el_out"]
        return -Inf
    end

    if (unit.output_interfaces[unit.m_el_out].target.sys_function == sf_transformer
        &&
        is_max_energy_nothing(unit.output_interfaces[unit.m_el_out].max_energy))
        # end of condition
        return -Inf
    end

    exchanges = balance_on(unit.output_interfaces[unit.m_el_out],
                           unit.output_interfaces[unit.m_el_out].target)
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
function check_h2_out(unit::Electrolyser,
                      sim_params::Dict{String,Any})
    if !unit.controller.parameters["consider_m_h2_out"]
        return -Inf
    end

    if (unit.output_interfaces[unit.m_h2_out].target.sys_function == sf_transformer
        &&
        is_max_energy_nothing(unit.output_interfaces[unit.m_h2_out].max_energy))
        # end of condition
        return -Inf
    end

    exchanges = balance_on(unit.output_interfaces[unit.m_h2_out],
                           unit.output_interfaces[unit.m_h2_out].target)
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
function check_o2_out(unit::Electrolyser,
                      sim_params::Dict{String,Any})
    if !unit.controller.parameters["consider_m_o2_out"]
        return -Inf
    end

    if (unit.output_interfaces[unit.m_o2_out].target.sys_function == sf_transformer
        &&
        is_max_energy_nothing(unit.output_interfaces[unit.m_o2_out].max_energy))
        # end of condition
        return -Inf
    end

    exchanges = balance_on(unit.output_interfaces[unit.m_o2_out],
                           unit.output_interfaces[unit.m_o2_out].target)
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
which can be more than one exchange, and checks if the output temperature of the component
falls into the minimum and maximum temperature range of the exchange, if any is given at all.

# Arguments
- `unit::Union{CHPP,Electrolyser,FuelBoiler}`: The component
- `sim_params::Dict{String,Any}`: Simulation parameters
# Returns
- `Floathing`: The available energy on the interface. If the value is nothing, that means
    no energy is available on this interface. The value can be `-Inf`, which is a special
    floating point value signifying an infinite value
"""
function check_heat_out_impl(unit::Union{CHPP,Electrolyser,FuelBoiler},
                             interface_name::String,
                             medium::Symbol,
                             output_temperature::Floathing,
                             sim_params::Dict{String,Any})
    if !unit.controller.parameters["consider_m_" * interface_name]
        return -Inf
    end

    if (unit.output_interfaces[medium].target.sys_function == sf_transformer
        &&
        is_max_energy_nothing(unit.output_interfaces[medium].max_energy))
        # end of condition
        return -Inf
    end

    exchanges = balance_on(unit.output_interfaces[medium],
                           unit.output_interfaces[medium].target)

    # if we get the exchanges from a bus, the temperature check has already been performed
    if unit.output_interfaces[medium].target.sys_function == sf_bus
        potential_energy_heat_out = balance(exchanges) + energy_potential(exchanges)
    else # check temperature for 1-to-1 connections and sum up energy
        potential_energy_heat_out, _, _, _ = check_temperatures_source(exchanges, output_temperature, Inf)
        potential_energy_heat_out = sum(potential_energy_heat_out; init=0.0)
    end

    if potential_energy_heat_out >= -sim_params["epsilon"]
        return nothing
    end
    return potential_energy_heat_out
end

"""
    check_heat_out_layered(unit, interface, medium, output_temperature, sim_params)

Checks the available energy on the output high temperature heat interface of the electrolyser.
Includes the required check for return temperatures for the control module "limit_cooling_input_temperature".
Therefore, it returns not only the energies but also the corresponding target_uac.

# Arguments
- `unit::Electrolyser`: The component
- `interface_name::String`: The name of the interface, e.g. "heat_ht_out"
- `medium::Symbol`: The medium of the interface, e.g. `:m_heat_ht_out`
- `output_temperature::Floathing`: The output temperature of the electrolyser for 1-to-1 temperature checks
- `sim_params::Dict{String,Any}`: Simulation parameters
# Returns
- `Vector{Floathing}`: The requested energy on the interface as one layer per target. The
    value can be `-Inf`, which is a special floating point value signifying a negative
    infinite value.
- `Vector{Stringing}`: The UACs of the targets on the interface corresponding to the energies.
"""
function check_heat_out_layered(unit::Electrolyser,
                                interface_name::String,
                                medium::Symbol,
                                output_temperature::Floathing,
                                sim_params::Dict{String,Any})
    if !unit.controller.parameters["consider_m_" * interface_name]
        # ignore the energy of the heat output, but target_uac is still required
        return ([-Inf for _ in exchanges],
                [e.purpose_uac for e in exchanges])
    end

    if (unit.output_interfaces[medium].target.sys_function == sf_transformer
        &&
        is_max_energy_nothing(unit.output_interfaces[medium].max_energy))
        # direct connection to transformer that has not had its potential
        # no temperature check required
        energy = [-Inf]
        purpose_uac = [unit.output_interfaces[medium].target.uac]
    else
        exchanges = balance_on(unit.output_interfaces[medium],
                               unit.output_interfaces[medium].target)

        if unit.output_interfaces[medium].target.sys_function == sf_bus
            # if we get the exchanges from a bus, the temperature check has already been performed
            energy = [e.balance + e.energy_potential for e in exchanges]
            purpose_uac = [e.purpose_uac for e in exchanges]
        else
            # check temperature for 1-to-1 connections
            energy, _, _, purpose_uac = check_temperatures_source(exchanges, output_temperature, Inf)
        end
    end

    # If the source is an Electrolyser, call the ControlModule limit_cooling_input_temperature
    # for every target to check if the energy transfer is allowed or not. 
    # If the temperature input limit is exceeded (control module returns true), remove the entry from returns.
    idx_to_remove = Int[]
    for idx in eachindex(purpose_uac)
        if cooling_input_temperature_exceeded(unit.controller, purpose_uac[idx], sim_params)
            push!(idx_to_remove, idx)
        end
    end

    if !isempty(idx_to_remove)
        deleteat!(energy, idx_to_remove)
        deleteat!(purpose_uac, idx_to_remove)
    end

    return (energy, purpose_uac)
end

"""
Alias to check_heat_out_impl for a heat output interface called heat_out.
"""
function check_heat_out(unit::Union{CHPP,FuelBoiler},
                        sim_params::Dict{String,Any})
    return check_heat_out_impl(unit,
                               "heat_out",
                               unit.m_heat_out,
                               unit.output_temperature,
                               sim_params)
end

"""
Alias to check_heat_out_impl for a heat output interface called heat_ht_out.
"""
function check_heat_ht_out(unit::Electrolyser,
                           sim_params::Dict{String,Any})
    return check_heat_out_layered(unit,
                                  "heat_ht_out",
                                  unit.m_heat_ht_out,
                                  unit.output_temperature_ht,
                                  sim_params)
end

"""
Alias to check_heat_out_impl for a heat output interface called heat_lt_out.
"""
function check_heat_lt_out(unit::Electrolyser,
                           sim_params::Dict{String,Any})
    return check_heat_out_impl(unit,
                               "heat_lt_out",
                               unit.m_heat_lt_out,
                               unit.output_temperature_lt,
                               sim_params)
end

"""
    check_heat_out_layered(unit::HeatPump, sim_params)

Checks the available energy on the output heat interface.

This specific implementation also considers potentially available secondary output interfaces on the heat pump. 
The function calls balance_on() on both output interfaces and merges the returned exchanges, depending on the
actual input priority of the primary and secondary input interface at the target bus.

# Arguments
- `unit::HeatPump`: The component
- `sim_params::Dict{String,Any}`: Simulation parameters
# Returns
- `Vector{Floathing}`: The requested energy on the interface as one layer per target. The
    value can be `-Inf`, which is a special floating point value signifying a negative
    infinite value.
- `Vector{Temperature}`: The minimum temperatures on the interface as one layer per target.
- `Vector{Temperature}`: The maximum temperatures on the interface as one layer per target.
- `Vector{Stringing}`: The UACs of the targets on the interface.
"""
function check_heat_out_layered(unit::HeatPump, sim_params::Dict{String,Any})
    if (unit.output_interfaces[unit.m_heat_out].target.sys_function == sf_transformer
        &&
        is_max_energy_nothing(unit.output_interfaces[unit.m_heat_out].max_energy)) ||
       (unit.has_secondary_interface &&
        unit.output_interfaces[unit.m_heat_out_secondary].target.sys_function == sf_transformer &&
        is_max_energy_nothing(unit.output_interfaces[unit.m_heat_out_secondary].max_energy))
        # direct connection to transformer that has not had its potential
        return ([-Inf],
                [unit.output_interfaces[unit.m_heat_out].max_energy.temperature_min[1]],
                [unit.output_interfaces[unit.m_heat_out].max_energy.temperature_max[1]],
                [unit.output_interfaces[unit.m_heat_out].target.uac])
    else
        exchanges = balance_on(unit.output_interfaces[unit.m_heat_out],
                               unit.output_interfaces[unit.m_heat_out].target)
        if unit.controller.parameters["consider_m_heat_out"]
            if unit.has_secondary_interface
                # If a secondary interface exists, call balance_on on this interface as well...
                exchanges_secondary = balance_on(unit.output_interfaces[unit.m_heat_out_secondary],
                                                 unit.output_interfaces[unit.m_heat_out_secondary].target)
                # ...and merge them together depending on the dynamic order of them at the target bus
                if unit.output_interfaces[unit.m_heat_out_secondary].target.sys_function === sf_bus &&
                   unit.output_interfaces[unit.m_heat_out].target.sys_function === sf_bus
                    secondary_prio = unit.output_interfaces[unit.m_heat_out_secondary].target.balance_table_inputs[adjust_name_if_secondary(unit.uac,
                                                                                                                                            true)].priority
                    regular_prio = unit.output_interfaces[unit.m_heat_out].target.balance_table_inputs[unit.uac].priority
                else
                    # if no bus is present, set secondary interface to second priority
                    secondary_prio = 2
                    regular_prio = 1
                end
                if secondary_prio < regular_prio
                    # secondary interface has higher priority (smaller number)
                    exchanges_first = exchanges_secondary
                    exchanges_second = exchanges
                else
                    # secondary interface has lower priority (higher number)                    
                    exchanges_first = exchanges
                    exchanges_second = exchanges_secondary
                end

                for exchange_second in exchanges_second
                    duplicate = false
                    for exchange_first in exchanges_first
                        duplicate = exchange_second.balance == exchange_first.balance &&
                                    exchange_second.energy_potential == exchange_first.energy_potential &&
                                    exchange_second.purpose_uac == exchange_first.purpose_uac &&
                                    exchange_second.temperature_min == exchange_first.temperature_min &&
                                    exchange_second.temperature_min == exchange_first.temperature_min
                        if duplicate
                            break
                        end
                    end
                    if !duplicate
                        push!(exchanges_first, exchange_second)
                    end
                end
            else
                exchanges_first = exchanges
            end
            return ([e.balance + e.energy_potential for e in exchanges_first],
                    temp_min_all(exchanges_first),
                    temp_max_all(exchanges_first),
                    [e.purpose_uac for e in exchanges_first])
        else # ignore the energy of the heat output, but temperature is still required
            return ([Inf for _ in exchanges],
                    temp_min_all(exchanges),
                    temp_max_all(exchanges),
                    [e.purpose_uac for e in exchanges])
        end
    end
end
