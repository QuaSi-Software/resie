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
