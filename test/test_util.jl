"""
    pwc_units_astr(expected, actual)

Piece-wise comparison of components noting differences in the return string.
"""
function pwc_units_astr(expected::Vector, actual::Vector)
    max_length = max(length(expected), length(actual))
    result = ""

    for i in 1:max_length
        if i > length(expected)
            result = result * "#$i: nothing != $(actual[i].uac) \n"
        elseif i > length(actual)
            result = result * "#$i: $(expected[i].uac) != nothing\n"
        elseif expected[i].uac != actual[i].uac
            result = result * "#$i: $(expected[i].uac) != $(actual[i].uac)\n"
        end
    end

    return result
end

"""
    pwc_steps_astr(expected, actual)

Piece-wise comparison of step instructions noting differences in the return string.
This function takes vecors of tuples (number, string).
"""
function pwc_steps_astr(expected::Vector, actual::Vector)
    max_length = max(length(expected), length(actual))
    result = ""

    for i in 1:max_length
        if i > length(expected)
            result = result * "#$i: nothing !=  $(actual[i][1])-$(actual[i][2])\n"
        elseif i > length(actual)
            result = result * "#$i: $(expected[i][1])-$(expected[i][2]) != nothing\n"
        elseif expected[i][1] != actual[i][1] || expected[i][2] != actual[i][2]
            result = result * "#$i: $(expected[i][1])-$(expected[i][2]) != $(actual[i][1])-$(actual[i][2])\n"
        end
    end

    return result
end

"""
    pwc_ooo_astr(expected, actual)

Piece-wise comparison of ooo noting differences in the return string.
This function takes vecors of strings.
"""
function pwc_ooo_astr(expected::Vector, actual::Vector)
    max_length = max(length(expected), length(actual))
    result = ""

    for i in 1:max_length
        if i > length(expected)
            result = result * "#$i: nothing !=  $(actual[i])\n"
        elseif i > length(actual)
            result = result * "#$i: $(expected[i]) != nothing\n"
        elseif expected[i] != actual[i]
            result = result * "#$i: $(expected[i]) != $(actual[i])\n"
        end
    end

    return result
end

"""
    get_default_sim_params()

Returns default simulation parameters useful for tests.

# Returns
- `Dict{String,Any}`: Default simulation parameters
"""
function get_default_sim_params()::Dict{String,Any}
    return Dict{String,Any}(
        "time_step_seconds" => 900,
        "time" => 0,
        "epsilon" => 1e-9,
        "is_first_timestep" => true,
    )
end
