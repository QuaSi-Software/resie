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
