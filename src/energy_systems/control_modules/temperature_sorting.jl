"""
Control module for reordering the inputs and outputs of a component to be most beneficial
in terms of temperature levels. Usually this means that inputs are sorted by descending
maximum temperatures and outputs by ascending minimum temperatures. This is also the default
configuration of the module and used for heat pumps because a minimum temperature difference
is optimal.
"""
mutable struct CM_Temperature_Sorting <: ControlModule
    name::String
    parameters::Dict{String,Any}

    function CM_Temperature_Sorting(parameters::Dict{String,Any},
                                    components::Grouping,
                                    sim_params::Dict{String,Any})
        default_parameters = Dict{String,Any}(
            "name" => "temperature_sorting",
            "input_temps" => "max",
            "input_order" => "desc",
            "output_temps" => "min",
            "output_order" => "asc",
        )
        params = Base.merge(default_parameters, parameters)

        return new("temperature_sorting", params)
    end
end

function has_method_for(mod::CM_Temperature_Sorting, func::ControlModuleFunction)::Bool
    return func == cmf_reorder_inputs || func == cmf_reorder_outputs
end

function update(mod::CM_Temperature_Sorting)
    # nothing to do
end

"""
Returns if a is strictly greater than b and considers values of nothing.

# Arguments
- `a::Temperature`: The first temperature
- `b::Temperature`: The second temperature
# Returns
- `Bool`: True if a is not nothing and (a > b or b is nothing)
"""
function coalesce_greater(a::Temperature, b::Temperature)::Bool
    if a === nothing
        return false
    elseif b === nothing
        return true
    else
        return a > b
    end
end

function reorder_inputs(mod::CM_Temperature_Sorting,
                        temps_min::Vector{<:Temperature},
                        temps_max::Vector{<:Temperature})::Vector{Integer}
    temps = mod.parameters["input_temps"] == "max" ? temps_max : temps_min
    do_reverse = mod.parameters["input_order"] == "asc"
    return sortperm(temps; by=x -> x, lt=coalesce_greater, rev=do_reverse)
end

function reorder_outputs(mod::CM_Temperature_Sorting,
                         temps_min::Vector{<:Temperature},
                         temps_max::Vector{<:Temperature})::Vector{Integer}
    temps = mod.parameters["output_temps"] == "max" ? temps_max : temps_min
    do_reverse = mod.parameters["output_order"] == "asc"
    return sortperm(temps; by=x -> x, lt=coalesce_greater, rev=do_reverse)
end
