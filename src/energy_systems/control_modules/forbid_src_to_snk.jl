"""
Control module for disallowing a specific source to be used to supply a specific sink.
This can be used, for example, to avoid a heat pump raising a heat storage component's
temperature using electricity, while still allowing the heat pump to draw from the storage
and putting excess heat into the storage, as long as there other sinks and sources involved.
"""
mutable struct CM_Forbid_Src_To_Snk <: ControlModule
    name::String
    parameters::Dict{String,Any}

    function CM_Forbid_Src_To_Snk(parameters::Dict{String,Any},
                                  components::Grouping,
                                  sim_params::Dict{String,Any},
                                  unit_uac::String)
        default_parameters = Dict{String,Any}(
            "name" => "forbid_src_to_snk",
            "src_uac" => nothing,
            "snk_uac" => nothing,
        )
        params = Base.merge(default_parameters, parameters)

        return new("forbid_src_to_snk", params)
    end
end

function has_method_for(mod::CM_Forbid_Src_To_Snk, func::ControlModuleFunction)::Bool
    return func == cmf_check_src_to_snk
end

function update(mod::CM_Forbid_Src_To_Snk)
    # nothing to do
end

function check_src_to_snk(mod::CM_Forbid_Src_To_Snk,
                          in_uac::Stringing,
                          out_uac::Stringing)::Bool
    return in_uac != mod.parameters["src_uac"] || out_uac != mod.parameters["snk_uac"]
end
