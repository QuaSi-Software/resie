"""

"""
mutable struct CM_ProfileLimited <: ControlModule
    name::String
    parameters::Dict{String,Any}
    profile::Profile

    function CM_ProfileLimited(
        parameters::Dict{String,Any},
        components::Grouping,
        sim_params::Dict{String,Any}
    )
        default_parameters=Dict{String,Any}(
            "name" => "default",
            "profile_path" => nothing
        )
        params = Base.merge(default_parameters, parameters)

        if params["profile_path"] === nothing
            @error "Required profile path for control module profile_limited not given"
        end
        profile = Profile(params["profile_path"], sim_params)

        return new("default", params, profile)
    end
end

function update(mod::CM_ProfileLimited)
    # nothing to do
end

function upper_plr_limit(mod::CM_ProfileLimited, sim_params::Dict{String,Any})
    return value_at_time(mod.profile, sim_params["time"])
end