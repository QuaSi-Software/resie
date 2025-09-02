using Dates

"""
Control module for setting limits to the PLR of a component according to a profile.
"""
mutable struct CM_ProfileLimited <: ControlModule
    name::String
    parameters::Dict{String,Any}
    profile::Profile

    function CM_ProfileLimited(parameters::Dict{String,Any},
                               components::Grouping,
                               sim_params::Dict{String,Any},
                               unit_uac::String)
        default_parameters = Dict{String,Any}(
            "name" => "profile_limited",
            "profile_path" => nothing,
        )
        params = Base.merge(default_parameters, parameters)

        if params["profile_path"] === nothing
            @error "Required profile path for control module profile_limited not given"
        end
        profile = Profile(params["profile_path"], sim_params)

        return new("profile_limited", params, profile)
    end

    function CM_ProfileLimited()
        dummy_dates = DateTime[]
        push!(dummy_dates, DateTime(2024, 1, 1, 0, 0, 0))
        dummy_values = Float64[]
        push!(dummy_values, 0.0)
        dummy_profile = Profile("dummy",
                                Dict{String,Any}(
                                    "time_step_seconds" => 1,
                                    "start_date" => DateTime(2024, 1, 1, 0, 0, 0),
                                    "end_date" => DateTime(2024, 1, 1, 0, 0, 0),
                                    "force_profiles_to_repeat" => false,
                                );
                                given_profile_values=dummy_values,
                                given_timestamps=dummy_dates,
                                given_time_step=Second(1),
                                given_data_type="intensive")
        return new("profile_limited", Dict{String,Any}(), dummy_profile)
    end
end

function has_method_for(mod::CM_ProfileLimited, func::ControlModuleFunction)::Bool
    return func == cmf_upper_plr_limit
end

function update(mod::CM_ProfileLimited)
    # nothing to do
end

function upper_plr_limit(mod::CM_ProfileLimited, sim_params::Dict{String,Any})::Float64
    return value_at_time(mod.profile, sim_params)
end
