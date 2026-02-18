using Dates

"""
Control module for setting limits to the PLR of a component according to a profile.
"""
mutable struct CM_ChargeControl <: ControlModule
    name::String
    parameters::Dict{String,Any}
    profile::Profile

    function CM_ChargeControl(parameters::Dict{String,Any},
                              components::Grouping,
                              sim_params::Dict{String,Any},
                              unit_uac::String)
        default_parameters = Dict{String,Any}(
            "name" => "charge_control",
            "charge_is_allowed" => true,
            "discharge_is_allowed" => true
        )
        params = Base.merge(default_parameters, parameters)

        return new(params["name"], params)
    end
end

# method for control module name on type-level
control_module_name(x::Type{CM_ChargeControl})::String = "charge_control"

function has_method_for(mod::CM_ChargeControl, func::ControlModuleFunction)::Bool
    return func == cmf_charge_is_allowed || func == cmf_discharge_is_allowed
end

function update(mod::CM_ChargeControl)
    # nothing to do
end

function charge_is_allowed(mod::CM_ChargeControl, sim_params::Dict{String,Any})::Bool
    return mod.parameters["charge_is_allowed"]
end

function discharge_is_allowed(mod::CM_ChargeControl, sim_params::Dict{String,Any})::Bool
    return mod.parameters["discharge_is_allowed"]
end
