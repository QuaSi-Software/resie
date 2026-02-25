"""
Control module Limit Return Temperature

This control module determines if the return flow temperature of the cooling 
(return flow from the target to the source into the energy output of the source) is low 
enough to satisfy the cooling demand of the source.
This is particularly relevant for the connection between electrolyser and seasonal thermal
storage, as the electrolyser may require a relatively large temperature spread.

Currently implemented for
    Sources: Electrolyser
    Targets: SeasonalThermalStorage

Source components need to have the following functions:
    - the need of a maximum input temperature for energy output

Target components need to have the following functions:
    - a given variable "current_energy_input_return_temperature" in the target unit that specifies 
      the temperature that is returned to the source during the energy transfer from source to target

"""

mutable struct CM_LimitCoolingInputTemperature <: ControlModule
    name::String
    parameters::Dict{String,Any}

    function CM_LimitCoolingInputTemperature(parameters::Dict{String,Any},
                                             components::Grouping,
                                             sim_params::Dict{String,Any},
                                             source_uac::String)
        default_parameters = Dict{String,Any}(
            "name" => "limit_cooling_input_temperature",
            "source_uac" => source_uac,
            "target_uac" => nothing,
            "temperature_limit" => nothing,
        )
        params = Base.merge(default_parameters, parameters)

        # check if "temperature_limit" is given
        if params["temperature_limit"] === nothing
            @error "The control module limit_cooling_input_temperature of component $(params["source_uac"]) " *
                   "has no temperature_limit given!"
            throw(InputError())
        end

        # check if "target_uac" is given
        if params["target_uac"] === nothing
            @error "The control module limit_cooling_input_temperature of component $(params["source_uac"]) " *
                   "has no target_uac given!"
            throw(InputError())
        end

        # check if a valid target uac is specified
        if !(params["target_uac"] !== nothing
             && params["target_uac"] in keys(components)
             && components[params["target_uac"]] isa LimitCoolingInputTemperatureTarget)
            @error "The target of the control module limit_cooling_input_temperature of component $(params["source_uac"]) is " *
                   "not a valid component for this control module!"
            throw(InputError())
        end

        # check if a valid source uac is specified
        if !(params["source_uac"] !== nothing
             && params["source_uac"] in keys(components)
             && components[params["source_uac"]] isa LimitCoolingInputTemperatureSource)
            @error "The source of the control module limit_cooling_input_temperature of component $(params["source_uac"]) " *
                   "is not a valid source for this control module"
            throw(InputError())
        end

        return new("limit_cooling_input_temperature", params)
    end
end

# method for control module name on type-level
control_module_name(x::Type{CM_LimitCoolingInputTemperature})::String = "limit_cooling_input_temperature"

function has_method_for(mod::CM_LimitCoolingInputTemperature, func::ControlModuleFunction)::Bool
    return func == cmf_limit_cooling_input_temperature
end

function update(mod::CM_LimitCoolingInputTemperature)
    # nothing to do
end

"""
    cooling_input_temperature_exceeded(mod::CM_LimitCoolingInputTemperature,
                                     components::Grouping,
                                     target_uac::String)

This function checks if the "temperature_limit" specified in the control module for the 
current input temperature of the cooling energy output of the source is exceeded by the
actual return temperature provided by the input interface of the target.

"""
function cooling_input_temperature_exceeded(mod::CM_LimitCoolingInputTemperature,
                                            components::Grouping,
                                            target_uac::Stringing)::Bool
    if mod.parameters["temperature_limit"] < components[target_uac].current_energy_input_return_temperature
        # temperature exceeded
        return true
    else
        return false
    end
end
