"""
Implementation of an energy system that models the demand consumers in a building require.

As the simulation does not encompass demand calculations, this is usually taken from other
tools that calculate the demand before a simulation of the energy systems is done. These
profiles usually are normalized to some degree, therefore Demand instances require a scaling
factor to turn the relative values to absolute values of required energy.
"""
mutable struct Demand <: ControlledSystem
    uac :: String
    controller :: Controller
    sys_function :: SystemFunction
    medium :: MediumCategory

    input_interfaces :: InterfaceMap
    output_interfaces :: InterfaceMap

    profile :: Profile
    scaling_factor :: Float64

    function Demand(uac :: String, config :: Dict{String, Any})
        profile = Profile(config["profile_file_path"])
        medium = getproperty(EnergySystems, Symbol(config["medium"]))

        return new(
            uac, # uac
            controller_for_strategy( # controller
                config["strategy"]["name"], config["strategy"]
            ),
            sf_fixed_sink, # sys_function
            medium, # medium
            InterfaceMap( # input_interfaces
                medium => nothing
            ),
            InterfaceMap( # output_interfaces
                medium => nothing
            ),
            profile, # profile
            config["scale"], # scaling_factor
        )
    end
end

function output_values(unit :: Demand) :: Vector{String}
    return ["IN", "Load"]
end

function output_value(unit :: Demand, key :: OutputKey) :: Float64
    if key.value_key == "IN"
        return unit.input_interfaces[key.medium].sum_abs_change * 0.5
    elseif key.value_key == "Load"
        return unit.scaling_factor # @TODO: Save the last calculated values and return it
    end
    raise(KeyError(key.value_key))
end

function produce(unit :: Demand, parameters :: Dict{String, Any}, watt_to_wh :: Function)
    inface = unit.input_interfaces[unit.medium]
    sub!(inface, unit.scaling_factor * work_at_time(unit.profile, parameters["time"]))
end

export Demand