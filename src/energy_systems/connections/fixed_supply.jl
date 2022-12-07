"""
Implementation of an energy system modeling an abstract fixed supply of some medium.

This is particularly useful for testing, but can also be used to model any energy system
or other equipment unit that produces energy in a medium, all of which has to be consumed
as the energy system cannot be dispatched like a grid connection can.
Note that "fixed" in this context means that the amount of energy the unit produces is
fixed within a timestep, but can vary over multiple timesteps. No calculation other than
scaling of profile values is performed in each timestep.
"""
mutable struct FixedSupply <: ControlledSystem
    uac :: String
    controller :: Controller
    sys_function :: SystemFunction
    medium :: MediumCategory

    input_interfaces :: InterfaceMap
    output_interfaces :: InterfaceMap

    profile :: Profile
    scaling_factor :: Float64

    function FixedSupply(uac :: String, config :: Dict{String, Any})
        profile = Profile(config["profile_file_path"])
        medium = getproperty(EnergySystems, Symbol(config["medium"]))

        return new(
            uac, # uac
            controller_for_strategy( # controller
                config["strategy"]["name"], config["strategy"]
            ),
            sf_fixed_source, # sys_function
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

function output_values(unit :: FixedSupply) :: Vector{String}
    return ["OUT"]
end

function output_value(unit :: FixedSupply, key :: OutputKey) :: Float64
    if key.value_key == "OUT"
        return unit.output_interfaces[key.medium].sum_abs_change * 0.5
    end
    raise(KeyError(key.value_key))
end

function produce(unit :: FixedSupply, parameters :: Dict{String, Any}, watt_to_wh :: Function)
    outface = unit.output_interfaces[unit.medium]
    add!(outface, unit.scaling_factor * work_at_time(unit.profile, parameters["time"]))
end

export FixedSupply