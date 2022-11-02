"""
Implementation of an energy system that models the demand consumers in a building require.

As the simulation does not encompass demand calculations, this is usually taken from other
tools that calculate the demand before a simulation of the energy systems is done. For the
moment this remains a simple implementation that does not consider any external profile
and instead approximates a demand profile with peaks in the morning and evening. The load
parameter is a scaling factor, but does not correspond to an average demand load.
"""
Base.@kwdef mutable struct Demand <: ControlledSystem
    controller :: StateMachine
    sys_function :: SystemFunction
    medium :: MediumCategory

    input_interfaces :: InterfaceMap
    output_interfaces :: InterfaceMap

    load :: Float64
end

function make_Demand(medium :: MediumCategory, load :: Float64) :: Demand
    return Demand(
        StateMachine(), # controller
        limited_sink, # sys_function
        medium, # medium
        InterfaceMap( # input_interfaces
            medium => nothing
        ),
        InterfaceMap( # output_interfaces
            medium => nothing
        ),
        load, # load
    )
end

function output_values(unit :: Demand) :: Vector{String}
    return ["IN", "Load"]
end

function output_value(unit :: Demand, key :: OutputKey) :: Float64
    if key.key_value == "IN"
        return unit.input_interfaces[key.medium].sum_abs_change * 0.5
    elseif key.key_value == "Load"
        return unit.load # @TODO: Save the last calculated values and return it
    end
    raise(KeyError(key.key_value))
end

function specific_values(unit :: Demand, time :: Int) :: Vector{Tuple}
    return []
end

function produce(unit :: Demand, parameters :: Dict{String, Any}, watt_to_wh :: Function)
    inface = unit.input_interfaces[unit.medium]
    sub!(inface, watt_to_wh(load_at_time(unit, parameters["time"])))
end

function load_at_time(unit :: Demand, time :: Int)
    seconds_in_day = 60 * 60 * 24
    time_of_day = Float64(time % seconds_in_day) / seconds_in_day

    if unit.medium == m_e_ac_230v
        if time_of_day < 0.25 || time_of_day >= 0.9
            return unit.load * 0.2
        elseif time_of_day >= 0.25 && time_of_day < 0.292
            return unit.load * 1.2
        elseif time_of_day >= 0.75 && time_of_day < 0.9
            return unit.load * 1.5
        else
            return unit.load * 0.5
        end
    end

    if time_of_day < 0.25 || time_of_day >= 0.8333
        return unit.load * 0.3
    elseif time_of_day >= 0.25 && time_of_day < 0.292
        return unit.load * 3.0
    elseif time_of_day >= 0.75 && time_of_day < 0.8333
        return unit.load * 2.0
    else
        return unit.load * 0.6
    end
end

export Demand, specific_values, load_at_time, make_Demand