Base.@kwdef mutable struct Demand <: ControlledSystem
    controller :: StateMachine
    is_storage :: Bool
    medium :: MediumCategory

    input_interfaces :: InterfaceMap
    output_interfaces :: InterfaceMap

    load :: Float64
end

function make_Demand(medium :: MediumCategory, load :: Float64) :: Demand
    return Demand(
        StateMachine(), # controller
        false, # is_storage
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

function specific_values(unit :: Demand, time :: Int) :: Vector{Tuple}
    return [("Load", "$(Wh(load_at_time(unit, time)))")]
end

function produce(unit :: Demand, parameters :: Dict{String, Any}, watt_to_wh :: Function)
    inface = unit.input_interfaces[unit.medium]
    sub!(inface, watt_to_wh(load_at_time(unit, parameters["time"])))
end

function load_at_time(unit :: Demand, time :: Int)
    if unit.medium == m_e_ac_230v
        return unit.load
    end

    seconds_in_day = 60 * 60 * 24
    time_of_day = Float64(time % seconds_in_day) / seconds_in_day

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