Base.@kwdef mutable struct PVPlant <: ControlledSystem
    controller :: StateMachine
    sys_function :: SystemFunction

    input_interfaces :: InterfaceMap
    output_interfaces :: InterfaceMap

    amplitude :: Float64
end

function make_PVPlant(amplitude :: Float64) :: PVPlant
    return PVPlant(
        StateMachine(), # controller
        limited_source, # sys_function
        InterfaceMap(), # input_interfaces
        InterfaceMap( # output_interfaces
            m_e_ac_230v => nothing
        ),
        amplitude, # amplitude
    )
end

function specific_values(unit :: PVPlant, time :: Int) :: Vector{Tuple}
    return []
end

function produce(unit :: PVPlant, parameters :: Dict{String, Any}, watt_to_wh :: Function)
    outface = unit.output_interfaces[m_e_ac_230v]
    add!(outface, watt_to_wh(power_at_time(unit, parameters["time"])))
end

function power_at_time(plant :: PVPlant, time :: Int) :: Float64
    seconds_in_day = 60 * 60 * 24
    base_sine = Base.Math.sin(Base.MathConstants.pi * (time % seconds_in_day) / seconds_in_day)
    return Base.Math.max(0.0, Base.Math.min(
        plant.amplitude,
        1.4 * plant.amplitude * base_sine * base_sine * base_sine - 0.2 * plant.amplitude
    ))
end

export PVPlant, specific_values, make_PVPlant