Base.@kwdef mutable struct PVPlant <: ControlledSystem
    controller :: StateMachine
    is_storage :: Bool

    input_interfaces :: InterfaceMap
    output_interfaces :: InterfaceMap

    last_produced_e :: Float64

    amplitude :: Float64
end

function make_PVPlant(amplitude :: Float64) :: PVPlant
    return PVPlant(
        StateMachine(), # controller
        false, # is_storage
        InterfaceMap(), # input_interfaces
        InterfaceMap( # output_interfaces
            m_e_ac_230v => nothing
        ),
        0.0, # last_produced_e
        amplitude, # amplitude
    )
end

function specific_values(unit :: PVPlant, time :: Int) :: Vector{Tuple}
    return [("Production E", "$(unit.last_produced_e)")]
end

function produce(unit :: PVPlant, parameters :: Dict{String, Any}, watt_to_wh :: Function)
    outface = unit.output_interfaces[m_e_ac_230v]
    outface.balance += watt_to_wh(power_at_time(unit, parameters["time"]))
    unit.last_produced_e = watt_to_wh(power_at_time(unit, parameters["time"]))
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