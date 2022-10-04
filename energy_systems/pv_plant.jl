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
            m_h_w_60c => nothing
        ),
        0.0, # last_produced_e
        amplitude, # amplitude
    )
end

function specific_values(unit :: PVPlant, time :: Int) :: Vector{Tuple}
    return [("Production E", "$(unit.last_produced_e)")]
end

export PVPlant, specific_values, make_PVPlant