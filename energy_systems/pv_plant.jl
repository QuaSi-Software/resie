Base.@kwdef mutable struct PVPlant <: ControlledSystem
    controller :: StateMachine
    inputs :: Dict{MediumCategory, ControlledSystem}
    outputs :: Dict{MediumCategory, ControlledSystem}
    accepted_inputs :: Vector{MediumCategory}
    accepted_outputs :: Vector{MediumCategory}

    input_interfaces :: Dict{MediumCategory, SystemInterface}
    output_interfaces :: Dict{MediumCategory, SystemInterface}

    last_produced_e :: Float64

    amplitude :: Float64
end

function make_PVPlant(amplitude :: Float64) :: PVPlant
    return PVPlant(
        StateMachine(), # controller
        Dict{MediumCategory, ControlledSystem}(), # inputs
        Dict{MediumCategory, ControlledSystem}(), # outputs
        [], # accepted_inputs
        [m_e_ac_230v], # accepted_outputs
        Dict{MediumCategory, SystemInterface}(), # input_interfaces
        Dict{MediumCategory, SystemInterface}( # output_interfaces
            m_h_w_60c => SystemInterface()
        ),
        0.0, # last_produced_e
        amplitude, # amplitude
    )
end

function specific_values(unit :: PVPlant, time :: Int) :: Vector{Tuple}
    return [("Production E", "$(unit.last_produced_e)")]
end

export PVPlant, specific_values, make_PVPlant