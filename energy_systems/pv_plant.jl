Base.@kwdef mutable struct PVPlant <: ControlledSystem
    controller :: StateMachine = StateMachine()
    inputs :: Dict{MediumCategory, ControlledSystem} = Dict{MediumCategory, ControlledSystem}()
    outputs :: Dict{MediumCategory, ControlledSystem} = Dict{MediumCategory, ControlledSystem}()
    accepted_inputs :: Vector{MediumCategory} = []
    accepted_outputs :: Vector{MediumCategory} = [m_e_ac_230v]

    last_produced_e :: Float64 = 0.0

    amplitude :: Float64
end

function specific_values(unit :: PVPlant, time :: Int) :: Vector{Tuple}
    return [("Production E", "$(unit.last_produced_e)")]
end

export PVPlant, specific_values