Base.@kwdef mutable struct PVPlant <: ControlledSystem
    controller :: StateMachine = StateMachine()

    last_produced_e :: Float64 = 0.0

    amplitude :: Float64
end

function specific_values(unit :: PVPlant, time :: Int) :: Vector{Tuple}
    return [("Production E", "$(unit.last_produced_e)")]
end

export PVPlant, specific_values