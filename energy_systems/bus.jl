Base.@kwdef mutable struct Bus <: ControlledSystem
    controller :: StateMachine = StateMachine()
    medium :: MediumCategory
    inputs :: Vector{ControlledSystem} = Vector{ControlledSystem}()
    outputs :: Vector{ControlledSystem} = Vector{ControlledSystem}()
    accepted_inputs :: Vector{MediumCategory}
    accepted_outputs :: Vector{MediumCategory}

    balance :: Float64 = 0.0
end

function specific_values(unit :: Bus, time :: Int) :: Vector{Tuple}
    return [("Balance", "$(unit.balance)")]
end

export Bus, specific_values