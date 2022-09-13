Base.@kwdef mutable struct Bus <: ControlledSystem
    controller :: StateMachine = StateMachine()
    medium :: MediumCategory

    balance :: Float64 = 0.0
end

function specific_values(unit :: Bus, time :: Int) :: Vector{Tuple}
    return [("Balance", "$(unit.balance)")]
end

export Bus, specific_values