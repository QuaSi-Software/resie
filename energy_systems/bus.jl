Base.@kwdef mutable struct Bus <: ControlledSystem
    controller :: StateMachine
    medium :: MediumCategory

    input_interfaces :: Vector{SystemInterface}
    output_interfaces :: Vector{SystemInterface}

    balance :: Float64
end

function make_Bus(medium :: MediumCategory) :: Bus
    return Bus(
        StateMachine(), # controller
        medium, # medium
        [], # input_interfaces
        [], # output_interfaces
        0.0 # balance
    )
end

function specific_values(unit :: Bus, time :: Int) :: Vector{Tuple}
    return [("Balance", "$(unit.balance)")]
end

export Bus, specific_values, make_Bus