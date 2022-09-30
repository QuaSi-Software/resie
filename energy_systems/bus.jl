Base.@kwdef mutable struct Bus <: ControlledSystem
    controller :: StateMachine
    medium :: MediumCategory
    inputs :: Vector{ControlledSystem}
    outputs :: Vector{ControlledSystem}
    accepted_inputs :: Vector{MediumCategory}
    accepted_outputs :: Vector{MediumCategory}

    input_interfaces :: Dict{MediumCategory, SystemInterface}
    output_interfaces :: Dict{MediumCategory, SystemInterface}

    balance :: Float64
end

function make_Bus(medium :: MediumCategory) :: Bus
    return Bus(
        StateMachine(), # controller
        medium, # medium
        [], # inputs
        [], # outputs
        [medium], # accepted_inputs
        [medium], # accepted_outputs
        Dict{MediumCategory, SystemInterface}( # input_interfaces
            medium => SystemInterface()
        ),
        Dict{MediumCategory, SystemInterface}( # output_interfaces
            medium => SystemInterface()
        ),
        0.0, # balance
    )
end

function specific_values(unit :: Bus, time :: Int) :: Vector{Tuple}
    return [("Balance", "$(unit.balance)")]
end

export Bus, specific_values, make_Bus