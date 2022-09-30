Base.@kwdef mutable struct GridConnection <: ControlledSystem
    controller :: StateMachine
    medium :: MediumCategory
    inputs :: Dict{MediumCategory, ControlledSystem}
    outputs :: Dict{MediumCategory, ControlledSystem}
    accepted_inputs :: Vector{MediumCategory}
    accepted_outputs :: Vector{MediumCategory}

    input_interfaces :: Dict{MediumCategory, SystemInterface}
    output_interfaces :: Dict{MediumCategory, SystemInterface}

    draw_sum :: Float64
    load_sum :: Float64
end

function make_GridConnection(medium :: MediumCategory) :: GridConnection
    return GridConnection(
        StateMachine(), # controller
        medium, # medium
        Dict{MediumCategory, ControlledSystem}(), # inputs
        Dict{MediumCategory, ControlledSystem}(), # outputs
        [medium], # accepted_inputs
        [medium], # accepted_outputs
        Dict{MediumCategory, SystemInterface}( # input_interfaces
            medium => SystemInterface()
        ),
        Dict{MediumCategory, SystemInterface}( # output_interfaces
            medium => SystemInterface()
        ),
        0.0, # draw_sum,
        0.0, # load_sum
    )
end

function specific_values(unit :: GridConnection, time :: Int) :: Vector{Tuple}
    return [
        ("Draw sum", "$(unit.draw_sum)"),
        ("Load sum", "$(unit.load_sum)")
    ]
end

export GridConnection, specific_values, make_GridConnection