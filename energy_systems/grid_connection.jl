Base.@kwdef mutable struct GridConnection <: ControlledSystem
    controller :: StateMachine
    medium :: MediumCategory

    input_interfaces :: InterfaceMap
    output_interfaces :: InterfaceMap

    draw_sum :: Float64
    load_sum :: Float64
end

function make_GridConnection(medium :: MediumCategory) :: GridConnection
    return GridConnection(
        StateMachine(), # controller
        medium, # medium
        InterfaceMap( # input_interfaces
            medium => nothing
        ),
        InterfaceMap( # output_interfaces
            medium => nothing
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