Base.@kwdef mutable struct GridConnection <: ControlledSystem
    controller :: StateMachine = StateMachine()
    medium :: MediumCategory
    inputs :: Dict{MediumCategory, ControlledSystem} = Dict{MediumCategory, ControlledSystem}()
    outputs :: Dict{MediumCategory, ControlledSystem} = Dict{MediumCategory, ControlledSystem}()
    accepted_inputs :: Vector{MediumCategory}
    accepted_outputs :: Vector{MediumCategory}

    draw_sum :: Float64 = 0.0
    load_sum :: Float64 = 0.0
end

function specific_values(unit :: GridConnection, time :: Int) :: Vector{Tuple}
    return [
        ("Draw sum", "$(unit.draw_sum)"),
        ("Load sum", "$(unit.load_sum)")
    ]
end

export GridConnection, specific_values