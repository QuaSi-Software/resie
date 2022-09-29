Base.@kwdef mutable struct BufferTank <: ControlledSystem
    controller :: StateMachine = StateMachine()
    medium = m_h_w_60c
    inputs :: Dict{MediumCategory, ControlledSystem} = Dict{MediumCategory, ControlledSystem}()
    outputs :: Dict{MediumCategory, ControlledSystem} = Dict{MediumCategory, ControlledSystem}()
    accepted_inputs :: Vector{MediumCategory} = [m_h_w_60c]
    accepted_outputs :: Vector{MediumCategory} = [m_h_w_60c]

    capacity :: Float64
    load :: Float64
end

function specific_values(unit :: BufferTank, time :: Int) :: Vector{Tuple}
    return [
        ("Load", "$(unit.load)"),
        ("Capacity", "$(unit.capacity)")
    ]
end

export BufferTank, specific_values