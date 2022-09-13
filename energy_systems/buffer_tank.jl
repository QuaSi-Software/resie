Base.@kwdef mutable struct BufferTank <: ControlledSystem
    controller :: StateMachine = StateMachine()

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