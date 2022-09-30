Base.@kwdef mutable struct BufferTank <: ControlledSystem
    controller :: StateMachine

    input_interfaces :: Dict{MediumCategory, SystemInterface}
    output_interfaces :: Dict{MediumCategory, SystemInterface}

    capacity :: Float64
    load :: Float64
end

function make_BufferTank(capacity :: Float64, load :: Float64) :: BufferTank
    return BufferTank(
        StateMachine(), # controller
        Dict{MediumCategory, SystemInterface}( # input_interfaces
            m_h_w_60c => SystemInterface()
        ),
        Dict{MediumCategory, SystemInterface}( # output_interfaces
            m_h_w_60c => SystemInterface()
        ),
        capacity, # capacity
        load # load
    )
end

function specific_values(unit :: BufferTank, time :: Int) :: Vector{Tuple}
    return [
        ("Load", "$(unit.load)"),
        ("Capacity", "$(unit.capacity)")
    ]
end

export BufferTank, specific_values, make_BufferTank