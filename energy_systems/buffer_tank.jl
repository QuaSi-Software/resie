Base.@kwdef mutable struct BufferTank <: ControlledSystem
    controller :: StateMachine
    is_storage :: Bool

    input_interfaces :: InterfaceMap
    output_interfaces :: InterfaceMap

    capacity :: Float64
    load :: Float64
end

function make_BufferTank(capacity :: Float64, load :: Float64) :: BufferTank
    return BufferTank(
        StateMachine(), # controller
        true, # is_storage
        InterfaceMap( # input_interfaces
            m_h_w_60c => nothing
        ),
        InterfaceMap( # output_interfaces
            m_h_w_60c => nothing
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