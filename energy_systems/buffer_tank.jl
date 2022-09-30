Base.@kwdef mutable struct BufferTank <: ControlledSystem
    controller :: StateMachine
    inputs :: Dict{MediumCategory, ControlledSystem}
    outputs :: Dict{MediumCategory, ControlledSystem}
    accepted_inputs :: Vector{MediumCategory}
    accepted_outputs :: Vector{MediumCategory}

    input_interfaces :: Dict{MediumCategory, SystemInterface}
    output_interfaces :: Dict{MediumCategory, SystemInterface}

    capacity :: Float64
    load :: Float64
end

function make_BufferTank(capacity :: Float64, load :: Float64) :: BufferTank
    return BufferTank(
        StateMachine(), # controller
        Dict{MediumCategory, ControlledSystem}(), # inputs
        Dict{MediumCategory, ControlledSystem}(), # outputs
        [m_h_w_60c], # accepted_inputs
        [m_h_w_60c], # accepted_outputs
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