Base.@kwdef mutable struct BufferTank <: ControlledSystem
    controller :: StateMachine
    sys_function :: SystemFunction

    input_interfaces :: InterfaceMap
    output_interfaces :: InterfaceMap

    capacity :: Float64
    load :: Float64
end

function make_BufferTank(capacity :: Float64, load :: Float64) :: BufferTank
    return BufferTank(
        StateMachine(), # controller
        storage, # sys_function
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

function produce(unit :: BufferTank, parameters :: Dict{String, Any}, watt_to_wh :: Function)
    outface = unit.output_interfaces[m_h_w_60c]
    balance, _ = balance_on(outface, outface.target)

    if balance >= 0.0
        return # produce is only concerned with moving energy to the target
    end

    if unit.load > abs(balance)
        unit.load += balance
        add!(outface, abs(balance))
    else
        add!(outface, unit.load)
        unit.load = 0.0
    end
end

function load(unit :: BufferTank, parameters :: Dict{String, Any}, watt_to_wh :: Function)
    inface = unit.input_interfaces[m_h_w_60c]
    balance, _ = balance_on(inface, inface.source)

    if balance <= 0.0
        return # load is only concerned with receiving energy from the target
    end

    diff = unit.capacity - unit.load
    if diff > balance
        unit.load += balance
        sub!(inface, balance)
    else
        unit.load = unit.capacity
        sub!(inface, diff)
    end
end

function specific_values(unit :: BufferTank, time :: Int) :: Vector{Tuple}
    return [
        ("Load", "$(unit.load)"),
        ("Capacity", "$(unit.capacity)")
    ]
end

export BufferTank, specific_values, make_BufferTank