"""
Implementation of a buffer tank holding hot water for heating or DHW purposes.

For the moment this remains a very simple implementation that only has a load (how much
energy is stored) and a capacity, with no temperatures being considered. Given how the
simulation engine works, there will likely always be the need to deal with energy being
transfered with water temperature being secondary input variables.
"""
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

function output_values(unit :: BufferTank) :: Vector{String}
    return ["IN", "OUT", "Load", "Capacity"]
end

function output_value(unit :: BufferTank, key :: OutputKey) :: Float64
    if key.key_value == "IN"
        return unit.input_interfaces[key.medium].sum_abs_change * 0.5
    elseif key.key_value == "OUT"
        return unit.output_interfaces[key.medium].sum_abs_change * 0.5
    elseif key.key_value == "Load"
        return unit.load
    elseif key.key_value == "Capacity"
        return unit.capacity
    end
    raise(KeyError(key.key_value))
end

export BufferTank, make_BufferTank, output_values, output_value