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

function produce(unit :: BufferTank, parameters :: Dict{String, Any}, watt_to_wh :: Function)
    outface = unit.output_interfaces[m_h_w_60c]
    balance = gather_from_all(outface, outface.target)

    if balance >= 0.0
        return # produce is only concerned with moving energy to the target
    end

    # first we only check if there is a balance that produce needs to handle
    # without already gathering energy in the output interface. otherwise,
    # when the balance is positive, no energy is moved but the act of calling
    # gather_from_all! has incorrectly recorded a move of energy
    gather_from_all!(outface, outface.target)

    if unit.load > outface.balance
        unit.load += outface.balance
        set!(outface, 0.0)
    else
        add!(outface, unit.load)
        unit.load = 0.0
    end
end

function load(unit :: BufferTank, parameters :: Dict{String, Any}, watt_to_wh :: Function)
    inface = unit.input_interfaces[m_h_w_60c]
    gather_from_all!(inface, inface.source)

    if inface.balance <= 0.0
        return # load is only concerned with receiving energy from the target
    end

    unit.load += inface.balance # @TODO: check if loading exceeds capacity
    set!(inface, 0.0)
end

function specific_values(unit :: BufferTank, time :: Int) :: Vector{Tuple}
    return [
        ("Load", "$(unit.load)"),
        ("Capacity", "$(unit.capacity)")
    ]
end

export BufferTank, specific_values, make_BufferTank