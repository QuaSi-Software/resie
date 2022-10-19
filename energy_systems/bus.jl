Base.@kwdef mutable struct Bus <: ControlledSystem
    controller :: StateMachine
    sys_function :: SystemFunction
    medium :: MediumCategory

    input_interfaces :: Vector{SystemInterface}
    output_interfaces :: Vector{SystemInterface}

    storage_space :: Float64
end

function make_Bus(medium :: MediumCategory) :: Bus
    return Bus(
        StateMachine(), # controller
        bus, # sys_function
        medium, # medium
        [], # input_interfaces
        [], # output_interfaces,
        0.0 # storage_space
    )
end

function gather_from_all!(interface :: SystemInterface, unit :: Bus)
    balance = 0.0

    for inface in unit.input_interfaces
        balance += inface.balance
        set!(inface, 0.0)
    end

    for outface in unit.output_interfaces
        balance += outface.balance
        set!(outface, 0.0)
    end

    set!(interface, balance)
end

function reset(unit :: Bus)
    for inface in unit.input_interfaces
        reset!(inface)
    end
    for outface in unit.output_interfaces
        reset!(outface)
    end
    unit.storage_space = 0.0
end

function produce(unit :: Bus, parameters :: Dict{String, Any}, watt_to_wh :: Function)
    for inface in unit.input_interfaces
        if inface.source.sys_function === storage
            unit.storage_space += inface.source.capacity - inface.source.load
        end
    end
end

function balance_on(interface :: SystemInterface, unit :: Bus) :: Float64
    balance = 0.0

    for inface in unit.input_interfaces
        balance += inface.balance
    end

    for outface in unit.output_interfaces
        balance += outface.balance
    end

    return balance, unit.storage_space
end

function specific_values(unit :: Bus, time :: Int) :: Vector{Tuple}
    return [("Balance", "$(balance(unit))")]
end

export Bus, specific_values, make_Bus, gather_from_all!, reset, balance_on