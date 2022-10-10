Base.@kwdef mutable struct Bus <: ControlledSystem
    controller :: StateMachine
    is_storage :: Bool
    medium :: MediumCategory

    input_interfaces :: Vector{SystemInterface}
    output_interfaces :: Vector{SystemInterface}
end

function make_Bus(medium :: MediumCategory) :: Bus
    return Bus(
        StateMachine(), # controller
        false, # is_storage
        medium, # medium
        [], # input_interfaces
        [], # output_interfaces
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
end

function produce(unit :: Bus, parameters :: Dict{String, Any}, watt_to_wh :: Function)
    # nothing to do
end

function check_balance(unit :: Bus) :: Float64
    balance = 0.0

    for inface in unit.input_interfaces
        balance += inface.balance
    end

    for outface in unit.output_interfaces
        balance += outface.balance
    end

    return balance
end

function specific_values(unit :: Bus, time :: Int) :: Vector{Tuple}
    return [("Balance", "$(check_balance(unit))")]
end

export Bus, specific_values, make_Bus, gather_from_all!, reset