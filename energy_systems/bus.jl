Base.@kwdef mutable struct Bus <: ControlledSystem
    controller :: StateMachine
    is_storage :: Bool
    medium :: MediumCategory

    input_interfaces :: Vector{SystemInterface}
    output_interfaces :: Vector{SystemInterface}

    balance :: Float64
end

function make_Bus(medium :: MediumCategory) :: Bus
    return Bus(
        StateMachine(), # controller
        false, # is_storage
        medium, # medium
        [], # input_interfaces
        [], # output_interfaces
        0.0 # balance
    )
end

function gather_from_all!(interface :: SystemInterface, unit :: Bus)
    balance = 0.0

    for inface in unit.input_interfaces
        balance += inface.balance
        inface.balance = 0.0
    end

    for outface in unit.output_interfaces
        balance += outface.balance
        outface.balance = 0.0
    end

    interface.balance = balance
end

function reset(unit :: Bus)
    for inface in unit.input_interfaces
        inface.balance = 0.0
    end
    for outface in unit.output_interfaces
        outface.balance = 0.0
    end
end

function produce(unit :: Bus, parameters :: Dict{String, Any}, watt_to_wh :: Function)
    # nothing to do
end

function specific_values(unit :: Bus, time :: Int) :: Vector{Tuple}
    return [("Balance", "$(unit.balance)")]
end

export Bus, specific_values, make_Bus, gather_from_all!, reset