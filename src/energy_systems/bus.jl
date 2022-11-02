"""
Imnplementation of a bus energy system for balancing multiple inputs and outputs.

This energy system is both a possible real system (mostly for electricity) as well as a
necessary abstraction of the model. The basic idea is that one or more energy systems feed
energy of the same medium into a bus and one or more energy systems draw that energy from
the bus. A bus with only one input and only one output can be replaced with a direct
connection between both systems.

The function and purpose is described in more detail in the accompanying documentation.
"""
Base.@kwdef mutable struct Bus <: ControlledSystem
    controller :: StateMachine
    sys_function :: SystemFunction
    medium :: MediumCategory

    input_interfaces :: Vector{SystemInterface}
    output_interfaces :: Vector{SystemInterface}

    storage_space :: Float64
    remainder :: Float64
end

function make_Bus(medium :: MediumCategory) :: Bus
    return Bus(
        StateMachine(), # controller
        bus, # sys_function
        medium, # medium
        [], # input_interfaces
        [], # output_interfaces,
        0.0, # storage_space
        0.0 # remainder
    )
end

function reset(unit :: Bus)
    for inface in unit.input_interfaces
        reset!(inface)
    end
    for outface in unit.output_interfaces
        reset!(outface)
    end
    unit.storage_space = 0.0
    unit.remainder = 0.0
end

"""
    produce(unit, parameters, watt_to_wh)

Bus-specific implementation of produce.

The production of a bus does not generate or consume energy, but calculates the potential
of storage systems connected to the bus and saves that value, which is later required to
distinguish between the actual demand of consuming energy systems and the potential to load
storage with excess energy.

See also [`produce`](@ref)
"""
function produce(unit :: Bus, parameters :: Dict{String, Any}, watt_to_wh :: Function)
    for inface in unit.input_interfaces
        if inface.source.sys_function === storage
            unit.storage_space += inface.source.capacity - inface.source.load
        end
    end
end

function balance(unit :: Bus) :: Float64
    balance = 0.0

    for inface in unit.input_interfaces
        balance += inface.balance
    end

    for outface in unit.output_interfaces
        balance += outface.balance
    end

    return balance + unit.remainder
end

function balance_on(
    interface :: SystemInterface,
    unit :: Bus
) :: Tuple{Float64, Float64}
    return balance(unit), -unit.storage_space #  negative is demand
end

"""
    distribute!(unit)

Bus-specific implementation of distribute!.

This implementation checks the balance of the bus (which includes the so-called remainder),
then sets all interfaces back to zero and saves the balance as the current remainder. This
has the intention that calling this method will distribute the energy without changing the
energy balance of the bus. That way the method can be called at any point after the balance
might have changed and should always consider the correct energy balance regardless of when
or how often it was called.
"""
function distribute!(unit :: Bus)
    remainder = balance(unit)

    for inface in unit.input_interfaces
        set!(inface, 0.0)
    end

    for outface in unit.output_interfaces
        set!(outface, 0.0)
    end

    unit.remainder = remainder
end

function output_values(unit :: Bus) :: Vector{String}
    return ["Balance"]
end

function output_value(unit :: Bus, key :: OutputKey) :: Float64
    if key.key_value == "Balance"
        return balance(unit)
    end
    raise(KeyError(key.key_value))
end

function specific_values(unit :: Bus, time :: Int) :: Vector{Tuple}
    return [("Balance", "$(balance(unit))")]
end

export Bus, specific_values, make_Bus, reset, balance, balance_on, distribute!, produce