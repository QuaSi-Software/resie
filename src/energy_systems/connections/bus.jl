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
    uac::String
    controller::Controller
    sys_function::SystemFunction
    medium::Symbol

    input_interfaces::Vector{SystemInterface}
    output_interfaces::Vector{SystemInterface}

    input_priorities::Vector{String}
    output_priorities::Vector{String}

    remainder::Float64

    function Bus(uac::String, config::Dict{String,Any})
        medium = Symbol(config["medium"])
        return new(
            uac, # uac
            controller_for_strategy( # controller
                config["strategy"]["name"], config["strategy"]
            ),
            sf_bus, # sys_function
            medium, # medium
            [], # input_interfaces
            [], # output_interfaces,
            config["input_priorities"],
            config["production_refs"],
            0.0 # remainder
        )
    end
end

function reset(unit::Bus)
    for inface in unit.input_interfaces
        reset!(inface)
    end
    for outface in unit.output_interfaces
        reset!(outface)
    end
    unit.remainder = 0.0
end

"""
    balance_nr(unit, caller)

Variant of [`balance`](@ref) that includes other connected bus systems and their energy
balance, but does so in a non-recursive manner such that any bus in the chain of connected
bus systems is only considered once.
"""
function balance_nr(unit::Bus, caller::Bus)
    balance = 0.0

    for inface in unit.input_interfaces
        if inface.source == caller
            continue
        end

        if isa(inface.source, Bus)
            supply = balance_nr(inface.source, unit)
            if supply < 0.0
                continue
            end
        else
            supply, _, _ = balance_on(inface, inface.source)
        end
        balance += supply
    end

    for outface in unit.output_interfaces
        if outface.target == caller
            continue
        end

        if isa(outface.target, Bus)
            demand = balance_nr(outface.target, unit)
            if demand > 0.0
                continue
            end
        else
            demand, _, _ = balance_on(outface, outface.target)
        end
        balance += demand
    end

    return balance + unit.remainder
end

"""
    balance_direct(unit)

Energy balance on a bus system without considering any other connected bus systems.
"""
function balance_direct(unit::Bus)::Float64
    balance = 0.0

    for inface in unit.input_interfaces
        if isa(inface.source, Bus)
            continue
        else
            supply, _, _ = balance_on(inface, inface.source)
            balance += supply
        end
    end

    for outface in unit.output_interfaces
        if isa(outface.target, Bus)
            continue
        else
            demand, _, _ = balance_on(outface, outface.target)
            balance += demand
        end
    end

    return balance + unit.remainder
end

function balance(unit::Bus)::Float64
    # we can use the non-recursive version of the method as a bus will never
    # be connected to itself... right?
    return balance_nr(unit, unit)
end

function balance_on(
    interface::SystemInterface,
    unit::Bus
)::Tuple{Float64,Float64,Temperature}
    highest_demand_temp = -1e9
    storage_space = 0.0

    for outface in unit.output_interfaces
        if outface.target.sys_function === sf_bus
            balance, potential, temperature = balance_on(outface, outface.target)
        else
            balance = outface.balance
            temperature = outface.temperature
            if outface.target.sys_function === sf_storage
                _, potential, _ = balance_on(outface, outface.target)
            else
                potential = 0.0
            end
        end

        if temperature !== nothing && balance < 0
            highest_demand_temp = (
                temperature > highest_demand_temp ? temperature : highest_demand_temp
            )
        end

        storage_space += potential
    end

    return (
        balance(unit),
        storage_space,
        highest_demand_temp <= -1e9 ? nothing : highest_demand_temp
    )
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
function distribute!(unit::Bus)
    remainder = balance_direct(unit)

    # reset all non-bus input interfaces
    for inface in unit.input_interfaces
        if inface.source.sys_function !== sf_bus
            set!(inface, 0.0, inface.temperature)
        end
    end

    for output_uac in unit.output_priorities # for every output UAC (to enshure the correct order)...
        for outface in unit.output_interfaces  # ...seach corresponding output inferface by...
            if outface.target.sys_function === outface.source.sys_function  # ...making sure the output interface is from type bus...
                if outface.target.uac === output_uac # ...and is corresponding to the uac given in the output_priority
                    # calculate energy balance on current bus (outface.source) in order to write
                    # sum_abs_change into interface to following bus (outface.taget) 
                    sum_abs_change = 0
                    for input_interface in outface.source.input_interfaces
                        sum_abs_change += input_interface.sum_abs_change
                    end
                    for output_interface in outface.source.output_interfaces
                        sum_abs_change -= output_interface.sum_abs_change
                    end
                    # check if sum_abs_change is bigger than need in outface.target to handel multiple busses on unit
                    if (sum_abs_change / 2) > abs(outface.target.remainder)
                        # order of distribution into multiple busses is equal to order of simulations for now (last buss will get least energy)
                        sum_abs_change = 2 * abs(outface.target.remainder)
                    end

                    # write energy flow from source to target into sum_abs_change
                    outface.sum_abs_change = sum_abs_change

                    # adjust remainder on outface target bus
                    outface.target.remainder += outface.sum_abs_change / 2

                    # update remainder of current bus
                    remainder = balance_direct(unit)

                    # reset output interfaces
                    # set!(outface, 0.0, outface.temperature)
                end
            end
        end
    end

    # reset all non-bus output interfaces
    for outface in unit.output_interfaces
        if outface.target.sys_function !== sf_bus
            set!(outface, 0.0, outface.temperature)
        end
    end

    unit.remainder = remainder
end

function output_values(unit::Bus)::Vector{String}
    return ["Balance"]
end

function output_value(unit::Bus, key::OutputKey)::Float64
    if key.value_key == "Balance"
        return balance(unit)
    end
    throw(KeyError(key.value_key))
end

export Bus