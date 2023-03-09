using ResumableFunctions

"""
Utility struct to contain the connections, input/output priorities and other related data
for bus systems.
"""
Base.@kwdef mutable struct ConnectionMatrix
    input_order::Union{Nothing,Vector{String}}
    output_order::Vector{String}
    storage_loading::Union{Nothing,Vector{Vector{Bool}}}

    function ConnectionMatrix(config::Dict{String,Any})
        input_order = nothing
        output_order = [String(u) for u in config["production_refs"]]
        storage_loading = nothing

        if "connection_matrix" in keys(config)
            if "input_order" in keys(config["connection_matrix"])
                input_order = [String(u) for u in config["connection_matrix"]["input_order"]]
            end

            if "output_order" in keys(config["connection_matrix"])
                output_order = [String(u) for u in config["connection_matrix"]["output_order"]]
            end

            if "storage_loading" in keys(config["connection_matrix"])
                storage_loading = []
                for row in config["connection_matrix"]["storage_loading"]
                    vec = [Bool(v) for v in row]
                    push!(storage_loading, vec)
                end
            end
        end

        return new(
            input_order,
            output_order,
            storage_loading,
        )
    end
end

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
    connectivity::ConnectionMatrix

    remainder::Float64

    function Bus(uac::String, config::Dict{String,Any})
        medium = Symbol(config["medium"])
        register_media([medium])

        return new(
            uac, # uac
            controller_for_strategy( # controller
                config["strategy"]["name"], config["strategy"]
            ),
            sf_bus, # sys_function
            medium, # medium
            [], # input_interfaces
            [], # output_interfaces,
            ConnectionMatrix(config),
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
            supply = max(balance_nr(inface.source, unit), inface.balance)
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
            demand = min(balance_nr(outface.target, unit), outface.balance)
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

# """
#     for x in bus_infaces(bus)

# Iterator over the input interfaces that connect the given bus to other busses.
# """
@resumable function bus_infaces(unit::Bus)
    # for every input UAC (to ensure the correct order)...
    for input_uac in unit.connectivity.input_order
        # ...seach corresponding input inferface by...
        for inface in unit.input_interfaces
            # ...making sure the input interface is of type bus...
            if inface.source.sys_function === sf_bus
                # ...and the source's UAC matches the one in the input_priority.
                if inface.source.uac === input_uac
                    @yield inface
                    break # we found the match, so we can break out of the inner loop.
                end
            end
        end
    end
end

# """
#     for x in bus_outfaces(bus)

# Iterator over the output interfaces that connect the given bus to other busses.
# """
@resumable function bus_outfaces(unit::Bus)
    # for every output UAC (to ensure the correct order)...
    for output_uac in unit.connectivity.output_order
        # ...seach corresponding output inferface by...
        for outface in unit.output_interfaces
            # ...making sure the output interface is of type bus...
            if outface.target.sys_function === sf_bus
                # ...and the target's UAC matches the one in the output_priority.
                if outface.target.uac === output_uac
                    @yield outface
                    break # we found the match, so we can break out of the inner loop.
                end
            end
        end
    end
end

"""
    distribute!(unit)

Bus-specific implementation of distribute!.

This moves the energy from connected energy system from supply to demand systems both
on the bus directly as well as taking other bus systems into account. This allows busses
to be connected in chains (but not loops) and "communicate" the energy across. The method
implicitly requires that each bus on the chain is called with distribute!() in a specific
order, which is explained in more detail in the documentation. Essentially it starts from
the leaves of the chain and progresses to the roots.
"""
function distribute!(unit::Bus)
    balance = balance_direct(unit)

    # reset all non-bus input interfaces
    for inface in unit.input_interfaces
        if inface.source.sys_function !== sf_bus
            set!(inface, 0.0, inface.temperature)
        end
    end

    # reset all non-bus output interfaces
    for outface in unit.output_interfaces
        if outface.target.sys_function !== sf_bus
            set!(outface, 0.0, outface.temperature)
        end
    end

    # distribute to outgoing busses according to output priority
    if balance > 0.0
        for outface in bus_outfaces(unit)
            if balance > abs(outface.balance)
                balance += outface.balance
                set!(outface, 0.0, outface.temperature)
            else
                add!(outface, balance, outface.temperature)
                balance = 0.0
            end
        end
    end

    # write any remaining demand into input bus interfaces (if any) according to input
    # priority, however as available supply is not considered (as this happens implicitly
    # through output priorities of the input bus), this effectively writes all the
    # remaining demand into the first input according to the priority.
    if balance < 0.0
        for inface in bus_infaces(unit)
            add!(inface, balance)
            balance = 0.0
        end
    end

    # if there is a balance unequal zero remaining, this happens either because there is
    # no input bus or the balance was positive and is thus not communicated "backwards" to
    # the input bus. the balance is saved in the remainder so it is available for further
    # balance calculations
    unit.remainder = balance
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