"""
Utility struct to contain the connections, input/output priorities and other related data
for bus components.
"""
Base.@kwdef mutable struct ConnectionMatrix
    input_order::Vector{String}
    output_order::Vector{String}
    storage_loading::Union{Nothing,Vector{Vector{Bool}}}

    function ConnectionMatrix(config::Dict{String,Any})
        input_order = []
        output_order = [String(u) for u in config["output_refs"]]
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
Imnplementation of a bus component for balancing multiple inputs and outputs.

This component is both a possible real energy system component (mostly for electricity) as well as a
necessary abstraction of the model. The basic idea is that one or more components feed
energy of the same medium into a bus and one or more components draw that energy from
the bus. A bus with only one input and only one output can be replaced with a direct
connection between both components.

The function and purpose is described in more detail in the accompanying documentation.
"""
Base.@kwdef mutable struct Bus <: Component
    uac::String
    controller::Controller
    sys_function::SystemFunction
    medium::Symbol

    input_interfaces::Vector{SystemInterface}
    output_interfaces::Vector{SystemInterface}
    connectivity::ConnectionMatrix

    remainder::Float64

    function Bus(uac::String, config::Dict{String,Any}, sim_params::Dict{String,Any})
        medium = Symbol(config["medium"])
        register_media([medium])

        return new(
            uac, # uac
            controller_for_strategy( # controller
                config["strategy"]["name"], config["strategy"], sim_params
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

Variant of [`balance`](@ref) that includes other connected bus components and their energy
balance, but does so in a non-recursive manner such that any bus in the chain of connected
bus components is only considered once.
"""
function balance_nr(unit::Bus, caller::Bus)::Float64
    blnc = 0.0

    for inface in unit.input_interfaces # supply
        if inface.source == caller
            continue
        end

        if isa(inface.source, Bus)
            other_bus_balance = balance_nr(inface.source, unit)
            balance_supply = max(other_bus_balance, inface.balance)
            if balance_supply < 0.0
                continue
            end
        else
            balance_supply = balance(balance_on(inface, inface.source))
        end
        blnc += balance_supply
    end

    for outface in unit.output_interfaces # demand
        if outface.target == caller
            continue
        end

        if isa(outface.target, Bus)
            other_bus_balance = balance_nr(outface.target, unit)
            balance_demand = min(other_bus_balance, outface.balance)
            if balance_demand > 0.0
                continue
            end
        else
            balance_demand = balance(balance_on(outface, outface.target))
        end
        blnc += balance_demand
    end

    return blnc + unit.remainder
end

"""
    balance_direct(unit)

Energy balance on a bus component without considering any other connected bus components.
"""
function balance_direct(unit::Bus)::Float64
    blnc = 0.0

    for inface in unit.input_interfaces # supply
        if isa(inface.source, Bus)
            continue
        else
            blnc += balance(balance_on(inface, inface.source))
        end
    end

    for outface in unit.output_interfaces # demand
        if isa(outface.target, Bus)
            continue
        else
            blnc += balance(balance_on(outface, outface.target))
        end
    end

    return blnc + unit.remainder
end

function balance(unit::Bus)::Float64
    # we can use the non-recursive version of the method as a bus will never
    # be connected to itself... right?
    return balance_nr(unit, unit)
end

"""
    energy_flow_is_allowed(bus, input_idx, output_idx)

Checks the connectivity matrix of the bus as returns if energy is allowed to flow from the
input with the given index to the output with the given index.

Args:
    `unit::Bus`: The bus to check
    `input_idx::Integer`: Input index
    `output_idx::Integer`: Output index
Returns:
    `Bool`: True if the flow is allowed, false otherwise
"""
function energy_flow_is_allowed(unit::Bus, input_idx::Integer, output_idx::Integer)::Bool
    return (
        unit.connectivity.storage_loading === nothing ||
        unit.connectivity.storage_loading[input_idx][output_idx]
    )
end

function balance_on(
    interface::SystemInterface,
    unit::Bus
)::Vector{EnergyExchange}
    input_index = nothing
    output_index = nothing
    caller_is_input = false # if interface is input of unit (caller puts energy in unit)
    caller_is_output = false # if interface is output of unit (caller gets energy from unit)

    # determine if the calling component is an input or output to the bus and remember the
    # index within the list of input/output interfaces for later
    for (idx, input_interface) in pairs(unit.input_interfaces)
        if input_interface.source.uac == interface.source.uac
            input_index = idx
            caller_is_input = true
            break
        end
    end

    for (idx, output_interface) in pairs(unit.output_interfaces)
        if output_interface.target.uac == interface.target.uac
            output_index = idx
            caller_is_output = true
            break
        end
    end

    # sanity check, as this situation should not happen
    if (caller_is_input && caller_is_output) || (!caller_is_input && !caller_is_output)
        throw(ArgumentError(
            "Error in connnection of components on bus \"$(unit.uac)\". " * 
            "Caller must be input XOR output."
        ))
    end

    return_exchanges = []

    for (idx, outface) in pairs(unit.output_interfaces)
        if caller_is_input && !energy_flow_is_allowed(unit, input_index, idx)
            continue
        end

        if isa(outface.target, Bus)
            # don't recurse back into the bus which called balance_on
            if caller_is_output && idx == output_index; continue end

            exchanges = balance_on(outface, outface.target)
            if caller_is_output
                # for other outputs on the bus, the entire chain of components and busses
                # behind the output interface we are currently checking is only relevant for
                # balance calculations. therefore we can replace it with a single exchange
                # that is capped up at zero (since energy can't flow back into the
                # current bus)
                push!(return_exchanges, EnEx(
                    balance=min(0.0, balance(exchanges)),
                    uac=outface.target.uac,
                    energy_potential=0.0,
                    storage_potential=0.0,
                    temperature=nothing,
                    pressure=nothing,
                    voltage=nothing
                ))
            else
                # for inputs on the bus, we can add all exchanges of the outgoing bus to
                # the list of exchanges. the potentials are filtered within the balance_on
                # calculation of the outgoing bus
                append!(return_exchanges, exchanges)
            end
        else
            exchanges = balance_on(outface, outface.target)

            # check storage potential only for outgoing storages and make sure storages
            # don't load themselves. also the storage potential is only
            # considered if no energy was transfered over the interface yet
            if (
                caller_is_input &&
                outface.target.sys_function === sf_storage &&
                outface.target.uac !== interface.source.uac &&
                outface.sum_abs_change == 0.0 && interface.sum_abs_change == 0.0
            )
                storage_pot = storage_potential(exchanges)
            else
                storage_pot = 0.0
            end

            # if energy was already transfered over the interface or no information is
            # available or we're checking other outputs, set the energy potential
            # to zero
            if (
                caller_is_output
                || outface.max_energy === nothing
                || outface.sum_abs_change > 0.0
                || interface.sum_abs_change > 0.0
            )
                energy_pot = 0.0
            else
                energy_pot = energy_potential(exchanges)
            end

            temperature = outface.temperature

            push!(return_exchanges, EnEx(
                balance=balance(exchanges),
                uac=exchanges[1].uac,
                energy_potential=energy_pot,
                storage_potential=storage_pot,
                temperature=temperature,
                pressure=nothing,
                voltage=nothing
            ))
        end
    end

    for (idx, inface) in pairs(unit.input_interfaces)
        if caller_is_output && !energy_flow_is_allowed(unit, idx, output_index)
            continue
        end

        if isa(inface.source, Bus)
            # don't recurse back into the bus which called balance_on
            if caller_is_input && idx == input_index; continue end

            exchanges = balance_on(inface, inface.source)
            if caller_is_input
                # for other inputs on the bus, the entire chain of components and busses
                # behind the input interface we are currently checking is only relevant for
                # balance calculations. therefore we can replace it with a single exchange
                # that is capped down at zero (since energy can't flow back into the
                # incoming bus)
                push!(return_exchanges, EnEx(
                    balance=max(0.0, balance(exchanges)),
                    uac=inface.source.uac,
                    energy_potential=0.0,
                    storage_potential=0.0,
                    temperature=nothing,
                    pressure=nothing,
                    voltage=nothing
                ))
            else
                # for outputs on the bus, we can add all exchanges of the incoming bus to
                # the list of exchanges. the potentials are filtered within the balance_on
                # calculation of the incoming bus
                append!(return_exchanges, exchanges)
            end
        else
            exchanges = balance_on(inface, inface.source)

            # check storage potential only for incoming storages and make sure storages
            # don't load themselves. also the storage potential is only
            # considered if no energy was transfered over the interface yet
            if (
                caller_is_output &&
                inface.source.sys_function === sf_storage &&
                inface.source.uac !== interface.target.uac &&
                inface.sum_abs_change == 0.0 && interface.sum_abs_change == 0.0
            )
                storage_pot = storage_potential(exchanges)
            else
                storage_pot = 0.0
            end

            # if energy was already transfered over the interface or no information is
            # available or we're checking other inputs, set the energy potential
            # to zero
            if (
                caller_is_input
                || inface.max_energy === nothing
                || inface.sum_abs_change > 0.0
                || interface.sum_abs_change > 0.0
            )
                energy_pot = 0.0
            else
                energy_pot = energy_potential(exchanges)
            end

            temperature = inface.temperature

            push!(return_exchanges, EnEx(
                balance=balance(exchanges),
                uac=exchanges[1].uac,
                energy_potential=energy_pot,
                storage_potential=storage_pot,
                temperature=temperature,
                pressure=nothing,
                voltage=nothing
            ))
        end
    end

    return return_exchanges
end

"""
    filter_inputs(unit, condition, inclusive)

Filters the input interfaces of the given bus based on a condition on the system function
of the components on the source side of the interfaces. The condition can be negated with
the argument `inclusive`.

Args:
    `unit::Bus`: The bus to check
    `condition::SystemFunction`: Determines which components to filter in/out
    `inclusive::Bool`: If true, the return list will include only interfaces to which the
        condition applies. If false, the return list will include only interfaces to which
        condition does not apply.
Returns:
    `Vector{SystemInterface}`: The filtered list of input interfaces of the bus
"""
function filter_inputs(unit::Bus, condition::SystemFunction, inclusive::Bool)
    return [f for f in unit.input_interfaces
        if (inclusive && f.source.sys_function == condition)
            || (!inclusive && f.source.sys_function != condition)
    ]
end

"""
    filter_outputs(unit, condition, inclusive)

Filters the output interfaces of the given bus based on a condition on the system function
of the components on the target side of the interfaces. The condition can be negated with
the argument `inclusive`.

Args:
    `unit::Bus`: The bus to check
    `condition::SystemFunction`: Determines which components to filter in/out
    `inclusive::Bool`: If true, the return list will include only interfaces to which the
        condition applies. If false, the return list will include only interfaces to which
        condition does not apply.
Returns:
    `Vector{SystemInterface}`: The filtered list of output interfaces of the bus
"""
function filter_outputs(unit::Bus, condition::SystemFunction, inclusive::Bool)
    return [f for f in unit.output_interfaces
        if (inclusive && f.target.sys_function == condition)
            || (!inclusive && f.target.sys_function != condition)
    ]
end

"""
    distribute!(unit)

Bus-specific implementation of distribute!.

This moves the energy from connected component from supply to demand components both
on the bus directly as well as taking other bus components into account. This allows busses
to be connected in chains (but not loops) and "communicate" the energy across. The method
implicitly requires that each bus on the chain is called with distribute!() in a specific
order, which is explained in more detail in the documentation. Essentially it starts from
the leaves of the chain and progresses to the roots.
"""
function distribute!(unit::Bus)
    balance = balance_direct(unit)

    # reset all non-bus input interfaces
    for inface in filter_inputs(unit, sf_bus, false)
        set!(inface, 0.0, inface.temperature)
    end

    # reset all non-bus output interfaces
    for outface in filter_outputs(unit, sf_bus, false)
        set!(outface, 0.0, outface.temperature)
    end

    # distribute to outgoing busses according to output priority
    if balance > 0.0
        for outface in filter_outputs(unit, sf_bus, true)
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
        for inface in filter_inputs(unit, sf_bus, true)
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