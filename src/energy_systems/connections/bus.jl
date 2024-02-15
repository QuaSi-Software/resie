"""
Utility struct to contain the connections, input/output priorities and other related data
for bus components.
"""
Base.@kwdef mutable struct ConnectionMatrix
    input_order::Vector{String}
    output_order::Vector{String}
    energy_flow::Union{Nothing,Vector{Vector{Bool}}}
end

function ConnectionMatrix(config::Dict{String,Any})::ConnectionMatrix
    if !haskey(config, "connections")
        return new([], [], nothing)
    end

    input_order = [String(u) for u in config["connections"]["input_order"]]
    output_order = [String(u) for u in config["connections"]["output_order"]]

    energy_flow = nothing
    if haskey(config["connections"], "energy_flow")
        energy_flow = []
        for row in config["connections"]["energy_flow"]
            vec = [Bool(v) for v in row]
            push!(energy_flow, vec)
        end
    end

    return ConnectionMatrix(
        input_order,
        output_order,
        energy_flow,
    )
end

Base.@kwdef mutable struct BTInputRow
    source::Component
    priority::Integer
    input_index::Integer
    do_storage_transfer::Bool
    energy_potential::Floathing = nothing
    energy_pool::Floathing = nothing
    storage_potential::Floathing = nothing
    storage_pool::Floathing = nothing
    temperature_min::Temperature = nothing
    temperature_max::Temperature = nothing
end

function reset!(row::BTInputRow)
    row.energy_potential = nothing
    row.energy_pool = nothing
    row.storage_potential = nothing
    row.storage_pool = nothing
    row.temperature_min = nothing
    row.temperature_max = nothing
end

Base.@kwdef mutable struct BTOutputRow
    target::Component
    priority::Integer
    output_index::Integer
    do_storage_transfer::Bool
    energy_potential::Floathing = nothing
    energy_pool::Floathing = nothing
    storage_potential::Floathing = nothing
    storage_pool::Floathing = nothing
    temperature_min::Temperature = nothing
    temperature_max::Temperature = nothing
end

function reset!(row::BTOutputRow)
    row.energy_potential = nothing
    row.energy_pool = nothing
    row.storage_potential = nothing
    row.storage_pool = nothing
    row.temperature_min = nothing
    row.temperature_max = nothing
end

function is_empty(row::Union{BTInputRow, BTOutputRow})
    return (
        row.energy_potential === nothing
        && row.energy_pool === nothing
        && row.storage_potential === nothing
        && row.storage_pool === nothing
    )
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

    balance_table_inputs::Dict{String,BTInputRow}
    balance_table_outputs::Dict{String,BTOutputRow}
    balance_table::Array{Union{Nothing, Float64}, 2}
    proxy::Union{Nothing,Bus}

    epsilon::Float64
end

function Bus(uac::String, config::Dict{String,Any}, sim_params::Dict{String,Any})::Bus
    medium = Symbol(config["medium"])
    register_media([medium])

    return Bus(
        uac, # uac
        controller_for_strategy( # controller
            config["strategy"]["name"], config["strategy"], sim_params
        ),
        sf_bus, # sys_function
        medium, # medium
        [], # input_interfaces
        [], # output_interfaces,
        ConnectionMatrix(config), # connectivity
        0.0, # remainder
        Dict{String,BTInputRow}(), # balance_table_inputs
        Dict{String,BTOutputRow}(), # balance_table_outputs
        Array{Union{Nothing, Float64}, 2}(undef, 0, 0), # balance_table, filled in reset()
        nothing, # proxy
        sim_params["epsilon"] # system-wide epsilon for easy access within the bus functions
    )
end

function Bus(
    uac::String,
    medium::Symbol,
    epsilon::Float64
)
    return Bus(
        uac,
        controller_for_strategy("default", Dict{String,Any}(), Dict{String,Any}()),
        sf_bus,
        medium,
        [], # input_interfaces
        [], # output_interfaces
        ConnectionMatrix([], [], nothing), # connectivity
        0.0, # remainder
        Dict{String,BTInputRow}(), # balance_table_inputs
        Dict{String,BTOutputRow}(), # balance_table_outputs
        Array{Union{Nothing, Float64}, 2}(undef, 0, 0), # balance_table
        nothing, # proxy
        epsilon
    )
end

function initialise!(unit::Bus, sim_params::Dict{String,Any})
    p = 1
    for (idx, inface) in pairs(unit.input_interfaces)
        unit.balance_table_inputs[inface.source.uac] = BTInputRow(
            source=inface.source,
            priority=p,
            input_index=idx,
            do_storage_transfer = inface.do_storage_transfer
        )
        p += 1
    end

    p = 1
    for (idx, outface) in pairs(unit.output_interfaces)
        unit.balance_table_outputs[outface.target.uac] = BTOutputRow(
            target=outface.target,
            priority=p,
            output_index=idx,
            do_storage_transfer = outface.do_storage_transfer
        )
        p += 1
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

    for row in values(unit.balance_table_inputs)
        reset!(row)
    end
    for row in values(unit.balance_table_outputs)
        reset!(row)
    end
   
    unit.balance_table = fill(0.0, (length(unit.balance_table_inputs), 2*length(unit.balance_table_outputs)))
    for i in 1:length(unit.balance_table_inputs)
        for j in 2:2:(2*length(unit.balance_table_outputs))
            unit.balance_table[i, j] = nothing
        end
    end
    
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
            principal = inface.source.output_interfaces[unit.medium]
            blnc += balance(balance_on(principal, principal.source))
        end
    end

    for outface in unit.output_interfaces # demand
        if isa(outface.target, Bus)
            continue
        else
            principal = outface.target.input_interfaces[unit.medium]
            blnc += balance(balance_on(principal, principal.target))
        end
    end

    return blnc + unit.remainder
end

function balance(unit::Bus)::Float64
    if unit.proxy === nothing
        # we can use the non-recursive version of the method as a bus will never
        # be connected to itself... right?
        return balance_nr(unit, unit)
    else
        return balance_direct(unit.proxy)
    end
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
        unit.connectivity.energy_flow === nothing ||
        unit.connectivity.energy_flow[input_idx][output_idx]
    )
end

function _sub(first::Float64, second::Float64) return first - second end
function _sub(first::Nothing, second::Float64) return -second end
function _sub(first::Float64, second::Nothing) return first end
function _sub(first::Nothing, second::Nothing) return nothing end

function _add(first::Float64, second::Float64) return first + second end
function _add(first::Nothing, second::Float64) return second end
function _add(first::Float64, second::Nothing) return first end
function _add(first::Nothing, second::Nothing) return nothing end

function _sum(vector::Vector{Union{Float64, Nothing}})
    sum = nothing
    for entry in vector
        sum = _add(sum, entry)
    end
    return sum
end

function set_max_energy!(unit::Bus, comp::Component, is_input::Bool, value::Float64)
    bus = unit.proxy === nothing ? unit : unit.proxy
    if is_input
        bus.balance_table_inputs[comp.uac].energy_potential = abs(value)
    else
        bus.balance_table_outputs[comp.uac].energy_potential = abs(value)
    end
end

function set_storage_potential!(unit::Bus, comp::Component, is_input::Bool, value::Float64)
    bus = unit.proxy === nothing ? unit : unit.proxy
    if is_input
        bus.balance_table_inputs[comp.uac].storage_potential = abs(value)
    else
        bus.balance_table_outputs[comp.uac].storage_potential = abs(value)
    end
end

function add_balance!(unit::Bus, comp::Component, is_input::Bool, value::Float64)
    bus = unit.proxy === nothing ? unit : unit.proxy
    if is_input
        if comp.sys_function == sf_storage
            bus.balance_table_inputs[comp.uac].storage_pool =
                _add(bus.balance_table_inputs[comp.uac].storage_pool, abs(value))
            bus.balance_table_inputs[comp.uac].storage_potential = nothing
        else
            bus.balance_table_inputs[comp.uac].energy_pool =
                _add(bus.balance_table_inputs[comp.uac].energy_pool, abs(value))
            bus.balance_table_inputs[comp.uac].energy_potential = nothing
        end
    else
        if comp.sys_function == sf_storage
            bus.balance_table_outputs[comp.uac].storage_pool =
                _add(bus.balance_table_outputs[comp.uac].storage_pool, abs(value))
            bus.balance_table_outputs[comp.uac].storage_potential = nothing
        else
            bus.balance_table_outputs[comp.uac].energy_pool =
                _add(bus.balance_table_outputs[comp.uac].energy_pool, abs(value))
            bus.balance_table_outputs[comp.uac].energy_potential = nothing
        end
    end
end

function sub_balance!(unit::Bus, comp::Component, is_input::Bool, value::Float64)
    bus = unit.proxy === nothing ? unit : unit.proxy
    if is_input
        if comp.sys_function == sf_storage
            bus.balance_table_inputs[comp.uac].storage_pool =
                _add(bus.balance_table_inputs[comp.uac].storage_pool, abs(value))
            bus.balance_table_inputs[comp.uac].storage_potential = nothing
        else
            bus.balance_table_inputs[comp.uac].energy_pool =
                _add(bus.balance_table_inputs[comp.uac].energy_pool, abs(value))
            bus.balance_table_inputs[comp.uac].energy_potential = nothing
        end
    else
        if comp.sys_function == sf_storage
            bus.balance_table_outputs[comp.uac].storage_pool =
                _add(bus.balance_table_outputs[comp.uac].storage_pool, abs(value))
            bus.balance_table_outputs[comp.uac].storage_potential = nothing
        else
            bus.balance_table_outputs[comp.uac].energy_pool =
                _add(bus.balance_table_outputs[comp.uac].energy_pool, abs(value))
            bus.balance_table_outputs[comp.uac].energy_potential = nothing
        end
    end
end

function set_balance!(unit::Bus, comp::Component, is_input::Bool, value::Float64)
    bus = unit.proxy === nothing ? unit : unit.proxy
    if is_input
        if comp.sys_function == sf_storage
            bus.balance_table_inputs[comp.uac].storage_pool = abs(value)
            bus.balance_table_inputs[comp.uac].storage_potential = nothing
        else
            bus.balance_table_inputs[comp.uac].energy_pool = abs(value)
            bus.balance_table_inputs[comp.uac].energy_potential = nothing
        end
    else
        if comp.sys_function == sf_storage
            bus.balance_table_outputs[comp.uac].storage_pool = abs(value)
            bus.balance_table_outputs[comp.uac].storage_potential = nothing
        else
            bus.balance_table_outputs[comp.uac].energy_pool = abs(value)
            bus.balance_table_outputs[comp.uac].energy_potential = nothing
        end
    end
end

function set_temperatures!(unit::Bus, comp::Component, is_input::Bool, value_min::Temperature, value_max::Temperature)
    bus = unit.proxy === nothing ? unit : unit.proxy
    if is_input
        bus.balance_table_inputs[comp.uac].temperature_min = value_min
        bus.balance_table_inputs[comp.uac].temperature_max = value_max
    else
        bus.balance_table_outputs[comp.uac].temperature_min = value_min
        bus.balance_table_outputs[comp.uac].temperature_max = value_max
    end
end

function find_interface_on_proxy(proxy::Bus, needle::SystemInterface)
    for inface in proxy.input_interfaces
        if inface.source == needle.source return inface end
    end
    for outface in proxy.output_interfaces
        if outface.target == needle.target return outface end
    end
end

function balance_on(
    interface::SystemInterface,
    unit::Bus
)::Vector{EnergyExchange}
    if unit.proxy !== nothing
        return balance_on(find_interface_on_proxy(unit.proxy, interface), unit.proxy)
    end

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

    inner_distribute!(unit)

    return_exchanges = []

    if caller_is_input
        input_row = [row for row in values(unit.balance_table_inputs) if row.source.uac == interface.source.uac][1]
        for output_row in sort(collect(values(unit.balance_table_outputs)), by=x->x.priority)
            if !energy_flow_is_allowed(unit, input_index, output_row.output_index)
                continue
            end

            is_storage = output_row.target.sys_function == sf_storage

           # if is_storage && !input_row.do_storage_transfer
            #     continue
            # end
            # then it is not neccessary to check this in each component! ToDo

            if interface.max_energy === nothing
                energy_pot = -(_sub(_add(output_row.energy_pool, output_row.energy_potential),
                    (is_storage ? 0.0 : _sum(unit.balance_table[:, output_row.priority*2-1]))))
                storage_pot = -(_sub(_add(output_row.storage_pool, output_row.storage_potential),
                    (is_storage ? _sum(unit.balance_table[:, output_row.priority*2-1]) : 0.0)))
            else
                energy_pot = -(is_storage ? 0.0 : unit.balance_table[input_row.priority, output_row.priority*2-1])
                storage_pot = -(is_storage ? unit.balance_table[input_row.priority, output_row.priority*2-1] : 0.0)
            end

            push!(return_exchanges, EnEx(
                balance=0.0,
                uac=output_row.target.uac,
                energy_potential=energy_pot,
                storage_potential=storage_pot,
                temperature_min=output_row.temperature_min,
                temperature_max=output_row.temperature_max,
                pressure=nothing,
                voltage=nothing
            ))
        end
    else
        output_row = [row for row in values(unit.balance_table_outputs) if row.target.uac == interface.target.uac][1]
        for input_row in sort(collect(values(unit.balance_table_inputs)), by=x->x.priority)
            if !energy_flow_is_allowed(unit, input_row.input_index, output_index)
                continue
            end

            is_storage = input_row.source.sys_function == sf_storage

            # if is_storage && !output_row.do_storage_transfer
            #     continue
            # end
            # then it is not neccessary to check this in each component!  ToDo

            if interface.max_energy === nothing
                energy_pot = _sub(_add(input_row.energy_pool, input_row.energy_potential),
                    (is_storage ? 0.0 : _sum(unit.balance_table[input_row.priority, 1:2:end])))
                storage_pot = _sub(_add(input_row.storage_pool, input_row.storage_potential),
                    (is_storage ? _sum(unit.balance_table[input_row.priority, 1:2:end]) : 0.0))
            else
                energy_pot = is_storage ? 0.0 : unit.balance_table[input_row.priority, output_row.priority*2-1]
                storage_pot = (is_storage ? unit.balance_table[input_row.priority, output_row.priority*2-1] : 0.0)
            end

            push!(return_exchanges, EnEx(
                balance=0.0,
                uac=input_row.source.uac,
                energy_potential=energy_pot,
                storage_potential=storage_pot,
                temperature_min=input_row.temperature_min,
                temperature_max=input_row.temperature_max,
                pressure=nothing,
                voltage=nothing
            ))
        end
    end

    return return_exchanges
end

function inner_distribute!(unit::Bus)
    continue_iteration = true

    for input_row in sort(collect(values(unit.balance_table_inputs)), by=x->x.priority)
        continue_iteration = continue_iteration && !is_empty(input_row)
        if !continue_iteration
            break
        end

        for output_row in sort(collect(values(unit.balance_table_outputs)), by=x->x.priority)
            continue_iteration = continue_iteration && !is_empty(output_row)
            if !continue_iteration
                break
            end

            if !energy_flow_is_allowed(unit, input_row.input_index, output_row.output_index)
                continue
            end

            max_min = highest(input_row.temperature_min, output_row.temperature_min)
            min_max = lowest(input_row.temperature_max, output_row.temperature_max)
            if max_min !== nothing && min_max !== nothing && max_min > min_max
                continue
            end

            if ((input_row.source.sys_function == sf_storage && !output_row.do_storage_transfer) ||
               (output_row.target.sys_function == sf_storage && !input_row.do_storage_transfer))
               continue
            end

            bt_input_row_sum = _sum(unit.balance_table[input_row.priority, 1:2:end])

            available_energy = _sub(_add(_add(_add(
                input_row.energy_potential,
                input_row.energy_pool),
                input_row.storage_potential),
                input_row.storage_pool
            ), bt_input_row_sum)

            bt_output_row_sum = _sum(unit.balance_table[:, output_row.priority*2-1]) 

            target_energy = _sub(_add(_add(_add(
                output_row.energy_potential,
                output_row.energy_pool),
                output_row.storage_potential),
                output_row.storage_pool
            ), bt_output_row_sum)               

            if available_energy < -unit.epsilon || target_energy < -unit.epsilon
                reset_balance_table!(unit::Bus)
                continue_iteration = false
                break
            end

            unit.balance_table[input_row.priority, output_row.priority*2-1] += min(target_energy, available_energy) 
            unit.balance_table[input_row.priority, output_row.priority*2] = max_min
        end
    end
end

function reset_balance_table!(unit::Bus)

    unit.balance_table = fill(0.0, (length(unit.balance_table_inputs), 2*length(unit.balance_table_outputs)))
    for i in 1:length(unit.balance_table_inputs)
        for j in 2:2:(2*length(unit.balance_table_outputs))
            unit.balance_table[i, j] = nothing
        end
    end

    inner_distribute!(unit)
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
    if unit.proxy !== nothing
        return # the proxy has its own distribute step
    end

    inner_distribute!(unit::Bus)
    balance = balance_direct(unit)

    # reset all non-bus input interfaces
    for inface in filter_inputs(unit, sf_bus, false)
        principal = inface.source.output_interfaces[unit.medium]
        set!(principal, 0.0)
    end

    # reset all non-bus output interfaces
    for outface in filter_outputs(unit, sf_bus, false)
        principal = outface.target.input_interfaces[unit.medium]
        set!(principal, 0.0)
    end

    # distribute to outgoing busses according to output priority
    if balance > 0.0
        for outface in filter_outputs(unit, sf_bus, true)
            if balance > abs(outface.balance)
                balance -= abs(outface.balance)
                set!(outface, 0.0)
            else
                set!(outface, outface.balance - balance)
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
            set!(inface, balance + inface.balance)
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

# merge functionality is in an extra file
include("./bus_merging.jl")

export Bus