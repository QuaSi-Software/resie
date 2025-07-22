using ..Resie: get_run
using UUIDs

"""
Utility struct to contain the connections, input/output priorities and other related data
for bus components.
"""
Base.@kwdef mutable struct ConnectionMatrix
    input_order::Vector{String}
    output_order::Vector{String}
    energy_flow::Union{Nothing,Vector{Vector{Bool}}}
end

"""
    ConnectionMatrix(config)

Constructor for ConnectionMatrix from a config in the user-input project file.

# Arguments
`config::Dict{String,Any}`: The component config

# Returns
`ConnectionMatrix`: The constructed instance
"""
function ConnectionMatrix(config::Dict{String,Any})::ConnectionMatrix
    if !haskey(config, "connections")
        return ConnectionMatrix([], [], nothing)
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

    return ConnectionMatrix(input_order,
                            output_order,
                            energy_flow)
end

"""
Container struct for information on an input of a bus relevant to the balance table.
"""
Base.@kwdef mutable struct BTInputRow
    source::Component
    priority::Integer
    input_index::Integer
    do_storage_transfer::Bool
    energy_potential::MaxEnergy = MaxEnergy()
    energy_potential_temp::MaxEnergy = MaxEnergy()
    energy_pool::MaxEnergy = MaxEnergy()
    energy_pool_temp::MaxEnergy = MaxEnergy()
end

"""
    reset!(bt_input_row)

Resets the volatile fields of the balance table input row back to zero/nothing.

# Arguments
`row::BTInputRow`: The row to reset
"""
function reset!(row::BTInputRow)
    row.energy_potential = MaxEnergy()
    row.energy_potential_temp = MaxEnergy()
    row.energy_pool = MaxEnergy()
    row.energy_pool_temp = MaxEnergy()
end

"""
Container struct for information on an output of a bus relevant to the balance table.
"""
Base.@kwdef mutable struct BTOutputRow
    target::Component
    priority::Integer
    output_index::Integer
    do_storage_transfer::Bool
    energy_potential::MaxEnergy = MaxEnergy()
    energy_potential_temp::MaxEnergy = MaxEnergy()
    energy_pool::MaxEnergy = MaxEnergy()
    energy_pool_temp::MaxEnergy = MaxEnergy()
end

"""
    reset!(bt_output_row)

Resets the volatile fields of the balance table output row back to zero/nothing.

# Arguments
`row::BTOutputRow`: The row to reset
"""
function reset!(row::BTOutputRow)
    row.energy_potential = MaxEnergy()
    row.energy_potential_temp = MaxEnergy()
    row.energy_pool = MaxEnergy()
    row.energy_pool_temp = MaxEnergy()
end

"""
    is_empty(balance_table_row)

Checks if the given input/output row is "empty".

Empty in this regard means if an energy potential or pool (utilized energy) has been
written.

# Arguments
`row::Union{BTInputRow, BTOutputRow}`: The row to check

# Returns
`Bool`: If the row is empty or not
"""
function is_empty(row::Union{BTInputRow,BTOutputRow})::Bool
    return (is_max_energy_nothing(row.energy_potential)
            &&
            is_max_energy_nothing(row.energy_pool))
end

"""
Implementation of a bus component for balancing multiple inputs and outputs.

This component is both a possible real energy system component (mostly for electricity) as
well as a necessary abstraction of the model. The basic idea is that one or more components
feed energy of the same medium into a bus and one or more components draw that energy from
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
    balance_table::Array{Union{Nothing,Float64},2}
    proxy::Union{Nothing,Bus}

    run_id::UUID

    epsilon::Float64
end

"""
    Bus(uac, config, sim_params)

Config-constructor for a Bus.

# Arguments
`uac::String`: The UAC of the new bus
`config::Dict{String,Any}`: The config from the project file
`sim_params::Dict{String,Any}`: Simulation parameters

# Returns
`Bus`: The constructed bus
"""
function Bus(uac::String, config::Dict{String,Any}, sim_params::Dict{String,Any})::Bus
    medium = Symbol(config["medium"])
    register_media([medium])

    return Bus(uac,                                          # uac
               Controller(default(config, "control_parameters", nothing)),
               sf_bus,                                       # sys_function
               medium,                                       # medium
               [],                                           # input_interfaces
               [],                                           # output_interfaces,
               ConnectionMatrix(config),                     # connectivity
               0.0,                                          # remainder
               Dict{String,BTInputRow}(),                    # balance_table_inputs
               Dict{String,BTOutputRow}(),                   # balance_table_outputs
               Array{Union{Nothing,Float64},2}(undef, 0, 0), # balance_table, filled in reset()
               nothing,                                      # proxy
               uuid1(),                                      # holds the current run ID later
               sim_params["epsilon"])                        # system-wide epsilon for easy access within the bus functions
end

"""
    Bus(uac, medium, epsilon)

Constructor for a Bus that creates a mostly empty bus with the minimal parameters.

# Arguments
`uac::String`: The UAC of the new bus
`medium::Symbol`: The medium of the bus
`epsilon::Float64`: Simulation parameter for epsilon. Can also just be a small number,
    e.g. a value of 1e-9

# Returns
`Bus`: The constructed bus
"""
function Bus(uac::String,
             medium::Symbol,
             run_id::UUID,
             epsilon::Float64)
    return Bus(uac,
               Controller(nothing),
               sf_bus,
               medium,
               [],                                           # input_interfaces
               [],                                           # output_interfaces
               ConnectionMatrix([], [], nothing),            # connectivity
               0.0,                                          # remainder
               Dict{String,BTInputRow}(),                    # balance_table_inputs
               Dict{String,BTOutputRow}(),                   # balance_table_outputs
               Array{Union{Nothing,Float64},2}(undef, 0, 0), # balance_table
               nothing,                                      # proxy
               run_id,                                       # run ID
               epsilon)
end

function initialise!(unit::Bus, sim_params::Dict{String,Any})
    p = 1
    for (idx, inface) in pairs(unit.input_interfaces)
        unit.balance_table_inputs[inface.source.uac] = BTInputRow(; source=inface.source,
                                                                  priority=p,
                                                                  input_index=idx,
                                                                  do_storage_transfer=inface.do_storage_transfer)
        p += 1
    end

    p = 1
    for (idx, outface) in pairs(unit.output_interfaces)
        unit.balance_table_outputs[outface.target.uac] = BTOutputRow(; target=outface.target,
                                                                     priority=p,
                                                                     output_index=idx,
                                                                     do_storage_transfer=outface.do_storage_transfer)
        p += 1
    end

    unit.run_id = sim_params["run_ID"]
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

    reset_balance_table!(unit::Bus, false)
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
    # if there is a proxy bus, the balance is only correct on its calculation, which happens
    # separately. if there is no proxy, there are also no bus components connected to the
    # given bus
    return balance_direct(unit)
end

"""
energy_flow_is_denied(bus, input_row, output_row)

Checks the connectivity matrix of the bus and returns if energy is allowed to flow from the
input_row to the output_row. Also checks the storage loading control specified by the component
in the interface and denies self-feeding of components and grids feeding into other grids.

Args:
    `unit::Bus`: The bus to check
    `input_row::BTInputRow`: Input row
    `output_row::BTOutputRow`: Output row
Returns:
    `Bool`: True if the flow is denied, false otherwise
"""
function energy_flow_is_denied(unit::Bus, input_row::BTInputRow, output_row::BTOutputRow)::Bool
    return (!(unit.connectivity.energy_flow === nothing ||                                       # check energy_flow matrix
              unit.connectivity.energy_flow[input_row.input_index][output_row.output_index]) ||  # check energy_flow matrix
            (output_row.target.sys_function == sf_storage && !input_row.do_storage_transfer) ||  # check storage loading control
            (input_row.source.sys_function == sf_storage && !output_row.do_storage_transfer) ||  # check storage unloading control
            output_row.target.uac == input_row.source.uac ||                                     # do not allow self-feeding of any component
            (get_max_energy(output_row.energy_potential, input_row.source.uac) == Inf &&         # do not allow grids to feed into grids
             get_max_energy(input_row.energy_potential, output_row.target.uac) == Inf))          # do not allow grids to feed into grids
end

"""
    set_max_energy!(bus, component, is_input, energy, purpose_uac, has_calculated_all_maxima, recalculate_max_energy)

Communicates the max_energy of an input/output on a bus.

This is required for the balance calculations on the bus chain and is typically called from
the set_max_energy! function on an interface, if one side is a bus.

# Arguments
`unit::Bus`: The bus
`comp::Component`: The component that is an input/output
`is_input::Bool`: If the component is an input
`energy::Float64`: The energy of max_energy
"""
function set_max_energy!(unit::Bus,
                         comp::Component,
                         is_input::Bool,
                         energy::Union{Floathing,Vector{<:Floathing}},
                         temperature_min::Union{Temperature,Vector{<:Temperature}},
                         temperature_max::Union{Temperature,Vector{<:Temperature}},
                         purpose_uac::Union{Stringing,Vector{<:Stringing}},
                         has_calculated_all_maxima::Bool,
                         recalculate_max_energy::Bool)
    bus = unit.proxy === nothing ? unit : unit.proxy

    if is_input
        set_max_energy!(bus.balance_table_inputs[comp.uac].energy_potential,
                        _abs(energy),
                        temperature_min,
                        temperature_max,
                        purpose_uac,
                        has_calculated_all_maxima,
                        recalculate_max_energy)
    else
        set_max_energy!(bus.balance_table_outputs[comp.uac].energy_potential,
                        _abs(energy),
                        temperature_min,
                        temperature_max,
                        purpose_uac,
                        has_calculated_all_maxima,
                        recalculate_max_energy)
    end

    if unit.proxy !== nothing
        proxy_interface = is_input ?
                          bus.input_interfaces[bus.balance_table_inputs[comp.uac].input_index] :
                          bus.output_interfaces[bus.balance_table_outputs[comp.uac].output_index]
        set_max_energy!(proxy_interface.max_energy,
                        energy,
                        temperature_min,
                        temperature_max,
                        purpose_uac,
                        has_calculated_all_maxima,
                        recalculate_max_energy)
    end
end

"""
    add_balance!(unit, component, is_input, energy, temperature, purpose_uac)

Communicates a balance addition of an input/output on a bus.

This is required for the balance calculations on the bus chain and is typically called from
the add! function on an interface, if one side is a bus.

# Arguments
`unit::Bus`: The bus
`comp::Component`: The component that is an input/output
`is_input::Bool`: If the component is an input (=true)
`energy::Union{Floathing, Vector{<:Floathing}}`: The energy to add to the balance
`temperature::Union{Temperature,Vector{<:Temperature}}`: The temperature of the energy, if present
`purpose_uac::Union{Stringing, Vector{<:Stringing}}`: The purpose uac, if given
"""
function add_balance!(unit::Bus,
                      comp::Component,
                      is_input::Bool,
                      energy::Union{Floathing,Vector{<:Floathing}},
                      temperature_min::Union{Temperature,Vector{<:Temperature}},
                      temperature_max::Union{Temperature,Vector{<:Temperature}},
                      purpose_uac::Union{Stringing,Vector{<:Stringing}}=nothing)
    bus = unit.proxy === nothing ? unit : unit.proxy
    if is_input
        increase_max_energy!(bus.balance_table_inputs[comp.uac].energy_pool,
                             energy,
                             temperature_min,
                             temperature_max,
                             purpose_uac)
        bus.balance_table_inputs[comp.uac].energy_potential = MaxEnergy()
    else
        increase_max_energy!(bus.balance_table_outputs[comp.uac].energy_pool,
                             energy,
                             temperature_min,
                             temperature_max,
                             purpose_uac)
        bus.balance_table_outputs[comp.uac].energy_potential = MaxEnergy()
    end
end

"""
    sub_balance!(unit, component, is_input, energy, temperature, purpose_uac)

Communicates a balance subtraction of an input/output on a bus.

This is required for the balance calculations on the bus chain and is typically called from
the sub! function on an interface, if one side is a bus.

# Arguments
`unit::Bus`: The bus
`comp::Component`: The component that is an input/output
`is_input::Bool`: If the component is an input (=true)
`energy::Union{Float64, Vector{Float64}}`: The energy to add to the balance
`temperature::Union{Temperature,Vector{<:Temperature}}`: The temperature of the energy, if present
`purpose_uac::Union{Stringing, Vector{Stringing}}`: The purpose uac if given
"""
function sub_balance!(unit::Bus,
                      comp::Component,
                      is_input::Bool,
                      energy::Union{Floathing,Vector{<:Floathing}},
                      temperature_min::Union{Temperature,Vector{<:Temperature}},
                      temperature_max::Union{Temperature,Vector{<:Temperature}},
                      purpose_uac::Union{Stringing,Vector{Stringing}}=nothing)
    bus = unit.proxy === nothing ? unit : unit.proxy
    if is_input
        increase_max_energy!(bus.balance_table_inputs[comp.uac].energy_pool,
                             abs.(energy),
                             temperature_min,
                             temperature_max,
                             purpose_uac)
        bus.balance_table_inputs[comp.uac].energy_potential = MaxEnergy()
    else
        increase_max_energy!(bus.balance_table_outputs[comp.uac].energy_pool,
                             abs.(energy),
                             temperature_min,
                             temperature_max,
                             purpose_uac)
        bus.balance_table_outputs[comp.uac].energy_potential = MaxEnergy()
    end
end

"""
    find_interface_on_proxy(proxy, interface)

Finds the interface on the given proxy bus that is the same as the given interface.

Interfaces are considered to be the same if the source/target is the same and the proxy
bus is on the corresponding other side. This mechanism is required because proxy busses
create new interfaces that point to the non-bus components of the principal busses.

# Arguments
`proxy::Bus`: The proxy bus
`needle::SystemInterface`: The interface to find

# Returns
`Union{Nothing,SystemInterface}`: The corresponding interface or nothing, if it can't be
    found.
"""
function find_interface_on_proxy(proxy::Bus,
                                 needle::SystemInterface)::Union{Nothing,SystemInterface}
    for inface in proxy.input_interfaces
        if inface.source == needle.source
            return inface
        end
    end
    for outface in proxy.output_interfaces
        if outface.target == needle.target
            return outface
        end
    end
    return nothing
end

"""
    balance_on(interface, bus)

Returns the energy exchange information on the given interface.

This is the same function as the generic balance_on for all components, but functions
vastly differently due to busses having a special role in the energy system simulation. A
bus keeps an internal store of how much energy is available and requested and updates this
store every time balance_on is called. This makes this function expensive to call, but also
very important for the correct energy flow calculations.

Because components, on a bus, communicate via balance_on with other components, the internal
balance calculation in some sense replaces the individual calculations of each component.

If a chain of busses has been set up with a proxy, the principal busses relay the call to
balance_on to the proxy bus and return the results. From the perspective of the non-bus
components this behaves the same, but allows for communication and energy flow across all
busses of the chain.

# Arguments
`interface::SystemInterface`: The interface connecting the component and the bus
`unit::Bus`: The bus

# Returns
`Vector{EnergyExchange}`: A list of energy exchanges, each of which encode one potential
    source or target for the component requesting a balance calculation.
"""
function balance_on(interface::SystemInterface, unit::Bus)::Vector{EnergyExchange}
    if unit.proxy !== nothing
        proxy_interface = find_interface_on_proxy(unit.proxy, interface)
        return proxy_interface === nothing ? [] : balance_on(proxy_interface, unit.proxy)
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
        throw(ArgumentError("Error in connection of components on bus \"$(unit.uac)\". " *
                            "Caller must be input XOR output."))
    end

    inner_distribute!(unit)

    return_exchanges = []

    if caller_is_input
        input_row = [row for row in values(unit.balance_table_inputs) if row.source.uac == interface.source.uac][1]
        for output_row in sort(collect(values(unit.balance_table_outputs)); by=x -> x.priority)
            if energy_flow_is_denied(unit, input_row, output_row)
                continue
            end

            #   target is transformer that has not been calculated its potential or process...
            if (output_row.target.sys_function === EnergySystems.sf_transformer
                &&
                is_max_energy_nothing(unit.balance_table_outputs[output_row.target.uac].energy_potential)
                &&
                is_max_energy_nothing(unit.balance_table_outputs[output_row.target.uac].energy_pool)
                # ...or has Inf written in its interface
                || get_max_energy(output_row.energy_pool, input_row.source.uac) == Inf
                || get_max_energy(output_row.energy_potential, input_row.source.uac) == Inf)
                # end of condition
                energy_pot = -Inf
                temperature_min = get_min_temperature(output_row.energy_potential, output_row.energy_pool,
                                                      input_row.source.uac)
                temperature_max = get_max_temperature(output_row.energy_potential, output_row.energy_pool,
                                                      input_row.source.uac)
            else
                if is_max_energy_nothing(interface.max_energy)  # the caller has not performed a potential
                    energy_pot = -_add(get_max_energy(output_row.energy_pool_temp, input_row.source.uac),
                                       get_max_energy(output_row.energy_potential_temp, input_row.source.uac))
                    temperature_min = get_min_temperature(output_row.energy_potential_temp, output_row.energy_pool_temp,
                                                          input_row.source.uac)
                    temperature_max = get_max_temperature(output_row.energy_potential_temp, output_row.energy_pool_temp,
                                                          input_row.source.uac)
                else  # the caller has performed a potential and has already written a max_energy itself
                    energy_pot = -(unit.balance_table[input_row.priority, output_row.priority * 2 - 1])
                    temperature_min = unit.balance_table[input_row.priority, output_row.priority * 2]
                    temperature_max = temperature_min
                end
            end

            if energy_pot < 0.0
                push!(return_exchanges,
                      EnEx(; balance=0.0,
                           energy_potential=energy_pot,
                           purpose_uac=output_row.target.uac,
                           temperature_min=temperature_min,
                           temperature_max=temperature_max,
                           pressure=nothing,
                           voltage=nothing))
            end
        end
    else
        output_row = [row for row in values(unit.balance_table_outputs) if row.target.uac == interface.target.uac][1]
        for input_row in sort(collect(values(unit.balance_table_inputs)); by=x -> x.priority)
            if energy_flow_is_denied(unit, input_row, output_row)
                continue
            end

            #   source is transformer that has not been calculated its potential or process...
            if (input_row.source.sys_function === EnergySystems.sf_transformer
                &&
                is_max_energy_nothing(unit.balance_table_inputs[input_row.source.uac].energy_potential)
                &&
                is_max_energy_nothing(unit.balance_table_inputs[input_row.source.uac].energy_pool)
                # ...or has Inf written in its interface
                || get_max_energy(input_row.energy_pool, output_row.target.uac) == Inf
                || get_max_energy(input_row.energy_potential, output_row.target.uac) == Inf)
                # end of condition
                energy_pot = Inf
                temperature_min = get_min_temperature(input_row.energy_potential, input_row.energy_pool,
                                                      output_row.target.uac)
                temperature_max = get_max_temperature(input_row.energy_pool, input_row.energy_potential,
                                                      output_row.target.uac)
            else
                if is_max_energy_nothing(interface.max_energy) # no max energy is written --> get infos from input_row
                    energy_pot = _add(get_max_energy(input_row.energy_pool_temp, output_row.target.uac),
                                      get_max_energy(input_row.energy_potential_temp, output_row.target.uac))
                    temperature_min = get_min_temperature(input_row.energy_potential_temp, input_row.energy_pool_temp,
                                                          output_row.target.uac)
                    temperature_max = get_max_temperature(input_row.energy_potential_temp, input_row.energy_pool_temp,
                                                          output_row.target.uac)
                else # max energy is written and was distributed by the bus --> get infos from balance_table
                    energy_pot = unit.balance_table[input_row.priority, output_row.priority * 2 - 1]
                    temperature_min = unit.balance_table[input_row.priority, output_row.priority * 2]
                    temperature_max = temperature_min
                end
            end

            if energy_pot > 0.0
                push!(return_exchanges,
                      EnEx(; balance=0.0,
                           energy_potential=energy_pot,
                           purpose_uac=input_row.source.uac,
                           temperature_min=temperature_min,
                           temperature_max=temperature_max,
                           pressure=nothing,
                           voltage=nothing))
            end
        end
    end

    if isempty(return_exchanges)
        push!(return_exchanges,
              EnEx(; balance=0.0,
                   energy_potential=0.0,
                   purpose_uac=nothing,
                   temperature_min=nothing,
                   temperature_max=nothing,
                   pressure=nothing,
                   voltage=nothing))
    end

    return return_exchanges
end

"""
    inner_distribute!(bus)

Perform energy distribution calculation on a bus.

This function should only be called within the balance_on function for a bus. It serves no
other purpose.

# Arguments
`unit::Bus`: The bus for which to calculate energy distribution
"""
function inner_distribute!(unit::Bus)
    reset_balance_table!(unit::Bus, false)
    continue_iteration = true

    for input_row in sort(collect(values(unit.balance_table_inputs)); by=x -> x.priority)
        continue_iteration = continue_iteration && !is_empty(input_row)
        if !continue_iteration
            break
        end

        for output_row in sort(collect(values(unit.balance_table_outputs)); by=x -> x.priority)
            if is_empty(output_row)
                break
            end

            if energy_flow_is_denied(unit, input_row, output_row)
                continue
            end

            temperature_highest_min = highest(get_min_temperature(input_row.energy_potential_temp,
                                                                  input_row.energy_pool_temp, output_row.target.uac),
                                              get_min_temperature(output_row.energy_potential_temp,
                                                                  output_row.energy_pool_temp, input_row.source.uac))
            temperature_lowest_max = lowest(get_max_temperature(input_row.energy_potential_temp,
                                                                input_row.energy_pool_temp, output_row.target.uac),
                                            get_max_temperature(output_row.energy_potential_temp,
                                                                output_row.energy_pool_temp, input_row.source.uac))

            if temperature_highest_min !== nothing && temperature_lowest_max !== nothing &&
               temperature_highest_min > temperature_lowest_max
                continue
            end

            available_energy = _add(get_max_energy(input_row.energy_potential_temp, output_row.target.uac),
                                    get_max_energy(input_row.energy_pool_temp, output_row.target.uac))

            target_energy = _add(get_max_energy(output_row.energy_potential_temp, input_row.source.uac),
                                 get_max_energy(output_row.energy_pool_temp, input_row.source.uac))

            if available_energy < -unit.epsilon || target_energy < -unit.epsilon
                reset_balance_table!(unit::Bus, true)
                continue_iteration = false
                break
            end

            energy_flow = min(target_energy, available_energy)
            unit.balance_table[input_row.priority, output_row.priority * 2 - 1] += energy_flow
            # if both min and max temperature are given and differing (can currently only happen during HP bypass), 
            # the lowest temperature is set:
            unit.balance_table[input_row.priority, output_row.priority * 2] = lowest(temperature_highest_min,
                                                                                     temperature_lowest_max)

            # extract the already distributed energy from balance table
            already_distributed_input_energies = Float64.(unit.balance_table[input_row.priority,
                                                                             1:2:(output_row.priority * 2 - 1)])
            already_distributed_input_temperatures = unit.balance_table[input_row.priority,
                                                                        2:2:(output_row.priority * 2)]

            already_distributed_output_energies = Float64.(unit.balance_table[1:(input_row.priority),
                                                                              output_row.priority * 2 - 1])
            already_distributed_output_temperatures = unit.balance_table[1:(input_row.priority),
                                                                         output_row.priority * 2]

            if energy_flow !== 0.0
                if !is_max_energy_nothing(input_row.energy_potential_temp)
                    reduce_max_energy!(input_row.energy_potential_temp,
                                       already_distributed_input_energies,
                                       already_distributed_input_temperatures,
                                       output_row.target.uac,
                                       input_row.source.uac,
                                       unit.run_id)
                end
                if !is_max_energy_nothing(input_row.energy_pool_temp)
                    reduce_max_energy!(input_row.energy_pool_temp,
                                       already_distributed_input_energies,
                                       already_distributed_input_temperatures,
                                       output_row.target.uac,
                                       input_row.source.uac,
                                       unit.run_id)
                end
                if !is_max_energy_nothing(output_row.energy_potential_temp)
                    reduce_max_energy!(output_row.energy_potential_temp,
                                       already_distributed_output_energies,
                                       already_distributed_output_temperatures,
                                       input_row.source.uac,
                                       output_row.target.uac,
                                       unit.run_id)
                end
                if !is_max_energy_nothing(output_row.energy_pool_temp)
                    reduce_max_energy!(output_row.energy_pool_temp,
                                       already_distributed_output_energies,
                                       already_distributed_output_temperatures,
                                       input_row.source.uac,
                                       output_row.target.uac,
                                       unit.run_id)
                end
            end
        end
    end
end

"""
    reset_balance_table(bus, true)

Resets the balance table of the bus.

This is typically called within inner_distribute! or when resetting a bus.

# Arguments
`unit::Bus`: The bus containing the balance table
`call_inner_distribute::Bool`: If true, calls inner_distribute! after resetting the balance
    table. Please note that this can lead to infinite recursion if used correctly.
"""
function reset_balance_table!(unit::Bus, call_inner_distribute::Bool)
    unit.balance_table = fill(0.0, (length(unit.balance_table_inputs), 2 * length(unit.balance_table_outputs)))
    for i in 1:length(unit.balance_table_inputs)
        for j in 2:2:(2 * length(unit.balance_table_outputs))
            unit.balance_table[i, j] = nothing
        end
    end

    for input_row in collect(values(unit.balance_table_inputs))
        input_row.energy_potential_temp = copy(input_row.energy_potential)
        input_row.energy_pool_temp = copy(input_row.energy_pool)
    end

    for output_row in collect(values(unit.balance_table_outputs))
        output_row.energy_potential_temp = copy(output_row.energy_potential)
        output_row.energy_pool_temp = copy(output_row.energy_pool)
    end

    if call_inner_distribute
        inner_distribute!(unit)
    end
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
    return [f
            for f in unit.input_interfaces
            if (inclusive && f.source.sys_function == condition)
               ||
               (!inclusive && f.source.sys_function != condition)]
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
    return [f
            for f in unit.output_interfaces
            if (inclusive && f.target.sys_function == condition)
               ||
               (!inclusive && f.target.sys_function != condition)]
end

"""
    inputs_recursive(bus)

List of all inputs to the given bus, including all inputs from incoming busses, recursively.

The components are ordered according to the input priorities of each bus, where the list
of any bus in inserted in place of the input priority of that bus in the succeeding bus.

# Arguments
`unit::Bus`: The bus from which to start the recursive ascent

# Returns
`Vector{Component}`: A list of input components that are "reachable" from the starting bus
"""
function inputs_recursive(unit::Bus)::Vector{Component}
    inputs = []
    for inface in unit.input_interfaces
        if inface.source.sys_function == sf_bus
            append!(inputs, inputs_recursive(inface.source))
        else
            push!(inputs, inface.source)
        end
    end
    return inputs
end

"""
    outputs_recursive(bus)

List of all outputs to the given bus, including all outputs from outgoing busses,
recursively.

The components are ordered according to the output priorities of each bus, where the list
of any bus in inserted in place of the output priority of that bus in the preceding bus.

# Arguments
`unit::Bus`: The bus from which to start the recursive descent

# Returns
`Vector{Component}`: A list of output components that are "reachable" from the starting bus
"""
function outputs_recursive(unit::Bus)::Vector{Component}
    outputs = []
    for outface in unit.output_interfaces
        if outface.target.sys_function == sf_bus
            append!(outputs, outputs_recursive(outface.target))
        else
            push!(outputs, outface.target)
        end
    end
    return outputs
end

"""
    bus_transfer_sum(proxy, left_bus, right_bus)

Sum of energy that was transferred from the left bus to the right bus.

# Arguments
`proxy::Bus`: The proxy bus to both busses
`left::Bus`: The bus providing energy
`right::Bus`: The bus receiving energy

# Returns
`Float64`: The sum of energy transferred
"""
function bus_transfer_sum(proxy::Bus, left::Bus, right::Bus)::Float64
    transfer_sum = 0.0

    for input in inputs_recursive(left)
        input_row = proxy.balance_table_inputs[input.uac]
        input_sum = 0.0

        for output in outputs_recursive(right)
            output_row = proxy.balance_table_outputs[output.uac]
            input_sum += proxy.balance_table[input_row.input_index,
                                             output_row.output_index * 2 - 1]
        end

        transfer_sum += input_sum
    end

    return transfer_sum
end

"""
    distribute!(unit)

Bus-specific implementation of distribute!.

This balances the busses in the chain and sets the energy transferred between busses as the
value on the connecting interfaces. The actual balancing is done in the balance table
calculations of balance_on, so this serves as the last call to distribution functions at the
end of a timestep. The function is designed to work both for single busses and busses in a
chain. However it makes it necessary that proxy busses get distributed before any of their
principals.
"""
function distribute!(unit::Bus)
    # the proxy bus has its own distribute step, such that the principal busses only
    # have to handle the interfaces between busses
    if unit.proxy !== nothing
        # set energy between busses by requesting the value from the proxy. we do this for
        # outgoing busses only so the connections are not counted twice
        for outface in filter_outputs(unit, sf_bus, true)
            val = bus_transfer_sum(unit.proxy, unit, outface.target)
            # set the value, then set the interface (back to) zero, so the balance is
            # upheld, but the transfer of energy is registered
            set!(outface, val)
            set!(outface, 0.0)
        end

        return
    end

    # this is not always necessary as calls to balance_on usually call the inner distribute
    # function, however some components might have changed the balance without calling
    # balance_on
    inner_distribute!(unit::Bus)

    # if there is a balance unequal zero remaining, this means the energy balance across
    # the chain of buses was not upheld. we save the balance in the remainder such that
    # subsequent calls to balance consider the missing/extra energy accordingly
    unit.remainder = balance_direct(unit)

    # reset all principal non-bus input interfaces from the original busses
    for inface in filter_inputs(unit, sf_bus, false)
        principal = inface.source.output_interfaces[unit.medium]
        set!(principal, 0.0)
    end

    # reset all principal non-bus output interfaces from the original busses
    for outface in filter_outputs(unit, sf_bus, false)
        principal = outface.target.input_interfaces[unit.medium]
        set!(principal, 0.0)
    end
end

function output_values(unit::Bus)::Vector{String}
    # dynamic output channels
    outputs = ["Transfer->" * outface.target.uac
               for outface in unit.output_interfaces
               if outface.target.sys_function == sf_bus]

    # add to static output channels
    return append!(["Balance"], outputs)
end

function output_value(unit::Bus, key::OutputKey)::Float64
    if key.value_key == "Balance"
        return balance(unit)

    elseif startswith(key.value_key, "Transfer")
        out_uac = last(split(key.value_key, "->"))
        outface = first([f for f in unit.output_interfaces
                         if f.target.uac == out_uac])
        return calculate_energy_flow(outface)
    end

    throw(KeyError(key.value_key))
end

# merge functionality is in an extra file
include("./bus_merging.jl")

export Bus
