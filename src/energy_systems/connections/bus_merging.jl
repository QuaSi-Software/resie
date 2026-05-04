Base.@kwdef mutable struct ComponentNode
    uac::String
    sys_function::SystemFunction
    do_storage_transfer::Bool
    is_secondary_interface::Bool = false
end

Base.@kwdef mutable struct BusNode
    uac::String
    sys_function::SystemFunction
    inputs::Vector{Union{ComponentNode,BusNode}}
    outputs::Vector{Union{ComponentNode,BusNode}}
    energy_flow::Array{Int,2}
end

Base.@kwdef struct PriorityConstraint
    before::String
    after::String
    kind::String          # "input" or "output"
    before_bus::String
    before_order::Int
    after_bus::String
    after_order::Int
end

function bus_inputs(node::BusNode)
    return [input for input in node.inputs if input.sys_function == sf_bus]
end

function bus_outputs(node::BusNode)
    return [output for output in node.outputs if output.sys_function == sf_bus]
end

function input_key(node::ComponentNode)::String
    return adjust_name_if_secondary(node.uac, node.is_secondary_interface)
end

function output_key(node::ComponentNode)::String
    return node.uac
end

function is_internal_bus(node, internal_uacs::Set{String})::Bool
    return node.sys_function == sf_bus && node.uac in internal_uacs
end

function proxy_uac_for(nodes::Dict{String,BusNode})::String
    return "Proxy-" * join(sort(collect(keys(nodes))), "|")
end

function lexless(a::Vector{Int}, b::Vector{Int})::Bool
    n = min(length(a), length(b))

    for idx in 1:n
        if a[idx] < b[idx]
            return true
        elseif a[idx] > b[idx]
            return false
        end
    end

    return length(a) < length(b)
end

function nodes_from_components(bus_components::Grouping)::Dict{String,BusNode}
    nodes = Dict{String,BusNode}()

    for uac in keys(bus_components)
        nodes[uac] = BusNode(; uac=uac,
                             sys_function=sf_bus,
                             inputs=Vector{Union{ComponentNode,BusNode}}(),
                             outputs=Vector{Union{ComponentNode,BusNode}}(),
                             energy_flow=Array{Int,2}(undef, 0, 0))
    end

    for (uac, bus) in pairs(bus_components)
        node = nodes[uac]

        for inface in values(bus.input_interfaces)
            inface === nothing && continue

            input = inface.source

            if input.sys_function == sf_bus && haskey(nodes, input.uac)
                push!(node.inputs, nodes[input.uac])
            else
                push!(node.inputs,
                      ComponentNode(; uac=input.uac,
                                    sys_function=input.sys_function,
                                    do_storage_transfer=inface.do_storage_transfer,
                                    is_secondary_interface=inface.is_secondary_interface))
            end
        end

        for outface in values(bus.output_interfaces)
            outface === nothing && continue

            output = outface.target

            if output.sys_function == sf_bus && haskey(nodes, output.uac)
                push!(node.outputs, nodes[output.uac])
            else
                push!(node.outputs,
                      ComponentNode(; uac=output.uac,
                                    sys_function=output.sys_function,
                                    do_storage_transfer=outface.do_storage_transfer,
                                    is_secondary_interface=outface.is_secondary_interface))
            end
        end

        node.energy_flow = fill(1, (length(node.inputs), length(node.outputs)))

        if bus.connectivity.energy_flow !== nothing
            for (row_idx, row) in pairs(bus.connectivity.energy_flow)
                for (col_idx, col_val) in pairs(row)
                    node.energy_flow[row_idx, col_idx] = Int(col_val)
                end
            end
        end
    end

    return nodes
end

function input_index(node::BusNode, uac::String)::Int
    indices = indexin([uac],
                      [adjust_name_if_secondary(n.uac,
                                                hasproperty(n, :is_secondary_interface) ?
                                                n.is_secondary_interface :
                                                false)
                       for n in node.inputs])

    return length(indices) != 1 || indices[1] === nothing ? 0 : indices[1]
end

function output_index(node::BusNode, uac::String)::Int
    indices = indexin([uac], [n.uac for n in node.outputs])
    return length(indices) != 1 || indices[1] === nothing ? 0 : indices[1]
end

function internal_bus_inputs(node::BusNode, nodes::Dict{String,BusNode})
    return [input for input in node.inputs if input isa BusNode && haskey(nodes, input.uac)]
end

function internal_bus_outputs(node::BusNode, nodes::Dict{String,BusNode})
    return [output for output in node.outputs if output isa BusNode && haskey(nodes, output.uac)]
end

function bus_discovery_rank(nodes::Dict{String,BusNode})::Dict{String,Int}
    roots = BusNode[]

    for bus in values(nodes)
        if isempty(internal_bus_inputs(bus, nodes))
            push!(roots, bus)
        end
    end

    if isempty(roots)
        roots = collect(values(nodes))
    end

    sort!(roots; by=bus -> bus.uac)

    ordered = String[]
    seen = Set{String}()
    queue = copy(roots)

    while !isempty(queue)
        bus = popfirst!(queue)

        if bus.uac in seen
            continue
        end

        push!(seen, bus.uac)
        push!(ordered, bus.uac)

        for child in internal_bus_outputs(bus, nodes)
            if !(child.uac in seen)
                push!(queue, child)
            end
        end
    end

    for uac in sort(collect(keys(nodes)))
        if !(uac in seen)
            push!(ordered, uac)
        end
    end

    return Dict(uac => idx for (idx, uac) in pairs(ordered))
end

function format_constraint_line(c::PriorityConstraint)::String
    if c.before_bus == c.after_bus
        return "$(c.kind): bus=$(c.before_bus), orders $(c.before_order) < $(c.after_order): $(c.before) before $(c.after)"
    else
        return "$(c.kind): before at bus=$(c.before_bus) order=$(c.before_order), after at bus=$(c.after_bus) order=$(c.after_order): $(c.before) before $(c.after)"
    end
end

function log_reversed_constraints(constraints::Vector{PriorityConstraint}, label::String)
    by_pair = Dict{Tuple{String,String},Vector{PriorityConstraint}}()

    for c in constraints
        push!(get!(by_pair, (c.before, c.after), PriorityConstraint[]), c)
    end

    already_logged = Set{Tuple{String,String}}()

    for ((a, b), forward_constraints) in by_pair
        reverse_pair = (b, a)

        if !haskey(by_pair, reverse_pair)
            continue
        end

        canonical_pair = a < b ? (a, b) : (b, a)

        if canonical_pair in already_logged
            continue
        end

        push!(already_logged, canonical_pair)

        reverse_constraints = by_pair[reverse_pair]

        for c1 in forward_constraints
            for c2 in reverse_constraints
                @warn "Direct contradictory priority constraint while sorting $label\n$(format_constraint_line(c1))\n$(format_constraint_line(c2))"
            end
        end
    end
end

function canonical_cycle_key(cycle::Vector{String})::String
    # cycle is expected to contain the start node again as the final element
    entries = cycle[1:(end - 1)]

    if isempty(entries)
        return ""
    end

    rotations = String[]

    for idx in eachindex(entries)
        rotated = vcat(entries[idx:end], entries[1:(idx - 1)])
        push!(rotations, join(rotated, " -> "))
    end

    return minimum(rotations)
end

function log_priority_cycles(priority_constraints::Vector{PriorityConstraint},
                             unresolved::Vector{String},
                             label::String;
                             max_cycles::Int=20)
    unresolved_set = Set(unresolved)

    if isempty(unresolved_set)
        return
    end

    adjacency = Dict{String,Set{String}}()
    by_pair = Dict{Tuple{String,String},Vector{PriorityConstraint}}()

    for c in priority_constraints
        if !(c.before in unresolved_set) || !(c.after in unresolved_set)
            continue
        end

        push!(get!(adjacency, c.before, Set{String}()), c.after)
        push!(get!(by_pair, (c.before, c.after), PriorityConstraint[]), c)
    end

    cycles = Vector{Vector{String}}()
    seen_cycle_keys = Set{String}()

    function dfs(start::String, current::String, path::Vector{String})
        if length(cycles) >= max_cycles
            return
        end

        if length(path) > length(unresolved_set)
            return
        end

        for next in sort(collect(get(adjacency, current, Set{String}())))
            if next == start && length(path) >= 2
                cycle = vcat(path, [start])
                key = canonical_cycle_key(cycle)

                if !(key in seen_cycle_keys)
                    push!(seen_cycle_keys, key)
                    push!(cycles, cycle)
                end

                if length(cycles) >= max_cycles
                    return
                end
            elseif next in unresolved_set && !(next in path)
                dfs(start, next, vcat(path, [next]))

                if length(cycles) >= max_cycles
                    return
                end
            end
        end
    end

    for start in sort(collect(unresolved_set))
        dfs(start, start, [start])

        if length(cycles) >= max_cycles
            break
        end
    end

    if isempty(cycles)
        @warn "No explicit priority cycle could be reconstructed while sorting $label, although unresolved entries exist." unresolved_entries = join(unresolved,
                                                                                                                                                     "\n")
        return
    end

    for (cycle_idx, cycle) in pairs(cycles)
        lines = String[]

        push!(lines, "Priority cycle $cycle_idx while sorting $label:")
        push!(lines, "cycle: " * join(cycle, " -> "))

        for idx in 1:(length(cycle) - 1)
            before = cycle[idx]
            after = cycle[idx + 1]

            constraints = get(by_pair, (before, after), PriorityConstraint[])

            if isempty(constraints)
                push!(lines, "$(idx). $(before) before $(after) [no provenance found]")
            else
                c = constraints[1]
                push!(lines, "$(idx). " * format_constraint_line(c))
            end
        end

        @warn join(lines, "\n")
    end
end

function topo_sort_keys(keys_in::Vector{String},
                        priority_constraints_in::Vector{PriorityConstraint},
                        fallback_rank::Dict{String,Tuple{Int,Int,Int,String}},
                        label::String)::Vector{String}
    keys_unique = unique(keys_in)
    key_set = Set(keys_unique)

    priority_constraints = PriorityConstraint[]
    pair_constraints = Tuple{String,String}[]

    for c in priority_constraints_in
        if c.before == c.after
            continue
        end

        if c.before in key_set && c.after in key_set
            push!(priority_constraints, c)
            push!(pair_constraints, (c.before, c.after))
        end
    end

    priority_constraints = unique(priority_constraints)
    pair_constraints = unique(pair_constraints)

    incoming = Dict(key => Set{String}() for key in keys_unique)
    outgoing = Dict(key => String[] for key in keys_unique)

    for (before, after) in pair_constraints
        push!(incoming[after], before)
        push!(outgoing[before], after)
    end

    fallback(key) = get(fallback_rank,
                        key,
                        (typemax(Int), typemax(Int), typemax(Int), key))

    ready = sort([key for key in keys_unique if isempty(incoming[key])];
                 by=fallback)

    result = String[]

    while !isempty(ready)
        key = popfirst!(ready)
        push!(result, key)

        for other in outgoing[key]
            delete!(incoming[other], key)

            if isempty(incoming[other]) && !(other in result) && !(other in ready)
                push!(ready, other)
                sort!(ready; by=fallback)
            end
        end
    end

    if length(result) != length(keys_unique)
        unresolved = [key for key in keys_unique if !(key in result)]
        @warn "Contradictory priority constraints while sorting $label. Falling back for unresolved entries. \n" unresolved_entries = join(unresolved,
                                                                                                                                           " \n ")

        # log_priority_cycles(priority_constraints,
        #                     unresolved,
        #                     label;
        #                     max_cycles=20)
        # TODO

        append!(result, sort(unresolved; by=fallback))

        log_reversed_constraints(priority_constraints, label)
    end

    return result
end

function reachable_bus_input_states_from_input(nodes::Dict{String,BusNode},
                                               start_bus_uac::String,
                                               start_input_idx::Int)::Set{Tuple{String,Int}}
    internal_uacs = Set(collect(keys(nodes)))

    reachable = Set{Tuple{String,Int}}()
    queue = Tuple{String,Int}[(start_bus_uac, start_input_idx)]

    while !isempty(queue)
        bus_uac, input_idx = popfirst!(queue)
        state = (bus_uac, input_idx)

        if state in reachable
            continue
        end

        push!(reachable, state)

        bus = nodes[bus_uac]

        if input_idx <= 0 || input_idx > size(bus.energy_flow, 1)
            continue
        end

        for (output_idx, output) in pairs(bus.outputs)
            if output_idx > size(bus.energy_flow, 2)
                continue
            end

            if bus.energy_flow[input_idx, output_idx] == 0
                continue
            end

            if is_internal_bus(output, internal_uacs)
                child = nodes[output.uac]
                child_input_idx = input_index(child, bus.uac)

                if child_input_idx == 0
                    @warn "Internal bus input priority could not be followed" parent = bus.uac child = child.uac
                    continue
                end

                push!(queue, (child.uac, child_input_idx))
            end
        end
    end

    return reachable
end

function input_priority_constraints(nodes::Dict{String,BusNode},
                                    input_ports::Vector{Tuple{String,Int,String}})::Vector{PriorityConstraint}
    bus_to_input_priority = Dict{String,Dict{String,Int}}()

    for (bus_uac, input_idx, key) in input_ports
        states = reachable_bus_input_states_from_input(nodes, bus_uac, input_idx)

        for (state_bus_uac, state_input_idx) in states
            entry = get!(bus_to_input_priority,
                         state_bus_uac,
                         Dict{String,Int}())

            entry[key] = min(get(entry, key, typemax(Int)), state_input_idx)
        end
    end

    constraints = PriorityConstraint[]

    for (bus_uac, priority_by_key) in bus_to_input_priority
        keys_here = collect(keys(priority_by_key))

        for before in keys_here
            for after in keys_here
                if before == after
                    continue
                end

                before_order = priority_by_key[before]
                after_order = priority_by_key[after]

                if before_order < after_order
                    push!(constraints,
                          PriorityConstraint(; before=before,
                                             after=after,
                                             kind="input",
                                             before_bus=bus_uac,
                                             before_order=before_order,
                                             after_bus=bus_uac,
                                             after_order=after_order))
                end
            end
        end
    end

    return unique(constraints)
end

function collect_proxy_inputs(nodes::Dict{String,BusNode})
    internal_uacs = Set(collect(keys(nodes)))
    bus_rank = bus_discovery_rank(nodes)

    by_key = Dict{String,ComponentNode}()
    input_ports = Tuple{String,Int,String}[]
    fallback_rank = Dict{String,Tuple{Int,Int,Int,String}}()

    occurrence = 0

    for bus_uac in sort(collect(keys(nodes)); by=uac -> bus_rank[uac])
        bus = nodes[bus_uac]

        for (input_idx, input) in pairs(bus.inputs)
            if is_internal_bus(input, internal_uacs)
                continue
            end

            if !(input isa ComponentNode)
                error("Unexpected non-component external input $(input.uac)")
            end

            occurrence += 1
            key = input_key(input)

            push!(input_ports, (bus_uac, input_idx, key))

            if !haskey(by_key, key)
                by_key[key] = input
                fallback_rank[key] = (bus_rank[bus_uac], input_idx, occurrence, key)
            else
                if by_key[key].do_storage_transfer != input.do_storage_transfer ||
                   by_key[key].is_secondary_interface != input.is_secondary_interface
                    @warn "Duplicate proxy input with differing interface metadata. Keeping first metadata." key
                end

                fallback_rank[key] = min(fallback_rank[key],
                                         (bus_rank[bus_uac], input_idx, occurrence, key))
            end
        end
    end

    constraints = input_priority_constraints(nodes, input_ports)

    ordered_keys = topo_sort_keys(collect(keys(by_key)),
                                  constraints,
                                  fallback_rank,
                                  "proxy inputs")

    proxy_inputs = ComponentNode[by_key[key] for key in ordered_keys]

    return proxy_inputs, input_ports
end

function path_order_indices(path::Vector{Tuple{String,Int}})::Vector{Int}
    return [idx for (_, idx) in path]
end

function path_lexless(a::Vector{Tuple{String,Int}},
                      b::Vector{Tuple{String,Int}})::Bool
    return lexless(path_order_indices(a), path_order_indices(b))
end

function prefer_path_priority!(priorities::Dict{String,Vector{Tuple{String,Int}}},
                               key::String,
                               priority::Vector{Tuple{String,Int}})
    if !haskey(priorities, key) || path_lexless(priority, priorities[key])
        priorities[key] = priority
    end
end

function reachable_outputs_with_priority(nodes::Dict{String,BusNode},
                                         start_bus_uac::String,
                                         start_input_idx::Int)::Dict{String,Vector{Tuple{String,Int}}}
    internal_uacs = Set(collect(keys(nodes)))

    reachable = Dict{String,Vector{Tuple{String,Int}}}()

    queue = Tuple{String,Int,Vector{Tuple{String,Int}}}[(start_bus_uac, start_input_idx, Tuple{String,Int}[])]

    best_state_priority = Dict{Tuple{String,Int},Vector{Tuple{String,Int}}}()

    while !isempty(queue)
        bus_uac, input_idx, path_priority = popfirst!(queue)
        state = (bus_uac, input_idx)

        previous = get(best_state_priority, state, nothing)

        if previous !== nothing && !path_lexless(path_priority, previous)
            continue
        end

        best_state_priority[state] = path_priority

        bus = nodes[bus_uac]

        if input_idx <= 0 || input_idx > size(bus.energy_flow, 1)
            continue
        end

        for (output_idx, output) in pairs(bus.outputs)
            if output_idx > size(bus.energy_flow, 2)
                continue
            end

            if bus.energy_flow[input_idx, output_idx] == 0
                continue
            end

            new_priority = copy(path_priority)
            push!(new_priority, (bus_uac, output_idx))

            if is_internal_bus(output, internal_uacs)
                child = nodes[output.uac]
                child_input_idx = input_index(child, bus.uac)

                if child_input_idx == 0
                    @warn "Internal bus output priority could not be followed" parent = bus.uac child = child.uac
                    continue
                end

                push!(queue, (child.uac, child_input_idx, new_priority))
            else
                if !(output isa ComponentNode)
                    error("Unexpected reachable non-component output $(output.uac)")
                end

                prefer_path_priority!(reachable, output_key(output), new_priority)
            end
        end
    end

    return reachable
end

function output_priority_constraint_from_paths(before::String,
                                               after::String,
                                               before_path::Vector{Tuple{String,Int}},
                                               after_path::Vector{Tuple{String,Int}})::PriorityConstraint
    n = min(length(before_path), length(after_path))

    for idx in 1:n
        before_bus, before_order = before_path[idx]
        after_bus, after_order = after_path[idx]

        if before_bus != after_bus || before_order != after_order
            return PriorityConstraint(; before=before,
                                      after=after,
                                      kind="output",
                                      before_bus=before_bus,
                                      before_order=before_order,
                                      after_bus=after_bus,
                                      after_order=after_order)
        end
    end

    if length(before_path) < length(after_path)
        before_bus, before_order = isempty(before_path) ? ("<start>", 0) : before_path[end]
        after_bus, after_order = after_path[length(before_path) + 1]
    else
        before_bus, before_order = before_path[length(after_path) + 1]
        after_bus, after_order = isempty(after_path) ? ("<start>", 0) : after_path[end]
    end

    return PriorityConstraint(; before=before,
                              after=after,
                              kind="output",
                              before_bus=before_bus,
                              before_order=before_order,
                              after_bus=after_bus,
                              after_order=after_order)
end

function output_priority_constraints(nodes::Dict{String,BusNode},
                                     input_ports::Vector{Tuple{String,Int,String}})::Vector{PriorityConstraint}
    constraints = PriorityConstraint[]

    for (bus_uac, input_idx, _) in input_ports
        reachable = reachable_outputs_with_priority(nodes, bus_uac, input_idx)
        keys_here = collect(keys(reachable))

        for before in keys_here
            for after in keys_here
                if before == after
                    continue
                end

                before_path = reachable[before]
                after_path = reachable[after]

                if path_lexless(before_path, after_path)
                    push!(constraints,
                          output_priority_constraint_from_paths(before,
                                                                after,
                                                                before_path,
                                                                after_path))
                end
            end
        end
    end

    return unique(constraints)
end

function collect_proxy_outputs(nodes::Dict{String,BusNode},
                               input_ports::Vector{Tuple{String,Int,String}})
    internal_uacs = Set(collect(keys(nodes)))
    bus_rank = bus_discovery_rank(nodes)

    by_key = Dict{String,ComponentNode}()
    fallback_rank = Dict{String,Tuple{Int,Int,Int,String}}()

    occurrence = 0

    for bus_uac in sort(collect(keys(nodes)); by=uac -> bus_rank[uac])
        bus = nodes[bus_uac]

        for (output_idx, output) in pairs(bus.outputs)
            if is_internal_bus(output, internal_uacs)
                continue
            end

            if !(output isa ComponentNode)
                error("Unexpected non-component external output $(output.uac)")
            end

            occurrence += 1
            key = output_key(output)

            if !haskey(by_key, key)
                by_key[key] = output
                fallback_rank[key] = (bus_rank[bus_uac], output_idx, occurrence, key)
            else
                if by_key[key].do_storage_transfer != output.do_storage_transfer
                    @warn "Duplicate proxy output with differing interface metadata. Keeping first metadata." key
                end

                fallback_rank[key] = min(fallback_rank[key],
                                         (bus_rank[bus_uac], output_idx, occurrence, key))
            end
        end
    end

    constraints = output_priority_constraints(nodes, input_ports)

    ordered_keys = topo_sort_keys(collect(keys(by_key)),
                                  constraints,
                                  fallback_rank,
                                  "proxy outputs")

    proxy_outputs = ComponentNode[by_key[key] for key in ordered_keys]

    return proxy_outputs
end

function bus_from_node(node::BusNode, template::Bus, components::Grouping)::Bus
    bus = Bus(node.uac, template.medium, template.run_id, template.epsilon)

    for (idx, input) in pairs(node.inputs)
        if !(input isa ComponentNode)
            error("Cannot build proxy bus $(node.uac): unresolved BusNode input $(input.uac)")
        end

        if !haskey(components, input.uac)
            error("Cannot build proxy bus $(node.uac): input component $(input.uac) not found in components")
        end

        push!(bus.input_interfaces,
              SystemInterface(; source=components[input.uac],
                              target=bus,
                              do_storage_transfer=input.do_storage_transfer,
                              is_secondary_interface=input.is_secondary_interface))

        input_uac = input_key(input)

        push!(bus.connectivity.input_order, input_uac)

        bus.balance_table_inputs[input_uac] = BTInputRow(; source=components[input.uac],
                                                         priority=idx,
                                                         do_storage_transfer=input.do_storage_transfer,
                                                         is_secondary_interface=input.is_secondary_interface)
    end

    for (idx, output) in pairs(node.outputs)
        if !(output isa ComponentNode)
            error("Cannot build proxy bus $(node.uac): unresolved BusNode output $(output.uac)")
        end

        if !haskey(components, output.uac)
            error("Cannot build proxy bus $(node.uac): output component $(output.uac) not found in components")
        end

        push!(bus.output_interfaces,
              SystemInterface(; source=bus,
                              target=components[output.uac],
                              do_storage_transfer=output.do_storage_transfer))

        output_uac = output_key(output)

        push!(bus.connectivity.output_order, output_uac)

        bus.balance_table_outputs[output_uac] = BTOutputRow(; target=components[output.uac],
                                                            priority=idx,
                                                            do_storage_transfer=output.do_storage_transfer)
    end

    bus.connectivity.energy_flow = []

    for row_idx in keys(bus.connectivity.input_order)
        row = []

        for col_idx in keys(bus.connectivity.output_order)
            push!(row, node.energy_flow[row_idx, col_idx])
        end

        push!(bus.connectivity.energy_flow, row)
    end

    reset(bus)

    bus.input_output_rows_iteration, bus.has_custom_order = iterate_balance_table(bus)

    return bus
end

function merge_busses(busses_to_merge::Grouping, components::Grouping)::Union{Nothing,Bus}
    if isempty(busses_to_merge)
        return nothing
    end

    nodes = nodes_from_components(busses_to_merge)

    proxy_inputs, input_ports = collect_proxy_inputs(nodes)
    proxy_outputs = collect_proxy_outputs(nodes, input_ports)

    if isempty(proxy_inputs)
        error("Cannot merge bus chain into proxy bus: no external inputs found.")
    end

    if isempty(proxy_outputs)
        error("Cannot merge bus chain into proxy bus: no external outputs found.")
    end

    new_node = BusNode(; uac=proxy_uac_for(nodes),
                       sys_function=sf_bus,
                       inputs=Union{ComponentNode,BusNode}[proxy_inputs...],
                       outputs=Union{ComponentNode,BusNode}[proxy_outputs...],
                       energy_flow=zeros(Int, length(proxy_inputs), length(proxy_outputs)))

    input_row = Dict{String,Int}()

    for (idx, input) in pairs(proxy_inputs)
        input_row[input_key(input)] = idx
    end

    output_col = Dict{String,Int}()

    for (idx, output) in pairs(proxy_outputs)
        output_col[output_key(output)] = idx
    end

    for (bus_uac, input_idx, key) in input_ports
        row_idx = input_row[key]
        reachable = reachable_outputs_with_priority(nodes, bus_uac, input_idx)

        for out_key in keys(reachable)
            col_idx = get(output_col, out_key, 0)

            if col_idx > 0
                new_node.energy_flow[row_idx, col_idx] = 1
            end
        end
    end

    return bus_from_node(new_node,
                         first(values(busses_to_merge)),
                         components)
end
