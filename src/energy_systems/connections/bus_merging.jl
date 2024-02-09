function merge(first::Bus, second::Bus, uac::String)::Bus
    new_bus = deepcopy(first)
    new_bus.uac = uac

    # insert outputs of second to outputs of first in place of bus->bus interface
    output_idx = [
        idx for (idx,f) in pairs(first.output_interfaces)
            if occursin(f.target.uac, second.uac)
    ][1]
    splice!(new_bus.output_interfaces, output_idx, second.output_interfaces)
    splice!(
        new_bus.connectivity.output_order,
        output_idx,
        second.connectivity.output_order
    )

    # insert inputs of first to inputs of second in place of bus->bus interface
    input_idx = [
        idx for (idx,f) in pairs(second.input_interfaces)
            if occursin(f.source.uac, first.uac)
    ][1]
    new_bus.input_interfaces = Base.deepcopy(second.input_interfaces)
    splice!(new_bus.input_interfaces, input_idx, first.input_interfaces)
    new_bus.connectivity.input_order = Base.deepcopy(second.connectivity.input_order)
    splice!(
        new_bus.connectivity.input_order,
        input_idx,
        first.connectivity.input_order
    )

    # append balance table inputs and rewrite priorities
    merge!(new_bus.balance_table_inputs, Base.deepcopy(second.balance_table_inputs))
    delete!(new_bus.balance_table_inputs, first.uac)
    for (idx,inface) in pairs(new_bus.input_interfaces)
        new_bus.balance_table_inputs[inface.source.uac].priority = idx
        new_bus.balance_table_inputs[inface.source.uac].input_index = idx
    end

    # append balance table outputs and rewrite priorities
    merge!(new_bus.balance_table_outputs, Base.deepcopy(second.balance_table_outputs))
    delete!(new_bus.balance_table_outputs, second.uac)
    for (idx,outface) in pairs(new_bus.output_interfaces)
        new_bus.balance_table_outputs[outface.target.uac].priority = idx
        new_bus.balance_table_outputs[outface.target.uac].output_index = idx
    end

    # merge energy_flow
    new_bus.connectivity.energy_flow = []

    for input_row in sort(collect(values(new_bus.balance_table_inputs)), by=x->x.priority)
        is_allowed_flags = []
        is_input_from_first = input_row.source.uac in keys(first.balance_table_inputs)

        for output_row in sort(collect(values(new_bus.balance_table_outputs)), by=x->x.priority)
            is_output_from_first = output_row.target.uac in keys(first.balance_table_outputs)

            if is_input_from_first && is_output_from_first
                input_index = first.balance_table_inputs[input_row.source.uac].input_index
                output_index = first.balance_table_outputs[output_row.target.uac].output_index
                is_allowed = first.connectivity.energy_flow[input_index][output_index]

            elseif is_input_from_first && !is_output_from_first
                input_index = first.balance_table_inputs[input_row.source.uac].input_index
                output_index = first.balance_table_outputs[second.uac].output_index
                is_allowed = first.connectivity.energy_flow[input_index][output_index]
                input_index = second.balance_table_inputs[first.uac].input_index
                output_index = second.balance_table_outputs[output_row.target.uac].output_index
                is_allowed = is_allowed && second.connectivity.energy_flow[input_index][output_index]

            elseif !is_input_from_first && is_output_from_first
                # inputs from the second bus can't see outputs of the first
                is_allowed = false

            else
                input_index = second.balance_table_inputs[input_row.source.uac].input_index
                output_index = second.balance_table_outputs[output_row.target.uac].output_index
                is_allowed = second.connectivity.energy_flow[input_index][output_index]
            end

            push!(is_allowed_flags, is_allowed)
        end

        push!(new_bus.connectivity.energy_flow, is_allowed_flags)
    end

    # reset recreates the balance table
    reset(new_bus)

    return new_bus
end

Base.@kwdef mutable struct BusNode
    uac::String
    inputs::Vector{BusNode}
    outputs::Vector{BusNode}
end

function is_linear_connection(parent::BusNode, child::BusNode)::Bool
    nr_children = length(parent.outputs)
    nr_parents = length(child.inputs)
    return nr_children == 1 && nr_parents == 1
end

function find_next_linear_connection(
    chain::Dict{String,BusNode}
)::Tuple{Union{Nothing,BusNode},Union{Nothing,BusNode}}
    for bus in values(chain)
        for child in bus.outputs
            if is_linear_connection(bus, child)
                return (bus, child)
            end
        end
        for parent in bus.inputs
            if is_linear_connection(parent, bus)
                return (parent, bus)
            end
        end
    end
    return (nothing, nothing)
end

function nodes_from_components(components::Grouping)::Dict{String,BusNode}
    nodes = Dict{String,BusNode}()
    for uac in keys(components)
        nodes[uac] = BusNode(
            uac=uac,
            inputs=Vector{BusNode}(),
            outputs=Vector{BusNode}()
        )
    end
    for (uac, bus) in pairs(components)
        node = nodes[uac]
        node.inputs = [nodes[i.source.uac] for i in filter_inputs(bus, sf_bus, true)]
        node.outputs = [nodes[o.target.uac] for o in filter_outputs(bus, sf_bus, true)]
    end
    return nodes
end

function merge(parent::BusNode, child::BusNode, new_uac::String)::BusNode
    new_node = BusNode(
        uac=new_uac,
        inputs=[],
        outputs=[]
    )
    append!(new_node.inputs, parent.inputs)
    append!(new_node.inputs, [n for n in child.inputs if n.uac != parent.uac])
    append!(new_node.outputs, child.outputs)
    append!(new_node.outputs, [n for n in parent.outputs if n.uac != child.uac])
    return new_node
end

function update_nodes!(
    nodes::Dict{String,BusNode},
    parent_uac::String,
    child_uac::String,
    new_node::BusNode
)
    for node in values(nodes)
        for (idx,output) in pairs(node.outputs)
            if output.uac == parent_uac || output.uac == child_uac
                splice!(node.outputs, idx, [new_node])
            end
        end
        for (idx,input) in pairs(node.inputs)
            if input.uac == parent_uac || input.uac == child_uac
                splice!(node.inputs, idx, [new_node])
            end
        end
    end
end

function priority_of_child(parent::BusNode, child::BusNode)::Integer
    indices = indexin([child], parent.outputs)
    return indices[1] === nothing ? 0 : indices[1]
end

function priority_of_parent(parent::BusNode, child::BusNode)::Integer
    indices = indexin([parent], child.inputs)
    return indices[1] === nothing ? 0 : indices[1]
end

function find_single_parent_leaves(
    nodes::Dict{String,BusNode}
)::Vector{Tuple{Integer,BusNode,BusNode}}
    leaves = []
    for node in values(nodes)
        if length(node.inputs) != 1 || length(node.outputs) > 0
            continue
        end
        parent = node.inputs[1]
        push!(leaves, (priority_of_child(parent, node), parent, node))
    end
    return leaves
end

function find_parents_of_leaves(
    nodes::Dict{String,BusNode}
)::Vector{Tuple{Integer,BusNode,BusNode}}
    parents = []
    for node in values(nodes)
        if length(node.outputs) > 0
            continue
        end
        for parent in node.inputs
            push!(parents, (priority_of_parent(parent, node), parent, node))
        end
    end
    return parents
end

function merge_busses(busses_to_merge::Grouping)::Union{Nothing,Bus}
    nodes = nodes_from_components(busses_to_merge)
    components = Grouping()
    for (key,val) in pairs(busses_to_merge)
        components[key] = val
    end

    # first remove all linear connections
    parent, child = find_next_linear_connection(nodes)
    while parent !== nothing
        new_uac = parent.uac * child.uac
        new_bus = merge(components[parent.uac], components[child.uac], new_uac)
        components[new_uac] = new_bus

        new_node = merge(parent, child, new_uac)
        nodes[new_uac] = new_node
        update_nodes!(nodes, parent.uac, child.uac, new_node)
        delete!(nodes, parent.uac)
        delete!(nodes, child.uac)

        parent, child = find_next_linear_connection(nodes)
    end

    # now keep merging until we have only one node left
    while length(nodes) > 1
        # prioritise leaves with only one parent and order by descending priority
        leaves = find_single_parent_leaves(nodes)
        if length(leaves) > 0
            next = sort(leaves, by=x->x[1], rev=true)[1]
            parent = next[2]
            child = next[3]

            new_uac = parent.uac * child.uac
            new_bus = merge(components[parent.uac], components[child.uac], new_uac)
            components[new_uac] = new_bus

            new_node = merge(parent, child, new_uac)
            nodes[new_uac] = new_node
            update_nodes!(nodes, parent.uac, child.uac, new_node)
            delete!(nodes, parent.uac)
            delete!(nodes, child.uac)

            continue
        end

        # for leaves with multiple parents, merge into those with the last priority
        parents = find_parents_of_leaves(nodes)
        if length(parents) > 0
            next = sort(parents, by=x->x[1], rev=true)[1]
            parent = next[2]
            child = next[3]

            new_uac = parent.uac * child.uac
            new_bus = merge(components[parent.uac], components[child.uac], new_uac)
            components[new_uac] = new_bus

            new_node = merge(parent, child, new_uac)
            nodes[new_uac] = new_node
            update_nodes!(nodes, parent.uac, child.uac, new_node)
            delete!(nodes, parent.uac)
            delete!(nodes, child.uac)

            continue
        end

        # if we got here, something went wrong and would result in an endless loop
        @error("Could not resolve merging of busses")
        break
    end

    return components[first(values(nodes)).uac]
end
