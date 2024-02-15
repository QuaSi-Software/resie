Base.@kwdef mutable struct ComponentNode
    uac::String
    sys_function::SystemFunction
    do_storage_transfer::Bool
end

Base.@kwdef mutable struct BusNode
    uac::String
    sys_function::SystemFunction
    inputs::Vector{Union{ComponentNode, BusNode}}
    outputs::Vector{Union{ComponentNode, BusNode}}
    energy_flow::Array{Int,2}
end

function bus_inputs(node::BusNode)
    return [input for input in node.inputs if input.sys_function == sf_bus]
end

function bus_outputs(node::BusNode)
    return [output for output in node.outputs if output.sys_function == sf_bus]
end

function is_linear_connection(parent::BusNode, child::BusNode)::Bool
    nr_children = length(bus_outputs(parent))
    nr_parents = length(bus_inputs(child))
    return nr_children == 1 && nr_parents == 1
end

function find_linear_connections(
    nodes::Dict{String,BusNode}
)::Vector{Tuple{Integer,BusNode,BusNode}}
    connections = []

    for bus in values(nodes)
        for child in bus_outputs(bus)
            if is_linear_connection(bus, child)
                push!(connections, (1, bus, child))
            end
        end

        for parent in bus_inputs(bus)
            if is_linear_connection(parent, bus)
                push!(connections, (1, parent, bus))
            end
        end
    end

    return connections
end

function nodes_from_components(bus_components::Grouping)::Dict{String,BusNode}
    nodes = Dict{String,BusNode}()

    for uac in keys(bus_components)
        nodes[uac] = BusNode(
            uac=uac,
            sys_function=sf_bus,
            inputs=Vector{Union{ComponentNode, BusNode}}(),
            outputs=Vector{Union{ComponentNode, BusNode}}(),
            energy_flow=Array{Int,2}(undef, 0, 0)
        )
    end

    for (uac, bus) in pairs(bus_components)
        node = nodes[uac]

        for inface in bus.input_interfaces
            input = inface.source
            if input.sys_function == sf_bus
                push!(node.inputs, nodes[input.uac])
            else
                push!(node.inputs, ComponentNode(
                    uac=input.uac,
                    sys_function=input.sys_function,
                    do_storage_transfer=inface.do_storage_transfer
                ))
            end
        end

        for outface in bus.output_interfaces
            output = outface.target
            if output.sys_function == sf_bus
                push!(node.outputs, nodes[output.uac])
            else
                push!(node.outputs, ComponentNode(
                    uac=output.uac,
                    sys_function=output.sys_function,
                    do_storage_transfer=outface.do_storage_transfer
                ))
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
    indices = indexin([uac], [n.uac for n in node.inputs])
    return length(indices) != 1 || indices[1] === nothing ? 0 : indices[1]
end

function output_index(node::BusNode, uac::String)::Int
    indices = indexin([uac], [n.uac for n in node.outputs])
    return length(indices) != 1 || indices[1] === nothing ? 0 : indices[1]
end

function merge(parent::BusNode, child::BusNode, new_uac::String)::BusNode
    new_node = BusNode(
        uac=new_uac,
        sys_function=sf_bus,
        inputs=Vector{Union{ComponentNode, BusNode}}(),
        outputs=Vector{Union{ComponentNode, BusNode}}(),
        energy_flow=Array{Int,2}(undef, 0, 0)
    )

    # insert outputs of child in outputs of parent in place of the bus->bus connection
    output_idx = output_index(parent, child.uac)
    new_node.outputs = [n for n in parent.outputs]
    splice!(new_node.outputs, output_idx, child.outputs)

    # insert inputs of parent to inputs of child in place of bus->bus connection
    input_idx = input_index(child, parent.uac)
    new_node.inputs = [n for n in child.inputs]
    splice!(new_node.inputs, input_idx, parent.inputs)

    # merge energy_flow matrices
    new_node.energy_flow = fill(1, (length(new_node.inputs), length(new_node.outputs)))

    for (new_in_idx,input) in pairs(new_node.inputs)
        is_input_from_first = input_index(parent, input.uac) > 0

        for (new_out_idx,output) in pairs(new_node.outputs)
            is_output_from_first = output_index(parent, output.uac) > 0

            if is_input_from_first && is_output_from_first
                in_idx = input_index(parent, input.uac)
                out_idx = output_index(parent, output.uac)
                is_allowed = parent.energy_flow[in_idx,out_idx]

            elseif is_input_from_first && !is_output_from_first
                in_idx = input_index(parent, input.uac)
                out_idx = output_index(parent, child.uac)
                is_allowed = parent.energy_flow[in_idx,out_idx]
                in_idx = input_index(child, parent.uac)
                out_idx = output_index(child, output.uac)
                is_allowed += child.energy_flow[in_idx,out_idx]
                is_allowed = is_allowed == 2 ? 1 : 0

            elseif !is_input_from_first && is_output_from_first
                # inputs from the second bus can't "see" outputs of the first
                is_allowed = 0

            else
                in_idx = input_index(child, input.uac)
                out_idx = output_index(child, output.uac)
                is_allowed = child.energy_flow[in_idx,out_idx]
            end

            new_node.energy_flow[new_in_idx,new_out_idx] = is_allowed
        end
    end

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
    indices = indexin([child], bus_outputs(parent))
    return length(indices) != 1 || indices[1] === nothing ? 0 : indices[1]
end

function priority_of_parent(parent::BusNode, child::BusNode)::Integer
    indices = indexin([parent], bus_inputs(child))
    return length(indices) != 1 || indices[1] === nothing ? 0 : indices[1]
end

function find_single_parent_sinks(
    nodes::Dict{String,BusNode}
)::Vector{Tuple{Integer,BusNode,BusNode}}
    sinks = []
    for node in values(nodes)
        if length(bus_inputs(node)) != 1 || length(bus_outputs(node)) > 0
            continue
        end
        parent = bus_inputs(node)[1]
        push!(sinks, (priority_of_child(parent, node), parent, node))
    end
    return sinks
end

function find_parents_of_sinks(
    nodes::Dict{String,BusNode}
)::Vector{Tuple{Integer,BusNode,BusNode}}
    parents = []
    for node in values(nodes)
        if length(bus_outputs(node)) > 0
            continue
        end
        for parent in bus_inputs(node)
            push!(parents, (priority_of_parent(parent, node), parent, node))
        end
    end
    return parents
end

function bus_from_node(node::BusNode, template::Bus, components::Grouping)::Bus
    bus = Bus(node.uac, template.medium, template.epsilon)

    # create inputs
    for (idx,input) in pairs(node.inputs)
        push!(bus.input_interfaces, SystemInterface(
            source=components[input.uac],
            target=bus,
            do_storage_transfer=input.do_storage_transfer
        ))
        push!(bus.connectivity.input_order, input.uac)
        bus.balance_table_inputs[input.uac] = BTInputRow(
            source=components[input.uac],
            priority=idx,
            input_index=idx,
            do_storage_transfer=input.do_storage_transfer
        )
    end

    # create outputs
    for (idx,output) in pairs(node.outputs)
        push!(bus.output_interfaces, SystemInterface(
            source=bus,
            target=components[output.uac],
            do_storage_transfer=output.do_storage_transfer
        ))
        push!(bus.connectivity.output_order, output.uac)
        bus.balance_table_outputs[output.uac] = BTOutputRow(
            target=components[output.uac],
            priority=idx,
            output_index=idx,
            do_storage_transfer=output.do_storage_transfer
        )
    end

    # for now the energy flow matrices have the same content, but different structure,
    # hence we need to set values individually
    bus.connectivity.energy_flow = []
    for row_idx in keys(bus.connectivity.input_order)
        row = []
        for col_idx in keys(bus.connectivity.output_order)
            push!(row, node.energy_flow[row_idx, col_idx])
        end
        push!(bus.connectivity.energy_flow, row)
    end

    # reset rebuilds the balance table
    reset(bus)

    return bus
end

function merge_busses(busses_to_merge::Grouping, components::Grouping)::Union{Nothing,Bus}
    nodes = nodes_from_components(busses_to_merge)

    # keep merging until we have only one node left
    while length(nodes) > 1
        next = nothing

        # first remove all linear connections
        candidates = find_linear_connections(nodes)
        if length(candidates) > 0
            next = candidates[1] # no particular order
        end

        # then prioritise sinks with only one parent and order by descending priority
        if next === nothing
            candidates = find_single_parent_sinks(nodes)
            if length(candidates) > 0
                next = sort(candidates, by=x->x[1], rev=true)[1]
            end
        end

        # finally, for sinks with multiple parents, merge into those with the last priority
        if next === nothing
            candidates = find_parents_of_sinks(nodes)
            if length(candidates) > 0
                next = sort(candidates, by=x->x[1], rev=true)[1]
            end
        end

        # now merge the candidate parent->child pair and update nodes
        if next !== nothing
            parent = next[2]
            child = next[3]

            new_uac = parent.uac * child.uac
            new_node = merge(parent, child, new_uac)
            nodes[new_uac] = new_node

            update_nodes!(nodes, parent.uac, child.uac, new_node)
            delete!(nodes, parent.uac)
            delete!(nodes, child.uac)
        else
            # if we got here, something went wrong and would result in an endless loop
            @error("Could not resolve merging of busses")
            break
        end
    end

    return bus_from_node(
        first(values(nodes)),
        first(values(busses_to_merge)),
        components
    )
end
