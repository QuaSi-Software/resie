# this file contains functionality pertaining to loading a project's metadata and the
# energy system components from the project config file, as well as constructing certain
# helpful information data structures from the inputs in the config

import JSON

"""
    read_JSON(filepath)

Read and parse the JSON-encoded Dict in the given file.
"""
function read_JSON(filepath::String)::Dict{AbstractString,Any}
    open(filepath, "r") do file_handle
        content = read(file_handle, String)
        return JSON.parse(content)
    end
end

"""
load_components(config, sim_params)

Construct instances of components from the given config.

The config must have the structure:
```
{
"UAC key": {
    "type": "PVPlant",
    ...
},
...
}
```

The required config to construct a component from one entry in the config must match what is
required for the particular component. The `type` parameter must be present and must match
the symbol of the component class exactly. The structure is described in more detail in the
accompanying documentation on the project file.
"""
function load_components(config::Dict{String,Any}, sim_params::Dict{String,Any})::Grouping
    components = Grouping()
    for (unit_key, entry) in pairs(config)
        default_dict = Dict{String,Any}(
            "strategy" => Dict{String,Any}("name" => "default")
        )
        unit_config = merge(default_dict, entry)

        symbol = Symbol(String(unit_config["type"]))
        unit_class = getproperty(EnergySystems, symbol)
        if unit_class <: EnergySystems.Component
            instance = unit_class(unit_key, unit_config, sim_params)
            components[unit_key] = instance
        end
    end

    for (unit_key, entry) in pairs(config)
        if length(entry["control_refs"]) > 0
            others = Grouping(key => components[key] for key in entry["control_refs"])
            link_control_with(components[unit_key], others)
        end

        if (
            String(entry["type"]) != "Bus"
            && haskey(entry, "output_refs")
            && length(entry["output_refs"]) > 0
        )
            others = Grouping(key => components[key] for key in entry["output_refs"])
            link_output_with(components[unit_key], others)

        elseif (
            String(entry["type"]) == "Bus"
            && haskey(entry, "connections")
            && length(entry["connections"]) > 0
        )
            others = Grouping(
                key => components[key]
                for key in entry["connections"]["output_order"]
            )
            link_output_with(components[unit_key], others)
        end
    end

    # the input/output interfaces of busses are constructed in the order of appearance in
    # the config, so after all components are loaded they need to be reordered to match
    # the input/output priorities
    components = reorder_interfaces_of_busses(components)

    EnergySystems.initialise_components(components, sim_params)

    return components
end

"""
reorder_interfaces_of_busses(components)

Reorder the input and output interfaces of busses according to their input and output
priorities given in the connectivity matrix.
"""
function reorder_interfaces_of_busses(components::Grouping)::Grouping
    for unit in each(components)
        if unit.sys_function == EnergySystems.sf_bus
            # get correct order according to connectivity matrix
            output_order = unit.connectivity.output_order
            input_order = unit.connectivity.input_order

            # check for misconfigured bus (it should have at least one input and at least
            # one output)
            if length(input_order) == 0 || length(output_order) == 0
                continue
            end

            # Create a dictionary to map 'uac' to its correct position
            output_order_dict = Dict(uac => idx for (idx, uac) in enumerate(output_order))
            input_order_dict = Dict(uac => idx for (idx, uac) in enumerate(input_order))

            # Get the permutation indices that would sort the 'source'/'target' field by
            # 'uac' order
            output_perm_indices = sortperm([
                output_order_dict[unit.output_interfaces[i].target.uac]
                for i in 1:length(unit.output_interfaces)
            ])
            input_perm_indices = sortperm([
                input_order_dict[unit.input_interfaces[i].source.uac]
                for i in 1:length(unit.input_interfaces)
            ])

            # Reorder the input and output interfaces using the permutation indices
            unit.output_interfaces = unit.output_interfaces[output_perm_indices]
            unit.input_interfaces = unit.input_interfaces[input_perm_indices]
        end
    end
    return components
end

"""
    categorize_by_function(components)

Sort the given components into buckets by their system function.
"""
function categorize_by_function(components)
    return [
        [unit for unit in each(components)
         if unit.sys_function == EnergySystems.sf_fixed_source],
        [unit for unit in each(components)
         if unit.sys_function == EnergySystems.sf_fixed_sink],
        [unit for unit in each(components)
         if unit.sys_function == EnergySystems.sf_bus],
        [unit for unit in each(components)
         if unit.sys_function == EnergySystems.sf_transformer],
        [unit for unit in each(components)
         if unit.sys_function == EnergySystems.sf_storage],
        [unit for unit in each(components)
         if unit.sys_function == EnergySystems.sf_bounded_source],
        [unit for unit in each(components)
         if unit.sys_function == EnergySystems.sf_bounded_sink],
    ]
end

"""
    base_order(components_by_function)

Calculate the base order for the simulation steps.

This is determined by the system functions having a certain "natural" order as well as the
simulation steps having a natural order as well.
"""
function base_order(components_by_function)
    simulation_order = []
    initial_nr = sum([length(bucket) for bucket in components_by_function]) * 100

    # reset all components, order doesn't matter
    for sf_order = 1:7
        for unit in values(components_by_function[sf_order])
            push!(simulation_order, [initial_nr, (unit.uac, EnergySystems.s_reset)])
            initial_nr -= 1
        end
    end


    # calculate control of all components. the order corresponds to the general order of
    # system functions
    for sf_order = 1:7
        for unit in values(components_by_function[sf_order])
            push!(simulation_order, [initial_nr, (unit.uac, EnergySystems.s_control)])
            initial_nr -= 1
        end
    end

    # process fixed sources/sinks and busses.
    for sf_order = 1:3
        for unit in values(components_by_function[sf_order])
            push!(simulation_order, [initial_nr, (unit.uac, EnergySystems.s_process)])
            initial_nr -= 1
        end
    end

    # place steps potential and process for transformers in order by "chains"
    chains = find_chains(components_by_function[4], EnergySystems.sf_transformer)
    for chain in chains
        nr_of_remaining_units_in_chain = length(chain)
        if nr_of_remaining_units_in_chain > 1
            for unit in iterate_chain(chain, EnergySystems.sf_transformer, reverse=false)
                if nr_of_remaining_units_in_chain > 1  # skip potential of last transformer in the chain
                    push!(simulation_order, [initial_nr, (unit.uac, EnergySystems.s_potential)])
                    initial_nr -= 1
                    nr_of_remaining_units_in_chain -= 1
                end
            end
        end
        for unit in iterate_chain(chain, EnergySystems.sf_transformer, reverse=true)
            push!(simulation_order, [initial_nr, (unit.uac, EnergySystems.s_process)])
            initial_nr -= 1
        end
    end

    # process, then load storages
    for unit in values(components_by_function[5])
        push!(simulation_order, [initial_nr, (unit.uac, EnergySystems.s_process)])
        initial_nr -= 1
    end
    for unit in values(components_by_function[5])
        push!(simulation_order, [initial_nr, (unit.uac, EnergySystems.s_load)])
        initial_nr -= 1
    end

    # process bounded sources/sinks
    for sf_order = 6:7
        for unit in values(components_by_function[sf_order])
            push!(simulation_order, [initial_nr, (unit.uac, EnergySystems.s_process)])
            initial_nr -= 1
        end
    end

    # distribute busses
    for unit in values(components_by_function[3])
        push!(simulation_order, [initial_nr, (unit.uac, EnergySystems.s_distribute)])
        initial_nr -= 1
    end

    return simulation_order
end

"""
    idx_of(order, uac, step)

Helper function to find the index of a given combination of step and unit (by its UAC).
"""
function idx_of(order, uac, step)
    for idx in eachindex(order)
        if order[idx][2][1].uac == uac && order[idx][2][2] == step
            return idx
        end
    end
    return 0
end

"""
    uac_is_bus(component, uac)

Helper function to check if the given UAC corresponds to a bus in the outputs of the
given component.
"""
function uac_is_bus(component, uac)
    for output_interface in component.output_interfaces
        if (output_interface.target.uac === uac && output_interface.target.sys_function === EnergySystems.sf_bus)
            return true
        end
    end
    return false
end

"""
    add_non_recursive!(node_set, unit, caller_uac)

Add connected units of the same system function to the node set in a non-recursive manner.
"""
function add_non_recursive!(node_set, unit, caller_uac)
    push!(node_set, unit)

    for inface in values(unit.input_interfaces)
        if inface.source.uac == caller_uac
            continue
        elseif inface.source.sys_function === unit.sys_function
            add_non_recursive!(node_set, inface.source, unit.uac)
        end
    end

    for outface in values(unit.output_interfaces)
        if outface.target.uac == caller_uac
            continue
        elseif outface.target.sys_function === unit.sys_function
            add_non_recursive!(node_set, outface.target, unit.uac)
        end
    end
end

"""
    find_chains(components, sys_function)

Find all chains of the given system function in the given collection of components.

A chain is a subgraph of the graph spanned by all connections of the given components,
which is a directed graph. The subgraph is defined by all connected components of the given
system function.
"""
function find_chains(components, sys_function)
    chains = []

    for unit in components
        if unit.sys_function !== sys_function
            continue
        end

        already_recorded = false
        for chain in chains
            if unit in chain
                already_recorded = true
            end
        end

        if !already_recorded
            chain = Set()
            add_non_recursive!(chain, unit, "")
            push!(chains, chain)
        end
    end

    return chains
end

"""
    distance_to_sink(node, sys_function)

Calculate the distance of the given node to the sinks of the chain.

A sink is defined as a node with no successors of the same system function. For the sinks
this distance is 0. For all other nodes it is the maximum over the distances of its
successors plus one.
"""
function distance_to_sink(node, sys_function)
    is_leaf = function(node)
        for outface in values(node.output_interfaces)
            if outface.target.sys_function === sys_function
                return false
            end
        end
        return true
    end

    if is_leaf(node)
        return 0
    else
        max_distance = 0
        for outface in values(node.output_interfaces)
            if outface.target.sys_function == sys_function
                distance = distance_to_sink(outface.target, sys_function)
                max_distance = distance > max_distance ? distance : max_distance
            end
        end
        return max_distance + 1
    end
end

"""
    iterate_chain(chain, sys_function, reverse)

Iteration order over a chain of units of the same system function.

This order is determined by the distance of a node to the sinks of the subgraph (the chain),
where the distance is maxed over the successors of a node. The order is then an ascending
ordering over the distances of nodes.

# Arguments
-`chain::Set`: A chain as a set of components
-`sys_function`: The system function of the components in the chain
-`reverse::Bool`: If true, the ordering will be in descending order
"""
function iterate_chain(chain, sys_function; reverse=false)
    distances = []
    for node in chain
        push!(distances, (distance_to_sink(node, sys_function), node))
    end
    fn_first = function (entry)
        return entry[1]
    end
    return [u[2] for u in sort(distances, by=fn_first, rev=reverse)]
end

"""
calculate_order_of_operations(components)

Calculate the order of steps that need to be performed to simulate the given components.

This function works by an algorithm described in more detail in the accompanying
documentation. The result of this are step-by-step instructions telling the simulation
engine in which order the system functions are performed for each unit. This algorithm is
not trivial and might not work for each possible grouping of components.

# Args
- `components::Grouping`: The components for which an order is required
# Returns
- `StepInstructions`: The order in the structure:
```
[
    ["UAC Key", s_step],
    ...
]
```
"""
function calculate_order_of_operations(components::Grouping)::StepInstructions
    components_by_function = categorize_by_function(components)
    simulation_order = base_order(components_by_function)

    reorder_for_input_priorities(simulation_order, components, components_by_function)
    reorder_distribution_of_busses(simulation_order, components, components_by_function)
    reorder_storage_loading(simulation_order, components, components_by_function)

    fn_first = function (entry)
        return entry[1]
    end
    return [(u[2][1], u[2][2]) for u in sort(simulation_order, by=fn_first, rev=true)]
end

"""
load_order_of_operations(order_of_operation_input, components)

Reads in the order of operation given in the input file.
The Vektor that is passed has to have one entry for each operation. All operations have to 
be in the following syntax:
[
    "UAC_Key s_step",
    ...
]
"UAC_Key s_step" can be "TST_01_HZG_02_DEM s_reset" for example.

# Args
- `order_of_operation_input::Vector{String}`: A Vektor of strings containing the operations 
                                              given in the input file 
- `components::Grouping`: All components of the input file
# Returns
- `step_instructions`: The order of opertaions given in the input file in the structure:
```
[
    ["UAC Key", s_step],
    ...
]
```
"""
function load_order_of_operations(order_of_operation_input, components::Grouping)::StepInstructions
    step_instructions = []
    all_components_uac = collect(keys(components)) # [unit.uac for unit in keys(components)]

    for entry in order_of_operation_input 
        uac = string(split(entry)[1])
        s_step = split(entry)[2]

        if s_step == "s_reset"
           s_step_component = EnergySystems.s_reset
        elseif s_step == "s_control"
            s_step_component = EnergySystems.s_control
        elseif s_step == "s_process"
            s_step_component = EnergySystems.s_process
        elseif s_step == "s_potential"
            s_step_component = EnergySystems.s_potential
        elseif s_step == "s_load"
            s_step_component = EnergySystems.s_load
        elseif s_step == "s_distribute"
            s_step_component = EnergySystems.s_distribute
        else 
            throw(ArgumentError("Unknown operation step \"$s_step\" in \"order_of_operation\" in input file. Must be one of: s_reset, s_control, s_process, s_potential, s_load, s_distribute"))
        end

        if !(uac in all_components_uac)
            throw(ArgumentError("Unknown component UAC \"$uac\" in \"order_of_operation\" given in input file. Each UAC must match the name of a component given in the input file."))
        end

        push!(step_instructions, [uac, s_step_component] )
    end
    return [(u[1], u[2]) for u in step_instructions]
end


"""
    find_indexes(steps, own, target)

Find the indexes of the specified two step instructions in the given collection.

# Arguments
-`steps::Vector{Tuple{Integer,Tuple{String,SystemFunction}}}`: The list of step instructions
-`own::Tuple{String,SystemFunction}`: The first step instruction to look for
-`target::Tuple{String,SystemFunction}`: The second step instruction to look for
# Returns
-`Integer`: Index in the collection of the first instruction
-`Integer`: Index in the collection of the second instruction
"""
function find_indexes(steps, own, target)
    own_idx = nothing
    target_idx = nothing

    for i in eachindex(steps)
        if steps[i][2][1] == own[1] && steps[i][2][2] == own[2]
            own_idx = i
        end
        if steps[i][2][1] == target[1] && steps[i][2][2] == target[2]
            target_idx = i
        end
    end

    return (own_idx, target_idx)
end

"""
    place_one!(steps, own, target, higher=true)

Give the target step instruction a priority one higher/lower than the given one.

# Arguments
-`steps::Vector{Tuple{Integer,Tuple{String,SystemFunction}}}`: The list of step instructions
-`own::Tuple{String,SystemFunction}`: The first step instruction to look for
-`target::Tuple{String,SystemFunction}`: The second step instruction to look for
-`higher::Bool`: If true places the target instruction one higher, otherwise one lower
-`force::Bool`: If true, forces the placement even if the target priority already is
    higher/lower than the first one. If false, no changes are performed if the priority
    already is higher/lower
"""
function place_one!(steps, own, target; higher=true, force=false)
    own_idx, target_idx = find_indexes(steps, own, target)
    if own_idx === nothing || target_idx === nothing
        return
    end

    # check if target priority already is higher/lower
    own_priority = steps[own_idx][1]
    target_priority = steps[target_idx][1]
    if !force && ((higher && target_priority > own_priority)
                  ||
                  (!higher && target_priority < own_priority))
        return
    end

    # shift other units with higher/lower priority up/down by one
    for i in eachindex(steps)
        if (higher && steps[i][1] > own_priority
            ||
            !higher && steps[i][1] < own_priority)

            steps[i][1] += higher ? 1 : -1
        end
    end

    # place target one step higher/lower than own
    steps[target_idx][1] = own_priority + (higher ? 1 : -1)
end

"""
    place_one_higher!(steps, own, target)

Alias to place_one! with the `higher` argument set to true.
"""
function place_one_higher!(steps, own, target; force=false)
    place_one!(steps, own, target, higher=true, force=force)
end

"""
    place_one_lower!(steps, own, target)

Alias to place_one! with the `higher` argument set to false.
"""
function place_one_lower!(steps, own, target; force=false)
    place_one!(steps, own, target, higher=false, force=force)
end

"""
    reorder_for_input_priorities(simulation_order, components, components_by_function)

Reorder components connected to a bus so they match the input priority defined on that bus.
"""
function reorder_for_input_priorities(simulation_order, components, components_by_function)
    # for every bus...
    for bus in values(components_by_function[3])
        # ...for each component in the bus' input priority...
        for own_idx = 1:length(bus.connectivity.input_order)
            # ...make sure every component following after...
            for other_idx = own_idx+1:length(bus.connectivity.input_order)
                # ...is of a lower priority in potential and process
                place_one_lower!(
                    simulation_order,
                    (bus.connectivity.input_order[own_idx], EnergySystems.s_potential),
                    (bus.connectivity.input_order[other_idx], EnergySystems.s_potential)
                )
                place_one_lower!(
                    simulation_order,
                    (bus.connectivity.input_order[own_idx], EnergySystems.s_process),
                    (bus.connectivity.input_order[other_idx], EnergySystems.s_process)
                )
            end
        end
    end
end

"""
    reorder_distribution_of_busses(simulation_order, components, components_by_function)

Reorder the distribution of busses so that any chain of busses connected to each other
have the "sink" busses before the "source" busses while also considering the output
priorities of two or more sink busses connected to the same source bus.

In the following, assume energy flows from left to right:

               ------------
               |   Bus 2  |   ------------
               ------------   |   Bus 5  |
             /              / ------------
------------   ------------
|   Bus 1  |---|   Bus 3  |
------------   -----------`
             `              ` ------------
               `-----------   |   Bus 6  |
               |   Bus 4  |   ------------
               ------------
With input priorities (4,3,2) on bus 1 and (6,5) on bus 3 this would result in an order
of: 4,6,5,3,2,1
"""
function reorder_distribution_of_busses(simulation_order, components, components_by_function)
    # for every bus chain...
    for bus_chain in find_chains(components_by_function[3], EnergySystems.sf_bus)
        # ...by sources-first ordering...
        for bus in iterate_chain(bus_chain, EnergySystems.sf_bus, reverse=true)
            output_order = length(bus.connectivity.output_order) > 0 ?
                           bus.connectivity.output_order :
                           [u.target.uac for u in bus.output_interfaces]
            # ...make sure every successor bus...
            for successor_uac in output_order
                # ...has a higher priority in the order they are appear in the list
                place_one_higher!(
                    simulation_order,
                    (bus.uac, EnergySystems.s_distribute),
                    (successor_uac, EnergySystems.s_distribute)
                )
            end
        end
    end
end


"""
    check_energy_flow(bus, target, source)

Checks if the energy is allowed to flow from the source to the target via the bus.

# Arguments
-`bus::Component`: The bus that is in the middle of source and storage
-`target::Component`: The target component to check
-`source::Union{Component,Nothing}`: The source component to check. If the source is
    nothing, the energy flow is allowed.
"""
function check_energy_flow(bus, target, source)::Bool
    if bus.connectivity.energy_flow === nothing || source === nothing
        # if no energy flow matrix or source is given, flow is assumed to be allowed
        return true
    end

    output_idx = nothing
    for (idx,output_uac) in pairs(bus.connectivity.output_order)
        if output_uac == target.uac
            output_idx = idx
        end
    end

    input_idx = nothing
    for (idx,input_uac) in pairs(bus.connectivity.input_order)
        if input_uac == source.uac
            input_idx = idx
        end
    end

    if input_idx !== nothing && output_idx !== nothing
        return bus.connectivity.energy_flow[input_idx][output_idx] == 1
    else
        return false
    end
end

"""
    find_storages_ordered(bus, components, source)

Finds storage components in the given bus and all successor busses ordered by the output priorities.

# Arguments
-`bus::Component`: The bus from which to start the search
-`components::Grouping`: All components in the energy system
-`source::Component`: The source from which bus was reached out. Can be nothing if bus is first bus.
-`reverse::Bool`: (Optional) If true, orders the storages in reverse order. Defaults to false.
"""
function find_storages_ordered(bus, components, source; reverse=false)
    storages = []
    limited = []
    pushing = reverse ? pushfirst! : push!

    for unit_uac in bus.connectivity.output_order
        unit = components[unit_uac]
        if unit.sys_function == EnergySystems.sf_storage
            pushing(storages, unit)
            if check_energy_flow(bus, unit, source)
                pushing(limited, true)
            else
                pushing(limited, false)
            end
        elseif unit.sys_function == EnergySystems.sf_bus
            storages_returned, limits = find_storages_ordered(unit, components, bus)
            for storage in storages_returned
                pushing(storages, storage)
            end
            for limit in limits
                pushing(limited, limit)
            end
        end
    end

    return storages, limited
end

"""
    reorder_storage_loading(simulation_order, components, components_by_function)

Reorder components such the loading (and unloading) of storages follows the priorities on
busses, including communication across connected busses.
"""
function reorder_storage_loading(simulation_order, components, components_by_function)
    for bus_chain in find_chains(components_by_function[3], EnergySystems.sf_bus)
        for bus in iterate_chain(bus_chain, EnergySystems.sf_bus)
            storages, limitations = find_storages_ordered(bus, components, nothing, reverse=true)
            if length(storages) < 2
                continue
            end

            # by continuosly placing lower than the last element, which has highest priority
            # due to the reverse ordering, the correct order is preserved
            for last_index in 1:length(storages)
                last_element = storages[length(storages)+1-last_index]
                for idx in 1:(length(storages)-last_index)
                    place_one_lower!(
                        simulation_order,
                        (last_element.uac, EnergySystems.s_process),
                        (storages[idx].uac, EnergySystems.s_process)
                    )
                    place_one_lower!(
                        simulation_order,
                        (last_element.uac, EnergySystems.s_load),
                        (storages[idx].uac, EnergySystems.s_load)
                    )
                end
            end

            # make sure that storages on the leaf-busses are loaded right after the processing of transformers on leaf-busses 
            # (only if storage-loading is allowed by the connectivity matrix of the bus)
            # to avoid that other storages in the system uses the not requested energy that was indented to load the storage, 
            # which can otherwise happen in the trunk bus or in other leaf busses as well. 
            for (_, storage) in pairs(storages)                                                         # iterate through all storages in the current bus_chain
                for (_, storage_input_interface) in pairs(storage.input_interfaces)                     # and get input interface of storage. Usually ony one interface is present in storage, but for-loop is used here to identify the interface without the need of the medium name
                    if storage_input_interface.source.sys_function == EnergySystems.sf_bus              # check if source of storage is a bus
                        for bus_input_interface in storage_input_interface.source.input_interfaces      # iterate through input interfaces of this bus
                            if bus_input_interface.source.sys_function == EnergySystems.sf_transformer  # and check if they are fed by transformers. 
                                if (                                                                    # if yes, check if the transfomer is allowed to load the storage
                                    check_energy_flow(storage_input_interface.source, storage_input_interface.target, bus_input_interface.source)
                                    &&                                                                  # and
                                    storage_input_interface.source !== bus                              # do not consider the trunk bus
                                )
                                    place_one_lower!(                                                   # and then make sure that the load step of the storage is right after the process step of the transformer
                                        simulation_order,
                                        (bus_input_interface.source.uac, EnergySystems.s_process),
                                        (storage.uac, EnergySystems.s_load),
                                        force=true
                                    )
                                end
                            end
                        end
                    end
                end
            end

            # sort storages into two list (storages with limited and unlimited loading)
            unimited_storages = []
            limited_storages = []
            for (idx, storage) in pairs(storages)
                if limitations[idx]
                    push!(unimited_storages, storage)
                else
                    push!(limited_storages, storage)
                end
            end 

            # make sure all limited storages are placed behind the unlimited ones
            for limited_storage in limited_storages
                for unimited_storage in unimited_storages
                    place_one_lower!(
                        simulation_order,
                        (unimited_storage.uac, EnergySystems.s_process),
                        (limited_storage.uac, EnergySystems.s_process)
                    )
                    place_one_lower!(
                        simulation_order,
                        (unimited_storage.uac, EnergySystems.s_load),
                        (limited_storage.uac, EnergySystems.s_load)
                    )
                end
            end
        end
    end
end

"""
get_timesteps(simulation_parameters)

Function to read in the time step information from the input file.
If no information is given in the input file, the following defaults 
will be set:
time_step = 900 s
start_timestamp = 0 s
end_timestamp = 900 s
"""
function get_timesteps(simulation_parameters::Dict{String, Any})
    time_step = 900
    if "time_step_seconds" in keys(simulation_parameters)
        time_step = UInt(simulation_parameters["time_step_seconds"])
    end

    start_timestamp = 0
    if "start" in keys(simulation_parameters)
        start_timestamp = Integer(simulation_parameters["start"])
    end

    end_timestamp = 900
    if "end" in keys(simulation_parameters)
        end_timestamp = Integer(simulation_parameters["end"])
    end

    return time_step, start_timestamp, end_timestamp
end