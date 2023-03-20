# this file contains functionality pertaining to loading a project's metadata and the
# energy systems from the project config file, as well as constructing certain helpful
# information data structures from the inputs in the config

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
load_systems(config)

Construct instances of energy systems from the given config.

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

The required parameters to construct an energy system from one entry in the config must
match what is required for the particular system. The `type` parameter must be present and
must match the symbol of the energy system class exactly. The structure is described in
more detail in the accompanying documentation on the project file.
"""
function load_systems(config::Dict{String,Any})::Grouping
    systems = Grouping()
    for (unit_key, entry) in pairs(config)
        default_dict = Dict{String,Any}(
            "strategy" => Dict{String,Any}("name" => "default")
        )
        unit_config = merge(default_dict, entry)

        symbol = Symbol(String(unit_config["type"]))
        unit_class = getproperty(EnergySystems, symbol)
        if unit_class <: EnergySystems.EnergySystem
            instance = unit_class(unit_key, unit_config)
            systems[unit_key] = instance
        end
    end

    for (unit_key, entry) in pairs(config)
        if length(entry["control_refs"]) > 0
            others = Grouping(key => systems[key] for key in entry["control_refs"])
            link_control_with(systems[unit_key], others)
        end

        if length(entry["production_refs"]) > 0
            others = Grouping(key => systems[key] for key in entry["production_refs"])
            link_production_with(systems[unit_key], others)
        end
    end

    return systems
end

"""
    categorize_by_function(systems)

Sort the given systems into buckets by their system functions.
"""
function categorize_by_function(systems)
    return [
        [unit for unit in each(systems)
         if unit.sys_function == EnergySystems.sf_fixed_source],
        [unit for unit in each(systems)
         if unit.sys_function == EnergySystems.sf_fixed_sink],
        [unit for unit in each(systems)
         if unit.sys_function == EnergySystems.sf_bus],
        [unit for unit in each(systems)
         if unit.sys_function == EnergySystems.sf_transformer],
        [unit for unit in each(systems)
         if unit.sys_function == EnergySystems.sf_storage],
        [unit for unit in each(systems)
         if unit.sys_function == EnergySystems.sf_dispatchable_source],
        [unit for unit in each(systems)
         if unit.sys_function == EnergySystems.sf_dispatchable_sink],
    ]
end

"""
    base_order(systems_by_function)

Calculate the base order for the simulation steps.

This is determined by the system functions having a certain "natural" order as well as the
simulation steps having a natural order as well.
"""
function base_order(systems_by_function)
    simulation_order = []
    initial_nr = sum([length(bucket) for bucket in systems_by_function]) * 100

    # reset all systems, order doesn't matter
    for sf_order = 1:7
        for unit in values(systems_by_function[sf_order])
            push!(simulation_order, [initial_nr, (unit, EnergySystems.s_reset)])
            initial_nr -= 1
        end
    end


    # calculate control of all systems. the order corresponds to the general order of
    # system functions
    for sf_order = 1:7
        for unit in values(systems_by_function[sf_order])
            push!(simulation_order, [initial_nr, (unit, EnergySystems.s_control)])
            initial_nr -= 1
        end
    end

    # produce fixed sources/sinks and busses.
    for sf_order = 1:3
        for unit in values(systems_by_function[sf_order])
            push!(simulation_order, [initial_nr, (unit, EnergySystems.s_produce)])
            initial_nr -= 1
        end
    end

    # place steps potential and produce for transformers in order by "chains"
    chains = find_chains(systems_by_function[4], EnergySystems.sf_transformer)
    for chain in chains
        for unit in iterate_chain(chain, EnergySystems.sf_transformer, reverse=true)
            push!(simulation_order, [initial_nr, (unit, EnergySystems.s_produce)])
            initial_nr -= 1
        end
    end

    # produce, then load storages
    for unit in values(systems_by_function[5])
        push!(simulation_order, [initial_nr, (unit, EnergySystems.s_produce)])
        initial_nr -= 1
    end
    for unit in values(systems_by_function[5])
        push!(simulation_order, [initial_nr, (unit, EnergySystems.s_load)])
        initial_nr -= 1
    end

    # produce dispatchable sources/sinks
    for sf_order = 6:7
        for unit in values(systems_by_function[sf_order])
            push!(simulation_order, [initial_nr, (unit, EnergySystems.s_produce)])
            initial_nr -= 1
        end
    end

    # distribute busses
    for unit in values(systems_by_function[3])
        push!(simulation_order, [initial_nr, (unit, EnergySystems.s_distribute)])
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
    uac_is_bus(energysystem, uac)

Helper function to check if the given UAC corresponds to a bus in the outputs of the
given energy system.
"""
function uac_is_bus(energysystem, uac)
    for output_interface in energysystem.output_interfaces
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
    find_chains(systems, sys_function)

Find all chains of the given system function in the given collection of systems.

A chain is a subgraph of the graph spanned by all connections of the given energy systems,
which is a directed graph. The subgraph is defined by all connected systems of the given
system function.
"""
function find_chains(systems, sys_function)
    chains = []

    for unit in systems
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
            distance = distance_to_sink(outface.target, sys_function)
            max_distance = distance > max_distance ? distance : max_distance
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
-`chain::Set`: A chain as a set of energy systems
-`sys_function`: The system function of the energy systems in the chain
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
order_of_operations(systems)

Calculate the order of steps that need to be performed to simulate the given systems.

This function works by an algorithm described in more detail in the accompanying
documentation. The result of this are step-by-step instructions telling the simulation
engine in which order the system functions are performed for each unit. This algorithm is
not trivial and might not work for each possible grouping of systems.

# Args
- `system::Grouping`: The systems for which an order is required
# Returns
- `StepInstructions`: The order in the structure:
```
[
    ["UAC Key", s_step],
    ...
]
```
"""
function order_of_operations(systems::Grouping)::StepInstructions
    systems_by_function = categorize_by_function(systems)
    simulation_order = base_order(systems_by_function)

    reorder_for_input_priorities(simulation_order, systems, systems_by_function)
    reorder_distribution_of_busses(simulation_order, systems, systems_by_function)
    reorder_storage_loading(simulation_order, systems, systems_by_function)
    reorder_for_control_dependencies(simulation_order, systems, systems_by_function)

    fn_first = function (entry)
        return entry[1]
    end
    return [(u[2][1].uac, u[2][2]) for u in sort(simulation_order, by=fn_first, rev=true)]
end

"""
    reorder_for_input_priorities(simulation_order, systems, systems_by_function)

Reorder systems connected to a bus so they match the input priority defined on that bus.
"""
function reorder_for_input_priorities(simulation_order, systems, systems_by_function)
    for bus in values(systems_by_function[3])
        # for each system in the bus' input priority...
        for own_idx = 1:length(bus.connectivity.input_order)
            own_uac = bus.connectivity.input_order[own_idx]
            own_ctrl_idx = idx_of(simulation_order, own_uac, EnergySystems.s_control)
            own_prod_idx = idx_of(simulation_order, own_uac, EnergySystems.s_produce)

            # ...make sure every system following after...
            for other_idx = own_idx:length(bus.connectivity.input_order)
                other_uac = bus.connectivity.input_order[other_idx]
                other_ctrl_idx = idx_of(simulation_order, other_uac, EnergySystems.s_control)
                other_prod_idx = idx_of(simulation_order, other_uac, EnergySystems.s_produce)

                # ...is of a lower priority. if not, swap the control steps...
                if simulation_order[own_ctrl_idx][1] < simulation_order[other_ctrl_idx][1]
                    tmp = simulation_order[own_ctrl_idx][1]
                    simulation_order[own_ctrl_idx][1] = simulation_order[other_ctrl_idx][1]
                    simulation_order[other_ctrl_idx][1] = tmp
                end

                # ...and swap the production too.
                if simulation_order[own_prod_idx][1] < simulation_order[other_prod_idx][1]
                    tmp = simulation_order[own_prod_idx][1]
                    simulation_order[own_prod_idx][1] = simulation_order[other_prod_idx][1]
                    simulation_order[other_prod_idx][1] = tmp
                end
            end
        end
    end
end

"""
    reorder_distribution_of_busses(simulation_order, systems, systems_by_function)

Reorder the distribution of busses so that any chain of busses connected to each other
have the "sink" busses before the "source" busses while also considering the output
priorities of two or more sink busses connected to the same source bus.

In the following, assume energy flow from left to right:

                ------------
                 Sink Bus 1
                ------------
              /
------------
Source Bus 1
------------
              \
                ------------
                Sink Bus 2
                ------------
"""
function reorder_distribution_of_busses(simulation_order, systems, systems_by_function)
    for bus in values(systems_by_function[3])
        # make sure that every following bus connected to a fist bus is calculated earlier than the first bus (only distribute here)
        own_uac = bus.uac
        own_idx = idx_of(simulation_order, own_uac, EnergySystems.s_distribute)
        # for every bus in the bus' output_interface...
        for other_unit in bus.output_interfaces
            if other_unit.target.sys_function == EnergySystems.sf_bus  # consider only busses
                other_uac = other_unit.target.uac
                other_idx = idx_of(simulation_order, other_uac, EnergySystems.s_distribute)

                # ...check if it is of a lower priority. If not, swap the distribute steps.
                if simulation_order[own_idx][1] > simulation_order[other_idx][1]
                    tmp = simulation_order[own_idx][1]
                    simulation_order[own_idx][1] = simulation_order[other_idx][1]
                    simulation_order[other_idx][1] = tmp
                end
            end
        end

        # make sure the order of distribution match the order according to the production_refs of the first bus
        # for bus in the first bus' production refs...
        for own_idx = 1:length(bus.connectivity.output_order)
            own_uac = bus.connectivity.output_order[own_idx]

            if uac_is_bus(bus, own_uac) # consider only busses
                own_dist_idx = idx_of(simulation_order, own_uac, EnergySystems.s_distribute)

                # ...make sure every bus following after...
                for other_idx = own_idx:length(bus.connectivity.output_order)
                    other_uac = bus.connectivity.output_order[other_idx]
                    if uac_is_bus(bus, other_uac) # consider only busses
                        other_dist_idx = idx_of(simulation_order, other_uac, EnergySystems.s_distribute)

                        # ...is of a lower priority. If not, swap the distribute steps.
                        if simulation_order[own_dist_idx][1] < simulation_order[other_dist_idx][1]
                            tmp = simulation_order[own_dist_idx][1]
                            simulation_order[own_dist_idx][1] = simulation_order[other_dist_idx][1]
                            simulation_order[other_dist_idx][1] = tmp
                        end
                    end
                end
            end
        end
    end
end

"""
    reorder_storage_loading(simulation_order, systems, systems_by_function)

Reorder systems such the loading (and unloading) of storages follows the priorities on
busses, including communication across connected busses.
"""
function reorder_storage_loading(simulation_order, systems, systems_by_function)
    for bus in values(systems_by_function[3])
        # make sure that the order of the load() function of the storages connected to the following busses have the same order than the busses.
        # for every bus in the first bus' production refs...
        for own_idx = 1:length(bus.connectivity.output_order)
            own_uac = bus.connectivity.output_order[own_idx]

            # predefine indexes for storages, set to zero
            own_storage_load_idx = 0
            own_storage_produce_idx = 0
            other_storage_load_idx = 0
            other_storage_produce_idx = 0

            if uac_is_bus(bus, own_uac) # consider only busses
                # check if a storage is connected to own bus
                for bus_output_interface in bus.output_interfaces # seach interfaces of original bus to find own bus
                    if bus_output_interface.target.uac == own_uac  # found correct interface
                        for own_bus_output_interfaces in bus_output_interface.target.output_interfaces # seach in interface.targets for interconnected storage
                            if own_bus_output_interfaces.target.sys_function == EnergySystems.sf_storage
                                own_storage_uac = own_bus_output_interfaces.target.uac
                                own_storage_produce_idx = idx_of(simulation_order, own_storage_uac, EnergySystems.s_produce)
                                own_storage_load_idx = idx_of(simulation_order, own_storage_uac, EnergySystems.s_load)
                            end
                        end
                    end
                end

                # ...make sure every storage connected to a bus following own_bus in the output priorities of the original bus...
                for other_idx = own_idx:length(bus.connectivity.output_order)
                    other_uac = bus.connectivity.output_order[other_idx]
                    if uac_is_bus(bus, other_uac) # consider only busses
                        # check if a storage is connected to other bus
                        for bus_output_interface in bus.output_interfaces # seach interfaces for other bus
                            if bus_output_interface.target.uac == other_uac  # found correct interface
                                for other_bus_output_interfaces in bus_output_interface.target.output_interfaces # seach interface.targets for storage
                                    if other_bus_output_interfaces.target.sys_function == EnergySystems.sf_storage
                                        other_storage_uac = other_bus_output_interfaces.target.uac
                                        other_storage_produce_idx = idx_of(simulation_order, other_storage_uac, EnergySystems.s_produce)
                                        other_storage_load_idx = idx_of(simulation_order, other_storage_uac, EnergySystems.s_load)
                                    end
                                end
                            end
                        end

                        #... is on lower index. So check if load() of storage of other bus is on lower index than own storage. If not, swap the load steps
                        if own_storage_load_idx > 0 && other_storage_load_idx > 0
                            if simulation_order[own_storage_load_idx][1] < simulation_order[other_storage_load_idx][1]
                                tmp = simulation_order[own_storage_load_idx][1]
                                simulation_order[own_storage_load_idx][1] = simulation_order[other_storage_load_idx][1]
                                simulation_order[other_storage_load_idx][1] = tmp
                            end
                        end
                        #... and also check if produce() of storage of other bus is on lower index than own storage. If not, swap the load steps
                        if own_storage_produce_idx > 0 && other_storage_produce_idx > 0
                            if simulation_order[own_storage_produce_idx][1] < simulation_order[other_storage_produce_idx][1]
                                tmp = simulation_order[own_storage_produce_idx][1]
                                simulation_order[own_storage_produce_idx][1] = simulation_order[other_storage_produce_idx][1]
                                simulation_order[other_storage_produce_idx][1] = tmp
                            end
                        end
                    end
                end
            end
        end
    end
end

"""
    reorder_for_control_dependencies(simulation_order, systems, systems_by_function)

Reorder systems such that all of their control dependencies have their control and
produce steps happen before the unit itself. Storage systems are excepted.
"""
function reorder_for_control_dependencies(simulation_order, systems, systems_by_function)
    for unit in values(systems)
        for other_uac in keys(unit.controller.linked_systems)
            other_unit = systems[other_uac]
            if other_unit.sys_function == EnergySystems.sf_storage
                continue
            end

            own_ctrl_idx = idx_of(simulation_order, unit.uac, EnergySystems.s_control)
            own_prod_idx = idx_of(simulation_order, unit.uac, EnergySystems.s_produce)
            other_ctrl_idx = idx_of(simulation_order, other_uac, EnergySystems.s_control)
            other_prod_idx = idx_of(simulation_order, other_uac, EnergySystems.s_produce)

            if simulation_order[own_ctrl_idx] > simulation_order[other_ctrl_idx]
                tmp = simulation_order[own_ctrl_idx][1]
                simulation_order[own_ctrl_idx][1] = simulation_order[other_ctrl_idx][1]
                simulation_order[other_ctrl_idx][1] = tmp
            end

            if simulation_order[own_prod_idx][1] > simulation_order[other_prod_idx][1]
                tmp = simulation_order[own_prod_idx][1]
                simulation_order[own_prod_idx][1] = simulation_order[other_prod_idx][1]
                simulation_order[other_prod_idx][1] = tmp
            end
        end
    end
end
