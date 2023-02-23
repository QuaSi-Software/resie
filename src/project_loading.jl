# this file contains functionality pertaining to loading a project's metadata and the
# energy systems from the project config file, as well as constructing certain helpful
# information data structures from the inputs in the config

import JSON

"""
    read_JSON(filepath)

Read and parse the JSON-encoded Dict in the given file.
"""
function read_JSON(filepath :: String) :: Dict{AbstractString, Any}
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
function load_systems(config :: Dict{String, Any}) :: Grouping
    systems = Grouping()

    for (unit_key, entry) in pairs(config)
        default_dict = Dict{String, Any}(
            "strategy" => Dict{String, Any}("name" => "default")
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
function order_of_operations(systems :: Grouping) :: StepInstructions
    systems_by_function = [
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

    simulation_order = []
    initial_nr = length(systems) * 100

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

    # produce all systems except dispatchable sources/sinks. the order corresponds
    # to the general order of system functions
    for sf_order = 1:5
        for unit in values(systems_by_function[sf_order])
            push!(simulation_order, [initial_nr, (unit, EnergySystems.s_produce)])
            initial_nr -= 1
        end
    end

    # sandwich loading of storage systems between the other systems and dispatchable ones
    for unit in values(systems_by_function[5])
        push!(simulation_order, [initial_nr, (unit, EnergySystems.s_load)])
        initial_nr -= 1
    end

    # handle dispatchable sources/sinks
    for sf_order = 6:7
        for unit in values(systems_by_function[sf_order])
            push!(simulation_order, [initial_nr, (unit, EnergySystems.s_produce)])
            initial_nr -= 1
        end
    end

    # finally, distribute bus systems
    for unit in values(systems_by_function[3])
        push!(simulation_order, [initial_nr, (unit, EnergySystems.s_distribute)])
        initial_nr -= 1
    end

    # helper function to find certain steps in the simulation order
    idx_of = function(order, uac, step)
        for idx in eachindex(order)
            if order[idx][2][1].uac == uac && order[idx][2][2] == step
                return idx
            end
        end
        return 0
    end

    # reorder systems connected to a bus so they match the input priority:
    for bus in values(systems_by_function[3])
        # for each system in the bus' input priority...
        for own_idx = 1:length(bus.input_priorities)
            own_uac = bus.input_priorities[own_idx]
            own_ctrl_idx = idx_of(simulation_order, own_uac, EnergySystems.s_control)
            own_prod_idx = idx_of(simulation_order, own_uac, EnergySystems.s_produce)

            # ...make sure every system following after...
            for other_idx = own_idx:length(bus.input_priorities)
                other_uac = bus.input_priorities[other_idx]
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

    # helper function to check if target system is bus
    uac_is_bus = function(energysystem,uac)
        bool = false
        for output_interface in energysystem.output_interfaces
            if output_interface.target.uac === uac
                if output_interface.target.sys_function === EnergySystems.sf_bus
                    bool = true
                end
            end
        end
        return bool    
    end

    # reorder distribution of busses and load of storages
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
        # also make sure that the order of the load() function of the storages connected to the following busses have the same order than the busses
        #
        # for bus in the first bus' production refs...
        for own_idx = 1:length(bus.output_priorities)
            own_uac = bus.output_priorities[own_idx]

            # predefine indexes for storages, set to zero
            own_storage_idx = 0
            other_storage_idx = 0

            if uac_is_bus(bus, own_uac) # consider only busses
                own_dist_idx = idx_of(simulation_order, own_uac, EnergySystems.s_distribute)

                # check if a storage is connected to own bus
                for bus_output_interface in bus.output_interfaces # seach interfaces for own bus
                    if bus_output_interface.target.uac == own_uac  # found correct interface
                        for own_bus_output_interfaces in bus_output_interface.target.output_interfaces # seach interface.targets for storage
                            if own_bus_output_interfaces.target.sys_function == EnergySystems.sf_storage
                                own_storage_uac = own_bus_output_interfaces.target.uac
                                own_storage_idx = idx_of(simulation_order, own_storage_uac, EnergySystems.s_load)
                            end
                        end
                    end
                end

                # ...make sure every bus following after...
                for other_idx = own_idx:length(bus.output_priorities)
                    other_uac = bus.output_priorities[other_idx]
                    if uac_is_bus(bus, other_uac) # consider only busses
                        other_dist_idx = idx_of(simulation_order, other_uac, EnergySystems.s_distribute)

                        # ...is of a lower priority. If not, swap the distribute steps.
                        if simulation_order[own_dist_idx][1] < simulation_order[other_dist_idx][1]
                            tmp = simulation_order[own_dist_idx][1]
                            simulation_order[own_dist_idx][1] = simulation_order[other_dist_idx][1]
                            simulation_order[other_dist_idx][1] = tmp
                        end

                        # check if a storage is connected to other bus
                        for bus_output_interface in bus.output_interfaces # seach interfaces for other bus
                            if bus_output_interface.target.uac == other_uac  # found correct interface
                                for other_bus_output_interfaces in bus_output_interface.target.output_interfaces # seach interface.targets for storage
                                    if other_bus_output_interfaces.target.sys_function == EnergySystems.sf_storage
                                        other_storage_uac = other_bus_output_interfaces.target.uac
                                        other_storage_idx = idx_of(simulation_order, other_storage_uac, EnergySystems.s_load)
                                    end
                                end
                            end
                        end
                                                
                        #... and check if load() of storage of other bus is on lower index than own storage. If not, swap the load steps
                        if own_storage_idx > 0 && other_storage_idx > 0
                            if simulation_order[own_storage_idx][1] < simulation_order[other_storage_idx][1]
                                tmp = simulation_order[own_storage_idx][1]
                                simulation_order[own_storage_idx][1] = simulation_order[other_storage_idx][1]
                                simulation_order[other_storage_idx][1] = tmp
                            end
                        end
                    end
                end
            end
        end      
    end

    # reorder systems such that their control dependencies are handled first, but only if
    # these are not storage systems (which are handled differently)
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

    fn_first = function(entry) return entry[1] end
    return [(u[2][1].uac, u[2][2]) for u in sort(simulation_order, by=fn_first, rev=true)]
end
