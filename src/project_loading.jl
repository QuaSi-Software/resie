# this file contains functionality pertaining to loading a project's metadata and the
# energy system components from the project config file, as well as constructing certain
# helpful information data structures from the inputs in the config

import JSON

export check_interface_for_transformer

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
        unit_config = Base.merge(default_dict, entry)

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

    chains = find_chains(values(components), EnergySystems.sf_bus)
    EnergySystems.merge_bus_chains(chains, components, sim_params)

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

    # determine and place potential and process steps for transformers
    simulation_order, initial_nr = wrapper_add_transformer_steps(components_by_function, simulation_order, initial_nr)

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
    wrapper_add_transformer_steps(components_by_function, simulation_order, initial_nr)

This is a wrapper for add_transformer_steps(). See the docstring of add_transformer_steps for detailed information!

The wrapper handles the different potential and process steps and ensures that all transformers are detected, 
especially if there are multiple non-connected energy systems (see definition of connection in the docstring of 
add_transformer_steps).
In case not all transformers are part of the simulation_order even after several attempts, the wrapper tries to set 
"reverse," which can help in energy systems with nested "middle busses."
If only one transformer is present in the current energy system, the wrapper just adds the process step of the single
transformer to the simulation_order without calling add_transformer_steps.

# Arguments
-`components_by_function`: All components of the current energy system sorted by their function
-`simulation_order`: The simulation order
-`initial_nr`: The current number of steps

# Returns simulation_order, initial_nr, checked_components
-`simulation_order`: The simulation order
-`initial_nr`: The current number of steps

"""
function wrapper_add_transformer_steps(components_by_function, simulation_order, initial_nr)
    transformers = components_by_function[4]
    if length(transformers) > 1
        transformers_and_busses = vcat(components_by_function[4], components_by_function[3])
        # detect parallel branches in the current energy system
        parallel_branches = find_parallels(transformers_and_busses)

        complete = false
        checked_components_pot = []
        checked_components_pro = []
        count = 1
        reverse = nothing
        while !complete
            simulation_order, initial_nr, checked_components_pot = add_transformer_steps(simulation_order,
                                                                            initial_nr,
                                                                            transformers_and_busses,
                                                                            parallel_branches,
                                                                            "potential",
                                                                            reverse=reverse,
                                                                            checked_components=checked_components_pot)
                                                                                    
            simulation_order, initial_nr, checked_components_pro = add_transformer_steps(simulation_order,
                                                                            initial_nr,
                                                                            transformers_and_busses,
                                                                            parallel_branches,
                                                                            "process",
                                                                            reverse=reverse,
                                                                            checked_components=checked_components_pro)
            
            # check success. If not all components can be found, try to preset a reverse as this might give a solution 
            # in case of multiple middle busses connected in circles.
            complete, transformers_and_busses = check_simulation_order_for_completeness(simulation_order,
                                                                                        transformers_and_busses)
            if count == 10
                reverse = false
            elseif count == 15
                reverse = true
            elseif count == 20
                @warn ("The order of operation is potentially wrong as the process-step of one or more transformers " * 
                       "are missig in the OoO. Check the input and the order in the aux_info file. May specify a custom OoO!")
                break
            end
            count += 1
        end
    elseif length(transformers) == 1
        # if only one transformer is present in the current energy system, we only need the process step
        push!(simulation_order, [initial_nr, (transformers[1].uac, EnergySystems.s_process)])
        initial_nr -= 1
    end
    
    return simulation_order, initial_nr
end

"""
    check_simulation_order_for_completeness(simulation_order, transformers_and_busses)

Checks if a process step for all transformers_and_busses is included in the simulation_order.

# Arguments
-`simulation_order`: The simulation order
-`transformers_and_busses::Array{Grouping}`: All transformers and busses of the current energy system

# Returns
-`complete::Bool`: False if there are process steps of transformers or busses are missing in the simulation_order, 
                   and true if they are all contained in the simulation_order.
-`not_included_components::Array{SystemComponets}`: List of all transformers and busses that are not part of the 
                                                    simulation_order
             
"""
function check_simulation_order_for_completeness(simulation_order, transformers_and_busses)
    not_included_components = []
    complete = true
    for component in transformers_and_busses
        if component.sys_function === EnergySystems.sf_transformer
            step_name = (component.uac, EnergySystems.s_process)
            if step_name in [pair[2] for pair in simulation_order] 
                continue
            else
                push!(not_included_components, component)
                complete = false
            end
        else
            push!(not_included_components, component)
        end
    end

    return complete, not_included_components
end

"""
    add_transformer_steps(simulation_order,
                          initial_nr,
                          current_components,
                          parallel_branches,
                          step_category,
                          reverse=nothing,
                          connecting_component=nothing,
                          checked_components=[])

This function determines the order of the potential or process calculation steps for transformers in a given EnergySystem.
Note that the algorithm may not determine the correct order of operation or find the most computationally efficient solution!

The following logic is used:
- At first, "middle transformers" (MT) are searched within the current_components. An MT is defined as a transformer with 
  at least one other transformer in each of at least two inputs and/or at least one transformer in each of at least 
  two outputs. If there is more than one MT in the current_components, the "first" MT is chosen to start with, 
  where "first" means that there are no other MTs in its input branches. 
- Of the chosen MT, the shorter branch is taken and the add_transformer_steps is called recursively with the 
  transformer of the current branch handed over as new current_components. The order within the branch is crucial,
  as the MT has to always be calculated between its branches. So the first branch starts towards the MT (reverse=false 
  for inputs and reverse=true for outputs), while the second branch is calculated away from the MT (reverse=true 
  for inputs and reverse=false for outputs).
- If there is no MT in the current_components, then "middle busses" (MB) are searched within the current_components.
  MBs are defined as busses with at least three interfaces, each containing at least one other transformer. 
  If there are several MBs in the current_components, the first MB is chosen that has no other MBs in its inputs 
  (reverse === nothing || reverse == false) or in its outputs (reverse == true).
- If the current reverse is false or nothing, first the input branches and then the output branches are given recursively 
  to add_transformer_steps(). If reverse is true, the output branches come first.
- If the current_components are part of a parallel_branch (PB), then they are merged together into one branch. PBs are 
  defined as branches that offer more than one way between two components, with each way containing at least one 
  other transformer (while considering the direction of energy flow). The starting or endpoint of a PB can be an MB or an MT.
  Note that MBs within the parallel branches cannot be considered correctly at the moment!
- If no MB or MT can be found in the current_components, the process or potential of them is determined by their order in 
  the energy system. The direction is set by reverse, which defaults to true (process step) or false (potential step). 
  If one of the components is part of a parallel branch, this is detected and the parallels are handled separately: 
  Depending on the output order (input order is not considered here, but logically they should be the same. Reordering
  for input is done later in calculate_order_of_operations, but should not do anything) of the MB, each of the 
  parallel branches performs its potentials in both ways to ensure that possible limitations reach both ends before 
  continuing with the next branch of the parallel branches.
- Generally: Connections between components are interrupted if:
    - there is a storage (no bypass directly through a storage!)
    - the connection matrix of a bus denies the connection
    - the input component into a bus is connected to an output grid out of the bus (considering the energy matrix)
    - an output component out of a bus is connected to an input grid in the bus (considering the energy matrix)
- To avoid double-counting components, all considered components are tracked in checked_components, which is updated 
  before each recursive function call. Note that in some cases, components can be considered multiple times if they are 
  part of multiple branches at one MB or MT! Therefore, after the order has been determined, only the first appearance of 
  the process is considered and all other process steps of the same component are deleted from the order of operation. 
  For potential steps, this is not done as in some cases, components have to have multiple potential steps, e.g., MTs.
- If an MB is connected to another MB, the first MB holds the second MB in one of its output branches. When determining
  the correct calculation order of the components to the rear MB, the connection branch has to perform its potentials
  after all of the other branches of the rear MB. To detect this, branches can get a connecting_component, defining the
  component that is connected to another MB (or MT) and implicitly also the direction.
  
# Arguments
- `simulation_order::Array[]`: The simulation order that should be extended by transformers' potentials and processes.
- `initial_nr::Int`: The current number of steps in the order of operation.
- `current_components::Array{Grouping}`: The currently considered components. Should be transformers and busses at
                                         initial call.
- `parallel_branches::Dict{Tuple, Array{Array}}`: A dictionary holding all parallels in the current energy system. 
                                                  Key is (starting-point-uac, end-point-uac).
- `step_category::String`: Can be either "potential" or "process", depending on the simulation step to find for 
                           current_components.
  Optional:
- `reverse::Bool`: Indicates if the current_components should be considered in reverse or not. Can also be nothing at 
                   the initial call and defaults to nothing.
- `connecting_component::Component`: Can be either nothing or a connecting component, indicates that the 
                                     current_components are part of a connecting branch. Defaults to nothing.
- `checked_components::Array{Grouping}`: A vector holding all checked components. Defaults to an empty vector. 

# Returns
- `simulation_order::Array[]`: The simulation order extended by transformer step_category.
- `initial_nr::Int`: The current number of steps in the order of operation.
- `checked_components::Array{Grouping}`: The components that have already been checked.
"""
function add_transformer_steps(simulation_order,
                               initial_nr,
                               current_components,
                               parallel_branches,
                               step_category;
                               reverse=nothing,
                               connecting_component=nothing,
                               checked_components=[])

    function add_component_to_checked_components(checked_components, branches)
        if branches !== nothing
            for branch in branches
                for component in branch
                    if !(component in checked_components)
                        push!(checked_components, component)
                    end
                end
            end
        end
        return checked_components
    end

    # detect the first "middle transformers" that has either at more than one input interface 
    # or at more than one output interface at least each one transformer. 
    # "First" means it has no other middle transformer in it's inputs.
    # inface_branches_with_transformers and outface_branches_with_transformers holds the branches 
    # of the input/output interfaces of first_middle_transformer that contain other transformers.
    # Considers only current components, but also possible connected upstream transformers in an 
    # output and the other way around.
    first_middle_transformer, 
            inface_branches_with_transformers,
            outface_branches_with_transformers = detect_first_middle_transformer(current_components, checked_components)

    checked_components = add_component_to_checked_components(checked_components, inface_branches_with_transformers)
    checked_components = add_component_to_checked_components(checked_components, outface_branches_with_transformers)
    if (first_middle_transformer !== nothing && !(first_middle_transformer in checked_components)) 
        push!(checked_components, first_middle_transformer) 
    end

    # Note: No parallels containing middle transformers are currently considered!

    # iterate through all branches of the middle transformer, if it exists, starting with the branch 
    # with the least amount of other transformers.
    if first_middle_transformer !== nothing
        # remove branches without transformers
        indices_to_remove = []
        for (idx,branch) in enumerate(inface_branches_with_transformers)
            if isempty([x for x in branch if x.sys_function === EnergySystems.sf_transformer])
                push!(indices_to_remove, idx)
            end
        end
        for idx in Iterators.reverse(sort(indices_to_remove))
            popat!(inface_branches_with_transformers, idx)
        end
        indices_to_remove = []
        for (idx,branch) in enumerate(outface_branches_with_transformers)
            if isempty([x for x in branch if x.sys_function === EnergySystems.sf_transformer])
                push!(indices_to_remove, idx)
            end
        end
        for idx in Iterators.reverse(sort(indices_to_remove))
            popat!(outface_branches_with_transformers, idx)
        end
        
        # sort the branches by the amount of "seen" transformers and classify them as input or output interface
        # Also reverse the order of the input branches so that they start with the first (nearest to source) element.
        combined_branches_with_transformers = vcat(
            [(collect(Iterators.reverse(inface)),
              length([x for x in inface if x.sys_function === EnergySystems.sf_transformer]),
              true) 
              for inface in inface_branches_with_transformers],
            [(outface,
              length([x for x in outface if x.sys_function === EnergySystems.sf_transformer]),
              false)
              for outface in outface_branches_with_transformers]
        )
        sorted_combined_branches_with_transformers = sort(combined_branches_with_transformers, by=x -> x[2])
        
        middle_transformer_branches = [x[1] for x in sorted_combined_branches_with_transformers]
        is_input = [x[3] for x in sorted_combined_branches_with_transformers]
        
        # move connecting branch to the end of the calculation
        is_connecting_branch = fill(false, length(is_input))
        if connecting_component !== nothing
            for (idx, middle_transformer_branch) in enumerate(middle_transformer_branches)
                for component in middle_transformer_branch
                    if component == connecting_component
                        is_connecting_branch[idx] = true
                        # should be only one.
                        break
                    end
                end
            end
        end
        if any(is_connecting_branch)
            push!(middle_transformer_branches, popat!(middle_transformer_branches, findfirst(is_connecting_branch)))
            push!(is_input, popat!(is_input, findfirst(is_connecting_branch)))
            push!(is_connecting_branch, popat!(is_connecting_branch, findfirst(is_connecting_branch)))
        end

        if step_category == "process"
            reverse!(middle_transformer_branches)
            reverse!(is_input)
            reverse!(is_connecting_branch)
        end

        for (idx, middle_transformer_branch) in enumerate(middle_transformer_branches)
            # change direction after the first branch while accounting for possible changes 
            # from input to output or vica versa.
            #
            # Should do the following:
            # Input --> Input --> Input: false --> true --> true
            # Input --> Input --> Output: false --> true --> false
            # Input --> Output --> Input: false --> false --> true
            # Input --> Output --> Output: false --> false --> false
            # Output --> Input --> Input: true --> true --> true
            # Output --> Input --> Output: true --> true --> false
            # Output --> Output --> Input: true --> false --> true
            # Output --> Output --> Output: true --> false --> false
            if idx == 1
                reverse = !is_input[1]
            else
                reverse = is_input[idx]
            end

            if is_input[idx]
                connecting_component = [x for x in middle_transformer_branch if x.sys_function === EnergySystems.sf_transformer][end]
            else
                connecting_component = [x for x in middle_transformer_branch if x.sys_function === EnergySystems.sf_transformer][1]
            end

            current_checked_components = setdiff(checked_components, middle_transformer_branch)

            simulation_order, initial_nr, _ = add_transformer_steps(simulation_order, 
                                                                 initial_nr,
                                                                 middle_transformer_branch,
                                                                 parallel_branches,
                                                                 step_category,
                                                                 reverse=reverse,
                                                                 connecting_component=connecting_component,
                                                                 checked_components=current_checked_components)

            if length(middle_transformer_branches) > 2 && idx > 1
                simulation_order, initial_nr, _ = add_transformer_steps(simulation_order, 
                                                                     initial_nr,
                                                                     middle_transformer_branch,
                                                                     parallel_branches,
                                                                     step_category,
                                                                     reverse=!reverse,
                                                                     connecting_component=connecting_component,
                                                                     checked_components=current_checked_components)
            end

            # add potentials of middle_transformer in between of the single branches except for the last one
            # and process bevore the last branch, if in process category            
            if idx !== length(middle_transformer_branches)
                if step_category == "process" && idx == length(middle_transformer_branches) - 1
                    push!(simulation_order, [initial_nr, (first_middle_transformer.uac, EnergySystems.s_process)])
                else
                    push!(simulation_order, [initial_nr, (first_middle_transformer.uac, EnergySystems.s_potential)])
                end
                initial_nr -= 1
            end
        end
    else
        # no middle transformer is present in the current components
        # search for middle_busses: busses with more than one input and/or output 
        # interface with as least one transformer each.
        # returns the first/last middle bus that has no other middle bus in its inputs/outputs,
        # depending if reverse is true (last) or false/nothing (first). 
        # Branches should be in the correct ascending input/output order of the bus,
        # starting with the highest priority.
        # Considers only current components that are not in checked_components,
        # but also possible connected upstream transformers in an output and the other way around.
        first_middle_bus, 
            inface_branches_with_transformers,
            outface_branches_with_transformers = detect_middle_bus(current_components, reverse, checked_components)
 
        checked_components = add_component_to_checked_components(checked_components, inface_branches_with_transformers)
        checked_components = add_component_to_checked_components(checked_components, outface_branches_with_transformers)
        if (first_middle_bus !== nothing && !(first_middle_bus in checked_components)) 
            push!(checked_components, first_middle_bus)
        end

        if first_middle_bus !== nothing
            if reverse === nothing || reverse == false
                combined_branches_with_transformers = vcat(
                    [(collect(Iterators.reverse(inface)), true) for inface in inface_branches_with_transformers],
                    [(outface, false) for outface in outface_branches_with_transformers]
                )
            else
                combined_branches_with_transformers = vcat(
                    [(outface, false) for outface in outface_branches_with_transformers],
                    [(collect(Iterators.reverse(inface)), true) for inface in inface_branches_with_transformers]
                )
            end
            
            middle_bus_branches = [x[1] for x in combined_branches_with_transformers]
            is_input = [x[2] for x in combined_branches_with_transformers]
    
            # detect connecting branch and move it to the end of the calculation if reverse == true
            is_connecting_branch = fill(false, length(is_input))
            if connecting_component !== nothing
                for (idx, middle_bus_branch) in enumerate(middle_bus_branches)
                    for component in middle_bus_branch
                        if component == connecting_component
                            is_connecting_branch[idx] = true
                            # should be only one.
                            break
                        end
                    end
                end
            end
            if any(is_connecting_branch) && (reverse === nothing ? false : reverse)
                push!(middle_bus_branches, popat!(middle_bus_branches, findfirst(is_connecting_branch)))
                push!(is_input, popat!(is_input, findfirst(is_connecting_branch)))
                push!(is_connecting_branch, popat!(is_connecting_branch, findfirst(is_connecting_branch)))
            end

            # detect parallel branches and merge them together
            parallel_branch_idx = fill(0, length(is_input))
            for (idx_middle_bus_branch, middle_bus_branch) in enumerate(middle_bus_branches)
                for component in middle_bus_branch
                    for (idx_parallel_branch, parallel_branch) in enumerate(parallel_branches)
                        for branch in parallel_branch[2]
                            if component in branch
                                parallel_branch_idx[idx_middle_bus_branch] = idx_parallel_branch
                            end
                        end
                    end
                end
            end
            merged_branches_dict = Dict{Int, Vector{Any}}()
            is_input_dict = Dict()
            is_connecting_branch_dict = Dict()
            order_dict = Dict{Int, Int}()
            # Iterate over the parallel_branch_idx and corresponding middle_bus_branches
            for (i, idx) in enumerate(parallel_branch_idx)
                if idx == 0
                    # Use the index i as key for branches that should be copied as is
                    merged_branches_dict[-i] = middle_bus_branches[i]
                    is_input_dict[-i] = is_input[i]
                    is_connecting_branch_dict[-i] = is_connecting_branch[i]
                    order_dict[-i] = i
                else
                    # If the index is not 0, merge the branches
                    if haskey(merged_branches_dict, idx)
                        if is_input_dict[idx] !== is_input[i]
                            # should actually not happen
                            @warn "The order of operation may be wrong..."
                        end
                        is_input_dict[idx] = copy(is_input[i])
                        is_connecting_branch_dict[idx] = is_connecting_branch_dict[idx] || is_connecting_branch[i]
                        if reverse === nothing 
                            temp_reverse = !is_input_dict[idx]
                            if is_connecting_branch_dict[idx]
                                temp_reverse = !temp_reverse
                            end
                        else
                            temp_reverse = reverse
                        end
                        merged_branches_dict[idx] = iterate_chain(vcat(merged_branches_dict[idx],
                                                                       middle_bus_branches[i]),
                                                                       EnergySystems.sf_transformer,
                                                                       reverse=temp_reverse)
                    else
                        merged_branches_dict[idx] = copy(middle_bus_branches[i])
                        is_input_dict[idx] = copy(is_input[i])
                        is_connecting_branch_dict[idx] = copy(is_connecting_branch[i])
                        order_dict[idx] = i  # Record the first occurrence for ordering
                    end
                end
            end
            # Sort keys by their first occurrence
            sorted_keys = sort(collect(keys(order_dict)), by=k -> order_dict[k])
            # Collect the results based on the sorted keys
            middle_bus_branches = [merged_branches_dict[k] for k in sorted_keys]
            is_input = [is_input_dict[k] for k in sorted_keys]
            is_connecting_branch = [is_connecting_branch_dict[k] for k in sorted_keys]

            for (idx,branch) in enumerate(middle_bus_branches)
                middle_bus_branches[idx] = unique(middle_bus_branches[idx])
            end

            # iterate over middle bus branches
            for (idx, middle_bus_branch) in enumerate(middle_bus_branches)
                if isempty(middle_bus_branch)
                    continue
                end
                if reverse === nothing || !is_connecting_branch[idx]
                    if step_category == "potential"
                        current_reverse = !is_input[idx]
                        if is_connecting_branch[idx]
                            current_reverse = !current_reverse
                        end
                    end
                    if step_category == "process"
                        current_reverse = is_input[idx]
                    end
                else
                    current_reverse = reverse
                end

                if is_input[idx]
                    connecting_component = [x for x in middle_bus_branch if x.sys_function === EnergySystems.sf_transformer][end]
                else
                    connecting_component = [x for x in middle_bus_branch if x.sys_function === EnergySystems.sf_transformer][1]
                end

                current_checked_components = setdiff(checked_components, middle_bus_branch)

                simulation_order, initial_nr, _ = add_transformer_steps(simulation_order,
                                                                     initial_nr,
                                                                     middle_bus_branch,
                                                                     parallel_branches,
                                                                     step_category,
                                                                     reverse=current_reverse,
                                                                     connecting_component=connecting_component,
                                                                     checked_components=current_checked_components)
            end

        else
            # no middle busses found
            # write steps according default or predefined order if given
            if reverse === nothing
                reverse = step_category == "process"
            end
            branch_finished = false
            if reverse
                current_components_ordered = collect(Iterators.reverse(current_components))
            else
                current_components_ordered = current_components
            end
            for unit in current_components_ordered
                if unit.sys_function === EnergySystems.sf_transformer && !branch_finished
                    # check for parallels. Assuming that if unit is in parallels, all current units are 
                    # in this parallel! May change this later...
                    is_parallel = false
                    current_parallel_branches = []
                    for parallel_branches in parallel_branches
                        for branch in parallel_branches[2]
                            if unit in branch
                                is_parallel = true
                                current_parallel_branches = parallel_branches[2]
                                break
                            end
                        end
                    end

                    if is_parallel
                        nr_parallel_branches = length(current_parallel_branches)
                        for (branch_idx, current_branch) in enumerate(current_parallel_branches)
                            # TODO call add_transformer_steps() recursively here to handle middle busses within parallels?
                            #      may do this in future versions as this is not straightforward...
                            current_branch = reverse ? collect(Iterators.reverse(current_branch)) : current_branch
                            current_branch_transformers = [x for x in current_branch if x.sys_function === EnergySystems.sf_transformer]
                            for component in current_branch_transformers
                                if step_category == "potential"
                                    push!(simulation_order, [initial_nr, (component.uac, EnergySystems.s_potential)])
                                else
                                    push!(simulation_order, [initial_nr, (component.uac, EnergySystems.s_process)])
                                end
                                initial_nr -= 1
                            end
                            if branch_idx !== nr_parallel_branches && step_category == "potential"
                                for component_rev in Iterators.reverse(current_branch_transformers[1:end-1])
                                    push!(simulation_order, [initial_nr, (component_rev.uac, EnergySystems.s_potential)])
                                    initial_nr -= 1
                                end
                            end
                        end
                        branch_finished = true
                    else
                        if step_category == "potential"
                            push!(simulation_order, [initial_nr, (unit.uac, EnergySystems.s_potential)])
                        else
                            push!(simulation_order, [initial_nr, (unit.uac, EnergySystems.s_process)])
                        end
                        initial_nr -= 1
                    end
                end
            end
        end
    end

    return simulation_order, initial_nr, checked_components
end

"""
    detect_middle_bus(current_components, reverse, checked_components)

This function searches for middle_busses in the current_components.
Middle_busses are busses with more than two interfaces, each with at least one transformer.
Returns the first/last middle bus that has no other middle bus in its inputs/outputs, depending on whether reverse is 
true (last) or false/nothing (first). 
Branches are returned in ascending input/output order of the bus, starting with the highest priority.
Considers only current components that are not in checked_components, but also possibly connected upstream transformers
in an output and the other way around.
The components in each branch are sorted by their distance to the middle_bus, meaning for input branches, the order
is reversed and for output branches, the order of energy flow is kept.

# Arguments
- `current_components::Array{Grouping}`: The currently considered components.
- `reverse::Bool`: Indicates if the current_components should be considered in reverse or not.
- `checked_components::Array{Grouping}`: A vector holding all checked components. 

# Returns
- `first_middle_bus::Component`: The first middle bus as a component.
- `inface_branches_with_transformers::Array{Array{Components}}`: An array of arrays holding the transformers and busses 
                                                                 of each input branch of the first_middle_bus.
- `outface_branches_with_transformers::Array{Array{Components}}`: An array of arrays holding the transformers and busses 
                                                                  of each output branch of the first_middle_bus.
"""
function detect_middle_bus(current_components, reverse, checked_components)
    middle_busses = Any[]
    transformers_in_infaces = Any[]
    transformers_in_outfaces = Any[]
    transformers_in_infaces_arr = Vector{Vector{Any}}()
    transformers_in_outfaces_arr = Vector{Vector{Any}}()

    for component in current_components
        outfaces_with_transformers = Any[]
        infaces_with_transformers = Any[]
        transformer_in_output_interface = Any[]
        transformer_in_input_interface = Any[]
        if component.sys_function === EnergySystems.sf_bus
            component = component.proxy === nothing ? component : component.proxy
            for inface in values(component.input_interfaces)
                if has_grid_output(component, inface.source.uac)
                    continue # skip all interfaces with connection to a grid
                end
                if check_interface_for_transformer(inface, "input")
                    push!(infaces_with_transformers, inface)
                    other_infaces = [x for x in values(component.input_interfaces) if x !== inface]
                    inface_transformers = Any[]
                    add_non_recursive_indirect_inputs!(inface_transformers,
                                                       component,
                                                       other_infaces,
                                                       [EnergySystems.sf_transformer, EnergySystems.sf_bus],
                                                       "",
                                                       true,
                                                       true)
                    push!(transformer_in_input_interface, inface_transformers)
                end
            end
            for outface in values(component.output_interfaces)
                if has_grid_input(component, outface.target.uac)
                    continue # skip all interfaces with connection to a grid
                end
                if check_interface_for_transformer(outface, "output")
                    push!(outfaces_with_transformers, outface)
                    other_outfaces = [x for x in values(component.output_interfaces) if x !== outface]
                    outface_transformers = Any[]
                    add_non_recursive_indirect_outputs!(outface_transformers,
                                                        component,
                                                        other_outfaces,
                                                        [EnergySystems.sf_transformer, EnergySystems.sf_bus],
                                                        "",
                                                        true,
                                                        true)
                    push!(transformer_in_output_interface, outface_transformers)
                end
            end 
        end
        if length(outfaces_with_transformers) > 1 || length(infaces_with_transformers) > 1
            push!(middle_busses, component)
            push!(transformers_in_infaces, transformer_in_input_interface)
            push!(transformers_in_outfaces, transformer_in_output_interface)
        end
    end

    # detect which middle bus is the first or the last, depending on reverse
    # default starts with the first bus
    has_middle_bus_in_interfaces = Any[]
    for (idx, middle_bus) in enumerate(middle_busses)
        bus_chain = Set()
        has_middle_bus_in_interfaces_temp = false

        if reverse === nothing || reverse == false
            add_non_recursive_indirect_inputs!(bus_chain, middle_bus, [], [EnergySystems.sf_bus], "", true, true)
        else
            add_non_recursive_indirect_outputs!(bus_chain, middle_bus, [], [EnergySystems.sf_bus], "", true, true)
        end
        for other_middle_bus in middle_busses
            if middle_bus == other_middle_bus
                continue
            elseif other_middle_bus in bus_chain
                has_middle_bus_in_interfaces_temp = true
            end
        end
        push!(has_middle_bus_in_interfaces, has_middle_bus_in_interfaces_temp)

        # deleat middle_bus from transformer sets and convert Set to array
        temp = []
        for interface in transformers_in_infaces[idx]
            temp_transformers = []
            for component in interface
                if component !== middle_bus && !(component in checked_components)
                    push!(temp_transformers, component)
                end
            end
            push!(temp, temp_transformers)
        end
        push!(transformers_in_infaces_arr, temp)

        temp = []
        for interface in transformers_in_outfaces[idx]
            temp_transformers = []
            for component in interface
                if component !== middle_bus && !(component in checked_components)
                    push!(temp_transformers, component)
                end
            end
            push!(temp, temp_transformers)
        end
        push!(transformers_in_outfaces_arr, temp)
    
    end

    false_indices = findall(x -> !x, has_middle_bus_in_interfaces)

    if length(false_indices) > 0
        # return the first one
        return middle_busses[false_indices[1]], 
               transformers_in_infaces_arr[false_indices[1]],
               transformers_in_outfaces_arr[false_indices[1]]
    else
        return nothing, nothing, nothing
    end
end

"""
    detect_first_middle_transformer(current_components, checked_components)

This function searches for a middle_transformer in the current_components.
Middle_transformers are transformers with either more than one input interface or more than one output interface, each 
with at least one transformer. Transformers with two input AND two output interfaces with one or more transformers each 
are not considered (currently, no component is implemented for which this could be the case).
Returns the first middle transformer that has no other middle transformer in its inputs.
Branches are returned in ascending input/output order of the transformer, starting with the highest priority.
Considers only current_components that are not in checked_components, but also possibly connected upstream transformers
in an output and the other way around.
The components in each branch are sorted by their distance to the middle_bus, meaning for input branches, the order
is reversed and for output branches, the order of energy flow is kept.

# Arguments
- `current_components::Array{Grouping}`: The currently considered components.
- `checked_components::Array{Grouping}`: A vector holding all checked components. 

# Returns
- `first_middle_transformer::Component`: The first middle transformer as a component.
- `inface_branches_with_transformers::Array{Array{Components}}`: An array of arrays holding the transformers and busses 
                                                                 of each input branch of the first_middle_transformer.
- `outface_branches_with_transformers::Array{Array{Components}}`: An array of arrays holding the transformers and busses
                                                                  of each output branch of the first_middle_transformer.
"""
function detect_first_middle_transformer(current_components, checked_components)
    middle_transformers = Any[]
    transformers_in_infaces = Any[]
    transformers_in_outfaces = Any[]
    transformers_in_infaces_arr = Vector{Vector{Any}}()
    transformers_in_outfaces_arr = Vector{Vector{Any}}()

    for component in current_components
        outfaces_with_transformers = Any[]
        infaces_with_transformers = Any[]
        transformer_in_output_interface = Any[]
        transformer_in_input_interface = Any[]
        if component.sys_function === EnergySystems.sf_transformer
            for inface in values(component.input_interfaces)
                if check_interface_for_transformer(inface, "input")
                    other_infaces = [x for x in values(component.input_interfaces) if x !== inface]
                    inface_transformers = Any[]
                    add_non_recursive_indirect_inputs!(inface_transformers,
                                                       component,
                                                       other_infaces,
                                                       [EnergySystems.sf_transformer, EnergySystems.sf_bus],
                                                       "", 
                                                       true, 
                                                       true)
                    if length([x for x in inface_transformers if x.sys_function === EnergySystems.sf_transformer]) > 1
                        push!(transformer_in_input_interface, inface_transformers)
                        push!(infaces_with_transformers, inface)
                    end
                end
            end
            for outface in values(component.output_interfaces)
                if check_interface_for_transformer(outface, "output")  
                    other_outfaces = [x for x in values(component.output_interfaces) if x !== outface]
                    outface_transformers = Any[]
                    add_non_recursive_indirect_outputs!(outface_transformers,
                                                        component,
                                                        other_outfaces,
                                                        [EnergySystems.sf_transformer, EnergySystems.sf_bus],
                                                        "",
                                                        true,
                                                        true)
                    if length([x for x in outface_transformers if x.sys_function === EnergySystems.sf_transformer]) > 1
                        push!(transformer_in_output_interface, outface_transformers)
                        push!(outfaces_with_transformers, outface)
                    end
                end
            end 
        end
        if length(outfaces_with_transformers) > 1 || length(infaces_with_transformers) > 1
            push!(middle_transformers, component)
            push!(transformers_in_infaces, transformer_in_input_interface)
            push!(transformers_in_outfaces, transformer_in_output_interface)
        end
    end

    has_middle_transformer_in_inputs = Any[]
    for (idx, middle_transformer) in enumerate(middle_transformers)
        transformer_input_chain = Set()
        has_middle_transformer_in_inputs_temp = false

        add_non_recursive_indirect_inputs!(transformer_input_chain, 
                                           middle_transformer,
                                           [],
                                           [EnergySystems.sf_transformer],
                                           "",
                                           true,
                                           true)
        for other_middle_transformer in middle_transformers
            if middle_transformer == other_middle_transformer
                continue
            elseif other_middle_transformer in transformer_input_chain
                has_middle_transformer_in_inputs_temp = true
            end
        end
        push!(has_middle_transformer_in_inputs, has_middle_transformer_in_inputs_temp)

        # deleat middle_transformer from transformer sets and convert Set to array
        temp = []
        for interface in transformers_in_infaces[idx]
            temp_transformers = []
            for transformer in interface
                if transformer !== middle_transformer && !((transformer in checked_components))
                    push!(temp_transformers, transformer)
                end
            end
            push!(temp, temp_transformers)
        end
        push!(transformers_in_infaces_arr, temp)

        temp = []
        for interface in transformers_in_outfaces[idx]
            temp_transformers = []
            for transformer in interface
                if transformer !== middle_transformer && !((transformer in checked_components))
                    push!(temp_transformers, transformer)
                end
            end
            push!(temp, temp_transformers)
        end
        push!(transformers_in_outfaces_arr, temp)
    
    end

    false_indices = findall(x -> !x, has_middle_transformer_in_inputs)

    if length(false_indices) > 0
        # return the first one
        return middle_transformers[false_indices[1]], 
               transformers_in_infaces_arr[false_indices[1]],
               transformers_in_outfaces_arr[false_indices[1]]
    else
        return nothing, nothing, nothing
    end
end

"""
    find_parallels(components)

This function searches for parallel branches in the components of the energy system.
Parallel branches are defined as at least two connections between two components, each containing at least one transformer.
Returns a dictionary with all parallel branches of the current energy system.
The components in each branch are sorted by the direction of the energy flow, starting with the start_point component.
The branches are sorted by the output order of the start_point component.
For busses with proxy busses, only the proxy busses are taken into account.

# Arguments
- `components::Array{Grouping}`: All components of the current energy system

# Returns
- `final_paths::Dict{Array{Array}}`: A dictionary containing parallel branches. The key is (start_point_uac, end_point_uac).
                                     Each key contains a vector of vectors where each one holds one path from start_point
                                     to the end_point. A path consists of the components, including the start_point and 
                                     end_point components.
"""
function find_parallels(components)
    function find_paths(unit, path, all_paths, old_uac)
        push!(path, unit)
        for outface in values(unit.output_interfaces)
            if (  outface === nothing 
               || outface.target in path 
               || outface.target.sys_function === EnergySystems.sf_storage
               || (unit.sys_function === EnergySystems.sf_bus 
                   && has_grid_input(unit, outface.target.uac)
                   && !(old_uac !== "" && !has_grid_output(unit, old_uac)) )
               || (old_uac !== "" && !connection_allowed(unit, old_uac, outface.target.uac))
               || (outface.target === EnergySystems.sf_transformer && !(outface.target in components) )
            )
                continue
            end
            new_path = copy(path)
            if outface.target.sys_function === EnergySystems.sf_bus && outface.target.proxy !== nothing
                find_paths(outface.target.proxy, new_path, all_paths, outface.source.uac)
            else
                find_paths(outface.target, new_path, all_paths, outface.source.uac)
            end
        end 
        push!(all_paths, copy(path))
        pop!(path)
    end
    
    all_paths = Vector{Vector}()
    paths_dict = Dict{Tuple, Vector{Vector}}()

    # find all possible paths between all components in current energy system
    # Paths are split by storages and by grid inputs and outputs. 
    # Connection matrixes of busses are taken into account.
    for unit in components
        if !(unit.sys_function === EnergySystems.sf_bus && unit.proxy !== nothing)
            find_paths(unit, [], all_paths, "")
        end
    end

    for path in all_paths
        start_unit = path[1]
        end_unit = path[end]
        # filter paths for...
        if (
            # starting with either a bus or transformer
            start_unit.sys_function in [EnergySystems.sf_bus, EnergySystems.sf_transformer]
            # ending with either a bus or transformer
            && end_unit.sys_function in [EnergySystems.sf_bus, EnergySystems.sf_transformer]
            # has a length > 2
            && length(path) > 2
            # path contains at least one transformer in between
            && EnergySystems.sf_transformer in [unit.sys_function for unit in path[2:end-1]]
        )
            # remove all non-transformers and non-busses from path
            path_reduced = Any[path[1]]
            for component in path[2:end-1]
                if component.sys_function in [EnergySystems.sf_transformer, EnergySystems.sf_bus]
                    push!(path_reduced, component)
                end
            end
            push!(path_reduced, path[end])

             # save path to paths_dict if it's not already included
             key = (start_unit.uac, end_unit.uac)
             if !haskey(paths_dict, key)
                 paths_dict[key] = Vector{Vector}()
             end
             if !any(path_reduced == p for p in paths_dict[key])
                 push!(paths_dict[key], path_reduced)
             end
        end
    end

    # filter for parallel paths that contain more than one variant from A to B
    parallel_paths = Dict(collect(filter(paths -> length(paths[2]) > 1, paths_dict)))

    # remove elements that are all the same for each path in each parallel path
    function values_are_identical(vectors, index)
        first_name = vectors[1][index]
        return all(v -> v[index] == first_name, vectors)
    end

    function last_are_identical(vectors)
        last_name = vectors[1][end]
        return all(v -> v[end] == last_name, vectors)
    end

    function second_last_are_identical(vectors)
        last_name = vectors[1][end-1]
        return all(v -> v[end-1] == last_name, vectors)
    end

    for (keys, paths) in parallel_paths
        while values_are_identical(paths, 1) && values_are_identical(paths, 2)
            for v in paths
                popfirst!(v)
            end
        end
        
        while last_are_identical(paths) && second_last_are_identical(paths)
            for v in paths
                pop!(v)
            end
        end
    end

    # build new dict with the updated unique paths
    final_paths = Dict()
    for paths in values(parallel_paths)
        for path in paths
            key = (path[1].uac, path[end].uac)
            if !haskey(final_paths, key)
                final_paths[key] = Vector{Vector}()
            end
            if !any(path == p for p in final_paths[key])
                push!(final_paths[key], path)
            end
        end
    end

    return final_paths
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
    add_non_recursive_indirect_outputs!(node_set, unit, checked_interfaces, sys_function)

Add connected units of the same system function to the node set in a non-recursive manner.
Here, "connected" does not necessarily mean that the components have to be connected directy 
to each other, they can also be connected via busses or other components within the directed graph.
This function only searches in the outputs of each component, starting with "unit".
Note: Storages are defined to interrupt chains!

# Arguments
-`node_set`: A globally defined Set() containing an (empty) chain of interconnected components (units)
-`unit`: The unit to start the search with
-`checked_interfaces::Array{SystemInterface}`: An array holding the already checked system 
                                               interfaces. Should be [] when calling externally.
-`sys_function`: The system function of the chain that should be found, as Vector to find multiple sys_function
-`last_unit_uac::String`: String to pass the uac of the last checked unit to the recursive call.
                          Should be an empty string ("") at the first external call.
"""
function add_non_recursive_indirect_outputs!(node_set,
                                             unit,
                                             checked_interfaces,
                                             sys_function,
                                             last_unit_uac,
                                             skip_storages=true,
                                             interrupt_at_grids=false)
    if unit.sys_function in sys_function && !(unit in node_set)
        push!(node_set, unit)
    end

    for outface in values(unit.output_interfaces)
        if outface !== nothing
            if outface in checked_interfaces
                continue
            elseif skip_storages && outface.target.sys_function === EnergySystems.sf_storage
                continue
            elseif interrupt_at_grids && unit.sys_function === EnergySystems.sf_bus && has_grid_input(unit, outface.target.uac)
                continue
            elseif last_unit_uac == "" || connection_allowed(unit, last_unit_uac, outface.target.uac)
                push!(checked_interfaces, outface)
                add_non_recursive_indirect_outputs!(node_set,
                                                    outface.target,
                                                    checked_interfaces,
                                                    sys_function,
                                                    unit.uac,
                                                    skip_storages,
                                                    interrupt_at_grids)
            end
        end
    end
end

"""
    add_non_recursive_indirect_inputs!(node_set, unit, checked_interfaces, sys_function)

Add connected units of the same system function to the node set in a non-recursive manner.
Here, "connected" does not necessarily mean that the components have to be connected directy 
to each other, they can also be connected via busses or other components within the directed graph.
This function only searches in the inputs of each component, starting with "unit".
Note: Storages are defined to interrupt chains!

# Arguments
-`node_set`: A globally defined Set() containing an (empty) chain of interconnected components (units)
-`unit`: The unit to start the search with
-`checked_interfaces::Array{SystemInterface}`: An array holding the already checked system 
                                               interfaces. Should be [] when calling externally.
-`sys_function`: The system function of the chain that should be found, as Vector to find multiple sys_function
-`last_unit_uac::String`: String to pass the uac of the last checked unit to the recursive call.
                          Should be an empty string ("") at the first external call.
"""
function add_non_recursive_indirect_inputs!(node_set,
                                            unit,
                                            checked_interfaces,
                                            sys_function,
                                            last_unit_uac,
                                            skip_storages=true,
                                            interrupt_at_grids=false)
    if unit.sys_function in sys_function && !(unit in node_set)
        push!(node_set, unit)
    end

    for inface in values(unit.input_interfaces)
        if inface !== nothing
            if inface in checked_interfaces
                continue
            elseif skip_storages && inface.source.sys_function === EnergySystems.sf_storage
                continue
            elseif interrupt_at_grids && unit.sys_function === EnergySystems.sf_bus && has_grid_output(unit, inface.source.uac)
                continue
            elseif last_unit_uac == "" || connection_allowed(unit, inface.source.uac, last_unit_uac)
                push!(checked_interfaces, inface)
                add_non_recursive_indirect_inputs!(node_set,
                                                   inface.source,
                                                   checked_interfaces,
                                                   sys_function,
                                                   unit.uac,
                                                   skip_storages,
                                                   interrupt_at_grids)
            end
        end
    end
end


"""
    find_chains(components, sys_function; direct_connection_only=true)

Find all chains of the given system function in the given collection of components.

A chain is a subgraph of the graph spanned by all connections of the given components,
which is a directed graph. The subgraph is defined by all connected components of the given
system function.
The optional flag direct_connection_only specifies if the chain should also be found
across busses and other component types. If the flag is set to false, the unique chains
that interconnect the components of sys_function are found. Components are only defined as
interconnected if they can be seen EITHER in the outputs or the inputs, but not via mixed 
input/output paths! 

Chains found with direct_connection_only=false are filtered to find unique chains. They are
merged to represent the longest possible chain and to avoid doublings or multible subset of
a bigger chain. Completely independent chains are returned as array of chains.

# Arguments
-`components`: All components of the current energy system
-`sys_function`: The system function for that the chains should be determined in components
-`direct_connection_only::Bool}`: Flag if direct (true) or indirect chains across other
                                  components should be determined (false)

# Returns
-`chain::Array{Set()}`: An array holding an unique collection of chains, each as Set() type
"""
function find_chains(components, sys_function; direct_connection_only=true)::Vector{Set{Component}}
    function merge_chains(original_chains)
        chains_unique = []
        for chain in original_chains 
            if chains_unique == []
                push!(chains_unique, chain)
            else
                chain_written = false
                for chain_unique in chains_unique
                    chain_uacs = [entry.uac for entry in chain]
                    chain_unique_uacs = [entry.uac for entry in chain_unique]

                    length_diff_chain = length(setdiff(chain_uacs, chain_unique_uacs))
                    elements_in_chain_but_not_in_unique = (length_diff_chain > 0 && length_diff_chain < length(chain_uacs) ) ? true : false

                    length_diff_unique = length(setdiff(chain_unique_uacs, chain_uacs))
                    elements_in_unique_but_not_in_chain = (length_diff_unique > 0 && length_diff_unique < length(chain_unique_uacs) ) ? true : false

                    if elements_in_chain_but_not_in_unique && elements_in_unique_but_not_in_chain
                        # chain is subset of chain_unique and vica versa 
                        # --> replace with union of chain and chains_unique
                        new_chain = union(chain, chain_unique)
                        replace!(chains_unique, chain_unique => new_chain)
                        chain_written = true
                        break
                    elseif !elements_in_chain_but_not_in_unique && elements_in_unique_but_not_in_chain
                        # chain is complete subset of chain_unique --> skip
                        chain_written = true
                        break
                    elseif elements_in_chain_but_not_in_unique && !elements_in_unique_but_not_in_chain
                        # chain_unique is complete subset of chain --> replace with chain
                        replace!(chains_unique, chain_unique => chain)
                        chain_written = true
                        break
                    elseif Set(chain_uacs) == Set(chain_unique_uacs)
                        # they are the same
                        chain_written = true
                        break
                    end

                end
                if !chain_written
                    # chain and all chains in chain_unique are completely different
                    push!(chains_unique, chain)
                end
            end
        end
        return chains_unique
    end

    chains = []

    for unit in components
        if unit.sys_function !== sys_function
            continue
        end

        if direct_connection_only
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
        else
            chain = Set()
            add_non_recursive_indirect_outputs!(chain, unit, [], [sys_function], "")
            add_non_recursive_indirect_inputs!(chain, unit, [], [sys_function], "")
            push!(chains, chain)
        end
    end

    # find unique chains for chains that has been found with direct_connection_only = false
    if !direct_connection_only
        check_chain_merging = function(merged_chains)
            # check if two chains are connected later through a third chain which is 
            # not recognized by merge_chains()
            merged_chains_again = merge_chains(merged_chains)
            if length(merged_chains_again) != length(merged_chains)
                merged_chains_again = check_chain_merging(merged_chains_again)
            end
            return merged_chains_again
        end

        chains = check_chain_merging(merge_chains(chains))
    end 
    return chains
end

"""
    distance_to_sink(node, sys_function, checked_interfaces)

Calculate the distance of the given node to the sinks of the chain.

A sink is defined as a node with no successors of the same system function. For the sinks
this distance is 0. For all other nodes it is the maximum over the distances of its
successors plus one. Only allowed connections are considered specified by the eneryg matrix 
of busses.
The parameter "checked_interfaces" is used to avoid loops when recursively calling the 
function. When calling distance_to_sink() from outside, the parameter checked_interfaces 
should be set as empty array ([]).

# Arguments
-`node`: A specific component (unit) for which the distance to the sink should be determined
-`sys_function`: The system function of the components in the chain
-`checked_interfaces::Array{SystemInterface}`: An array holding the already checked system 
                                               interfaces. Should be [] when calling externally.
-`last_unit_uac::String`: string of the uac of the last unit. Should be "" when calling externally.                                               

# Returns
The maximum distance to the farest sink as Int.
"""
function distance_to_sink(node, sys_function, checked_interfaces, last_unit_uac)
    is_leaf = function(current_node, checked_interfaces_leaf; is_leafe_result=true, last_uac="")
        for outface in values(current_node.output_interfaces)
            if outface !== nothing
                if outface in checked_interfaces_leaf || outface.target == node
                    continue
                elseif last_uac == "" || connection_allowed(current_node, last_uac, outface.target.uac)
                    push!(checked_interfaces_leaf, outface)
                    if outface.target.sys_function === sys_function
                        return false
                    else
                        is_leafe_result = is_leafe_result && 
                                          is_leaf(outface.target,
                                                  checked_interfaces_leaf,
                                                  is_leafe_result=is_leafe_result,
                                                  last_uac=current_node.uac)
                    end
                end
            end
        end
        return is_leafe_result
    end

    if is_leaf(node, [])
        return 0
    else
        max_distance = 0
        for outface in values(node.output_interfaces)
            if outface !== nothing
                if outface in checked_interfaces
                    continue
                elseif last_unit_uac == "" || connection_allowed(node, last_unit_uac, outface.target.uac)
                    push!(checked_interfaces, outface)
                    if outface.target.sys_function === sys_function || !is_leaf(outface.target, [])
                        distance = distance_to_sink(outface.target, sys_function, checked_interfaces, node.uac)
                        max_distance = distance > max_distance ? distance : max_distance
                    end
                end
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
        push!(distances, (distance_to_sink(node, sys_function, [], ""), node))
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
    simulation_order = remove_double_transformer_process_steps(simulation_order, components)
    # TODO: detect and remove unnecessary potentials in future versions?

    # Note that the input order at a bus has a higher priority compared to the output order! 
    # If there are contradictions, the input order applies. Currently, there is no warning 
    # message if this leads to missalignement with the output order.
    reorder_transformer_for_output_priorities(simulation_order, components, components_by_function)
    reorder_for_input_priorities(simulation_order, components, components_by_function)
    reorder_distribution_of_busses(simulation_order, components, components_by_function)
    reorder_storage_loading(simulation_order, components, components_by_function)

    fn_first = function (entry)
        return entry[1]
    end

    step_order = [(u[2][1], u[2][2]) for u in sort(simulation_order, by=fn_first, rev=true)]
    
    # remove transformers potential step if process is directly consecutive to a potential step
    while contains_double_potental_produce(step_order)
        step_order = remove_double_potental_produce(step_order)
    end

    return step_order
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
This does not apply, if there is a grid output from the bus and the connection is allowed from 
the inputs. This does not take into account any attributs like temperatures that deny the 
energy flow. Up to now, grids do not have any attributes, so that doesn't matter.
But, the energy flow matrix of the bus is taken into account.

# Arguments
-`simulation_order`: A global parameter holding the simulation order
-`components`: All components of the current energy system
-`components_by_function`: The mapping of component functions

"""
function reorder_for_input_priorities(simulation_order, components, components_by_function)
    # for every bus except for those with proxies...
    for bus in values(components_by_function[3])
        if bus.proxy !== nothing continue end
        # ...for each component in the bus' input priority...
        for own_idx = 1:length(bus.connectivity.input_order)
            # ...make sure every component following after...
            for other_idx = own_idx+1:length(bus.connectivity.input_order)
                #(...if there is a connected grid with an allowed connection to both components...)
                #(...then the order doesn't matter as the components can deliver their energy anyway)
                if has_grid_output(bus, bus.connectivity.input_order[own_idx]) &&
                   has_grid_output(bus, bus.connectivity.input_order[other_idx]) 
                   continue 
                end
                # ...is of a lower priority in process
                place_one_lower!(
                    simulation_order,
                    (bus.connectivity.input_order[own_idx], EnergySystems.s_process),
                    (bus.connectivity.input_order[other_idx], EnergySystems.s_process)
                )
            end
        end
    end
end

function has_grid_input(bus, output_interface_uac)
    for inface in values(bus.input_interfaces)
        if inface !== nothing && nameof(typeof(inface.source)) == :GridConnection
            input_idx = bus.balance_table_inputs[inface.source.uac].input_index
            output_idx = bus.balance_table_outputs[output_interface_uac].output_index
            if (bus.connectivity.energy_flow === nothing ||
                bus.connectivity.energy_flow[input_idx][output_idx])
                return true
            end
        end
    end
    return false
end

function has_grid_output(bus, input_interface_uac)
    for outface in values(bus.output_interfaces)
        if outface !== nothing && nameof(typeof(outface.target)) == :GridConnection
            input_idx = bus.balance_table_inputs[input_interface_uac].input_index
            output_idx = bus.balance_table_outputs[outface.target.uac].output_index
            if (bus.connectivity.energy_flow === nothing ||
                bus.connectivity.energy_flow[input_idx][output_idx])
                return true
            end
        end
    end
    return false
end

"""
    reorder_transformer_for_output_priorities(simulation_order, components, components_by_function)

Reorder transformers connected to a bus so they match the output priority defined on that bus.
This does not apply, if there is a grid input to the bus and the connection is allowed to 
the output transformers. This does not take into account any attributs like temperatures that deny the 
energy flow. Up to now, grids do not have any attributes, so that doesn't matter.
But, the energy flow matrix of the bus is taken into account.

# Arguments
-`simulation_order`: A global parameter holding the simulation order
-`components`: All components of the current energy system
-`components_by_function`: The mapping of component functions
"""
function reorder_transformer_for_output_priorities(simulation_order, components, components_by_function)
    # for every bus except for those with proxies...
    for bus in values(components_by_function[3])
        if bus.proxy !== nothing continue end
        # ...for each transformer in the bus' output priority...
        for own_idx = 1:length(bus.connectivity.output_order)
            own_uac = bus.connectivity.output_order[own_idx]
            if bus.balance_table_outputs[own_uac].target.sys_function !== EnergySystems.sf_transformer continue end
            # ...make sure every transformer following after...
            for other_idx = own_idx+1:length(bus.connectivity.output_order)
                other_uac = bus.connectivity.output_order[other_idx]
                if bus.balance_table_outputs[other_uac].target.sys_function !== EnergySystems.sf_transformer continue end
                #(...if there is a connected grid with an allowed connection to both transformers...)
                #(...then the order doesn't matter as the transformers can get their required energy anyway)
                if has_grid_input(bus, own_uac) &&
                   has_grid_input(bus, other_uac) 
                   continue 
                end
                # ...is of a lower priority in process
                place_one_lower!(
                    simulation_order,
                    (own_uac, EnergySystems.s_process),
                    (other_uac, EnergySystems.s_process)
                )
            end
        end
    end
end

"""
    reorder_distribution_of_busses(simulation_order, components, components_by_function)

Reorder the distribution of busses so that for chains of busses larger than one the proxy
bus for that chain is distributed first.
"""
function reorder_distribution_of_busses(simulation_order, components, components_by_function)
    non_proxy_busses = [c for c in components_by_function[3] if !startswith(c.uac, "Proxy")]
    for bus_chain in find_chains(non_proxy_busses, EnergySystems.sf_bus)
        # for a single bus we don't need to do anything. for chains we use the proxy
        if length(bus_chain) <= 1
            continue
        end
        bus = first(bus_chain)
        proxy_bus = bus.proxy !== nothing ? bus.proxy : bus

        # for principal busses in the chain, place their distribute step after that of the
        # proxy bus. the order of the principals within themselves doesn't matter, but it
        # should be deterministic for testing and debugging purposes, hence the sorting
        for principal in sort(collect(bus_chain), by=c->c.uac)
            place_one_lower!(
                simulation_order,
                (proxy_bus.uac, EnergySystems.s_distribute),
                (principal.uac, EnergySystems.s_distribute)
            )
        end
    end
end

"""
    reorder_storage_loading(simulation_order, components, components_by_function)

Reorder components such the loading (and unloading) of storages follows the priorities on
busses, including communication across connected busses.
"""
function reorder_storage_loading(simulation_order, components, components_by_function)
    non_proxy_busses = [c for c in components_by_function[3] if !startswith(c.uac, "Proxy")]
    for bus_chain in find_chains(non_proxy_busses, EnergySystems.sf_bus)
        # for a single bus we use the bus itself, for chains we use the proxy
        bus = first(bus_chain)
        if length(bus_chain) > 1 && bus.proxy !== nothing
            bus = bus.proxy
        end

        # identify first input and output storage
        output_storages = [
            i.target for i in bus.output_interfaces
            if i.target.sys_function == EnergySystems.sf_storage
        ]
        first_output_storage = length(output_storages) < 1 ? nothing : output_storages[1]

        input_storages = [
            i.source for i in bus.input_interfaces
            if i.source.sys_function == EnergySystems.sf_storage
        ]
        first_input_storage = length(input_storages) < 1 ? nothing : input_storages[1]
        last_input_storage = length(input_storages) < 1 ? nothing : input_storages[end]

        # for the storage with the highest output priority, place its load step after the
        # process step of the input storage with highest priority, so that the subsequent
        # insertions maintain the order of process steps first, then load steps
        if first_input_storage !== nothing && first_output_storage !== nothing
            place_one_lower!(
                simulation_order,
                (last_input_storage.uac, EnergySystems.s_process),
                (first_output_storage.uac, EnergySystems.s_load),
            )
        end

        # for the storage with the highest output priority, place its load step after the
        # process step of the transformer with the lowest input priority, so that the load
        # of storages happen after the last input transformer had its process.
        # storages act like a bounded source, they have already written a valid max_energy
        # in their control.
        input_transformer = [
            i.source for i in bus.input_interfaces
            if i.source.sys_function == EnergySystems.sf_transformer
        ]
        last_input_transformer = length(input_transformer) < 1 ? nothing : input_transformer[end]
        if first_input_storage !== nothing && last_input_transformer !== nothing
            place_one_lower!(
                simulation_order,
                (last_input_transformer.uac, EnergySystems.s_process),
                (first_output_storage.uac, EnergySystems.s_load),
            )
        end

        # reorder load for storages according to output priorities. this works by
        # continuously placing lower than the first storage in reverse order, which results
        # in the correct ordering
        if first_output_storage !== nothing
            for output_storage in reverse(output_storages)
                if output_storage.uac == first_output_storage.uac continue end
                place_one_lower!(
                    simulation_order,
                    (first_output_storage.uac, EnergySystems.s_load),
                    (output_storage.uac, EnergySystems.s_load),
                    force=true
                )
            end
        end

        # same as above, but for inputs and the process step
        # this is done in reorder_for_input_priorities()
    end
end

"""
    contains_double_potental_produce(step_order)

Checks the simulation step order for directly consecutive transformers potential and produce step and 
removes the potential step.
"""
function contains_double_potental_produce(step_order)
    last_unit = ""
    last_step = ""
    for entry in step_order
        if entry[1] == last_unit && entry[2] == last_step
            return true
        else
            last_unit = entry[1]
            last_step = entry[2]
        end
    end
    return false
end

"""
    remove_double_transformer_process_steps(simulation_order)

Checks the simulation order for multiple process steps of one transformer.
If there are multiple process steps for one tranformer, only the first one will be kept.
"""
function remove_double_transformer_process_steps(simulation_order, components)
    detected_transformer_processes = []
    for (idx,(_, step)) in enumerate(simulation_order)
        if step[2] == EnergySystems.s_process && components[step[1]].sys_function === EnergySystems.sf_transformer
            if step in detected_transformer_processes
                deleteat!(simulation_order, idx)
            else
                push!(detected_transformer_processes, step)
            end
        end
    end
    return simulation_order
end

"""
    remove_double_potental_produce(step_order)

Checks the simulation step order for directly consecutive transformers potential and produce step and 
removes the potential step.
"""
function remove_double_potental_produce(step_order)
    last_unit = ""
    last_step = ""
    to_remove = Int[]
    for (idx, entry) in enumerate(step_order)
        unit = entry[1]
        step = entry[2]
        
        if unit == last_unit 
            if (step == EnergySystems.s_process
                && last_step == EnergySystems.s_potential)
                push!(to_remove, idx-1)
            elseif (step == EnergySystems.s_potential
                    && last_step == EnergySystems.s_process)
                push!(to_remove, idx)
            elseif (step == last_step)
                push!(to_remove, idx)
            end
        end
        last_unit = unit
        last_step = step
    end

    for index in sort(to_remove, rev=true)
        deleteat!(step_order, index)
    end

    return step_order
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

"""
    check_interface_for_transformer(interface::SystemInterface, type::String)

Checks a given interface if there is a transformer in the following or previous chain.
The function search through the whole energy system.

# Arguments
-`interface::SystemInterface`: The interface that should be checked
-`type::String`: Can be either "input" or "output". Defines if the "interface" should 
                 be handeled as an input or an output interface.

# Returns
Returns either "true" if a transformer was found in the following or previous chain or 
"false" if no transformer has been found.
"""
function check_interface_for_transformer(interface, type)
    if type =="input"
        has_transformer_in_input = function(current_node, last_node_uac, checked_interfaces)
            if current_node.sys_function === EnergySystems.sf_transformer
                return true
            end
            for inface in values(current_node.input_interfaces)
                if inface !== nothing
                    if inface in checked_interfaces || inface.source == current_node
                        continue
                    elseif !connection_allowed(current_node, inface.source.uac, last_node_uac)
                        continue
                    else
                        push!(checked_interfaces, inface)
                        if inface.source.sys_function === EnergySystems.sf_transformer
                            return true
                        else
                            if has_transformer_in_input(inface.source, current_node.uac, checked_interfaces)
                                return true
                            end
                        end
                    end
                end
            end
            return false
        end

        return has_transformer_in_input(interface.source, interface.target.uac, [])

    elseif type == "output"
        has_transformer_in_output = function(current_node, last_node_uac, checked_interfaces)
            if current_node.sys_function === EnergySystems.sf_transformer
                return true
            end
            for outface in values(current_node.output_interfaces)
                if outface !== nothing
                    if outface in checked_interfaces || outface.target == current_node
                        continue
                    elseif !connection_allowed(current_node, last_node_uac, outface.target.uac)
                        continue
                    else
                        push!(checked_interfaces, outface)
                        if outface.target.sys_function === EnergySystems.sf_transformer
                            return true
                        else
                            if has_transformer_in_output(outface.target, current_node.uac, checked_interfaces)
                                return true
                            end
                        end
                    end
                end
            end
            return false
        end
        
        return has_transformer_in_output(interface.target, interface.source.uac, [])
    
    else
        @error "The function check_interface_for_transformer() was not able to detect if it is an output or an input interface. Check the function call."
        exit()
    end
end

"""
    connection_allowed(component::Component, input_uac::String, output_uac::String)

Checks a given connection defined by `input_uac` and `output_uac` is allowed. The `component`
is the commponent between the two other components. If `component` is not a bus, this function
will return true. If component is a bus, then the connection matrix of the bus is checked 
to determine if the connection from input_uac to output_uac is allowed. 

# Arguments
-`component::Component`: A component beween input_uac and output_uac
-`input_uac::String`: The uac of the input component
-`output_uac::String`: The uac of the output component


# Returns
Returns either "true" if `component` is a non-bus or if the connection is allowed, or "false" 
if `component` is a bus and the energy flow is denied by the energy flow matrix.
"""
function connection_allowed(component::Component, input_uac::String, output_uac::String)
    if component.sys_function === EnergySystems.sf_bus
        input_idx = component.balance_table_inputs[input_uac].input_index
        output_idx = component.balance_table_outputs[output_uac].output_index
        if (component.connectivity.energy_flow === nothing ||
            component.connectivity.energy_flow[input_idx][output_idx])
            return true
        else
            return false
        end
    else
        return true
    end
end
