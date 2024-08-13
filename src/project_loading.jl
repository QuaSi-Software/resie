# this file contains functionality pertaining to loading a project's metadata and the
# energy system components from the project config file, as well as constructing certain
# helpful information data structures from the inputs in the config
using JSON: JSON

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

    # create instances
    for (unit_key, entry) in pairs(config)
        default_dict = Dict{String,Any}()
        unit_config = Base.merge(default_dict, entry)

        symbol = Symbol(String(unit_config["type"]))
        unit_class = getproperty(EnergySystems, symbol)
        if unit_class <: EnergySystems.Component
            instance = unit_class(unit_key, unit_config, sim_params)
            components[unit_key] = instance
        end
    end

    # link inputs/outputs
    for (unit_key, entry) in pairs(config)
        if (String(entry["type"]) != "Bus"
            && haskey(entry, "output_refs")
            && length(entry["output_refs"]) > 0)
            # end of condition
            others = Grouping(key => components[key] for key in entry["output_refs"])
            link_output_with(components[unit_key], others)

        elseif (String(entry["type"]) == "Bus"
                && haskey(entry, "connections")
                && length(entry["connections"]) > 0)
            # end of condition
            others = Grouping(key => components[key]
                              for key in entry["connections"]["output_order"])
            link_output_with(components[unit_key], others)
        end
    end

    # add control modules to components
    for (unit_key, entry) in pairs(config)
        unit = components[unit_key]

        # TODO: rewrite this for automatic selection of modules so they don't need to be
        # registered here, compare automatic selection of component class
        for module_config in default(entry, "control_modules", [])
            if lowercase(module_config["name"]) === "economical_discharge"
                push!(unit.controller.modules,
                      EnergySystems.CM_EconomicalDischarge(module_config, components, sim_params))
            elseif lowercase(module_config["name"]) === "profile_limited"
                push!(unit.controller.modules,
                      EnergySystems.CM_ProfileLimited(module_config, components, sim_params))
            elseif lowercase(module_config["name"]) === "storage_driven"
                push!(unit.controller.modules,
                      EnergySystems.CM_StorageDriven(module_config, components, sim_params))
            end
        end
    end

    # the input/output interfaces of busses are constructed in the order of appearance in
    # the config, so after all components are loaded they need to be reordered to match
    # the input/output priorities
    components = reorder_interfaces_of_busses(components)

    # other type-specific initialisation
    EnergySystems.initialise_components(components, sim_params)

    # create proxy busses from bus chains
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
            output_perm_indices = sortperm([output_order_dict[unit.output_interfaces[i].target.uac]
                                            for i in 1:length(unit.output_interfaces)])
            input_perm_indices = sortperm([input_order_dict[unit.input_interfaces[i].source.uac]
                                           for i in 1:length(unit.input_interfaces)])

            # Reorder the input and output interfaces using the permutation indices
            unit.output_interfaces = unit.output_interfaces[output_perm_indices]
            unit.input_interfaces = unit.input_interfaces[input_perm_indices]
        end
    end
    return components
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
function get_timesteps(simulation_parameters::Dict{String,Any})
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

# calculation of the order of operations has its own include files due to its complexity
include("order_of_operations.jl")
