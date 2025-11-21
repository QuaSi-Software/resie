# this file contains functionality pertaining to loading a project's metadata and the
# energy system components from the project config file, as well as constructing certain
# helpful information data structures from the inputs in the config
using JSON: JSON
using OrderedCollections: OrderedDict

"""
    read_JSON(filepath)

Read and parse the JSON-encoded Dict in the given file.
"""
function read_JSON(filepath::String)::OrderedDict{AbstractString,Any}
    open(filepath, "r") do file_handle
        content = read(file_handle, String)
        return JSON.parse(content; dicttype=OrderedDict)
    end
end

"""
    load_control_module_class_mapping()

Loads the control modules' classes index by their name as used in the input file.

Returns:
-`Dict{String, Any}`: The mapping from name (String) to the module class, which is probably
    of type `Symbol`, however `getproperty` does not specify the return type. In any case
    the entry can be used for calling the constructor of the control module as a function.
"""
function load_control_module_class_mapping()::Dict{String,Any}
    mapping = Dict{String,Any}()

    for name in names(EnergySystems; all=true)
        if startswith(String(name), "CM_")
            symbol = Symbol(String(name))
            unit_class = getproperty(EnergySystems, symbol)

            if unit_class <: EnergySystems.ControlModule
                module_name = nothing

                # type-level accessor function implemented per control module as following:
                # control_module_name(::Type{CM_ModuleTypeName})::String = "module_type_name"
                try
                    module_name = EnergySystems.control_module_name(unit_class)
                    mapping[module_name] = unit_class
                catch
                    @error("Control module type $name does not have a method defined for " *
                           "function control_module_name.")
                end
            end
        end
    end

    return mapping
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
function load_components(config_ordered::AbstractDict{String,Any}, sim_params::Dict{String,Any})::Grouping
    # convert OrderedDict to normal Dict to have a normal dict in all components as they do not
    # require any sorting
    to_dict(x) = x
    to_dict(x::OrderedDict) = Dict{String,Any}(k => to_dict(v) for (k, v) in x)
    to_dict(x::AbstractVector) = map(to_dict, x)
    config = to_dict(config_ordered)

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
        if String(entry["type"]) != "Bus" && haskey(entry, "output_refs") && length(entry["output_refs"]) > 0
            if isa(entry["output_refs"], AbstractDict)
                # components with multiple outputs should enter the output_refs as Dict to achieve uniqueness 
                media_keys = collect(keys(entry["output_refs"]))
                target_components = [Grouping(uac => components[uac])
                                     for uac in [entry["output_refs"][key] for key in media_keys]]
                media_sym = Symbol[]
                for medium in media_keys
                    if hasproperty(components[unit_key], Symbol(medium))
                        push!(media_sym, getproperty(components[unit_key], Symbol(medium)))
                    else
                        @error "For component $unit_key, the key `$medium` in the `output_refs` could not be found!"
                        throw(InputError)
                    end
                end
                link_output_with(components[unit_key], target_components; given_media=media_sym)
            else
                if length(entry["output_refs"]) > 1
                    @warn "The component $unit_key has more than one output interface, but the `output_refs` are not " *
                          "specified explicitly! This can work, but it can also cause wrong interconnection between " *
                          "components! Consider using a mapping of the media to the target components "
                end
                target_components = Grouping(uac => components[uac] for uac in entry["output_refs"])
                link_output_with(components[unit_key], target_components)
            end
        elseif String(entry["type"]) == "Bus" && haskey(entry, "connections") && length(entry["connections"]) > 0
            target_components = Grouping(uac => components[uac] for uac in entry["connections"]["output_order"])
            link_output_with(components[unit_key], target_components)
        end
    end

    # add control modules to components
    mapping = load_control_module_class_mapping()
    for (unit_key, entry) in pairs(config)
        unit = components[unit_key]

        for module_config in default(entry, "control_modules", [])
            if !haskey(mapping, module_config["name"])
                @warn("Unknown control module type $(module_config["name"]) while loading " *
                      "unit $(unit.uac)")
                continue
            end
            module_class = mapping[module_config["name"]]
            push!(unit.controller.modules, module_class(module_config, components, sim_params, unit.uac))
        end
    end

    # the input/output interfaces of busses are constructed in the order of appearance in
    # the config, so after all components are loaded they need to be reordered to match
    # the input/output priorities
    components = reorder_interfaces_of_busses!(components)

    # other type-specific initialisation
    EnergySystems.initialise_components(components, sim_params)

    # create proxy busses from bus chains
    chains = find_chains(values(components), EnergySystems.sf_bus)
    EnergySystems.merge_bus_chains(chains, components, sim_params)

    return components
end

"""
    reorder_interfaces_of_busses!(components)

Calls reorder_interfaces_of_bus!() for all busses in the given grouping of components.

Args:
-`components::Grouping`: The components
Return:
-`Grouping`: The components with busses having their interfaces reordered
"""
function reorder_interfaces_of_busses!(components::Grouping)::Grouping
    for unit in each(components)
        if unit.sys_function == EnergySystems.sf_bus
            reorder_interfaces_of_bus!(unit)
        end
    end
    return components
end

"""
    reorder_interfaces_of_bus!(bus)

Reorder the input and output interfaces of busses according to their input and output
priorities given in the connectivity matrix.

Args:
-`bus::EnergySystems.Bus`: The bus for which to reorder interfaces
"""
function reorder_interfaces_of_bus!(bus::EnergySystems.Bus)
    # get correct order according to connectivity matrix
    output_order = bus.connectivity.output_order
    input_order = bus.connectivity.input_order

    # check for misconfigured bus (it should have at least one input and at least
    # one output)
    if length(input_order) == 0 || length(output_order) == 0
        return
    end

    # Create a dictionary to map 'uac' to its correct position
    output_order_dict = Dict(uac => idx for (idx, uac) in enumerate(output_order))
    input_order_dict = Dict(uac => idx for (idx, uac) in enumerate(input_order))

    # Get the permutation indices that would sort the 'source'/'target' field by
    # 'uac' order
    output_perm_indices = sortperm([output_order_dict[bus.output_interfaces[i].target.uac]
                                    for i in 1:length(bus.output_interfaces)])

    # Input side: include is_linked in the key
    component_key(uac::AbstractString, is_linked::Bool) = is_linked ? string(uac, "#linked") : uac
    input_perm_indices = sortperm([input_order_dict[component_key(bus.input_interfaces[i].source.uac,
                                                                  bus.input_interfaces[i].is_linked)]
                                   for i in eachindex(bus.input_interfaces)])

    # Reorder the input and output interfaces using the permutation indices
    bus.output_interfaces = bus.output_interfaces[output_perm_indices]
    bus.input_interfaces = bus.input_interfaces[input_perm_indices]
end

"""
get_timesteps(simulation_parameters)

Function to read in the time step information from the input file.
If no information is given in the input file, the following defaults 
will be set:
time_step = 900 s
"""
function get_timesteps(simulation_parameters::AbstractDict{String,Any})
    start_date = DateTime(0)
    start_date_output = DateTime(0)
    end_date = DateTime(0)
    try
        start_date = Dates.DateTime(simulation_parameters["start"], simulation_parameters["start_end_unit"])
        end_date = Dates.DateTime(simulation_parameters["end"], simulation_parameters["start_end_unit"])
        if haskey(simulation_parameters, "start_output")
            start_date_output = Dates.DateTime(simulation_parameters["start_output"],
                                               simulation_parameters["start_end_unit"])
        else
            start_date_output = start_date
        end
    catch e
        @error("Time given 'start_end_unit' of the simulation parameters does not fit to the data.\n" *
               "'start_end_unit' has to be a daytime format, e.g. 'dd-mm-yyyy HH:MM:SS'.\n" *
               "'start_end_unit' is `$(simulation_parameters["start_end_unit"])` which does not fit to the start" *
               "and end time given: `$(simulation_parameters["start"])` and `$(simulation_parameters["end"])`.\n" *
               "The following error occurred: $e")
        throw(InputError)
    end
    if start_date_output < start_date
        @error "The start date of the output can not be prior to the start date of the simulation!"
        throw(InputError)
    end
    if simulation_parameters["time_step_unit"] == "seconds"
        time_step = simulation_parameters["time_step"]
    elseif simulation_parameters["time_step_unit"] == "minutes"
        time_step = simulation_parameters["time_step"] * 60
    elseif simulation_parameters["time_step_unit"] == "hours"
        time_step = simulation_parameters["time_step"] * 60 * 60
    else
        time_step = 900
        @info("The simulation time step is set to 900 s as default, as it could not be found in the input file" *
              "(`time_step` and `time_step_unit` have to be given!).")
    end

    nr_of_steps = UInt(max(0, floor(Dates.value(Second(sub_ignoring_leap_days(end_date, start_date))) / time_step)) + 1)
    nr_of_steps_output = UInt(max(0,
                                  floor(Dates.value(Second(sub_ignoring_leap_days(end_date, start_date_output))) /
                                        time_step)) + 1)

    # set end_date to be integer dividable by the timestep
    end_date = add_ignoring_leap_days(start_date, (nr_of_steps - 1) * Second(time_step))

    if (month(start_date) == 2 && day(start_date) == 29) ||
       (month(end_date) == 2 && day(end_date) == 29) ||
       (month(start_date_output) == 2 && day(start_date_output) == 29)
        @error "The simulation start and end date and the start date of the output can not be at a leap day!"
        throw(InputError)
    end
    return UInt(time_step), start_date, start_date_output, end_date, nr_of_steps, nr_of_steps_output
end

# calculation of the order of operations has its own include files due to its complexity
include("order_of_operations.jl")
