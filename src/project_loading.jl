# this file contains functionality pertaining to loading a project's metadata and the
# energy system components from the project config file, as well as constructing certain
# helpful information data structures from the inputs in the config
using JSON: JSON
using OrderedCollections: OrderedDict

const HOURS_PER_SECOND::Float64 = 1.0 / 3600.0
const SECONDS_PER_HOUR::Float64 = 3600.0

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
    get_io_settings(project_config)

Constructs the dictionary of IO settings from the given config, considering default values.

# Arguments
-`project_config::Dict{AbstractString,Any}`: The project config
# Returns
-`Dict{String,Any}`: The IO settings dictionary
"""
function get_io_settings(project_config::AbstractDict{AbstractString,Any})::Dict{String,Any}
    io_settings = Dict{String,Any}(
        "plot_weather_data" => default(project_config["io_settings"], "plot_weather_data", false),
        "csv_output_weather" => default(project_config["io_settings"], "csv_output_weather", false),
        "write_csv_continuously" => default(project_config["io_settings"], "write_csv_continuously", false),
        "csv_output_file" => default(project_config["io_settings"], "csv_output_file", "./output/out.csv"),
        "csv_time_unit" => default(project_config["io_settings"], "csv_time_unit", "seconds"),
        "output_plot_file" => default(project_config["io_settings"], "output_plot_file", "./output/output_plot.html"),
        "output_plot_time_unit" => default(project_config["io_settings"], "output_plot_time_unit", "date"),
        "sankey_plot_file" => default(project_config["io_settings"], "sankey_plot_file", "./output/output_sankey.html"),
        "auxiliary_plots" => default(project_config["io_settings"], "auxiliary_plots", false),
        "show_detailed_errors" => default(project_config["io_settings"], "show_detailed_errors", false),
        "auxiliary_info" => default(project_config["io_settings"], "auxiliary_info", false),
        "auxiliary_info_file" => default(project_config["io_settings"], "auxiliary_info_file",
                                         "./output/auxiliary_info.md"),
        "auxiliary_plots_path" => default(project_config["io_settings"], "auxiliary_plots_path", "./output/"),
        "auxiliary_plots_formats" => default(project_config["io_settings"], "auxiliary_plots_formats", ["png"]),
        "economy_plot_file" => default(project_config["io_settings"], "economy_plot_file",
                                       "./output/economy_results.html"),
    )

    if !(io_settings["csv_time_unit"] in ["seconds", "minutes", "hours", "date"])
        @error "The IO setting `csv_time_unit` has to be one of: seconds, minutes, hours, date"
        throw(InputError())
    end

    if !(io_settings["output_plot_time_unit"] in ["seconds", "minutes", "hours", "date"])
        @info "The IO setting `output_plot_time_unit` has to be one of: seconds, minutes, " *
              "hours, date. It will be set to date as default."
        io_settings["output_plot_time_unit"] = "date"
    end

    if haskey(project_config["io_settings"], "sankey_plot")
        io_settings["sankey_plot"] = project_config["io_settings"]["sankey_plot"]
    end

    if haskey(project_config["io_settings"], "output_plot")
        io_settings["output_plot"] = project_config["io_settings"]["output_plot"]
    end

    if haskey(project_config["io_settings"], "step_info_interval")
        io_settings["step_info_interval"] = project_config["io_settings"]["step_info_interval"]
    end

    if haskey(project_config["io_settings"], "base_path")
        io_settings["base_path"] = abspath(project_config["io_settings"]["base_path"])
    else
        io_settings["base_path"] = abspath(joinpath(dirname(@__FILE__), ".."))
    end

    return io_settings
end

"""
    get_simulation_params(project_config, io_settings)

Constructs the dictionary of simulation parameters.

# Arguments
-`project_config::Dict{AbstractString,Any}`: The project config
-`io_settings::Dict{String,Any}`: IO settings, already extracted from the project config
# Returns
-`Dict{String,Any}`: The simulation parameter dictionary
"""
function get_simulation_params(project_config::AbstractDict{AbstractString,Any},
                               io_settings::Dict{String,Any})::Dict{String,Any}
    time_step,
    start_date,
    start_date_output,
    end_date,
    nr_of_steps,
    nr_of_steps_output = get_timesteps(project_config["simulation_parameters"])

    sim_params = Dict{String,Any}(
        "time" => 0,
        "time_since_output" => 0,
        "current_date" => start_date,
        "time_step_seconds" => time_step,
        "number_of_time_steps" => nr_of_steps,
        "number_of_time_steps_output" => nr_of_steps_output,
        "start_date" => start_date,
        "start_date_output" => start_date_output,
        "end_date" => end_date,
        "epsilon" => default(project_config["simulation_parameters"], "epsilon", 1e-9),
        "latitude" => default(project_config["simulation_parameters"], "latitude", nothing),
        "longitude" => default(project_config["simulation_parameters"], "longitude", nothing),
        "timezone" => default(project_config["simulation_parameters"], "time_zone", nothing),
        "step_info_interval" => default(io_settings, "step_info_interval", Integer(floor(nr_of_steps / 20))),
        "force_profiles_to_repeat" => default(project_config["simulation_parameters"], "force_profiles_to_repeat",
                                              false),
        "show_detailed_errors" => io_settings["show_detailed_errors"],
    )
    sim_params["economy_parameter"] = get_economy_parameter(project_config)
    sim_params["emissions_parameter"] = get_emission_parameter(project_config)

    # add helper functions to convert power to work and vice-versa. this uses the time step
    # of the simulation as the duration required for the conversion.
    sim_params["watt_to_wh"] = function (watts::Float64)
        return watts * time_step * HOURS_PER_SECOND
    end
    sim_params["wh_to_watts"] = function (wh::Float64)
        return wh * SECONDS_PER_HOUR / time_step
    end

    # add helper function for using paths, absolute or relative to the run base path
    sim_params["run_path"] = function (path)
        return isabspath(path) ? path : abspath(joinpath(io_settings["base_path"], path))
    end

    # load weather profiles accesible for all components
    weather_file_path = default(project_config["simulation_parameters"],
                                "weather_file_path",
                                nothing)

    if weather_file_path !== nothing
        weather_interpolation_type_solar = default(project_config["simulation_parameters"],
                                                   "weather_interpolation_type_solar", "linear_solar_radiation")
        weather_interpolation_type_general = default(project_config["simulation_parameters"],
                                                     "weather_interpolation_type_general", "linear_classic")
        # WeatherData() writes the lat and long to sim_params if they are not given in the input file
        sim_params["weather_data"] = WeatherData(sim_params["run_path"](weather_file_path),
                                                 sim_params,
                                                 guess_file_format(sim_params["run_path"](weather_file_path)),
                                                 weather_interpolation_type_solar,
                                                 weather_interpolation_type_general)
    end

    return sim_params
end

"""
    prepare_inputs(project_config)

Construct and prepare parameters, energy system components and the order of operation.

# Arguments
-`project_config::Dict{AbstractString,Any}`: The project config
# Returns
-`Dict{String,Any}`: Simulation parameters
-`Dict{String,Any}`: IO settings
-`Grouping`: The constructed energy system components
-`OrderOfOperations`: Order of operations
"""
function prepare_inputs(project_config::AbstractDict{AbstractString,Any}, run_ID::UUID)
    io_settings = get_io_settings(project_config)
    sim_params = get_simulation_params(project_config, io_settings)
    sim_params["run_ID"] = run_ID

    components = load_components(project_config["components"], sim_params)

    if haskey(project_config, "order_of_operation") && length(project_config["order_of_operation"]) > 0
        operations = load_order_of_operations(project_config["order_of_operation"], components)
        @info "The order of operations was successfully imported from the input file.\n" *
              "Note that the order of operations has a major impact on the simulation " *
              "result and should only be changed by experienced users!"
    else
        operations = calculate_order_of_operations(components)
    end

    return sim_params, io_settings, components, operations
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
                        @error "For component $unit_key, the given key `$medium` in the `output_refs` is not a valid key!"
                        throw(InputError())
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

    # Input side: include is_secondary_interface in the key
    component_key(uac::AbstractString, is_secondary_interface::Bool) = is_secondary_interface ?
                                                                       string(uac, "#secondary") : uac
    input_perm_indices = sortperm([input_order_dict[component_key(bus.input_interfaces[i].source.uac,
                                                                  bus.input_interfaces[i].is_secondary_interface)]
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
        throw(InputError())
    end
    if start_date_output < start_date
        @error "The start date of the output can not be prior to the start date of the simulation!"
        throw(InputError())
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
        throw(InputError())
    end
    return UInt(time_step), start_date, start_date_output, end_date, nr_of_steps, nr_of_steps_output
end

"""
    get_economy_parameter(project_config)

Extract economy parameters form input file.

Args:
-`project_config::AbstractDict{}`: The project config data
Return:
-`Dict{String,Any}`: The economy parameters from the input file. If non are given, calculate_economy will be set to false.
"""
function get_economy_parameter(project_config::AbstractDict{AbstractString,Any})::Dict{String,Any}
    if haskey(project_config, "economy_parameter")
        return Dict{String,Any}(
            "calculate_economy" => default(project_config["economy_parameter"], "calculate_economy", false),
            "observation_period_in_years" => default(project_config["economy_parameter"], "observation_period_in_years",
                                                     20.0),
            "interest_rate" => default(project_config["economy_parameter"], "interest_rate", 0.02),
            "labour_costs_per_hour" => default(project_config["economy_parameter"], "labour_costs_per_hour", 100.0),
            "labour_costs_price_change_rate_per_year" => default(project_config["economy_parameter"],
                                                                 "labour_costs_price_change_rate_per_year", 0.035),
        )
    else
        return Dict{String,Any}(
            "calculate_economy" => false,
        )
    end
end

"""
    get_emission_parameter(project_config)

Extract emission parameters form input file.

Args:
-`project_config::AbstractDict{}`: The project config data
Return:
-`Dict{String,Any}`: The emission parameters from the input file. If non are given, calculate_emissions will be set to false.
"""
function get_emission_parameter(project_config::AbstractDict{AbstractString,Any})::Dict{String,Any}
    if haskey(project_config, "emissions_parameter")
        return Dict{String,Any}(
            "calculate_emissions" => default(project_config["emissions_parameter"], "calculate_emissions", false),
            "observation_period_in_years" => default(project_config["emissions_parameter"],
                                                     "observation_period_in_years", 20.0),
            "include_embodied_emissions" => default(project_config["emissions_parameter"], "include_embodied_emissions",
                                                    true),
        )
    else
        return Dict{String,Any}(
            "calculate_emissions" => false,
        )
    end
end

# calculation of the order of operations has its own include files due to its complexity
include("order_of_operations.jl")
