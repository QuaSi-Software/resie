# this file contains functionality for writing output of the simulation to files.
using Dates
using Random
using CSV

"""
    get_output_keys(io_settings, components)

Determines output keys for:
  - lineplot
  - csv export
  - economic output (filtered to value_key containing "OUT" or "IN")

For each output channel:
  - if not requested, returns `nothing`
  - if requested, returns either:
      - all possible keys (in-/excluding flows), or
      - only requested keys from input file
"""
function get_output_keys(io_settings::AbstractDict{String,Any},
                         economic_parameters::Union{Nothing,AbstractDict{String,Any}},
                         emissions_parameters::Union{Nothing,AbstractDict{String,Any}},
                         optimiser_parameters::Union{Nothing,AbstractDict{String,Any}},
                         components::Grouping,
                         suppress_all_output::Bool)::Tuple{Union{Nothing,Vector{EnergySystems.OutputKey}},
                                                           Union{Nothing,Vector{EnergySystems.OutputKey}},
                                                           Union{Nothing,Vector{EnergySystems.OutputKey}},
                                                           Union{Nothing,Vector{EnergySystems.OutputKey}}}
    function parse_all_mode(io_settings, setting_name::String, suppress_all_output::Bool)
        if suppress_all_output
            return false, false, false
        end
        do_create = false
        do_all_excl = false
        do_all_incl = false

        if haskey(io_settings, setting_name)
            do_create = true
            v = io_settings[setting_name]
            if v == "all_excl_flows"
                do_all_excl = true
            elseif v == "all_incl_flows"
                do_all_incl = true
            elseif v == "nothing"
                do_create = false
            elseif v == "all"
                @error "For \"$(setting_name)\", the input \"all\" is no longer supported. Use \"all_incl_flows\" or \"all_excl_flows\"."
                throw(InputError())
            end
        end

        return do_create, do_all_excl, do_all_incl
    end

    # Sorting the keys
    function sort_by(output_key)
        any(startswith(output_key.value_key, p) for p in ("EnergyFlow", "TemperatureFlow")) ?
        # flows after others: primary = medium (always present), secondary
        lowercase("zzzzzzzzzzzzz" * string(output_key.medium) * output_key.value_key) :
        # default: primary = unit.uac, secondary = medium (if present), tertiary
        lowercase(string(output_key.unit.uac) * string(something(output_key.medium, "")) * output_key.value_key)
    end

    # Build "all outputs" once per include/exclude flows
    function collect_all_output_keys(components::Grouping; include_flows::Bool)::Vector{EnergySystems.OutputKey}
        all_keys = Vector{EnergySystems.OutputKey}()

        for unit in components
            output_vals = output_values(unit[2])
            temp_dict = Dict{String,Any}()

            for output_val in output_vals
                if startswith(output_val, "TemperatureFlow")
                    # do nothing, temperature are added in output_keys()
                    continue
                end

                if startswith(output_val, "EnergyFlow")
                    if !include_flows
                        continue
                    end

                    # handle secondary interfaces
                    if startswith(output_val, create_secondary_name("EnergyFlow"))
                        key = adjust_name_if_secondary(String(unit[2].medium), true)
                        nr_skip = 2 + length(create_secondary_name("EnergyFlow"))
                    else
                        key = String(unit[2].medium)
                        nr_skip = 12
                    end

                    if haskey(temp_dict, key)
                        push!(temp_dict[key], output_val[nr_skip:end])
                    else
                        temp_dict[key] = [output_val[nr_skip:end]]
                    end
                else
                    # Non-flow output => keyed by unit.uac
                    key = unit[2].uac
                    if haskey(temp_dict, key)
                        push!(temp_dict[key], output_val)
                    else
                        temp_dict[key] = [output_val]
                    end
                end
            end
            append!(all_keys, output_keys(components, temp_dict))
        end

        sort!(all_keys; by=sort_by)
        return all_keys
    end

    # Select keys for a channel, given its parsed mode and custom inputs
    #TODO not used right now
    function select_keys_for_channel(do_create::Bool,
                                     mode::Symbol,
                                     setting_name::String,
                                     all_excl::Union{Nothing,Vector{EnergySystems.OutputKey}},
                                     all_incl::Union{Nothing,Vector{EnergySystems.OutputKey}};
                                     custom_extractor::Function)
        if !do_create
            return nothing
        end

        if mode == :all_incl_flows
            return all_incl
        elseif mode == :all_excl_flows
            return all_excl
        else
            return custom_extractor(io_settings[setting_name])
        end
    end

    # Economy and emissions filter
    function is_economic_emissions_key(ok::EnergySystems.OutputKey)
        occursin("OUT", ok.value_key) ||
            occursin("IN", ok.value_key) ||
            occursin("Supply", ok.value_key) ||
            occursin("Demand", ok.value_key)
    end

    # get requirements
    do_create_plot, do_plot_all_excl, do_plot_all_incl = parse_all_mode(io_settings, "output_plot", suppress_all_output)
    do_write_CSV, do_csv_all_excl, do_csv_all_incl = parse_all_mode(io_settings, "csv_output", suppress_all_output)
    do_economy_emissions = economic_parameters["calculate_economy"] || emissions_parameters["calculate_emissions"]
    do_optimise = optimiser_parameters["run_optimisation"]
    do_matrix_plot = io_settings["matrix_plot"] == "custom"

    # Decide if we need all-keys lists
    need_all_excl = do_plot_all_excl || do_csv_all_excl || do_economy_emissions
    need_all_incl = do_plot_all_incl || do_csv_all_incl

    all_output_keys_excl_flows = need_all_excl ? collect_all_output_keys(components; include_flows=false) : nothing
    all_output_keys_incl_flows = need_all_incl ? collect_all_output_keys(components; include_flows=true) : nothing

    # Lineplot keys
    if do_create_plot
        if do_plot_all_incl
            output_keys_lineplot = all_output_keys_incl_flows
        elseif do_plot_all_excl
            output_keys_lineplot = all_output_keys_excl_flows
        else
            output_keys_lineplot = Vector{EnergySystems.OutputKey}()
            for (_, plot) in sort(collect(io_settings["output_plot_spec"]); by=p -> parse(Int, p[1]))
                append!(output_keys_lineplot, output_keys(components, plot["key"]))
            end
        end
    else
        output_keys_lineplot = nothing
    end

    # CSV keys
    if do_write_CSV
        if do_csv_all_incl
            output_keys_to_csv = all_output_keys_incl_flows
        elseif do_csv_all_excl
            output_keys_to_csv = all_output_keys_excl_flows
        else
            output_keys_to_csv = output_keys(components, io_settings["csv_output_keys"])
        end
    else
        output_keys_to_csv = nothing
    end

    # Economy or emissions keys
    if do_economy_emissions
        # Use excl_flows "all" as base
        output_keys_economic_emissions = copy(all_output_keys_excl_flows)
        filter!(is_economic_emissions_key, output_keys_economic_emissions)
        # keep stable ordering (all_output_keys_excl_flows is already sorted)
    else
        output_keys_economic_emissions = nothing
    end

    # optimiser keys
    if do_optimise
        output_keys_optimise = output_keys(components, optimiser_parameters["objective_keys_sum_mean"])
        if do_matrix_plot
            for (func, spec) in pairs(io_settings["matrix_plot_spec"])
                if func == "sum" || func == "mean"
                    matrix_keys = output_keys(components, spec)
                    append!(output_keys_optimise, matrix_keys)
                end
            end
        end
        output_keys_optimise = unique(output_keys_optimise)
    else
        output_keys_optimise = nothing
    end

    return output_keys_lineplot, output_keys_to_csv, output_keys_economic_emissions, output_keys_optimise
end

"""
get_interface_information(components)


Function to gather information for the sankey diagram.
Determines 
- the total number of present system interfaces [int]
- the corresponding medium in each interface [medium]
- the source component of each interface and [unit]
- the target component of each interface [unit]

The information is returned as single vectors with the indices matching together.
"""
function get_interface_information(components::Grouping)::Tuple{Int64,Vector{Any},Vector{Any},Vector{Any}}
    nr_of_interfaces = 0
    medium_of_interfaces = []
    output_sourcenames_sankey = []
    output_targetnames_sankey = []
    for each_component in components
        for each_outputinterface in each_component[2].output_interfaces
            medium = nothing  # reset medium
            if isa(each_outputinterface, Pair) # some output_interfaces are wrapped in a Tuple
                medium = each_outputinterface[1] # then, the medium is stored separately
                each_outputinterface = each_outputinterface[2]
            end

            if (isdefined(each_outputinterface, :target) &&
                !startswith(each_outputinterface.target.uac, "Proxy") &&
                !startswith(each_outputinterface.source.uac, "Proxy"))

                # count interface
                nr_of_interfaces += 1

                #get name of source and sink
                push!(output_sourcenames_sankey, each_outputinterface.source.uac)
                push!(output_targetnames_sankey, each_outputinterface.target.uac)

                # get name of medium
                if !(medium === nothing)
                    push!(medium_of_interfaces, medium)
                elseif isdefined(each_outputinterface.target, :medium)
                    push!(medium_of_interfaces, each_outputinterface.target.medium)
                elseif isdefined(each_outputinterface.source, :medium)
                    push!(medium_of_interfaces, each_outputinterface.source.medium)
                else
                    @warn "The name of the medium was not detected. This may lead to wrong colouring in Sankey plot."
                end

                # add "real" demands and sources
                if each_outputinterface.source.sys_function == EnergySystems.sf_fixed_source
                    push!(output_sourcenames_sankey, string(each_outputinterface.source.uac, "_total_supply"))
                    push!(output_targetnames_sankey, each_outputinterface.source.uac)
                    push!(medium_of_interfaces, "hide_medium")
                    nr_of_interfaces += 1
                end
                if each_outputinterface.target.sys_function == EnergySystems.sf_fixed_sink
                    push!(output_sourcenames_sankey, each_outputinterface.target.uac)
                    push!(output_targetnames_sankey, string(each_outputinterface.target.uac, "_total_demand"))
                    push!(medium_of_interfaces, "hide_medium")
                    nr_of_interfaces += 1
                end
            end
        end

        # add losses
        if "LossesGains" in output_values(each_component[2])
            push!(output_sourcenames_sankey, each_component[2].uac)
            push!(output_targetnames_sankey, "LossesGains")
            push!(medium_of_interfaces, "LossesGains")
            nr_of_interfaces += 1
        end
    end
    if length(medium_of_interfaces) !== nr_of_interfaces
        @error "Error in extracting information from input file for sankey plot."
    end

    return nr_of_interfaces, medium_of_interfaces, output_sourcenames_sankey, output_targetnames_sankey
end

"""
collect_interface_energies(components, nr_of_interfaces)

Collects and returns the energy that was transported through every interface.
If the balance of an interface was not zero, the actual energy that was flowing
is written to the outputs.
Attention: This can lead to overfilling of demands which is currently not visible
in the sankey diagram!
"""
function collect_interface_energies(components::Grouping, nr_of_interfaces::Int)
    n = 1
    energies = zeros(Float64, nr_of_interfaces)
    for each_component in components
        for each_outputinterface in each_component[2].output_interfaces
            if isa(each_outputinterface, Pair) # some output_interfaces are wrapped in a Tuple
                each_outputinterface = each_outputinterface[2]
            end
            if (isdefined(each_outputinterface, :target)
                && !startswith(each_outputinterface.target.uac, "Proxy")
                && !startswith(each_outputinterface.source.uac, "Proxy"))
                # end of condition
                energies[n] = calculate_energy_flow(each_outputinterface)
                n += 1

                # If source or target is fixed source or sink, gather also demand and supply
                if each_outputinterface.source.sys_function == EnergySystems.sf_fixed_source
                    energies[n] = each_outputinterface.source.supply
                    n += 1
                end

                if each_outputinterface.target.sys_function == EnergySystems.sf_fixed_sink
                    energies[n] = each_outputinterface.target.demand
                    n += 1
                end
            end
        end

        # add losses
        if "LossesGains" in output_values(each_component[2])
            energies[n] = output_value(each_component[2],
                                       EnergySystems.OutputKey(; unit=each_component[2],
                                                               medium=nothing,
                                                               value_key="LossesGains"))
            n += 1
        end
    end
    return energies
end

"""
output_keys(components, from_config)

Transform the output keys definition in the project config file into a list of OutputKey
items. This is done to speed up selection of values for the output in each time step,
as this transformation has to be done only once at the beginning.
"""
function output_keys(components::Grouping, from_config::AbstractDict{String,Any})::Vector{EnergySystems.OutputKey}
    outputs = Vector{EnergySystems.OutputKey}()

    all_current_media = []
    for component in values(components)
        if component.sys_function === EnergySystems.sf_bus
            push!(all_current_media, String(component.medium))
            for inface in component.input_interfaces
                if inface.is_secondary_interface
                    push!(all_current_media, adjust_name_if_secondary(String(component.medium), true))
                end
            end
        end
    end
    all_current_media = unique(all_current_media)

    for (key, _) in sort(collect(from_config); by=k -> k[1])
        if key in keys(components)
            unit_key = key
            unit = components[unit_key]

            for entry in sort(from_config[unit_key])
                splitted = split(String(entry), ":")
                if length(splitted) > 1
                    medium_key = splitted[1]
                    medium = Symbol(String(medium_key))
                    unit_fields = fieldnames(typeof(unit))
                    unit_media = [getfield(unit, f)
                                  for f in unit_fields if startswith(String(f), "m_") || String(f) == "medium"]
                    if medium in unit_media
                        value_key = splitted[2]
                    else
                        @error "In unit \"$(unit.uac)\", the given output key \"$entry\" could not be mapped to an " *
                               "output key. Make sure that the medium \"$(String(medium_key))\" exists in the current " *
                               "component and that you have used \":\" as separator without any extra spaces."
                        throw(InputError())
                    end
                else
                    medium = nothing
                    value_key = splitted[1]
                end

                push!(outputs, EnergySystems.OutputKey(; unit=unit,
                                                       medium=medium,
                                                       value_key=value_key))
            end
        elseif key in all_current_media
            media_key = key
            medium = Symbol(String(media_key))
            for value_key in sort(from_config[media_key])
                success = false
                in_uac, out_uac = split(value_key, "->")
                for bus in [unit for unit in values(components) if unit.sys_function === EnergySystems.sf_bus]
                    # consider only proxy busses or busses without proxies and busses with correct media
                    medium_trimmed, _ = trim_secondary_medium(medium)
                    if bus.proxy === nothing && bus.medium == medium_trimmed
                        # check if input and output exists
                        if in_uac in keys(bus.balance_table_inputs) && out_uac in keys(bus.balance_table_outputs)
                            push!(outputs,
                                  EnergySystems.OutputKey(; unit=bus,
                                                          medium=medium,
                                                          value_key="EnergyFlow " * value_key))
                            push!(outputs,
                                  EnergySystems.OutputKey(; unit=bus,
                                                          medium=medium,
                                                          value_key="TemperatureFlow " * value_key))
                            success = true
                            break
                        end
                    end
                end
                if !success
                    @error "The requested energy flow between components \"$(value_key)\" for medium \"$(media_key)\" could not " *
                           "be found for the CSV output or the plot output. Note that only connections between " *
                           "components with one or more busses but without any other component in between can be exported!"
                    throw(InputError())
                end
            end
        else
            @error "The key \"$(key)\" in the provided output keys for CSV output or plot output could not be found. " *
                   "It either has to be a medium or a component used in the current energy system."
            throw(InputError())
        end
    end

    return outputs
end

function parse_outkeys(output_keys::Vector{EnergySystems.OutputKey})::Array{String}
    keys = Array{String}(undef, length(output_keys))
    for (idx, outkey) in enumerate(output_keys)
        if outkey.medium === nothing
            keys[idx] = "$(outkey.unit.uac) $(outkey.value_key)"
        else
            if startswith(outkey.value_key, "EnergyFlow") || startswith(outkey.value_key, "TemperatureFlow")
                keys[idx] = "$(outkey.medium) $(outkey.value_key)"
            else
                keys[idx] = "$(outkey.unit.uac) $(outkey.medium) $(outkey.value_key)"
            end
        end
    end
    return keys
end

function parse_outkeys(output_keys::AbstractDict{String,Any})::Array{String}
    keys = Array{String}(undef, 0)
    for (uac, values) in pairs(output_keys)
        for entry in values
            splitted = split(String(entry), ":")
            if length(splitted) > 1
                medium = splitted[1]
                value_key = splitted[2]
                if startswith(value_key, "EnergyFlow") || startswith(value_key, "TemperatureFlow")
                    push!(keys, "$medium $value_key")
                else
                    push!(keys, "$uac $medium $value_key")
                end
            else
                medium = nothing
                value_key = splitted[1]
                push!(keys, "$uac $value_key")
            end
        end
    end

    return keys
end

"""
get_output_header(output_keys, weather_data_keys, csv_time_unit)

Get the output header for the given outputs to used in output file or dictionary.
"""
function get_output_header(output_keys::Union{Nothing,Vector{EnergySystems.OutputKey}},
                           weather_data_keys::Union{Nothing,Vector{String}},
                           csv_time_unit::String)
    header = Array{String}(undef, 0)
    if csv_time_unit !== nothing
        if csv_time_unit == "seconds"
            time_unit = "[s]"
        elseif csv_time_unit == "minutes"
            time_unit = "[min]"
        elseif csv_time_unit == "hours"
            time_unit = "[h]"
        elseif csv_time_unit == "date"
            time_unit = "[dd.mm.yyyy HH:MM:SS]"
        end

        push!(header, "Time $time_unit")
    end

    if output_keys !== nothing
        output_keys_names = parse_outkeys(output_keys)
        header = vcat(header, output_keys_names)
    end

    if weather_data_keys !== nothing
        for key in weather_data_keys
            push!(header, "Weather $key")
        end
    end

    return header
end

"""
get_output_row(output_keys, weather_data_keys, sim_params, csv_time_unit)

Create a row with values for given outputs to be written to file or dictionary.
"""
function get_output_row(output_keys::Union{Nothing,Vector{EnergySystems.OutputKey}},
                        weather_data_keys::Union{Nothing,Vector{String}},
                        sim_params::Dict{String,Any},
                        csv_time_unit::String,
                        io_settings::Dict{String,Any})
    row = Array{Union{Float64,String}}(undef, 0)
    if csv_time_unit == "seconds"
        time = sim_params["time_since_output"]
    elseif csv_time_unit == "minutes"
        time = sim_params["time_since_output"] / 60
    elseif csv_time_unit == "hours"
        time = sim_params["time_since_output"] / 60 / 60
    elseif csv_time_unit == "date"
        time = Dates.format(sim_params["current_date"], "dd.mm.yyyy HH:MM:SS")
    end

    push!(row, string(time))

    interpolator = v -> "$v"
    if io_settings["fixed_output_precision"] > 0
        interpolator = v -> "$(round(v; sigdigits=io_settings["fixed_output_precision"]))"
    end

    if output_keys !== nothing
        for outkey in output_keys
            value = output_value(outkey.unit, outkey)
            push!(row, interpolator(value))
        end
    end
    if weather_data_keys !== nothing
        for key in weather_data_keys
            value = Profiles.value_at_time(getfield(sim_params["weather_data"], Symbol(key)), sim_params)
            push!(row, interpolator(value))
        end
    end

    return row
end

"""
    listify_operations(operations)

Turns the given order of operations into a list with entries surrounded in quotation marks
and seperated by a comma and line feed.

Args:
-`operations::OrderOfOperations`: The operations to listify
Returns:
-`String`: The listified operations
"""
function listify_operations(operations::OrderOfOperations)::String
    list = ""
    for entry in operations
        comma = ","
        if entry == last(operations)
            comma = ""
        end
        list = list * "\"$(entry[1]):$(entry[2])\"$(comma)\n"
    end
    return list
end

"""
    dump_auxiliary_outputs(io_settings, components, order_of_operations, sim_params)

Dump a bunch of information to file that might be useful to explain the result of a run.

This is mostly used for debugging and development purposes, but might prove useful in
general to find out why the energy system behaves in the simulation as it does.
"""
function dump_auxiliary_outputs(io_settings::Dict{String,Any},
                                components::Grouping,
                                order_of_operations::OrderOfOperations,
                                sim_params::Dict{String,Any},
                                suppress_all_output::Bool)
    if suppress_all_output
        return
    end
    # export order of operations
    if io_settings["auxiliary_info"]
        aux_info_file_path = io_settings["auxiliary_info_file"]
        open(sim_params["run_path"](aux_info_file_path), "w") do file_handle
            # write base order (from input or calculated)
            write(file_handle, "# Order of operations\n")
            write(file_handle, listify_operations(order_of_operations))

            # look for any control modules that modify it and print the modified one
            for component in values(components)
                if component.sys_function == EnergySystems.sf_bus
                    for control_module in component.controller.modules
                        # this is very specific for the current implementation of this exact
                        # control module. @TODO make this more generalized once more control
                        # modules also modify the order of operations
                        if control_module.name == "economic_control"
                            for (state_id, order) in pairs(control_module.ooo_by_state)
                                write(file_handle, "\n# Order of operations $(component.uac) state #$(state_id)\n")
                                write(file_handle, listify_operations(order))
                            end
                        end
                    end
                end
            end
        end

        @info "Auxiliary info dumped to file $(sim_params["run_path"](aux_info_file_path))"
    end

    # plot additional figures potentially available from components after initialisation
    if io_settings["auxiliary_plots"]
        aux_plots_output_path = sim_params["run_path"](io_settings["auxiliary_plots_path"])
        aux_plots_formats = io_settings["auxiliary_plots_formats"]
        aux_plots_formats = Vector{String}(aux_plots_formats)
        component_list = []
        for component in components
            if plot_optional_figures_begin(component[2], aux_plots_output_path, aux_plots_formats, sim_params)
                push!(component_list, component[2].uac)
            end
        end
        if length(component_list) > 0
            @info "Auxiliary plots are saved to folder $(aux_plots_output_path) for the following components: $(join(component_list, ", "))"
        end
    end
end

"""
gather_output_data(output_keys, time)

returns a vector with the requested data in output_keys
"""
function gather_output_data(output_keys::Vector{EnergySystems.OutputKey}, time::Int)
    return_values = zeros(Union{Int,Float64}, length(output_keys) + 1)
    return_values[1] = time

    for (idx, outkey) in enumerate(output_keys)
        return_values[idx + 1] = output_value(outkey.unit, outkey)
    end

    return return_values
end

"""
    create_profile_line_plots(data, keys, weather_data, weather_keys, io_settings, sim_params)

Creates the line plots for the given output configuration.
"""
function create_profile_line_plots(outputs_plot_data::Union{Nothing,Matrix{Float64}},
                                   outputs_plot_keys::Union{Nothing,Vector{EnergySystems.OutputKey}},
                                   outputs_plot_weather::Union{Nothing,Matrix{Float64}},
                                   outputs_plot_weather_keys::Union{Nothing,Vector{String}},
                                   io_settings::Dict{String,Any},
                                   sim_params::Dict{String,Any})
    plot_all = isa(io_settings["output_plot"], String) &&
               io_settings["output_plot"][1:3] == "all"
    plot_data = outputs_plot_data !== nothing
    plot_weather = outputs_plot_weather !== nothing

    # set Axis, unit and scale factor if given
    if plot_all  # plot all outputs. Here no units or scaling factors are available.
        labels = parse_outkeys(outputs_plot_keys)
        if plot_weather
            for outkey in outputs_plot_weather_keys
                push!(labels, string("Weather $(outkey)"))
            end
        end
    else # plot only defined outputs. Here units and scaling factors are available.
        axis = String[]
        unit = String[]
        scale_fact = Float64[]
        if plot_data
            for (nr, plot) in sort(collect(io_settings["output_plot_spec"]); by=p -> parse(Int, p[1]))
                if occursin("->", first(plot["key"])[2][1])
                    # Here we are dealing with EnergyFlow and TemperatureFlow --> Two meta information required if 
                    # temperature should be plotted
                    if !isa(plot["axis"], AbstractVector) || length(plot["axis"]) !== 2 ||
                       !isa(plot["unit"], AbstractVector) || length(plot["unit"]) !== 2 ||
                       !isa(plot["scale_factor"], AbstractVector) || length(plot["scale_factor"]) !== 2
                        # only one set of meta information given. Use the given one both for energy and temperature flow
                        push!(axis, isa(plot["axis"], AbstractVector) ? plot["axis"][1] : plot["axis"])
                        push!(unit, isa(plot["unit"], AbstractVector) ? plot["unit"][1] : plot["unit"])
                        push!(scale_fact,
                              isa(plot["scale_factor"], AbstractVector) ? plot["scale_factor"][1] :
                              plot["scale_factor"])

                        push!(axis, "nothing")
                        push!(unit, "nothing")
                        push!(scale_fact, NaN)
                        @info "For the generation of the output plot, the meta information for entry $nr " *
                              "do not contain two values. Therefore, only energy values will be output. If you want " *
                              "to output also a corresponding temperature, provide two meta information: " *
                              "[EnergyFlow, TemperatureFlow] for axis, unit and scale_factor."
                    else
                        push!(axis, string(plot["axis"][1]))
                        push!(unit, string(plot["unit"][1]))
                        push!(scale_fact, plot["scale_factor"][1])

                        push!(axis, string(plot["axis"][2]))
                        push!(unit, string(plot["unit"][2]))
                        push!(scale_fact, plot["scale_factor"][2])
                    end
                else
                    push!(axis, string(plot["axis"]))
                    push!(unit, string(plot["unit"]))
                    push!(scale_fact, plot["scale_factor"])
                end
            end
        end
        if plot_weather
            for _ in outputs_plot_weather_keys
                push!(axis, "left")
                push!(scale_fact, 1.0)
            end
        end

        # create legend entries
        labels = String[]
        if plot_data
            labels = parse_outkeys(outputs_plot_keys) .* " [" .* unit .* "] (" .*
                     axis[1:length(outputs_plot_keys)] .* ")"
        end
        if plot_weather
            for outkey in outputs_plot_weather_keys
                push!(labels, string("Weather $(outkey) [($(axis[end]))"))
            end
        end
    end

    # create plot
    if plot_data
        time_x = outputs_plot_data[:, 1]
    else
        time_x = outputs_plot_weather[:, 1]
    end
    # filter NaN values
    if plot_data
        idxs_to_remove = Int[]
        for col in axes(outputs_plot_data, 2)
            if col == 1 # never remove date
                continue
            elseif all(isnan, outputs_plot_data[:, col])
                push!(idxs_to_remove, col)
            end
        end
        if !plot_all
            for (idx, ax) in enumerate(axis)
                if ax == "nothing"
                    push!(idxs_to_remove, idx + 1)
                end
            end
            idxs_to_remove = unique(idxs_to_remove)
        end
        if !isempty(idxs_to_remove)
            outputs_plot_data = outputs_plot_data[:, setdiff(axes(outputs_plot_data, 2), idxs_to_remove)]
            labels = labels[setdiff(1:length(labels), idxs_to_remove .- 1)]
            if !plot_all
                axis = axis[setdiff(1:length(axis), idxs_to_remove .- 1)]
                scale_fact = scale_fact[setdiff(1:length(scale_fact), idxs_to_remove .- 1)]
            end
        end
    end
    output_plot_time_unit = io_settings["output_plot_time_unit"]
    if output_plot_time_unit == "seconds"
        x = time_x
    elseif output_plot_time_unit == "minutes"
        x = time_x / 60
    elseif output_plot_time_unit == "hours"
        x = time_x / 60 / 60
    else
        if sim_params["start_date"] === nothing
            start_date = Dates.DateTime("2015/1/1 00:00:00", "yyyy/m/d HH:MM:SS")
            @info ("Date of first data point in output line plot is set to 01-01-2015 00:00:00, as the simulation start time is not given as date.")
        else
            start_date = sim_params["start_date_output"]
        end
        x = [add_ignoring_leap_days(start_date, Dates.Second(s)) for s in time_x]
    end

    y = hcat(plot_data ? outputs_plot_data[:, 2:end] : zeros(Float64, size(outputs_plot_weather, 1), 0),
             plot_weather ? outputs_plot_weather[:, 2:end] : zeros(Float64, size(outputs_plot_data, 1), 0))

    if io_settings["fixed_output_precision"] > 0
        y = round.(y; sigdigits=io_settings["fixed_output_precision"])
    end

    traces = GenericTrace[]
    for i in axes(y, 2)
        if plot_all
            trace = scatter(; x=x, y=y[:, i], mode="lines", name=labels[i])
        else
            trace = scatter(; x=x, y=scale_fact[i] * y[:, i], mode="lines", name=labels[i])
            if axis[i] == "right"
                trace.yaxis = "y2"
            else  # default is left axis
                trace.yaxis = "y1"
            end
        end
        push!(traces, trace)
    end

    if output_plot_time_unit == "date"
        leap_days_str = string.([Date(year, 2, 29)
                                 for year in
                                     Dates.value(Year(sim_params["start_date_output"])):Dates.value(Year(sim_params["end_date"]))
                                 if isleapyear(year)])

        layout = Layout(;
                        title_text="Plot of outputs as defined in the input-file. Attention: Energies are given within " *
                                   "the simulation time step of $(Int(sim_params["time_step_seconds"])) s",
                        xaxis_title_text="Time [$(output_plot_time_unit)]",
                        yaxis_title_text="",
                        yaxis2=attr(; title="", overlaying="y", side="right"),
                        xaxis=attr(; type="date",
                                   rangebreaks=[Dict("values" => leap_days_str)]))
    else
        layout = Layout(;
                        title_text="Plot of outputs as defined in the input-file. Attention: Energies are given within " *
                                   "the simulation time step of $(Int(sim_params["time_step_seconds"])) s",
                        xaxis_title_text="Time [$(output_plot_time_unit)]",
                        yaxis_title_text="",
                        yaxis2=attr(; title="", overlaying="y", side="right"))
    end

    p = plot(traces, layout)
    file_path = sim_params["run_path"](io_settings["output_plot_file_path"])
    savefig(p, file_path)
end

"""
    create_sankey(output_all_sourcenames, output_all_targetnames, output_all_values,
                  medium_of_interfaces, nr_of_interfaces, io_settings, sim_params)

create a sankey plot. 
Inputs:
output_all_sourcenames and *sinknames are vectors with names of the source and sink of each interface
output_all_values are logs with data from each timestep in the shape [timestep,interface].
medium_of_interface is a vector of the medium corresponding to each interface
nr_of_interfaces is the total number of interfaces in the current energy system
"""
function create_sankey(output_all_sourcenames::Vector{Any},
                       output_all_targetnames::Vector{Any},
                       output_all_values::Matrix{Float64},
                       medium_of_interfaces::Vector{Any},
                       nr_of_interfaces::Int64,
                       io_settings::Dict{String,Any},
                       sim_params::Dict{String,Any})

    # sum up data of each interface
    output_all_value_sum = zeros(Float64, nr_of_interfaces)
    for interface in 1:nr_of_interfaces
        output_all_value_sum[interface] = sum(output_all_values[:, interface])
    end

    # convert Losses into Gains if they are negative
    for idx in 1:nr_of_interfaces
        if medium_of_interfaces[idx] == "LossesGains"
            if output_all_value_sum[idx] < 0.0  # Losses
                medium_of_interfaces[idx] = "Losses"
                output_all_value_sum[idx] = -output_all_value_sum[idx]
                output_all_targetnames[idx] = "Losses"
            else # Gains
                medium_of_interfaces[idx] = "Gains"
                output_all_targetnames[idx] = output_all_sourcenames[idx]
                output_all_sourcenames[idx] = "Gains"
            end
        end
    end

    # remove data that should not be plotted in Sankey
    interface_new = 1
    for _ in 1:nr_of_interfaces
        if (
            # remove oxygen from data as the energy of oxygen is considered to be zero
            medium_of_interfaces[interface_new] == :m_c_g_o2
            # remove real sinks and sources if they match the delivered energy
            # to enable this to work, the "real" demand/supply needs to be always one entry
            # below the delivered/requested one in the array!
            ||
            (medium_of_interfaces[interface_new] == "hide_medium"
             &&
             (abs(output_all_value_sum[interface_new] - output_all_value_sum[interface_new - 1])) <
             sim_params["epsilon"]))
            # end of condition 
            deleteat!(output_all_sourcenames, interface_new)
            deleteat!(output_all_targetnames, interface_new)
            deleteat!(output_all_value_sum, interface_new)
            deleteat!(medium_of_interfaces, interface_new)
            interface_new -= 1
            nr_of_interfaces -= 1
        end
        interface_new += 1
    end
    interface_new -= 1

    # apply fixed precision before adding for non-zero as it may otherwise be rounded to zero again
    if io_settings["fixed_output_precision"] > 0
        output_all_value_sum = round.(output_all_value_sum; sigdigits=io_settings["fixed_output_precision"])
    end

    # add 0.000001 to all interfaces (except of losses and gains) to display interfaces that are zero
    output_all_value_sum += .![medium in ["Losses", "Gains"] for medium in medium_of_interfaces] * 0.000001

    # prepare data for sankey diagram and create sankey
    # set label of blocks
    block_labels = union(output_all_sourcenames, output_all_targetnames)
    block_labels_unique = Dict(blockname => i for (i, blockname) in enumerate(block_labels))
    output_all_source_num = [block_labels_unique[blockname] for blockname in output_all_sourcenames]
    output_all_target_num = [block_labels_unique[blockname] for blockname in output_all_targetnames]

    # set label and colour of interfaces
    medium_labels = [split(string(s), '.')[end] for s in medium_of_interfaces]
    unique_medium_labels = unique(medium_labels)

    if io_settings["sankey_plot"] == "default" || io_settings["sankey_plot"] == "custom"
        # in both cases of "default" or "custom", the setting sankey_plot_spec already
        # contains a color map definition, but it may not cover all media, thus we randomly
        # assign colors to the missing media
        color_map = Dict{String,Any}()
        for label in unique_medium_labels
            color = RGB(0, 0, 0)
            try
                if haskey(io_settings["sankey_plot_spec"], label)
                    color_entry = Colors.color_names[io_settings["sankey_plot_spec"][label]]
                    color = RGB(color_entry[1] / 255, color_entry[2] / 255, color_entry[3] / 255)
                else
                    # use deterministic random color based on label hash. colors are drawn
                    # from the roma color scheme
                    rng = Random.MersenneTwister(hash(label))
                    color = get(ColorSchemes.roma, rand(rng))
                end
            catch
                @error "The given color '$(io_settings["sankey_plot_spec"][label])' of " *
                       "medium '$label' for the sankey plot is not one of the available " *
                       "colors in `Colors.color_names`"
                throw(InputError())
            end
            color_map[label] = color
        end

        for medium in unique_medium_labels
            if medium in keys(color_map)
                continue
            elseif medium == "hide_medium"
                color_map[medium] = parse(RGBA, "rgba(0,0,0,0)")
            else
                @error "The color for the medium '$medium' for the sankey could not be " *
                       "found in the input file. Please add the medium and its color in " *
                       "IO setting 'sankey_plot_spec'."
                throw(InputError())
            end
        end
    end
    colors_for_medium = map(x -> color_map[x], medium_labels)

    do_hide_real_demands = true
    if do_hide_real_demands
        # hide sf_fixed_sink and sf_fixed_source interfaces
        colors_for_medium_RGBA = Array{Any}(nothing, interface_new)
        for (idx, medium) in pairs(medium_labels)
            if medium == "hide_medium"
                colors_for_medium_RGBA[idx] = parse(RGBA, "rgba(0,0,0,0)")
            else
                colors_for_medium_RGBA[idx] = colors_for_medium[idx]
            end
        end

        # hide blocks and set position of sf_fixed_sink and sf_fixed_source blocks
        block_colors = Array{Any}(nothing, length(block_labels))
        for (idx, block) in pairs(block_labels)
            if last(block, 12) == "total_demand" || last(block, 12) == "total_supply"
                block_labels[idx] = ""
                block_colors[idx] = parse(RGBA, "rgba(0,0,0,0)")
            else
                block_colors[idx] = "blue"
            end
        end
    else
        colors_for_medium_RGBA = colors_for_medium
    end

    # create plot
    p = plot(sankey(;
                    node=attr(; pad=25,
                              thickness=20,
                              line=do_hide_real_demands ? attr(; color="white", width=0.0) : nothing,
                              label=block_labels,
                              color=do_hide_real_demands ? block_colors : "blue"),
                    link=attr(; source=output_all_source_num .- 1, # indices correspond to block_labels starting from index 0
                              target=output_all_target_num .- 1, # indices correspond to block_labels starting from index 0
                              value=output_all_value_sum,
                              label=medium_labels,
                              color=colors_for_medium_RGBA)),
             Layout(; title_text="Sankey diagram of system topology and energy flows",
                    font_size=14))

    # save plot
    file_path = sim_params["run_path"](io_settings["sankey_plot_file_path"])
    savefig(p, file_path)
end

"""
    aggregate_csv(input_path, output_path, output_keys, weather_data_keys,
                  energy_terms, cumulative_terms, zero_as_missing_terms,
                  time_columns, separator, threshold, fixed_output_precision)

Create a summary CSV from a time-series CSV result file.

The input CSV is expected to contain one header row and one column per logged
parameter. Time columns are ignored. The classification of component and flow
outputs is based on the corresponding `OutputKey.value_key` entries in
`output_keys`, not on the CSV column names. This avoids incorrect classifications
caused by component names, media names or user-defined UACs.

Columns matching `cumulative_terms` are treated as cumulative time series and
represented by their last valid value. Columns matching `energy_terms` are
summed. Weather data and all other numeric columns are treated as intensive
quantities and averaged.

Empty values, non-numeric values, `NaN` and infinite values are ignored. For
columns matching `zero_as_missing_terms`, zero values are ignored as well, e.g.
for COP-like quantities where zero means inactive or undefined.

Returns `true` if the summary CSV was created successfully and `false` otherwise.
"""
function aggregate_csv(input_path::AbstractString,
                       output_path::AbstractString,
                       output_keys::Union{Nothing,Vector{EnergySystems.OutputKey}},
                       weather_data_keys::Union{Nothing,Vector{String}},
                       energy_terms::Vector{String},
                       cumulative_terms::Vector{String},
                       zero_as_missing_terms::Vector{String},
                       time_columns::Vector{String},
                       separator::Char,
                       threshold::Float64,
                       fixed_output_precision::Int)::Bool
    if !isfile(input_path)
        @info "No summary CSV could be created, as the CSV result file could not be found at: $(input_path)"
        return false
    end

    # identify columns that should not be aggregated
    function is_time_column(output_key::AbstractString)::Bool
        col = strip(output_key)
        return col in time_columns || startswith(lowercase(col), "time")
    end

    # identify cumulative time series for which the last value is used
    function is_cumulative_column(output_key::AbstractString)::Bool
        return any(term -> occursin(term, output_key), cumulative_terms)
    end

    # identify extensive energy-related quantities that should be summed
    function is_energy_column(output_key::AbstractString)::Bool
        return any(term -> occursin(term, output_key), energy_terms)
    end

    # identify intensive quantities where zero represents an invalid value
    function is_zero_as_missing_column(output_key::AbstractString)::Bool
        return any(term -> occursin(term, output_key), zero_as_missing_terms)
    end

    # identify weather quantities that should be summed
    function is_weather_energy_column(weather_key::AbstractString)::Bool
        # hardcoded at this point as they will probably not change...
        weather_energy_terms = ["beamHorIrr",
                                "difHorIrr",
                                "globHorIrr",
                                "longWaveIrr"]
        return any(term -> occursin(term, weather_key), weather_energy_terms)
    end

    # parse numbers with German decimal comma and tolerate simple thousands separators
    function parse_decimal_number(value::AbstractString)::Union{Float64,Missing}
        s = strip(value)

        if isempty(s)
            return missing
        end

        s = replace(s, "\ufeff" => "")
        s = replace(s, "\u00a0" => "")
        s = replace(s, " " => "")

        x = tryparse(Float64, replace(s, "," => "."))
        if x !== nothing
            return x
        end

        if occursin(",", s)
            s2 = replace(s, "." => "")
            s2 = replace(s2, "," => ".")

            x = tryparse(Float64, s2)
            if x !== nothing
                return x
            end
        end

        return missing
    end

    # avoid very small numerical residuals in the summary output
    function clean_small_value(value::Float64)::Float64
        return abs(value) < threshold ? 0.0 : value
    end

    # write numbers with decimal comma for consistency with the ReSiE CSV format
    function decimal_string(value::Real)::String
        # apply fixed precision to the summary output
        if fixed_output_precision > 0
            s = string(round(Float64(value); sigdigits=fixed_output_precision))
        else
            s = string(Float64(value))
        end
        return replace(s, "." => ",")
    end

    # escape fields only if required by the CSV format
    function csv_escape(value::AbstractString)::String
        s = String(value)

        if occursin(string(separator), s) || occursin("\"", s) ||
           occursin("\n", s) || occursin("\r", s)
            return "\"" * replace(s, "\"" => "\"\"") * "\""
        end

        return s
    end

    function write_csv_row(io, row::AbstractVector{<:AbstractString})::Nothing
        println(io, join(csv_escape.(row), string(separator)))
        return nothing
    end

    csv = CSV.File(input_path;
                   delim=separator,
                   header=1,
                   normalizenames=false,
                   stringtype=String,
                   types=String)

    if length(csv) < 1
        @info "No summary CSV could be created, as the CSV result file contains no data rows: $(input_path)"
        return false
    end

    # Get column names in CSV order.
    csv_column_names = collect(propertynames(first(csv)))

    # Clean header names for output and time-column checks.
    headers = [String(strip(replace(String(name), "\ufeff" => "")))
               for name in csv_column_names]

    # create classification keys in the same order as the CSV columns
    classification_keys = String["Time"]

    if output_keys !== nothing
        append!(classification_keys, [outkey.value_key for outkey in output_keys])
    end

    if weather_data_keys !== nothing
        append!(classification_keys, ["Weather " * key for key in weather_data_keys])
    end

    if length(classification_keys) != length(headers)
        @info "No summary CSV could be created, as the number of output keys does not match the number of CSV columns."
        return false
    end

    # ignore time columns and process all remaining columns alphabetically
    data_indices = [j for j in eachindex(headers)
                    if !is_time_column(headers[j])]

    sort!(data_indices; by=j -> lowercase(headers[j]))

    open(output_path, "w") do io
        write_csv_row(io,
                      ["Parameter",
                       "Type",
                       "Aggregation",
                       "Values used",
                       "Value",
                       "Min",
                       "Max",
                       "Value kWh",
                       "Value MWh"])

        for j in data_indices
            column_name = headers[j]
            values = Float64[]
            classification_key = classification_keys[j]

            column = csv[csv_column_names[j]]

            for raw_value in column
                if raw_value === missing
                    continue
                end

                parsed = parse_decimal_number(String(raw_value))

                if parsed !== missing && isfinite(parsed)
                    if is_zero_as_missing_column(classification_key) && parsed == 0.0
                        continue
                    end

                    push!(values, parsed)
                end
            end
            if isempty(values)
                continue
            end

            n = length(values)

            if is_cumulative_column(classification_key)
                # cumulative columns are already integrated over time
                value_wh = clean_small_value(values[end])
                min_value = clean_small_value(minimum(values))
                max_value = clean_small_value(maximum(values))

                write_csv_row(io,
                              [column_name,
                               "cumulative",
                               "last",
                               string(n),
                               decimal_string(value_wh),
                               decimal_string(min_value),
                               decimal_string(max_value),
                               decimal_string(value_wh / 1_000),
                               decimal_string(value_wh / 1_000_000)])

            elseif is_energy_column(classification_key) || is_weather_energy_column(classification_key)
                # energy columns contain timestep values and are therefore summed
                value_wh = clean_small_value(sum(values))
                min_value = clean_small_value(minimum(values))
                max_value = clean_small_value(maximum(values))

                write_csv_row(io,
                              [column_name,
                               "energy",
                               "sum",
                               string(n),
                               decimal_string(value_wh),
                               decimal_string(min_value),
                               decimal_string(max_value),
                               decimal_string(value_wh / 1_000),
                               decimal_string(value_wh / 1_000_000)])

            else
                # remaining numeric columns are interpreted as intensive quantities
                mean_value = sum(values) / n
                min_value = minimum(values)
                max_value = maximum(values)
                aggregation = is_zero_as_missing_column(classification_key) ? "mean_excluding_zero" : "mean"

                write_csv_row(io,
                              [column_name,
                               "intensive",
                               aggregation,
                               string(n),
                               decimal_string(mean_value),
                               decimal_string(min_value),
                               decimal_string(max_value),
                               "",
                               ""])
            end
        end
    end

    return true
end

# ------------------------------------------------------------------------------------------
# Optimisation plotting
# ------------------------------------------------------------------------------------------

"""
    create_matrix_plot(results, io_settings, sim_params)

Create a scatter-plot matrix of all optimisation parameters. The selected objective is used
as marker colour and the best objective value is highlighted.
"""
function create_matrix_plot(results::Vector{Any},
                            io_settings::Dict{String,Any},
                            sim_params::Dict{String,Any})
    if io_settings["matrix_plot"] == "default"
        obj_key = "objective"
    elseif io_settings["matrix_plot"] == "custom"
        func, spec = first(pairs(io_settings["matrix_plot_spec"]))
        if func == "sum" || func == "mean"
            obj_key = func * " " * parse_outkeys(spec)[1]
        elseif func == "economic" || func == "emissions"
            obj_key = func * " " * spec[1]
        end
    end
    param_names = sim_params["optimisation"]["optim_params_keys"]
    results_dict = optimisation_results_dict(results)

    if !haskey(results_dict, obj_key)
        @error("Cannot create matrix plot: objective key \"$obj_key\" not found in results.")
    end

    missing_params = filter(param -> !haskey(results_dict, param), param_names)
    if !isempty(missing_params)
        error("Cannot create matrix plot. Missing parameters: $(join(missing_params, ", ")).")
    end

    objective = Union{Missing,Float64}[is_finite_number(value) ? Float64(value) : missing
                                       for value in results_dict[obj_key]]
    valid_objective = Float64[value for value in skipmissing(objective)]
    if isempty(valid_objective)
        @error("Cannot create matrix plot: no finite objective values.")
    end

    cmin, cmax = optimisation_color_bounds(valid_objective)
    min_objective = minimum(valid_objective)
    line_width = [value !== missing && value == min_objective ? 2 : 0
                  for value in objective]

    trace = splom(;
                  dimensions=[attr(; label=param, values=results_dict[param])
                              for param in param_names],
                  marker=attr(; color=objective,
                              colorscale=reverse(ColorSchemes.viridis),
                              cmin=cmin,
                              cmax=cmax,
                              showscale=true,
                              line=attr(; color="red", width=line_width),
                              colorbar=attr(; title=obj_key)),
                  text=string.(objective),
                  hovertemplate="%{x}, %{y}<br>$obj_key = %{text}<extra></extra>")

    p = plot(trace)
    file_path = optimisation_plot_path(sim_params, io_settings, "matrix_plot")
    savefig(p, file_path)

    return file_path
end

"""
    optimisation_results_dict(results)

Collect every result column in a dictionary of vectors. Missing entries are represented by
`missing`.
"""
function optimisation_results_dict(results::Vector{Any})::Dict{String,Vector{Any}}
    if isempty(results)
        @error("Cannot create optimisation plots: results are empty.")
    end

    all_keys = unique([String(key) for result in results for key in keys(result)])

    return Dict(key => [haskey(result, key) ? result[key] : missing for result in results]
                for key in all_keys)
end

"""
    is_finite_number(value)

Return `true` for finite real values.
"""
is_finite_number(value)::Bool = value isa Real && isfinite(Float64(value))

"""
    optimisation_plot_path(sim_params, io_settings, suffix)

Create a plot path derived from `optim_plots_file_path`.
"""
function optimisation_plot_path(sim_params::Dict{String,Any},
                                io_settings::Dict{String,Any},
                                suffix::String)::String
    base_path = sim_params["run_path"](io_settings["optim_plots_file_path"])
    dir, filename = splitdir(base_path)
    root, _ = splitext(filename)
    ext = ".html"  # html only here
    return joinpath(dir, "$(root)_$(suffix)$(ext)")
end

"""
    safe_plot_name(name)

Convert a result key to a filename-safe string.
"""
safe_plot_name(name::AbstractString)::String = replace(name, r"[^A-Za-z0-9_]+" => "_")

"""
    optimisation_color_bounds(values)

Return robust colour limits. The upper colour limit is based on the lower quartile to retain
contrast around good objective values.
"""
function optimisation_color_bounds(values::Vector{Float64})::Tuple{Float64,Float64}
    if isempty(values)
        @error("Cannot determine colour bounds from an empty vector.")
    end

    sorted_values = sort(values)
    cmin = first(sorted_values)
    cmax = sorted_values[max(1, cld(length(sorted_values), 4))]

    if !isfinite(cmax) || cmax <= cmin
        cmax = last(sorted_values)
    end

    if cmax <= cmin
        cmax = cmin + max(abs(cmin) * 1.0e-9, 1.0e-9)
    end

    return cmin, cmax
end

"""
    create_objective_convergence_plot(results, io_settings, sim_params)

Create a plot containing the objective value of every run and the best value found so far.
"""
function create_objective_convergence_plot(results::Vector{Any},
                                           io_settings::Dict{String,Any},
                                           sim_params::Dict{String,Any})
    obj_key = "objective"
    results_dict = optimisation_results_dict(results)

    objective = Union{Missing,Float64}[is_finite_number(value) ? Float64(value) : missing
                                       for value in results_dict[obj_key]]
    valid_objective = Float64[value for value in skipmissing(objective)]
    if isempty(valid_objective)
        @error("Cannot create convergence plot: no finite objective values.")
    end

    best_so_far = Vector{Union{Missing,Float64}}(missing, length(objective))
    current_best = Inf

    for idx in eachindex(objective)
        if objective[idx] !== missing
            current_best = min(current_best, objective[idx])
        end
        if isfinite(current_best)
            best_so_far[idx] = current_best
        end
    end

    use_log_axis = all(value -> value > 0.0, valid_objective)

    objective_trace = scatter(; x=eachindex(objective),
                              y=objective,
                              mode="markers",
                              name=obj_key,
                              text=["run $idx<br>$obj_key = $(objective[idx])"
                                    for idx in eachindex(objective)],
                              hovertemplate="%{text}<extra></extra>")

    best_trace = scatter(; x=eachindex(best_so_far),
                         y=best_so_far,
                         mode="lines",
                         name="Best so far",
                         hovertemplate="run %{x}<br>best = %{y}<extra></extra>")

    layout = Layout(; title="Optimisation convergence",
                    xaxis=attr(; title="Run"),
                    yaxis=attr(; title=obj_key,
                               type=use_log_axis ? "log" : "linear"))

    p = plot([objective_trace, best_trace], layout)
    file_path = optimisation_plot_path(sim_params, io_settings, "convergence")
    savefig(p, file_path)

    return file_path
end

"""
    create_objective_parameter_plots(results, io_settings, sim_params)

Create one objective-versus-parameter plot for every optimisation parameter.
"""
function create_objective_parameter_plots(results::Vector{Any},
                                          io_settings::Dict{String,Any},
                                          sim_params::Dict{String,Any})
    obj_key = "objective"
    param_names = sim_params["optimisation"]["optim_params_keys"]
    results_dict = optimisation_results_dict(results)

    objective = Union{Missing,Float64}[is_finite_number(value) ? Float64(value) : missing
                                       for value in results_dict[obj_key]]
    valid_objective = Float64[value for value in skipmissing(objective)]
    if isempty(valid_objective)
        @error("Cannot create parameter plots: no finite objective values.")
    end

    use_log_axis = all(value -> value > 0.0, valid_objective)
    saved_paths = String[]

    for param_name in param_names
        if !haskey(results_dict, param_name)
            @warn "Skipping objective-versus-parameter plot. Parameter \"$param_name\" was not found in results."
            continue
        end

        param_values = results_dict[param_name]
        valid_idx = [idx
                     for idx in eachindex(objective)
                     if objective[idx] !== missing && is_finite_number(param_values[idx])]

        if isempty(valid_idx)
            @warn "Skipping objective-versus-parameter plot for \"$param_name\": no valid points."
            continue
        end

        x_values = Float64[Float64(param_values[idx]) for idx in valid_idx]
        obj_values = Float64[objective[idx] for idx in valid_idx]
        best_local_idx = argmin(obj_values)
        best_global_idx = valid_idx[best_local_idx]
        cmin, cmax = optimisation_color_bounds(obj_values)

        runs_trace = scatter(; x=x_values,
                             y=obj_values,
                             mode="markers",
                             name="Runs",
                             marker=attr(; color=obj_values,
                                         colorscale=reverse(ColorSchemes.viridis),
                                         cmin=cmin,
                                         cmax=cmax,
                                         showscale=true,
                                         colorbar=attr(; title=obj_key)),
                             text=["run $(valid_idx[idx])<br>$param_name = $(x_values[idx])" *
                                   "<br>$obj_key = $(obj_values[idx])"
                                   for idx in eachindex(valid_idx)],
                             hovertemplate="%{text}<extra></extra>")

        best_trace = scatter(; x=[x_values[best_local_idx]],
                             y=[obj_values[best_local_idx]],
                             mode="markers",
                             name="Best",
                             marker=attr(; symbol="x", size=14),
                             text=["best run $best_global_idx<br>$param_name = $(x_values[best_local_idx])" *
                                   "<br>$obj_key = $(obj_values[best_local_idx])"],
                             hovertemplate="%{text}<extra></extra>")

        layout = Layout(; title="Objective vs $param_name",
                        xaxis=attr(; title=param_name),
                        yaxis=attr(; title=obj_key,
                                   type=use_log_axis ? "log" : "linear"))

        p = plot([runs_trace, best_trace], layout)
        file_path = optimisation_plot_path(sim_params,
                                           io_settings,
                                           "objective_vs_$(safe_plot_name(param_name))")
        savefig(p, file_path)
        push!(saved_paths, file_path)
    end

    return saved_paths
end

"""
    create_parallel_coordinates_plot(results, io_settings, sim_params;
                                                objective_keys=nothing)

Create a filterable parallel-coordinates plot. Optimisation parameters and result axes are
visually separated. A floating control panel permits precise axis zooming.
"""
function create_parallel_coordinates_plot(results::Vector{Any},
                                          io_settings::Dict{String,Any},
                                          sim_params::Dict{String,Any};
                                          objective_keys=nothing)
    obj_key = "objective"
    param_names = sim_params["optimisation"]["optim_params_keys"]
    results_dict = optimisation_results_dict(results)

    all_objective_keys = optimisation_objective_axis_keys(results,
                                                          sim_params,
                                                          obj_key;
                                                          objective_keys=objective_keys)
    individual_objective_keys = filter(!=(obj_key), all_objective_keys)
    is_multiobjective = !isempty(individual_objective_keys)
    required_keys = unique(vcat(param_names, all_objective_keys))
    missing_keys = filter(key -> !haskey(results_dict, key), required_keys)

    if !isempty(missing_keys)
        @error("Cannot create parallel-coordinates plot. Missing result keys: " * join(missing_keys, ", "))
    end

    valid_idx = [idx
                 for idx in eachindex(results)
                 if all(key -> is_finite_number(results_dict[key][idx]), required_keys)]
    if isempty(valid_idx)
        @error("Cannot create parallel-coordinates plot: no complete numeric result rows.")
    end

    objective_values = Dict(
        key => Float64[Float64(results_dict[key][idx])
                       for idx in valid_idx]
        for key in all_objective_keys)
    primary_values = objective_values[obj_key]
    dimensions = Any[]

    # parameter axes
    for param_name in param_names
        values = Float64[Float64(results_dict[param_name][idx]) for idx in valid_idx]
        push!(dimensions, attr(; label=param_name, values=values))
    end

    if is_multiobjective
        # individual objective axes
        for key in individual_objective_keys
            push!(dimensions, attr(; label=key, values=objective_values[key]))
        end

        if !all(value -> value > 0.0, primary_values)
            @error("Cannot display log10($obj_key): the combined objective contains non-positive values.")
        end

        color_values = log10.(primary_values)
        color_title = "log10($obj_key)"
        push!(dimensions, attr(; label=color_title, values=color_values))
    else
        color_values = primary_values
        color_title = obj_key
        push!(dimensions, attr(; label=obj_key, values=primary_values))
    end

    cmin, cmax = optimisation_color_bounds(color_values)
    n_parameters = length(param_names)
    n_dimensions = length(dimensions)
    n_results = n_dimensions - n_parameters

    # keep full outer labels visible and leave room for group labels
    plot_x_min = 0.035
    plot_x_max = 0.955
    plot_y_max = 0.915
    group_boundary = if n_dimensions > 1
        plot_x_min +
        ((n_parameters - 0.5) / (n_dimensions - 1)) * (plot_x_max - plot_x_min)
    else
        plot_x_max
    end

    trace = parcoords(; ids=string.(valid_idx),
                      domain=attr(; x=[plot_x_min, plot_x_max], y=[0.0, plot_y_max]),
                      labelfont=attr(; size=11),
                      line=attr(; color=color_values,
                                colorscale=reverse(ColorSchemes.viridis),
                                cmin=cmin,
                                cmax=cmax,
                                showscale=true,
                                colorbar=attr(; title=color_title,
                                              len=0.82,
                                              y=0.45,
                                              x=0.985,
                                              xanchor="left",
                                              xpad=2,
                                              thickness=13)),
                      unselected=attr(; line=attr(; color="rgb(145,145,145)",
                                                  opacity=0.20)),
                      dimensions=dimensions)

    group_shapes, group_annotations = parallel_group_decorations(n_parameters,
                                                                 n_results,
                                                                 plot_x_min,
                                                                 group_boundary,
                                                                 plot_x_max,
                                                                 plot_y_max)
    layout = Layout(;
                    title=attr(; text="Interactive Optimisation Design Space",
                               x=0.5,
                               xanchor="center",
                               y=0.997,
                               yanchor="top",
                               font=attr(; size=18)),
                    margin=attr(; t=42, b=24, l=28, r=58),
                    shapes=group_shapes,
                    annotations=group_annotations)

    p = plot(trace, layout)
    file_path = optimisation_plot_path(sim_params,
                                       io_settings,
                                       "parallel_coordinates")
    savefig(p, file_path)
    inject_parallel_axis_zoom_controls!(file_path)

    return file_path
end

"""
    parallel_group_decorations(n_parameters, n_results, plot_x_min, group_boundary,
                               plot_x_max, plot_y_max)

Create background regions, group labels and the separator for a parallel-coordinates plot.
"""
function parallel_group_decorations(n_parameters::Int,
                                    n_results::Int,
                                    plot_x_min::Float64,
                                    group_boundary::Float64,
                                    plot_x_max::Float64,
                                    plot_y_max::Float64)
    shapes = Any[]
    annotations = Any[]
    groups = [(enabled=n_parameters > 0,
               label="Variable parameters",
               x0=plot_x_min,
               x1=group_boundary,
               fill="rgba(70,130,180,0.07)",
               text_color="rgb(55,90,120)"),
              (enabled=n_results > 0,
               label="Results",
               x0=group_boundary,
               x1=plot_x_max,
               fill="rgba(220,140,50,0.07)",
               text_color="rgb(140,85,30)")]

    for group in groups
        group.enabled || continue

        push!(shapes,
              attr(; type="rect",
                   xref="paper",
                   yref="paper",
                   x0=group.x0,
                   x1=group.x1,
                   y0=0.0,
                   y1=plot_y_max,
                   fillcolor=group.fill,
                   line=attr(; width=0),
                   layer="below"))

        push!(annotations,
              attr(; text="<b>$(group.label)</b>",
                   x=(group.x0 + group.x1) / 2,
                   y=0.992,
                   xref="paper",
                   yref="paper",
                   showarrow=false,
                   xanchor="center",
                   yanchor="top",
                   font=attr(; size=11, color=group.text_color)))
    end

    if n_parameters > 0 && n_results > 0
        push!(shapes,
              attr(; type="line",
                   xref="paper",
                   yref="paper",
                   x0=group_boundary,
                   x1=group_boundary,
                   y0=0.0,
                   y1=plot_y_max,
                   line=attr(; color="rgba(80,80,80,0.65)",
                             width=2,
                             dash="dot"),
                   layer="above"))
    end

    return shapes, annotations
end

"""
    optimisation_objective_axis_keys(results, sim_params, primary_obj_key;
                                     objective_keys=nothing)

Return the result columns displayed as objective or KPI axes.
"""
function optimisation_objective_axis_keys(results::Vector{Any},
                                          sim_params::Dict{String,Any},
                                          primary_obj_key::String;
                                          objective_keys=nothing)::Vector{String}
    if isempty(results)
        @error("Cannot determine objective axes: results are empty.")
    end

    param_names = sim_params["optimisation"]["optim_params_keys"]

    if objective_keys !== nothing
        selected_keys = unique(String.(objective_keys))
        primary_obj_key in selected_keys || pushfirst!(selected_keys, primary_obj_key)

        missing_keys = [key for key in selected_keys
                        if !any(result -> haskey(result, key), results)]
        if !isempty(missing_keys)
            @error("Cannot create objective axes. Missing result keys: " * join(missing_keys, ", "))
        end
        return selected_keys
    end

    excluded_keys = Set(vcat(param_names, ["error", "run", "run_id", "sample_id"]))
    all_result_keys = unique([String(key) for result in results for key in keys(result)])
    detected_keys = String[]

    for key in all_result_keys
        key in excluded_keys && continue
        values = [haskey(result, key) ? result[key] : missing for result in results]
        any(is_finite_number, values) && push!(detected_keys, key)
    end

    filter!(key -> key != primary_obj_key, detected_keys)
    pushfirst!(detected_keys, primary_obj_key)

    return unique(detected_keys)
end

"""
    inject_parallel_axis_zoom_controls!(file_path)

Inject a floating axis-zoom panel into a saved parallel-coordinates HTML plot.
"""
function inject_parallel_axis_zoom_controls!(file_path::String)
    html = read(file_path, String)

    injection = raw"""
<script>
(function () {
    function findPlotlyDiv() {
        const divs = document.querySelectorAll(".js-plotly-plot");
        return divs.length > 0 ? divs[0] : null;
    }

    function setupControls() {
        const gd = findPlotlyDiv();

        if (gd === null || !gd.data || gd.data.length === 0) {
            window.setTimeout(setupControls, 200);
            return;
        }

        if (document.getElementById("parallel-axis-zoom-controls")) {
            return;
        }

        const traceIndex = gd.data.findIndex(
            trace => trace.type === "parcoords"
        );

        if (traceIndex < 0) {
            return;
        }

        const dims = gd.data[traceIndex].dimensions;

        const originalRanges = dims.map(
            dimension =>
                Array.isArray(dimension.range)
                    ? dimension.range.slice()
                    : null
        );

        document.documentElement.style.height = "100%";
        document.body.style.height = "100%";
        document.body.style.margin = "0";
        document.body.style.overflow = "hidden";

        gd.style.width = "100vw";
        gd.style.height = "100vh";

        Plotly.Plots.resize(gd);

        const panel = document.createElement("details");
        panel.id = "parallel-axis-zoom-controls";
        panel.open = false;

        panel.style.position = "fixed";
        panel.style.top = "10px";
        panel.style.right = "12px";
        panel.style.zIndex = "10000";
        panel.style.fontFamily = "Arial, sans-serif";
        panel.style.fontSize = "13px";
        panel.style.background = "rgba(255,255,255,0.96)";
        panel.style.border = "1px solid #bbb";
        panel.style.borderRadius = "6px";
        panel.style.boxShadow = "0 2px 8px rgba(0,0,0,0.15)";
        panel.style.maxWidth = "360px";

        const summary = document.createElement("summary");
        summary.textContent = "Axis zoom";
        summary.style.cursor = "pointer";
        summary.style.fontWeight = "600";
        summary.style.padding = "7px 10px";
        summary.style.userSelect = "none";

        panel.appendChild(summary);

        const content = document.createElement("div");
        content.style.padding = "2px 10px 10px 10px";
        content.style.display = "grid";
        content.style.gridTemplateColumns = "auto 1fr";
        content.style.gap = "7px";
        content.style.alignItems = "center";

        panel.appendChild(content);

        function addLabel(text) {
            const label = document.createElement("span");
            label.textContent = text;
            content.appendChild(label);
        }

        const axisSelect = document.createElement("select");
        axisSelect.style.width = "220px";
        axisSelect.style.maxWidth = "220px";

        dims.forEach((dimension, index) => {
            const option = document.createElement("option");
            option.value = index;
            option.textContent =
                dimension.label || ("Axis " + (index + 1));

            axisSelect.appendChild(option);
        });

        addLabel("Axis");
        content.appendChild(axisSelect);

        const minInput = document.createElement("input");
        minInput.type = "number";
        minInput.step = "any";
        minInput.style.width = "130px";

        addLabel("Minimum");
        content.appendChild(minInput);

        const maxInput = document.createElement("input");
        maxInput.type = "number";
        maxInput.step = "any";
        maxInput.style.width = "130px";

        addLabel("Maximum");
        content.appendChild(maxInput);

        const buttonRow = document.createElement("div");
        buttonRow.style.gridColumn = "1 / span 2";
        buttonRow.style.display = "flex";
        buttonRow.style.flexWrap = "wrap";
        buttonRow.style.gap = "6px";
        buttonRow.style.marginTop = "3px";

        content.appendChild(buttonRow);

        const applyButton = document.createElement("button");
        applyButton.textContent = "Apply";
        buttonRow.appendChild(applyButton);

        const resetButton = document.createElement("button");
        resetButton.textContent = "Reset axis";
        buttonRow.appendChild(resetButton);

        const resetAllButton = document.createElement("button");
        resetAllButton.textContent = "Reset all";
        buttonRow.appendChild(resetAllButton);

        const hint = document.createElement("div");
        hint.textContent =
            "Zoom first, then drag on an axis to filter.";
        hint.style.gridColumn = "1 / span 2";
        hint.style.color = "#666";
        hint.style.fontSize = "11px";
        hint.style.marginTop = "2px";

        content.appendChild(hint);

        document.body.appendChild(panel);

        function currentDimensionIndex() {
            return Number(axisSelect.value);
        }

        function finiteValues(index) {
            return gd.data[traceIndex]
                .dimensions[index]
                .values
                .filter(
                    value =>
                        typeof value === "number" &&
                        Number.isFinite(value)
                );
        }

        function formatInputValue(value) {
            return Number(value.toPrecision(8)).toString();
        }

        function updateInputValues() {
            const index = currentDimensionIndex();
            const dimension =
                gd.data[traceIndex].dimensions[index];

            const values = finiteValues(index);

            if (values.length === 0) {
                minInput.value = "";
                maxInput.value = "";
                return;
            }

            if (
                Array.isArray(dimension.range) &&
                dimension.range.length === 2 &&
                Number.isFinite(dimension.range[0]) &&
                Number.isFinite(dimension.range[1])
            ) {
                minInput.value =
                    formatInputValue(dimension.range[0]);

                maxInput.value =
                    formatInputValue(dimension.range[1]);
            } else {
                minInput.value =
                    formatInputValue(Math.min(...values));

                maxInput.value =
                    formatInputValue(Math.max(...values));
            }
        }

        axisSelect.addEventListener(
            "change",
            updateInputValues
        );

        applyButton.addEventListener("click", function () {
            const index = currentDimensionIndex();
            const minimum = Number(minInput.value);
            const maximum = Number(maxInput.value);

            if (
                !Number.isFinite(minimum) ||
                !Number.isFinite(maximum) ||
                minimum >= maximum
            ) {
                alert("Please enter a valid minimum and maximum.");
                return;
            }

            const newDimensions =
                gd.data[traceIndex].dimensions.map(
                    (dimension, dimensionIndex) => {
                        const copied =
                            Object.assign({}, dimension);

                        if (dimensionIndex === index) {
                            copied.range = [
                                minimum,
                                maximum
                            ];

                            delete copied.constraintrange;
                        }

                        return copied;
                    }
                );

            Plotly.restyle(
                gd,
                {
                    dimensions: [newDimensions]
                },
                [traceIndex]
            ).then(updateInputValues);
        });

        resetButton.addEventListener("click", function () {
            const index = currentDimensionIndex();

            const newDimensions =
                gd.data[traceIndex].dimensions.map(
                    (dimension, dimensionIndex) => {
                        const copied =
                            Object.assign({}, dimension);

                        if (dimensionIndex === index) {
                            if (originalRanges[index] === null) {
                                delete copied.range;
                            } else {
                                copied.range =
                                    originalRanges[index].slice();
                            }

                            delete copied.constraintrange;
                        }

                        return copied;
                    }
                );

            Plotly.restyle(
                gd,
                {
                    dimensions: [newDimensions]
                },
                [traceIndex]
            ).then(updateInputValues);
        });

        resetAllButton.addEventListener("click", function () {
            const newDimensions =
                gd.data[traceIndex].dimensions.map(
                    (dimension, index) => {
                        const copied =
                            Object.assign({}, dimension);

                        if (originalRanges[index] === null) {
                            delete copied.range;
                        } else {
                            copied.range =
                                originalRanges[index].slice();
                        }

                        delete copied.constraintrange;

                        return copied;
                    }
                );

            Plotly.restyle(
                gd,
                {
                    dimensions: [newDimensions]
                },
                [traceIndex]
            ).then(updateInputValues);
        });

        window.addEventListener("resize", function () {
            Plotly.Plots.resize(gd);
        });

        updateInputValues();
    }

    if (document.readyState === "loading") {
        document.addEventListener(
            "DOMContentLoaded",
            setupControls
        );
    } else {
        setupControls();
    }
})();
</script>
"""

    if occursin("</body>", html)
        html = replace(html,
                       "</body>" => injection * "\n</body>")
    else
        html *= injection
    end

    write(file_path, html)

    return file_path
end

"""
    create_interactive_3d_optimisation_plot(results, io_settings, sim_params;
                                            objective_keys=nothing, x_key=nothing,
                                            y_key=nothing, z_key=nothing, color_key=nothing,
                                            objective_sense=:min)

Create an interactive 3D optimisation explorer with selectable axes, colour variable and
independent axis zoom controls.
"""
function create_interactive_3d_optimisation_plot(results::Vector{Any},
                                                 io_settings::Dict{String,Any},
                                                 sim_params::Dict{String,Any};
                                                 objective_keys=nothing,
                                                 x_key=nothing,
                                                 y_key=nothing,
                                                 z_key=nothing,
                                                 color_key=nothing,
                                                 objective_sense::Symbol=:min)
    if !(objective_sense in (:min, :max))
        @error("objective_sense must be :min or :max.")
    end

    obj_key = "objective"
    param_names = sim_params["optimisation"]["optim_params_keys"]
    results_dict = optimisation_results_dict(results)

    all_objective_keys = optimisation_objective_axis_keys(results,
                                                          sim_params,
                                                          obj_key;
                                                          objective_keys=objective_keys)
    individual_objective_keys = filter(!=(obj_key), all_objective_keys)
    is_multiobjective = !isempty(individual_objective_keys)
    source_keys = unique(vcat(param_names, all_objective_keys))
    missing_keys = filter(key -> !haskey(results_dict, key), source_keys)

    if !isempty(missing_keys)
        @error("Cannot create 3D plot. Missing result keys: $(join(missing_keys, ", ")).")
    end

    valid_idx = [idx for idx in eachindex(results)
                 if all(key -> is_finite_number(results_dict[key][idx]), source_keys)]
    if isempty(valid_idx)
        @error("Cannot create 3D plot: no complete numeric result rows.")
    end

    source_data = Dict(key => Float64[Float64(results_dict[key][idx])
                                      for idx in valid_idx]
                       for key in source_keys)
    primary_values = source_data[obj_key]
    axis_data = Dict{String,Vector{Float64}}()
    axis_keys = String[]

    # parameter axes
    for key in param_names
        axis_data[key] = source_data[key]
        push!(axis_keys, key)
    end

    if is_multiobjective
        # individual objectives and combined objective
        for key in individual_objective_keys
            axis_data[key] = source_data[key]
            push!(axis_keys, key)
        end

        if !all(value -> value > 0.0, primary_values)
            @error("Cannot create log10($obj_key): the objective contains non-positive values.")
        end

        default_objective_key = "log10($obj_key)"
        axis_data[default_objective_key] = log10.(primary_values)
        push!(axis_keys, default_objective_key)
    else
        default_objective_key = obj_key
        axis_data[obj_key] = primary_values
        push!(axis_keys, obj_key)
    end

    unique!(axis_keys)
    if length(axis_keys) < 3
        @error("The 3D plot requires at least three selectable numeric quantities.")
    end

    function select_axis_key(explicit_key,
                             preferred_keys::Vector{String},
                             used_keys::Vector{String})::String
        if explicit_key !== nothing
            selected = String(explicit_key)
            selected in axis_keys || @error("Axis key \"$selected\" is not available.")
            selected in used_keys &&
                @error("The x, y and z axes must use different keys.")
            return selected
        end

        for candidate in unique(vcat(preferred_keys, axis_keys))
            if candidate in axis_keys && !(candidate in used_keys)
                return candidate
            end
        end

        @error("Cannot find a distinct axis key.")
    end

    x_selected = select_axis_key(x_key, param_names, String[])
    y_preferred = length(param_names) >= 2 ? param_names[2:end] : axis_keys
    y_selected = select_axis_key(y_key, y_preferred, [x_selected])
    z_selected = select_axis_key(z_key,
                                 [default_objective_key],
                                 [x_selected, y_selected])
    color_selected = color_key === nothing ? default_objective_key : String(color_key)
    if !(color_selected in axis_keys)
        !error("Colour key \"$color_selected\" is not available.")
    end

    best_local_idx = objective_sense == :min ? argmin(primary_values) : argmax(primary_values)
    color_values = axis_data[color_selected]
    cmin = minimum(color_values)
    cmax = maximum(color_values)
    if cmax <= cmin
        cmax = cmin + max(abs(cmin), 1.0) * 1.0e-9
    end

    hover_text = ["run $(valid_idx[idx])" *
                  "<br>$x_selected = $(axis_data[x_selected][idx])" *
                  "<br>$y_selected = $(axis_data[y_selected][idx])" *
                  "<br>$z_selected = $(axis_data[z_selected][idx])" *
                  "<br>$color_selected = $(axis_data[color_selected][idx])"
                  for idx in eachindex(valid_idx)]

    runs_trace = PlotlyJS.scatter(; type="scatter3d",
                                  x=axis_data[x_selected],
                                  y=axis_data[y_selected],
                                  z=axis_data[z_selected],
                                  mode="markers",
                                  name="Optimisation runs",
                                  marker=attr(; size=5,
                                              opacity=0.75,
                                              color=color_values,
                                              colorscale=reverse(ColorSchemes.viridis),
                                              cmin=cmin,
                                              cmax=cmax,
                                              showscale=true,
                                              colorbar=attr(; title=color_selected,
                                                            thickness=16),
                                              line=attr(; width=0.3,
                                                        color="rgba(50,50,50,0.35)")),
                                  text=hover_text,
                                  hovertemplate="%{text}<extra></extra>")

    best_trace = PlotlyJS.scatter(; type="scatter3d",
                                  x=[axis_data[x_selected][best_local_idx]],
                                  y=[axis_data[y_selected][best_local_idx]],
                                  z=[axis_data[z_selected][best_local_idx]],
                                  mode="markers",
                                  name="Best objective",
                                  marker=attr(; size=9,
                                              symbol="diamond",
                                              color="black",
                                              line=attr(; color="white", width=1.5)),
                                  text=[hover_text[best_local_idx]],
                                  hovertemplate="%{text}<extra></extra>")

    layout = Layout(; title=attr(; text="Interactive 3D Optimisation Explorer",
                                 x=0.5,
                                 xanchor="center"),
                    scene=attr(; xaxis=attr(; title=x_selected, autorange=true),
                               yaxis=attr(; title=y_selected, autorange=true),
                               zaxis=attr(; title=z_selected, autorange=true),
                               aspectmode="cube",
                               camera=attr(; eye=attr(; x=1.35, y=1.35, z=1.10))),
                    margin=attr(; t=60, b=20, l=20, r=90),
                    height=720)

    p = plot([runs_trace, best_trace], layout)
    file_path = optimisation_plot_path(sim_params, io_settings, "interactive_3d")
    if lowercase(splitext(file_path)[2]) !== ".html"
        @error("The interactive 3D plot requires an HTML output path.")
    end

    mkpath(dirname(file_path))
    savefig(p, file_path)
    inject_3d_axis_selection_controls!(file_path,
                                       axis_data,
                                       axis_keys,
                                       valid_idx,
                                       best_local_idx,
                                       x_selected,
                                       y_selected,
                                       z_selected,
                                       color_selected)

    return file_path
end

"""
    inject_3d_axis_selection_controls!(file_path, axis_data, axis_keys, run_ids,
                                       best_local_idx, initial_x, initial_y, initial_z,
                                       initial_color)

Inject axis selection, colour selection, camera reset and axis-zoom controls into a saved
3D optimisation HTML plot.
"""
function inject_3d_axis_selection_controls!(file_path::String,
                                            axis_data::Dict{String,Vector{Float64}},
                                            axis_keys::Vector{String},
                                            run_ids::Vector{Int},
                                            best_local_idx::Int,
                                            initial_x::String,
                                            initial_y::String,
                                            initial_z::String,
                                            initial_color::String)
    html = read(file_path, String)

    json_for_html(value) = replace(JSON.json(value),
                                   "</" => "<\\/")

    data_json = json_for_html(axis_data)
    keys_json = json_for_html(axis_keys)
    runs_json = json_for_html(run_ids)

    x_json = json_for_html(initial_x)
    y_json = json_for_html(initial_y)
    z_json = json_for_html(initial_z)
    color_json = json_for_html(initial_color)

    best_index = best_local_idx - 1

    injection = """
<script>
(function () {
    const plotData = $data_json;
    const axisKeys = $keys_json;
    const runIds = $runs_json;
    const bestIndex = $best_index;

    const initialSelection = {
        x: $x_json,
        y: $y_json,
        z: $z_json,
        color: $color_json
    };

    function findPlotlyDiv() {
        const divs = document.querySelectorAll(".js-plotly-plot");
        return divs.length > 0 ? divs[0] : null;
    }

    function setupControls() {
        const gd = findPlotlyDiv();

        if (gd === null || !gd.data || gd.data.length === 0) {
            window.setTimeout(setupControls, 200);
            return;
        }

        if (document.getElementById("optimisation-3d-controls")) {
            return;
        }

        const runTraceIndex = gd.data.findIndex(
            trace =>
                trace.type === "scatter3d" &&
                trace.name === "Optimisation runs"
        );

        const bestTraceIndex = gd.data.findIndex(
            trace =>
                trace.type === "scatter3d" &&
                trace.name === "Best objective"
        );

        if (runTraceIndex < 0 || bestTraceIndex < 0) {
            return;
        }

        const initialCamera =
            gd.layout.scene && gd.layout.scene.camera
                ? JSON.parse(JSON.stringify(gd.layout.scene.camera))
                : {
                    eye: {
                        x: 1.35,
                        y: 1.35,
                        z: 1.10
                    }
                };

        const manualRanges = {
            x: null,
            y: null,
            z: null
        };

        const controls = document.createElement("div");
        controls.id = "optimisation-3d-controls";
        controls.style.fontFamily = "Arial, sans-serif";
        controls.style.fontSize = "13px";
        controls.style.margin = "8px 0 4px 0";
        controls.style.padding = "10px";
        controls.style.border = "1px solid #ccc";
        controls.style.borderRadius = "5px";
        controls.style.display = "flex";
        controls.style.alignItems = "center";
        controls.style.flexWrap = "wrap";
        controls.style.gap = "10px";

        function createSelect(labelText, values, initialValue) {
            const container = document.createElement("label");
            container.style.display = "flex";
            container.style.alignItems = "center";
            container.style.gap = "5px";

            const label = document.createElement("span");
            label.textContent = labelText;

            const select = document.createElement("select");
            select.style.maxWidth = "280px";

            values.forEach(value => {
                const option = document.createElement("option");
                option.value = value;
                option.textContent = value;
                select.appendChild(option);
            });

            select.value = initialValue;

            container.appendChild(label);
            container.appendChild(select);
            controls.appendChild(container);

            return select;
        }

        const xSelect = createSelect(
            "X axis:",
            axisKeys,
            initialSelection.x
        );

        const ySelect = createSelect(
            "Y axis:",
            axisKeys,
            initialSelection.y
        );

        const zSelect = createSelect(
            "Z axis:",
            axisKeys,
            initialSelection.z
        );

        const colorSelect = createSelect(
            "Color:",
            axisKeys,
            initialSelection.color
        );

        const resetSelectionButton = document.createElement("button");
        resetSelectionButton.textContent = "Reset selection";
        controls.appendChild(resetSelectionButton);

        const resetViewButton = document.createElement("button");
        resetViewButton.textContent = "Reset view";
        controls.appendChild(resetViewButton);

        const separator = document.createElement("span");
        separator.style.width = "1px";
        separator.style.height = "28px";
        separator.style.background = "#bbb";
        separator.style.margin = "0 4px";
        controls.appendChild(separator);

        const zoomLabel = document.createElement("span");
        zoomLabel.textContent = "Axis zoom:";
        controls.appendChild(zoomLabel);

        const zoomAxisSelect = document.createElement("select");

        [
            ["x", "X axis"],
            ["y", "Y axis"],
            ["z", "Z axis"]
        ].forEach(entry => {
            const option = document.createElement("option");
            option.value = entry[0];
            option.textContent = entry[1];
            zoomAxisSelect.appendChild(option);
        });

        controls.appendChild(zoomAxisSelect);

        const zoomMinInput = document.createElement("input");
        zoomMinInput.type = "number";
        zoomMinInput.step = "any";
        zoomMinInput.style.width = "120px";
        controls.appendChild(zoomMinInput);

        const zoomMaxInput = document.createElement("input");
        zoomMaxInput.type = "number";
        zoomMaxInput.step = "any";
        zoomMaxInput.style.width = "120px";
        controls.appendChild(zoomMaxInput);

        const applyZoomButton = document.createElement("button");
        applyZoomButton.textContent = "Apply zoom";
        controls.appendChild(applyZoomButton);

        const resetZoomButton = document.createElement("button");
        resetZoomButton.textContent = "Reset axis";
        controls.appendChild(resetZoomButton);

        gd.parentNode.insertBefore(controls, gd);

        function formatValue(value) {
            if (!Number.isFinite(value)) {
                return String(value);
            }

            const absolute = Math.abs(value);

            if (
                absolute !== 0 &&
                (absolute >= 100000 || absolute < 0.001)
            ) {
                return value.toExponential(4);
            }

            return Number(value.toPrecision(6)).toString();
        }

        function inputValue(value) {
            return Number(value.toPrecision(8)).toString();
        }

        function colorBounds(values) {
            const finite = values.filter(Number.isFinite);

            let minimum = Math.min(...finite);
            let maximum = Math.max(...finite);

            if (!(maximum > minimum)) {
                const delta =
                    Math.max(Math.abs(minimum), 1.0) * 1.0e-9;

                maximum = minimum + delta;
            }

            return [minimum, maximum];
        }

        function selectedKeyForAxis(axis) {
            if (axis === "x") {
                return xSelect.value;
            }

            if (axis === "y") {
                return ySelect.value;
            }

            return zSelect.value;
        }

        function dataRangeForAxis(axis) {
            const key = selectedKeyForAxis(axis);

            const values = plotData[key].filter(
                value =>
                    typeof value === "number" &&
                    Number.isFinite(value)
            );

            if (values.length === 0) {
                return null;
            }

            let minimum = Math.min(...values);
            let maximum = Math.max(...values);

            if (!(maximum > minimum)) {
                const delta =
                    Math.max(Math.abs(minimum), 1.0) * 1.0e-9;

                minimum -= delta;
                maximum += delta;
            }

            return [minimum, maximum];
        }

        function updateAxisZoomInputs() {
            const axis = zoomAxisSelect.value;

            const range =
                manualRanges[axis] !== null
                    ? manualRanges[axis]
                    : dataRangeForAxis(axis);

            if (range === null) {
                zoomMinInput.value = "";
                zoomMaxInput.value = "";
                return;
            }

            zoomMinInput.value = inputValue(range[0]);
            zoomMaxInput.value = inputValue(range[1]);
        }

        function hoverText(xKey, yKey, zKey, colorKey) {
            return plotData[xKey].map((value, index) =>
                "run " + runIds[index] +
                "<br>" + xKey + " = " +
                formatValue(plotData[xKey][index]) +
                "<br>" + yKey + " = " +
                formatValue(plotData[yKey][index]) +
                "<br>" + zKey + " = " +
                formatValue(plotData[zKey][index]) +
                "<br>" + colorKey + " = " +
                formatValue(plotData[colorKey][index])
            );
        }

        function updateDisabledOptions() {
            const selected = [
                xSelect.value,
                ySelect.value,
                zSelect.value
            ];

            [xSelect, ySelect, zSelect].forEach(select => {
                Array.from(select.options).forEach(option => {
                    option.disabled =
                        option.value !== select.value &&
                        selected.includes(option.value);
                });
            });
        }

        function updatePlot(resetAxisRanges) {
            const xKey = xSelect.value;
            const yKey = ySelect.value;
            const zKey = zSelect.value;
            const colorKey = colorSelect.value;

            const text = hoverText(
                xKey,
                yKey,
                zKey,
                colorKey
            );

            const bounds = colorBounds(plotData[colorKey]);

            const currentCamera =
                gd.layout.scene && gd.layout.scene.camera
                    ? JSON.parse(
                        JSON.stringify(gd.layout.scene.camera)
                    )
                    : initialCamera;

            const runUpdate = {
                x: [plotData[xKey]],
                y: [plotData[yKey]],
                z: [plotData[zKey]],
                text: [text],
                "marker.color": [plotData[colorKey]],
                "marker.cmin": bounds[0],
                "marker.cmax": bounds[1],
                "marker.colorbar.title.text": colorKey
            };

            const bestUpdate = {
                x: [[plotData[xKey][bestIndex]]],
                y: [[plotData[yKey][bestIndex]]],
                z: [[plotData[zKey][bestIndex]]],
                text: [[text[bestIndex]]]
            };

            if (resetAxisRanges) {
                manualRanges.x = null;
                manualRanges.y = null;
                manualRanges.z = null;
            }

            Promise.all([
                Plotly.restyle(
                    gd,
                    runUpdate,
                    [runTraceIndex]
                ),

                Plotly.restyle(
                    gd,
                    bestUpdate,
                    [bestTraceIndex]
                )
            ]).then(function () {
                const layoutUpdate = {
                    "scene.xaxis.title.text": xKey,
                    "scene.yaxis.title.text": yKey,
                    "scene.zaxis.title.text": zKey,
                    "scene.aspectmode": "cube",
                    "scene.camera": currentCamera
                };

                if (resetAxisRanges) {
                    layoutUpdate["scene.xaxis.range"] = null;
                    layoutUpdate["scene.yaxis.range"] = null;
                    layoutUpdate["scene.zaxis.range"] = null;

                    layoutUpdate["scene.xaxis.autorange"] = true;
                    layoutUpdate["scene.yaxis.autorange"] = true;
                    layoutUpdate["scene.zaxis.autorange"] = true;
                }

                return Plotly.relayout(
                    gd,
                    layoutUpdate
                );
            }).then(function () {
                updateAxisZoomInputs();
            });

            updateDisabledOptions();
        }

        xSelect.addEventListener("change", function () {
            updatePlot(true);
        });

        ySelect.addEventListener("change", function () {
            updatePlot(true);
        });

        zSelect.addEventListener("change", function () {
            updatePlot(true);
        });

        colorSelect.addEventListener("change", function () {
            updatePlot(false);
        });

        zoomAxisSelect.addEventListener(
            "change",
            updateAxisZoomInputs
        );

        applyZoomButton.addEventListener("click", function () {
            const axis = zoomAxisSelect.value;
            const minimum = Number(zoomMinInput.value);
            const maximum = Number(zoomMaxInput.value);

            if (
                !Number.isFinite(minimum) ||
                !Number.isFinite(maximum) ||
                minimum >= maximum
            ) {
                alert(
                    "Please enter a valid minimum and maximum."
                );
                return;
            }

            manualRanges[axis] = [
                minimum,
                maximum
            ];

            const update = {};

            update[
                "scene." + axis + "axis.range"
            ] = [
                minimum,
                maximum
            ];

            update[
                "scene." + axis + "axis.autorange"
            ] = false;

            Plotly.relayout(gd, update);
        });

        resetZoomButton.addEventListener("click", function () {
            const axis = zoomAxisSelect.value;

            manualRanges[axis] = null;

            const update = {};

            update[
                "scene." + axis + "axis.range"
            ] = null;

            update[
                "scene." + axis + "axis.autorange"
            ] = true;

            Plotly.relayout(
                gd,
                update
            ).then(function () {
                updateAxisZoomInputs();
            });
        });

        resetSelectionButton.addEventListener(
            "click",
            function () {
                xSelect.value = initialSelection.x;
                ySelect.value = initialSelection.y;
                zSelect.value = initialSelection.z;
                colorSelect.value = initialSelection.color;

                updatePlot(true);
            }
        );

        resetViewButton.addEventListener(
            "click",
            function () {
                Plotly.relayout(gd, {
                    "scene.camera": initialCamera
                });
            }
        );

        updateDisabledOptions();
        updateAxisZoomInputs();
    }

    if (document.readyState === "loading") {
        document.addEventListener(
            "DOMContentLoaded",
            setupControls
        );
    } else {
        setupControls();
    }
})();
</script>
"""

    if occursin("</body>", html)
        html = replace(html,
                       "</body>" => injection * "\n</body>")
    else
        html *= injection
    end

    write(file_path, html)

    return file_path
end

"""
    create_optimisation_diagnostic_plots(results, io_settings, sim_params)

Create all optimisation diagnostic plots and return their output paths.
"""
function create_optimisation_diagnostic_plots(results::Vector{Any},
                                              io_settings::Dict{String,Any},
                                              sim_params::Dict{String,Any})
    matrix = create_matrix_plot(results, io_settings, sim_params)
    convergence = create_objective_convergence_plot(results, io_settings, sim_params)
    parameter_plots = create_objective_parameter_plots(results, io_settings, sim_params)
    parallel_coordinates = create_parallel_coordinates_plot(results,
                                                            io_settings,
                                                            sim_params)
    interactive_3d = create_interactive_3d_optimisation_plot(results,
                                                             io_settings,
                                                             sim_params)

    return (; matrix,
            convergence,
            parameter_plots,
            parallel_coordinates,
            interactive_3d)
end

# TODO remove
function save_to_prf(timestamps::Array{Int,1}, values::Array{Float64,1}, filepath::String)
    header_variables = ["# data_type:", "# time_definition:", "# profile_start_date:",
                        "# profile_start_date_format:", "# timestamp_format:",
                        "# interpolation_type:"]
    header_values = ["intensive", "startdate_timestamp", "01.01.2024 00:00",
                     "dd.mm.yyyy HH:MM", "seconds", "stepwise"]
    open(filepath, "w") do file_handle
        for (var, val) in zip(header_variables, header_values)
            write(file_handle, var * "\t" * val * "\n")
        end
        for (ts, val) in zip(timestamps, values)
            write(file_handle, string(ts) * ";" * string(val) * "\n")
        end
    end

    @info "Profile file at $filepath created."
end

# TODO remove
function save_to_prf(dates::Array{DateTime,1}, values::Array{Float64,1}, filepath::String)
    header_variables = ["# data_type:", "# time_definition:", "# timestamp_format:",
                        "# time_zone:", "# interpolation_type:"]
    header_values = ["intensive", "datestamp", "dd.mm.yyyy HH:MM",
                     "Europe/Berlin", "stepwise"]
    open(filepath, "w") do file_handle
        for (var, val) in zip(header_variables, header_values)
            write(file_handle, var * "\t" * val * "\n")
        end
        for (ts, val) in zip(dates, values)
            write(file_handle, string(ts) * ";" * string(val) * "\n")
        end
    end

    @info "Profile file at $filepath created."
end
