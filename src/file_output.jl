# this file contains functionality for writing output of the simulation to files.
using Dates

"""
get_output_keys(config[io_settings], components)

This function determines both for the lineplot and the csv output if
they should be created:
    - if not, "nothing" will be returned as key list
    - if yes, the key lists of the outputs will be returned, either containing
        - all possible keys (in-/excluding flows) if this is requested in the input file or
        - only the requested keys as requested in the input file
"""
function get_output_keys(io_settings::AbstractDict{String,Any},
                         components::Grouping)::Tuple{Union{Nothing,Vector{EnergySystems.OutputKey}},
                                                      Union{Nothing,Vector{EnergySystems.OutputKey}}}
    # determine if lineplot and csv should be created
    do_plot_all_outputs_excl_flows = false
    do_plot_all_outputs_incl_flows = false
    do_create_plot = false
    if haskey(io_settings, "output_plot")
        do_create_plot = true  # if the key exists, set do_create_plot to true by default
        if io_settings["output_plot"] == "all_excl_flows"
            do_plot_all_outputs_excl_flows = true
        elseif io_settings["output_plot"] == "all_incl_flows"
            do_plot_all_outputs_incl_flows = true
        elseif io_settings["output_plot"] == "nothing"
            do_create_plot = false
        elseif io_settings["output_plot"] == "all"
            @error "For \"output_plot\", the input \"all\" is no longer supported. Use \"all_incl_flows\" or \"all_excl_flows\"."
            throw(InputError)
        end
    end

    do_write_all_CSV_outputs_excl_flows = false
    do_write_all_CSV_outputs_incl_flows = false
    do_write_CSV = false
    if haskey(io_settings, "csv_output_keys")
        do_write_CSV = true  # if the key exists, set do_create_plot to true by default
        if io_settings["csv_output_keys"] == "all_excl_flows"
            do_write_all_CSV_outputs_excl_flows = true
        elseif io_settings["csv_output_keys"] == "all_incl_flows"
            do_write_all_CSV_outputs_incl_flows = true
        elseif io_settings["csv_output_keys"] == "nothing"
            do_write_CSV = false
        elseif io_settings["csv_output_keys"] == "all"
            @error "For \"csv_output_keys\", the input \"all\" is no longer supported. Use \"all_incl_flows\" or \"all_excl_flows\"."
            throw(InputError)
        end
    end

    plot_all_excluding_flows = do_plot_all_outputs_excl_flows || do_write_all_CSV_outputs_excl_flows
    plot_all_including_flows = do_write_all_CSV_outputs_incl_flows || do_plot_all_outputs_incl_flows

    # collect all possible outputs of all units if needed
    if plot_all_excluding_flows || plot_all_including_flows
        if plot_all_excluding_flows
            all_output_keys_excl_flows = Vector{EnergySystems.OutputKey}()
        end
        if plot_all_including_flows
            all_output_keys_incl_flows = Vector{EnergySystems.OutputKey}()
        end
        for unit in components
            output_vals = output_values(unit[2])
            if plot_all_including_flows
                temp_dict_incl_flows = Dict{String,Any}()
                for output_val in output_vals
                    if startswith(output_val, "EnergyFlow")
                        key = String(unit[2].medium)
                        if haskey(temp_dict_incl_flows, key)
                            push!(temp_dict_incl_flows[key], output_val[12:end])
                        else
                            temp_dict_incl_flows[key] = [output_val[12:end]]
                        end
                    elseif startswith(output_val, "TemperatureFlow")
                        # do nothing, temperature are added in output_keys()
                        continue
                    else
                        key = unit[2].uac
                        if haskey(temp_dict_incl_flows, key)
                            push!(temp_dict_incl_flows[key], output_val)
                        else
                            temp_dict_incl_flows[key] = [output_val]
                        end
                    end
                end
                append!(all_output_keys_incl_flows, output_keys(components, temp_dict_incl_flows))
            end
            if plot_all_excluding_flows
                temp_dict_excl_flows = Dict{String,Any}()
                for output_val in output_vals
                    if startswith(output_val, "EnergyFlow")
                        continue
                    elseif startswith(output_val, "TemperatureFlow")
                        continue
                    else
                        key = unit[2].uac
                        if haskey(temp_dict_excl_flows, key)
                            push!(temp_dict_excl_flows[key], output_val)
                        else
                            temp_dict_excl_flows[key] = [output_val]
                        end
                    end
                end
                append!(all_output_keys_excl_flows, output_keys(components, temp_dict_excl_flows))
            end
        end
        # sort keys
        function sort_by(output_key)
            any(startswith(output_key.value_key, p) for p in ("EnergyFlow", "TemperatureFlow")) ?
            # flows after others: primary = medium (always present), secondary
            lowercase("zzzzzzzzzzzzz" * string(output_key.medium) * output_key.value_key) :
            # default: primary = unit.uac, secondary = medium (if present), tertiary
            lowercase(string(output_key.unit.uac) * string(something(output_key.medium, "")) * output_key.value_key)
        end
        if plot_all_excluding_flows
            sort!(all_output_keys_excl_flows; by=sort_by)
        end
        if plot_all_including_flows
            sort!(all_output_keys_incl_flows; by=sort_by)
        end
    end

    # collect output keys for lineplot and csv output
    if do_create_plot  # line plot
        if do_plot_all_outputs_incl_flows
            output_keys_lineplot = all_output_keys_incl_flows
        elseif do_plot_all_outputs_excl_flows
            output_keys_lineplot = all_output_keys_excl_flows
        else
            output_keys_lineplot = Vector{EnergySystems.OutputKey}()
            for plot in io_settings["output_plot"]
                key = plot[2]["key"]
                append!(output_keys_lineplot, output_keys(components, key))
            end
        end
    else
        output_keys_lineplot = nothing
    end

    if do_write_CSV  # csv export
        if do_write_all_CSV_outputs_incl_flows
            output_keys_to_csv = all_output_keys_incl_flows
        elseif do_write_all_CSV_outputs_excl_flows
            output_keys_to_csv = all_output_keys_excl_flows
        else  # get only requested output keys from input file for CSV-export
            output_keys_to_csv = output_keys(components, io_settings["csv_output_keys"])
        end
    else
        output_keys_to_csv = nothing
    end

    return output_keys_lineplot, output_keys_to_csv
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
                if isdefined(each_outputinterface.target, :medium)
                    push!(medium_of_interfaces, each_outputinterface.target.medium)
                elseif isdefined(each_outputinterface.source, :medium)
                    push!(medium_of_interfaces, each_outputinterface.source.medium)
                elseif !(medium === nothing)
                    push!(medium_of_interfaces, medium)
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

    all_current_media = String.(unique([bus.medium
                                        for bus in values(components) if bus.sys_function === EnergySystems.sf_bus]))

    for key in keys(from_config)
        if key in keys(components)
            unit_key = key
            unit = components[unit_key]

            for entry in from_config[unit_key]
                splitted = split(String(entry), ":")
                if length(splitted) > 1
                    medium_key = splitted[1]
                    medium = Symbol(String(medium_key))
                    if medium in EnergySystems.medium_categories
                        value_key = splitted[2]
                    else
                        @error "In unit \"$(unit.uac)\", the given output key \"$entry\" could not be mapped to an " *
                               "output key. Make sure that the medium \"$(String(medium_key))\" exists in the current " *
                               "component and that you have used \":\" as separator without any extra spaces."
                        throw(InputError)
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
            for value_key in from_config[media_key]
                success = false
                in_uac, out_uac = split(value_key, "->")
                for bus in [unit for unit in values(components) if unit.sys_function === EnergySystems.sf_bus]
                    # consider only proxy busses or busses without proxies and busses with correct media
                    if bus.proxy === nothing && bus.medium == medium
                        # check if input and output exists
                        if in_uac in keys(bus.balance_table_inputs) && out_uac in keys(bus.balance_table_outputs)
                            output_key = "EnergyFlow " * value_key
                            push!(outputs, EnergySystems.OutputKey(; unit=bus,
                                                                   medium=medium,
                                                                   value_key=output_key))
                            output_key = "TemperatureFlow " * value_key
                            push!(outputs, EnergySystems.OutputKey(; unit=bus,
                                                                   medium=medium,
                                                                   value_key=output_key))
                            success = true
                            break
                        end
                    end
                end
                if !success
                    @error "The requested energy flow between components \"$(value_key)\" for medium \"$(media_key)\" could not " *
                           "be found for the CSV output or the plot output. Note that only connections between " *
                           "components with one or more busses but without any other component in between can be exported!"
                    throw(InputError)
                end
            end
        else
            @error "The key \"$(key)\" in the provided output keys for CSV output or plot output could not be found. " *
                   "It either has to be a medium or a component used in the current energy system."
            throw(InputError)
        end
    end

    return outputs
end

"""
reset_file(filepath, output_keys)

Reset the output file and add headers for the given outputs.
"""
function reset_file(filepath::String,
                    output_keys::Union{Nothing,Vector{EnergySystems.OutputKey}},
                    weather_data_keys::Union{Nothing,Vector{String}},
                    csv_time_unit::String)
    open(abspath(filepath), "w") do file_handle
        if csv_time_unit == "seconds"
            time_unit = "[s]"
        elseif csv_time_unit == "minutes"
            time_unit = "[min]"
        elseif csv_time_unit == "hours"
            time_unit = "[h]"
        elseif csv_time_unit == "date"
            time_unit = "[dd.mm.yyyy HH:MM:SS]"
        end

        write(file_handle, "Time $time_unit")

        if output_keys !== nothing
            for outkey in output_keys
                if outkey.medium === nothing
                    header = "$(outkey.unit.uac) $(outkey.value_key)"
                else
                    if startswith(outkey.value_key, "EnergyFlow") || startswith(outkey.value_key, "TemperatureFlow")
                        header = "$(outkey.medium) $(outkey.value_key)"
                    else
                        header = "$(outkey.unit.uac) $(outkey.medium) $(outkey.value_key)"
                    end
                end
                write(file_handle, ";$header")
            end
        end
        if weather_data_keys !== nothing
            for key in weather_data_keys
                write(file_handle, ";Weather $key")
            end
        end

        write(file_handle, "\n")
    end
end

"""
write_to_file(filepath, output_keys, time)

Write the given outputs for the given time to file.
"""
function write_to_file(filepath::String,
                       output_keys::Union{Nothing,Vector{EnergySystems.OutputKey}},
                       weather_data_keys::Union{Nothing,Vector{String}},
                       sim_params::Dict{String,Any},
                       csv_time_unit::String)
    open(abspath(filepath), "a") do file_handle
        if csv_time_unit == "seconds"
            time = sim_params["time_since_output"]
        elseif csv_time_unit == "minutes"
            time = sim_params["time_since_output"] / 60
        elseif csv_time_unit == "hours"
            time = sim_params["time_since_output"] / 60 / 60
        elseif csv_time_unit == "date"
            time = Dates.format(sim_params["current_date"], "dd.mm.yyyy HH:MM:SS")
        end

        write(file_handle, "$time")

        if output_keys !== nothing
            for outkey in output_keys
                value = output_value(outkey.unit, outkey)
                value = replace("$value", "." => ",")
                write(file_handle, ";$value")
            end
        end
        if weather_data_keys !== nothing
            for key in weather_data_keys
                value = Profiles.value_at_time(getfield(sim_params["weather_data"], Symbol(key)), sim_params)
                value = replace("$value", "." => ",")
                write(file_handle, ";$value")
            end
        end

        write(file_handle, "\n")
    end
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
    dump_auxiliary_outputs(file_path, components, order_of_operations, sim_params)

Dump a bunch of information to file that might be useful to explain the result of a run.

This is mostly used for debugging and development purposes, but might prove useful in
general to find out why the energy system behaves in the simulation as it does.
"""
function dump_auxiliary_outputs(project_config::AbstractDict{AbstractString,Any},
                                components::Grouping,
                                order_of_operations::OrderOfOperations,
                                sim_params::Dict{String,Any})
    # export order of operations
    if default(project_config["io_settings"], "auxiliary_info", false)
        aux_info_file_path = default(project_config["io_settings"], "auxiliary_info_file", "./output/auxiliary_info.md")
        open(abspath(aux_info_file_path), "w") do file_handle
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

        @info "Auxiliary info dumped to file $(aux_info_file_path)"
    end

    # plot additional figures potentially available from components after initialisation
    if default(project_config["io_settings"], "auxiliary_plots", false)
        aux_plots_output_path = default(project_config["io_settings"], "auxiliary_plots_path", "./output/")
        aux_plots_formats = default(project_config["io_settings"], "auxiliary_plots_formats", ["png"])
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
create_profile_line_plots(data, keys, user_input)

create a line plot with data and label. user_input is dict from input file
"""
function create_profile_line_plots(outputs_plot_data::Union{Nothing,Matrix{Float64}},
                                   outputs_plot_keys::Union{Nothing,Vector{EnergySystems.OutputKey}},
                                   outputs_plot_weather::Union{Nothing,Matrix{Float64}},
                                   outputs_plot_weather_keys::Union{Nothing,Vector{String}},
                                   project_config::AbstractDict{AbstractString,Any},
                                   sim_params::Dict{String,Any})
    plot_all = isa(project_config["io_settings"]["output_plot"], String) &&
               project_config["io_settings"]["output_plot"][1:3] == "all"
    plot_data = outputs_plot_data !== nothing
    plot_weather = outputs_plot_weather !== nothing

    # set Axis, unit and scale factor if given
    if plot_all  # plot all outputs. Here no units or scaling factors are available.
        labels = String[]
        n = 1
        for outkey in outputs_plot_keys
            if outkey.medium === nothing
                push!(labels, string("$(outkey.unit.uac) $(outkey.value_key)"))
            elseif startswith(outkey.value_key, "EnergyFlow") || startswith(outkey.value_key, "TemperatureFlow")
                push!(labels, string("$(outkey.medium) $(outkey.value_key)"))
            else
                push!(labels, string("$(outkey.unit.uac) $(outkey.medium) $(outkey.value_key)"))
            end
            n += 1
        end
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
            for plot in project_config["io_settings"]["output_plot"]
                if occursin("->", first(plot[2]["key"])[2][1])
                    # Here we are dealing with EnergyFlow and TemperatureFlow --> Two meta information required if 
                    # temperature should be plotted
                    if !isa(plot[2]["axis"], AbstractVector) || length(plot[2]["axis"]) !== 2 ||
                       !isa(plot[2]["unit"], AbstractVector) || length(plot[2]["unit"]) !== 2 ||
                       !isa(plot[2]["scale_factor"], AbstractVector) || length(plot[2]["scale_factor"]) !== 2
                        # only one set of meta information given. Use the given one both for energy and temperature flow
                        push!(axis, isa(plot[2]["axis"], AbstractVector) ? plot[2]["axis"][1] : plot[2]["axis"])
                        push!(unit, isa(plot[2]["unit"], AbstractVector) ? plot[2]["unit"][1] : plot[2]["unit"])
                        push!(scale_fact,
                              isa(plot[2]["scale_factor"], AbstractVector) ? plot[2]["scale_factor"][1] :
                              plot[2]["scale_factor"])

                        push!(axis, "nothing")
                        push!(unit, "nothing")
                        push!(scale_fact, NaN)
                        @info "For the generation of the output_plot, the meta information for entry $(plot[1]) " *
                              "do not contain two values. Therefore, only energy values will be output. If you want " *
                              "to output also a corresponding temperature, provide two meta information: " *
                              "[EnergyFlow, TemperatureFlow] for axis, unit and scale_factor."
                    else
                        push!(axis, string(plot[2]["axis"][1]))
                        push!(unit, string(plot[2]["unit"][1]))
                        push!(scale_fact, plot[2]["scale_factor"][1])

                        push!(axis, string(plot[2]["axis"][2]))
                        push!(unit, string(plot[2]["unit"][2]))
                        push!(scale_fact, plot[2]["scale_factor"][2])
                    end
                else
                    push!(axis, string(plot[2]["axis"]))
                    push!(unit, string(plot[2]["unit"]))
                    push!(scale_fact, plot[2]["scale_factor"])
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
        n = 1
        if plot_data
            for outkey in outputs_plot_keys
                if outkey.medium === nothing
                    push!(labels, string("$(outkey.unit.uac) $(outkey.value_key) [$(unit[n])] ($(axis[n]))"))
                elseif startswith(outkey.value_key, "EnergyFlow") || startswith(outkey.value_key, "TemperatureFlow")
                    push!(labels, string("$(outkey.medium) $(outkey.value_key) [$(unit[n])] ($(axis[n]))"))
                else
                    push!(labels,
                          string("$(outkey.unit.uac) $(outkey.medium) $(outkey.value_key) [$(unit[n])] ($(axis[n]))"))
                end
                n += 1
            end
        end
        if plot_weather
            for outkey in outputs_plot_weather_keys
                push!(labels, string("Weather $(outkey) [($(axis[n]))"))
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
    output_plot_time_unit = haskey(project_config["io_settings"], "output_plot_time_unit") ?
                            project_config["io_settings"]["output_plot_time_unit"] : "date"
    if !(output_plot_time_unit in ["seconds", "minutes", "hours", "date"])
        @info "The `output_plot_time_unit` has to be one of `seconds`, `minutes`, `hours` or `date. " *
              "It will be set to `date` as default."
        output_plot_time_unit = "date"
    end
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

    file_path = default(project_config["io_settings"],
                        "output_plot_file",
                        "./output/output_plot.html")
    savefig(p, file_path)
end

"""
create_sankey(output_all_sourcenames, output_all_targetnames, output_all_values, medium_of_interfaces, nr_of_interfaces)

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
                       io_settings::AbstractDict{String,Any})

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
             output_all_value_sum[interface_new] == output_all_value_sum[interface_new - 1]))
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

    if io_settings["sankey_plot"] == "default" || !isa(io_settings["sankey_plot"], AbstractDict{})
        if length(unique_medium_labels) > 1
            colors = get(ColorSchemes.roma,
                         (0:(length(unique_medium_labels) - 1)) ./ (length(unique_medium_labels) - 1))
            color_map = Dict(zip(unique_medium_labels, colors))
        else
            color_map = Dict(unique_medium_labels[1] => get(ColorSchemes.roma, 0.7))
        end
    else
        color_map = Dict{String,Any}(io_settings["sankey_plot"])
        for (medium, color) in color_map
            try
                color_entry = Colors.color_names[color]
                color_map[medium] = RGB(color_entry[1] / 255, color_entry[2] / 255, color_entry[3] / 255)
            catch e
                @error "The color for the sankey '$color' of medium '$medium' could not be detected. " *
                       "The following error occurred: $e\n" *
                       "Color has to be one of: $(collect(keys(Colors.color_names)))"
                throw(InputError)
            end
        end
        for medium in unique_medium_labels
            if medium in keys(color_map)
                continue
            elseif medium == "hide_medium"
                color_map[medium] = parse(RGBA, "rgba(0,0,0,0)")
            else
                @error "The color for the medium '$medium' for the sankey could not be found in the input file. Please add the medium and its color in 'sankey_plot'."
                throw(InputError)
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
    file_path = default(io_settings, "sankey_plot_file", "./output/output_sankey.html")
    savefig(p, file_path)
end
