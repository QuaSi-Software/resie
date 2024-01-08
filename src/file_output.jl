# this file contains functionality for writing output of the simulation to files.


"""
get_output_keys(config[io_settings], components)

This function determines both for the lineplot and the csv output if
they should be created:
    - if not, "nothing" will be returned as key list
    - if yes, the key lists of the outputs will be returned, either containing
        - all possible keys if this is requested in the input file or
        - only the requested keys as requested in the input file
"""
function get_output_keys(io_settings::Dict{String, Any}, 
                         components::Grouping
                         )::Tuple{Union{Nothing, Vector{EnergySystems.OutputKey}}, Union{Nothing, Vector{EnergySystems.OutputKey}}}
    # determine if lineplot and csv should be created
    do_plot_all_outputs = false
    do_create_plot = false
    if haskey(io_settings, "output_plot")
        do_create_plot = true  # if the key exists, set do_create_plot to true by default
        if io_settings["output_plot"] == "all"
            do_plot_all_outputs = true
        elseif io_settings["output_plot"] == "nothing"
            do_create_plot = false
        end
    end

    do_write_all_CSV_outputs = false
    do_write_CSV = false
    if haskey(io_settings, "output_keys")
        do_write_CSV = true  # if the key exists, set do_create_plot to true by default
        if io_settings["output_keys"] == "all"
            do_write_all_CSV_outputs = true
        elseif io_settings["output_keys"] == "nothing"
            do_write_CSV = false
        end
    end

     # collect all possible outputs of all units if needed
    if do_plot_all_outputs || do_write_all_CSV_outputs
        all_output_keys = Vector{EnergySystems.OutputKey}()
        for unit in components
            temp_dict = Dict{String, Any}(unit[2].uac => output_values(unit[2]))
            append!(all_output_keys, output_keys(components, temp_dict))
        end
    end

    # collect output keys for lineplot and csv output
    if do_create_plot  # line plot
        output_keys_lineplot = do_plot_all_outputs ? all_output_keys : Vector{EnergySystems.OutputKey}()
        # Collect output keys if not plotting all outputs
        if !do_plot_all_outputs
            for plot in io_settings["output_plot"]
                key = plot[2]["key"]
                append!(output_keys_lineplot, output_keys(components, key))
            end
        end
    else
        output_keys_lineplot = nothing
    end

    if do_write_CSV  # csv export
        if do_write_all_CSV_outputs # gather all outputs
            output_keys_to_csv = all_output_keys
        else  # get only requested output keys from input file for CSV-export
            output_keys_to_csv = output_keys(components, io_settings["output_keys"])
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

The information is returned as single vectors with the indicees matching together.
"""
function get_interface_information(components::Grouping)::Tuple{Int64,Vector{Any},Vector{Any},Vector{Any}} 
    nr_of_interfaces = 0
    medium_of_interfaces = []
    output_sourcenames_sankey = []
    output_targetnames_sankey = []
    for each_component in components
        for each_outputinterface in each_component[2].output_interfaces
            medium = nothing  # reset medium
            if isa(each_outputinterface, Pair) # some output_interfaces are wrapped in a Touple
                medium = each_outputinterface[1] # then, the medium is stored separately
                each_outputinterface = each_outputinterface[2]
            end

            if isdefined(each_outputinterface, :target)
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
                    push!(output_sourcenames_sankey, string(each_outputinterface.source.uac,"_total_supply"))
                    push!(output_targetnames_sankey, each_outputinterface.source.uac)       
                    push!(medium_of_interfaces, "hide_medium")
                    nr_of_interfaces += 1
                end
                if each_outputinterface.target.sys_function == EnergySystems.sf_fixed_sink
                    push!(output_sourcenames_sankey, each_outputinterface.target.uac)
                    push!(output_targetnames_sankey, string(each_outputinterface.target.uac,"_total_demand"))       
                    push!(medium_of_interfaces, "hide_medium")
                    nr_of_interfaces += 1
                end 
            end
        end

        # add losses
        if hasfield(typeof(each_component[2]), Symbol("losses"))
            push!(output_sourcenames_sankey, each_component[2].uac)
            push!(output_targetnames_sankey, "Losses")       
            push!(medium_of_interfaces, "Losses")
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

Collects and returns the energy that was transportet through every interface.
If the balance of an interface was not zero, the actual energy that was flowing
is written to the outputs.
Attention: This can lead to overfilling of demands which is currenlty not visible
in the sankey diagram!
"""
function collect_interface_energies(components::Grouping, nr_of_interfaces::Int)
    n = 1
    energies = zeros(Float64, nr_of_interfaces)
    for each_component in components
        for each_outputinterface in each_component[2].output_interfaces
            if isa(each_outputinterface, Pair) # some output_interfaces are wrapped in a Touple
                each_outputinterface = each_outputinterface[2]
            end
            if isdefined(each_outputinterface, :target)
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
        if hasfield(typeof(each_component[2]), Symbol("losses"))
            energies[n] = each_component[2].losses
            n += 1
        end
    end
    return energies
end    

"""
output_keys(from_config)

Transform the output keys definition in the project config file into a list of OutputKey
items. This is done to speed up selection of values for the output in each time step,
as this transformation has to be done only once at the beginning.
"""
function output_keys(
    components::Grouping,
    from_config::Dict{String,Any}
)::Vector{EnergySystems.OutputKey}
    outputs = Vector{EnergySystems.OutputKey}()

    for unit_key in keys(from_config)
        unit = components[unit_key]

        for entry in from_config[unit_key]
            splitted = split(String(entry))
            if length(splitted) > 1
                medium_key = splitted[1]
                medium = Symbol(String(medium_key))
                value_key = splitted[2]
            else
                medium = nothing
                value_key = splitted[1]
            end

            push!(outputs, EnergySystems.OutputKey(
                unit=unit,
                medium=medium,
                value_key=value_key
            ))
        end
    end

    return outputs
end

"""
reset_file(filepath, output_keys)

Reset the output file and add headers for the given outputs.
"""
function reset_file(
    filepath::String,
    output_keys::Vector{EnergySystems.OutputKey}
)
    open(abspath(filepath), "w") do file_handle
        write(file_handle, "Time [s]")

        for outkey in output_keys
            if outkey.medium === nothing
                header = "$(outkey.unit.uac) $(outkey.value_key)"
            else
                header = "$(outkey.unit.uac) $(outkey.medium) $(outkey.value_key)"
            end
            write(file_handle, ";$header")
        end

        write(file_handle, "\n")
    end
end

"""
write_to_file(filepath, output_keys, time)

Write the given outputs for the given time to file.
"""
function write_to_file(
    filepath::String,
    output_keys::Vector{EnergySystems.OutputKey},
    time::Int
)
    open(abspath(filepath), "a") do file_handle
        write(file_handle, "$time")

        for outkey in output_keys
            value = output_value(outkey.unit, outkey)
            value = replace("$value", "." => ",")
            write(file_handle, ";$value")
        end

        write(file_handle, "\n")
    end
end

"""
dump_info(file_path, components, order_of_operations, sim_params)

Dump a bunch of information to file that might be useful to explain the result of a run.

This is mostly used for debugging and development purposes, but might prove useful in
general to find out why the energy system behaves in the simulation as it does.
"""
function dump_info(
    file_path::String,
    components::Grouping,
    order_of_operations::StepInstructions,
    sim_params::Dict{String,Any}
)
    open(abspath(file_path), "w") do file_handle
        write(file_handle, "# Simulation step order\n")

        for entry in order_of_operations
            for step in entry[2:lastindex(entry)]
                if entry == last(order_of_operations)
                    write(file_handle, "\"$(entry[1]) $(entry[2])\"\n")
                else
                    write(file_handle, "\"$(entry[1]) $(entry[2])\",\n")
                end
            end
        end
    end
end


"""
geather_output_data(output_keys, time)

returns a vektor with the requested data in output_keys
"""
function geather_output_data(
    output_keys::Vector{EnergySystems.OutputKey},
    time::Int
)
    return_values = Vector{Any}()
    append!(return_values, time)

    for outkey in output_keys
        append!(return_values, output_value(outkey.unit, outkey))
    end

    return return_values
end


"""
create_profile_line_plots(data, keys, user_input)

create a line plot with data and label. user_input is dict from input file
"""
function create_profile_line_plots(
    outputs_plot_data::Matrix{Float64},
    outputs_plot_keys::Vector{EnergySystems.OutputKey},
    project_config::Dict{AbstractString, Any},
    sim_params::Dict{String, Any}
)
    plot_all = project_config["io_settings"]["output_plot"] == "all"
    
    # set Axis, unit and scale factor if given
    if plot_all  # plot all outputs. Here no units or scaling factors are available.
        labels = String[]
        n = 1
        for outkey in outputs_plot_keys
            if outkey.medium === nothing
                push!(labels, string("$(outkey.unit.uac) $(outkey.value_key)"))
            else
                push!(labels, string("$(outkey.unit.uac) $(outkey.medium) $(outkey.value_key)"))
            end
            n += 1
        end
    else # plot only defined outputs. Here units and scaling factors are available.
        axis = String[]
        unit = String[]
        scale_fact = Float64[]
        for plot in project_config["io_settings"]["output_plot"]
            push!(axis, string(plot[2]["axis"]))
            push!(unit, string(plot[2]["unit"]))
            push!(scale_fact, plot[2]["scale_factor"])
        end
   
        # create legend entries
        labels = String[]
        n = 1
        for outkey in outputs_plot_keys
            if outkey.medium === nothing
                push!(labels, string("$(outkey.unit.uac) $(outkey.value_key) [$(unit[n])] ($(axis[n]))"))
            else
                push!(labels, string("$(outkey.unit.uac) $(outkey.medium) $(outkey.value_key) [$(unit[n])] ($(axis[n]))"))
            end
            n += 1
        end
    end

    # create plot
    x = outputs_plot_data[:, 1]
    y = outputs_plot_data[:, 2:end]
    traces = GenericTrace[]
    for i in axes(y, 2)
        if plot_all
            trace = scatter(x=x / 60 / 60, y=y[:, i], mode="lines", name=labels[i])            
        else
            trace = scatter(x=x / 60 / 60, y=scale_fact[i] * y[:, i], mode="lines", name=labels[i])
            if axis[i] == "right"
                trace.yaxis = "y2"
            else  # default is left axis
                trace.yaxis = "y1"
            end
        end
        push!(traces, trace)
    end

    layout = Layout(title_text="Plot of outputs as defined in the input-file. Attention: Energies are given within the current time step of $(Int(sim_params["time_step_seconds"]/60)) min!",
        xaxis_title_text="Time [hour]",
        yaxis_title_text="",
        yaxis2=attr(title="", overlaying="y", side="right"))

    p = plot(traces, layout)
    savefig(p, "output/output_plot.html")

end


"""
create_sankey(output_all_sourcenames, output_all_targetnames, output_all_values, medium_of_interfaces, nr_of_interfaces)

create a sankey plot. 
Inputs:
output_all_sourcenames and *sinknames are vectors with names of the source and sink of each interface
output_all_values are logs with data from each timestep in the shape [timestep,interface].
medium_of_interface is a vector of the medium corresponding to each interface
nr_of_interfaces is the total number of iterfaces in the current energy system
"""
function create_sankey(
    output_all_sourcenames::Vector{Any},
    output_all_targetnames::Vector{Any},
    output_all_values::Matrix{Float64},
    medium_of_interfaces::Vector{Any},
    nr_of_interfaces::Int64,
)

    # sum up data of each interface
    output_all_value_sum = zeros(Float64, nr_of_interfaces)
    for interface in 1:nr_of_interfaces
        output_all_value_sum[interface] = sum(output_all_values[:, interface])
    end

    # remove data that should not be plotted in Sankey
    interface_new = 1
    for _ in 1:nr_of_interfaces
        if ( 
            # remove oxygen from data as the energy of oxygen is considered to be zero
            medium_of_interfaces[interface_new] == :m_c_g_o2 ||  
            # remove real sinks and sources if they match the delivered energy
            # to enable this to work, the "real" demand/supply needs to be always one entry below the delivered/requested one in the array!
             (medium_of_interfaces[interface_new] == "hide_medium" && output_all_value_sum[interface_new] == output_all_value_sum[interface_new-1])
        )
            deleteat!(output_all_sourcenames, interface_new)
            deleteat!(output_all_targetnames, interface_new)
            deleteat!(output_all_value_sum, interface_new)
            deleteat!(medium_of_interfaces, interface_new)
            interface_new -= 1
            nr_of_interfaces -= 1
        end
        interface_new += 1
    end

    # add 0.000001 to all interfaces (except of losses) to display interfaces that are zero
    output_all_value_sum += (medium_of_interfaces .!= "Losses") .* 0.000001

    # prepare data for sankey diagram and create sankey
    # set label of blocks
    block_labels = union(output_all_sourcenames, output_all_targetnames)
    block_labels_unique = Dict(blockname => i for (i, blockname) in enumerate(block_labels))
    output_all_source_num = [block_labels_unique[blockname] for blockname in output_all_sourcenames]
    output_all_target_num = [block_labels_unique[blockname] for blockname in output_all_targetnames]

    # set label and colour of interfaces
    medium_labels = [split(string(s), '.')[end] for s in medium_of_interfaces]
    unique_medium_labels = unique(medium_labels)
    if length(unique_medium_labels) > 1
        colors = get(ColorSchemes.roma, (0:length(unique_medium_labels)-1) ./ (length(unique_medium_labels) - 1))
        color_map = Dict(zip(unique_medium_labels, colors))
        colors_for_medium = map(x -> color_map[x], medium_labels)
    else # account for cases with only one medium in the system topology
        colors_for_medium = get(ColorSchemes.roma, 0.5)
    end

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
    p = plot(sankey(
            node=attr(
                pad=25,
                thickness=20,
                line=do_hide_real_demands ? attr(color="white", width=0.0) : nothing,
                label=block_labels,
                color=do_hide_real_demands ? block_colors : "blue" ,
            ),
            link=attr(
                source=output_all_source_num .- 1, # indices correspond to block_labels starting from index 0
                target=output_all_target_num .- 1, # indices correspond to block_labels starting from index 0
                value=output_all_value_sum,
                label=medium_labels,
                color=colors_for_medium_RGBA
            )),
        Layout( title_text="Sankey diagram of system topology and energy flows", 
                font_size=14,
                )
    )

    # save plot
    savefig(p, "output/output_sankey.html")

end