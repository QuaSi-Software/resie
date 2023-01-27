# this file contains functionality for writing output of the simulation to files.

"""
output_keys(from_config)

Transform the output keys definition in the project config file into a list of OutputKey
items. This is done to speed up selection of values for the output in each time step,
as this transformation has to be done only once at the beginning.
"""
function output_keys(
    systems :: Grouping,
    from_config :: Dict{String, Any}
) :: Vector{EnergySystems.OutputKey}
    outputs = Vector{EnergySystems.OutputKey}()

    for unit_key in keys(from_config)
        unit = systems[unit_key]

        for entry in from_config[unit_key]
            splitted = split(String(entry))
            if length(splitted) > 1
                medium_key = splitted[1]
                medium = getproperty(EnergySystems, Symbol(String(medium_key)))
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
    filepath :: String,
    output_keys :: Vector{EnergySystems.OutputKey}
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
    filepath :: String,
    output_keys :: Vector{EnergySystems.OutputKey},
    time :: Int
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
dump_info(file_path, systems, order_of_operations, parameters)

Dump a bunch of information to file that might be useful to explain the result of a run.

This is mostly used for debugging and development purposes, but might prove useful in
general to find out why the systems behave in the simulation as they do.
"""
function dump_info(
    file_path :: String,
    systems :: Grouping,
    order_of_operations :: StepInstructions,
    parameters :: Dict{String, Any}
)
    open(abspath(file_path), "w") do file_handle
        write(file_handle, "# Simulation step order\n")

        for entry in order_of_operations
            for step in entry[2:lastindex(entry)]
                write(file_handle, "1. `$(entry[1]) $(entry[2])`\n")
            end
        end
    end
end


"""
geather_output_data(output_keys, time)

returns a vektor with the requested data in output_keys
"""
function geather_output_data(
    output_keys :: Vector{EnergySystems.OutputKey},
    time :: Int
)
    return_values = Vector{Any}()
    append!( return_values, time )

    for outkey in output_keys
        append!( return_values, output_value(outkey.unit, outkey) )
    end

    return return_values
end


"""
create_plots(data, keys, user_input)

create a line plot with data and label. user_input is dict from input file
"""
function create_plots(
    outputs_plot_data :: Matrix{Float64}, 
    outputs_plot_keys :: Vector{EnergySystems.OutputKey},
    user_input :: Dict{String, Any}
    )

    # set Axis, unit and scale factor
    axis = String[]
    unit = String[]
    scale_fact = Float64[]
    for plot in user_input
        push!(axis, string(plot[2]["axis"]))
        push!(unit, string(plot[2]["unit"]))
        push!(scale_fact, plot[2]["scale_factor"])
    end

    # create legend entries
    labels = String[] 
    n = 1
    for outkey in outputs_plot_keys
        if outkey.medium === nothing
            push!( labels, string("$(outkey.unit.uac) $(outkey.value_key) [$(unit[n])] ($(axis[n]))"))
        else
            push!( labels, string("$(outkey.unit.uac) $(outkey.medium) $(outkey.value_key) [$(unit[n])] ($(axis[n]))"))
        end
        n += 1
    end

    # create plot
    x = outputs_plot_data[:,1]
    y = outputs_plot_data[:,2:end]
    traces = GenericTrace[]
    for i in axes(y,2)
        trace = scatter(x = x/60, y = scale_fact[i] * y[:,i], mode = "lines", name = labels[i])
        if axis[i] == "right"
            trace.yaxis = "y2"
        else  # default is left axis
            trace.yaxis = "y1"
        end
        push!(traces, trace)
    end

    layout = Layout(title_text = "Plot of outputs as defined in the input-file", 
                    xaxis_title_text = "Time [min]", 
                    yaxis_title_text = "", 
                    yaxis2 = attr(title = "", overlaying = "y", side = "right"))

    p = plot(traces, layout)
    savefig(p, "output/output_plot.html") 

end    
