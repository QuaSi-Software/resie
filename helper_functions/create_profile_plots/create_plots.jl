# This is a standalone helper function to create plots from a text/csv file containing profile data.
# The textfile should have the following structure (or adjust the separator below...):
# Date             ; profile1 ; profile2 ...
# 01.01.2014 00:00 ; 5.5      ; 7.456    
# ....             ; ....     ; ....
#
# Several adjustments can be made like aggregation, selecting a time span to plot, calculate comparison
# metrics and others. More function will may be added later if needed. Currently only one line plot with
# several lines is created. Multi-plots may will follow...

using CSV
using DataFrames
using Plots
using Dates
using Statistics

function plot_data(x, profiles::Array, legends::Array{String,1}, colors::Array, title::String, output_name::String,
                   file_formats::Array, x_label::String, y_label::String, showplots::Bool)
    x_min, x_max = minimum(x), maximum(x)
    custom_y_min = 0 / 1.1  # set to zero if no custom ymin should be set as lower limit
    custom_y_max = 0 / 1.1  # set to zero if no custom ymax should be set as upper limit
    y_min, y_max = custom_y_min, custom_y_max
    for i in 1:size(profiles, 2)
        y_min, y_max = minimum([custom_y_min; y_min; profiles[:, i]]), maximum([custom_y_max; y_max; profiles[:, i]])
    end
    y_min, y_max = 1.1 * y_min, 1.1 * y_max

    plot() # Initialize an empty plot
    for i in 1:size(profiles, 2)
        plot!(x,
              profiles[:, i];
              color=colors[i],
              label=legends[i],
              title=title,
              xlabel=x_label,
              ylabel=y_label,
              xlims=(x_min, x_max),
              ylims=(y_min, y_max),
              legend=true,
              linewidth=8,
              gridlinewidth=1,
              size=(1800, 1200),
              titlefontsize=30,
              guidefontsize=24,
              tickfontsize=24,
              legendfontsize=24,
              grid=true,
              minorgrid=true,
              margin=15Plots.mm)
    end
    for file_format in file_formats
        savefig(output_name * "." * file_format)
    end
    if showplots
        gui()
        println("Press Enter to close program and all figures...")
        readline()
    end
end

"""
Aggregate profiles in a matrix over specified timesteps, including datetime aggregation.

# Arguments
- `datetime_vector`:        Vector of DateTime values corresponding to each timestep.
- `data_matrix`:            Matrix of profile data (rows: timesteps, columns: profiles).
- `timestep_aggregation`:   Number of timesteps to aggregate.
- `aggregation_methods`:    Array specifying aggregation method (`:mean` or `:sum`) for each profile.
- `time_aggregation`:       value specifying aggregation method (`:begin`, `:end`, `:mid`) for time.

# Returns
- Tuple containing the aggregated datetime vector and the aggregated data matrix .
"""
function aggregate_profiles_with_datetime(datetime_vector, data_matrix, timestep_aggregation, aggregation_methods,
                                          time_aggregation)
    matrix_size = size(data_matrix)
    if length(matrix_size) > 1
        nrows, ncols = matrix_size
    else
        nrows = matrix_size
        ncols = 1
    end

    aggregated_nrows = ceil(Int, nrows / timestep_aggregation)

    aggregated_matrix = Matrix{Float64}(undef, aggregated_nrows, ncols)
    aggregated_datetime = Vector{DateTime}(undef, aggregated_nrows)

    for col in 1:ncols
        method = aggregation_methods[col]
        for row in 1:aggregated_nrows
            start_index = (row - 1) * timestep_aggregation + 1
            end_index = min(row * timestep_aggregation, nrows)

            if method == :mean
                aggregated_matrix[row, col] = mean(data_matrix[start_index:end_index, col])
            elseif method == :sum
                aggregated_matrix[row, col] = sum(data_matrix[start_index:end_index, col])
            else
                error("Unsupported aggregation method for data: $method")
            end
        end
    end

    # Aggregate datetime
    for row in 1:aggregated_nrows
        start_index = (row - 1) * timestep_aggregation + 1
        end_index = min(row * timestep_aggregation, nrows)
        if time_aggregation == :begin
            aggregated_datetime[row] = datetime_vector[start_index]
        elseif time_aggregation == :end
            aggregated_datetime[row] = datetime_vector[end_index]
        elseif time_aggregation == :mid
            datetime_range = datetime_vector[start_index:end_index]
            aggregated_datetime[row] = datetime_range[1] + (datetime_range[end] - datetime_range[1]) ÷ 2
        else
            error("Unsupported aggregation method for time step: $time_aggregation")
        end
    end

    return aggregated_datetime, aggregated_matrix
end

"""
Calculates the mean and max difference of two profiles.

# Arguments
- `data`:        Matrix with profile data
- `compares`:    Array of Tuples containing the row number in data to be compared
- `labels`:      Array of Strings with the name of the columns in data

# Returns
- mean and maximum difference as array for each tuple in compares
"""
function mean_differences(data::Matrix, compares, labels)
    mean_diffs = []
    max_diffs = []
    for compar in compares
        col1 = compar[1]
        col2 = compar[2]

        mean_diff = mean(abs.(data[:, col1] .- data[:, col2]))
        max_diff = maximum(abs.(data[:, col1] .- data[:, col2]))
        push!(mean_diffs, mean_diff)
        push!(max_diffs, max_diff)
        print("Comparing $(labels[col1]) and $(labels[col2]) results in an absolute mean different of $mean_diff " *
              "and an abolute max difference of $max_diff.\n")
    end
    return mean_diffs, max_diffs
end

# enter link to file, first column is y axis in datetime format
path = "helper_create_plots/example_input_plot.txt"
data_matrix = CSV.Tables.matrix(CSV.File(path; delim=';'))

# strip data to specific time frame if needed
# data_matrix = data_matrix[14*7*24:15*7*24,:]

# may change to DateTime format of first column
x_vals_date = DateTime.(String.(strip.(data_matrix[:, 1])), "dd.mm.yyyy HH:MM")

# select columns containing the profiles to plot and set corresponding labels and colors
data = data_matrix[:, 2:4]
labels = ["Measurement", "EED", "ReSiE"]
colors = ["#004488", "#DDAA33", "#BB5566"]
#          blue        yellow     red          from Paul Tol's Colour Schemes: high contrast

# calculate profile metrix. The second argument is an array of tuples reffering to the column numbers 
# of data that should be compared.
mean_diffs, max_diffs = mean_differences(data, [(1, 3), (2, 3)], labels)

# set aggregation method for each column. Can be :mean or :sum.
aggregation_methods = [:mean, :mean, :mean]
# set number of datapoints to aggregate into one
number_of_datapoints_to_aggregate = 1
# set time stamp that should be used, can be :begin, :mid or :end
time_stamp_aggregation = :begin

# set figure title and axis labels
title = "Example\nTitle"
x_label = "date"
y_label = "temperature [°C]"
# set name and file formats to save figure
file_formats = ["png", "svg"]
filename = "helper_functions/create_profile_plots/comparison_plot"   # without file format ending!
# set if plot should be shown. Attention: Task has to be closed by pressing Enter if set to true.
showplots = false

x_vals_date, data = aggregate_profiles_with_datetime(x_vals_date, data, number_of_datapoints_to_aggregate,
                                                     aggregation_methods, time_stamp_aggregation)
plot_data(x_vals_date, data, labels, colors, title, filename, file_formats, x_label, y_label, showplots)
