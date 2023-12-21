# Pre-Processing for the geothermal probes model
# Get pre-calculated g-function values for big time steps out of json-file
# calculate g-function values for short timesteps
# overlap calculated g-function values to composite g-function array

# Packages
using SpecialFunctions  # To use "expint"-function; Pkg.add("SpecialFunctions").
using JSON              # To read g-function values out of open source library in JSON-format. 
using JSON2              # To read g-function values out of open source library in JSON-format. 
using Interpolations    # Pkg.add("SpecialFunctions")
using CSV               # For validation purposes of this PreProcessing file.             
using Plots             # For validation purposes of this PreProcessing file. 
using DataFrames        # For writing CSV

# Funktionen:

function infinite_line_source(borehole_radius, soil_diffusivity, g_time_vector)
    # function calculates g-function solution of an infinite line source/sink by Kelvin. 
    g_ILS = zeros(length(g_time_vector))
    
    for i in 1:length(g_time_vector)
        if g_time_vector[i] == 0    # if time vector starts at 0, g-function will be set to 0.
            g_ILS[i] = 0
        else
            g_ILS[i] = round.(0.5 * expint(borehole_radius^2 / (4 * soil_diffusivity * g_time_vector[i])), digits =4)
        end
    end

    return g_ILS
end

function ln_to_normal(library_grid_points_normalized_time, steady_state_time)
    # calculates absolut time values out of normalized time values from spitler/cook library.
    # time format of grid points: library_grid_points_normalized_time = ln(time/steady_state_time)
    # --> time = exp(library_grid_points_normalized_time)*steady_state_time

    n = length(library_grid_points_normalized_time)     # number of grid points of spitler/cook data set.
    time = zeros(n)

    for i in 1:n
        time[i] = round(exp(library_grid_points_normalized_time[i]) * steady_state_time)
    end

    return time
end

function round_to_simulation_step_width(grid_point,simulation_step_width)
# making sure, that time stamps from library fit into simulation step width.
grid_point_rounded = zeros(length(grid_point))
for i=1:length(grid_point)
    grid_point_rounded[i] = floor(Int(round(grid_point[i]/simulation_step_width))) * simulation_step_width
end 
    return grid_point_rounded
end

function get_library_g_values(borehole_depth, borehole_radius,borehole_spacing, m, n, t)
# Load library file using the JSON package - change path to library-file-path
# There are various sets of grid points refering to different borehole depths. Two sets of grid points have to be read out of .json library-file to do interpolation, as they are refered to specific borehole depth.
# check, which two sets of grid points (refered to borehole depth borehole_depth_library_lower and borehole depth borehole_depth_library_upper) will be read out of json file later.
# first, each set of grid points is refered to a default borehole radius (borehole_radius_library_lower, borehole_radius_library_upper). Later, the grid points will be corrected to the real borehole radius.

# borehole_depth will be set by user. Check, which two sets of grid-points will be read out of the library-file to be interpolated later.
if borehole_depth >= 24 && borehole_depth <= 48
    borehole_depth_library_lower = 24
    borehole_radius_library_lower = 0.075

    borehole_depth_library_upper = 48
    borehole_radius_library_upper = 0.075

elseif borehole_depth >= 48 && borehole_depth <= 96
    borehole_depth_library_lower = 48
    borehole_radius_library_lower = 0.075

    borehole_depth_library_upper = 96 
    borehole_radius_library_upper = 075

elseif borehole_depth >= 96 && borehole_depth <= 192
    borehole_depth_library_lower = 96
    borehole_radius_library_lower = 0.075
    
    borehole_depth_library_upper = 192
    borehole_radius_library_upper = 0.08

elseif borehole_depth >= 192 && borehole_depth <= 384
    borehole_depth_library_lower = 192
    borehole_radius_library_lower = 0.08

    borehole_depth_library_upper = 384 
    borehole_radius_library_upper = 0.0875

end

# Attention. At the moment, only Open rectangle configurations are considered. 
# More library files can be downloaded from the open source library. User-Input: Probe-field geometry.
Library = JSON.parsefile("C:/Users/vollmer/Documents/GitHub/resie/src/energy_systems/heat_sources/Open_configurations_5m_v1.0.json")

# Create keys to read gridpoint-sets out of library-file.
key1 = "$(m)_$(n)"
key2 = "$(t)"

# Get a specific configuration from the library
configuration = Library[key1][key2]

# Get values from the configuration
gVals = configuration["g"]

# Getting a Specific G-Function for the default configuration
borehole_spacing_library_default = 5   # Default borehole spacing for eacht configuration.
g_values_library_lower = gVals["$borehole_spacing_library_default._$borehole_depth_library_lower._$borehole_radius_library_lower"]
g_values_library_upper = gVals["$borehole_spacing_library_default._$borehole_depth_library_upper._$borehole_radius_library_upper"]

# Linear Intepolation between upper and lower default g_values on B/H-ratio. B: Borehole spacing. H: Borehole depth.
B_H_ratio_configuration = borehole_spacing / borehole_depth
B_H_ratio_library_default_lower = borehole_spacing_library_default / borehole_depth_library_lower
B_H_ratio_library_default_upper = borehole_spacing_library_default / borehole_depth_library_upper

# Substitution "Interpolation_Faktor" to save space. Describes, if B_H_ratio_configuration is closer to B_H_ratio_library_default_lower or B_H_ratio_library_default_upper
Interpolation_Faktor = (B_H_ratio_configuration-B_H_ratio_library_default_lower)/(B_H_ratio_library_default_upper-B_H_ratio_library_default_lower)

g_values_library_interpolated = g_values_library_upper-(1-Interpolation_Faktor)*(g_values_library_upper-g_values_library_lower)

# g_values_library_interpolated are refered to a not existing boreholeradius (borehole_radius_interpolated), because Interpolation is based on B/H.
# Within the interpolation based on B/H, the ratio r_b/H is interpolated, too. That's why a correction is needed to refer to the real borehole radius (which is set by user or default).
r_b_H_ratio_library_default_lower = borehole_radius_library_lower/borehole_depth_library_lower
r_b_H_ratio_library_default_upper = borehole_radius_library_upper/borehole_depth_library_upper

borehole_radius_interpolated = (r_b_H_ratio_library_default_upper+(1-Interpolation_Faktor)*
(r_b_H_ratio_library_default_lower-r_b_H_ratio_library_default_upper)) * borehole_depth 

g_values_library_corrected = g_values_library_interpolated .- log(borehole_radius/(borehole_radius_interpolated))

return g_values_library_corrected

end

# 
# Set Parameters
soil_density = 2000                     # [kg/m^3]
soil_specific_heat_capacity = 2400      # [J/(kg K)]
soil_heat_conductivity = 1.6            # [W/(m K)]
soil_diffusivity = soil_heat_conductivity / (soil_density * soil_specific_heat_capacity)    # [m^2/s]

borehole_radius = 0.08                 # [m]
borehole_depth = 150                    # [m]
borehole_spacing = 5                    # [m] distance between boreholes in the field, assumed to be constant. Set average spacing. 

steady_state_time = borehole_depth^2/(9* soil_diffusivity)  # [s]

simulation_time = 31536000              # [s] # 31536000 s = 1 year
simulation_step_width = 3600            # [s]

g_time_vector = 1:simulation_step_width:simulation_time

# User Input: Set n,m,t for Library geometrie properties. Currently, only open rectangle configuration is implemented.  
n = 10
m = 10
t = 1

# Get pre-caluclated g-function values out of open source library json file
# time stamps are the same for each set of grid points out of the library, as follows:
library_grid_points_normalized_time = [-8.5, -7.8, -7.2, -6.5, -5.9, -5.2, -4.5, 
    -3.963, -3.27, -2.864, -2.577, -2.171, -1.884, 
    -1.191, -0.497, -0.274, -0.051, 0.196, 0.419, 
    0.642, 0.873, 1.112, 1.335, 1.679, 2.028, 2.275, 3.003]

library_grid_points_g_values = get_library_g_values(borehole_depth,borehole_radius,borehole_spacing,m,n,t)
println(library_grid_points_g_values)

# Change normalized time values to absolut time and round to fit into simulation step width.
library_time_grid_points = round_to_simulation_step_width(ln_to_normal(library_grid_points_normalized_time,steady_state_time), simulation_step_width)

# Time-Interpolation
# Interpolation between grid points of library, to get g-function values for eacht time-step. Short time-step solution will be calculated later.
library_time_range = library_time_grid_points[1]:simulation_step_width:library_time_grid_points[end]  # time array, where g-values are available
library_time = collect(library_time_range) # make type array out of range

itp = LinearInterpolation(library_time_grid_points, library_grid_points_g_values,)
library_g_values_interpolated = itp(library_time)

library_time_slots_for_unknown_g_valulues = 0:simulation_step_width:(library_time[1] - 1)      # create time, where g-values are unknown
library_unknown_g_valulues = zeros(length(library_time_slots_for_unknown_g_valulues))   # create 0-values for unknown library_total_g_values (short time-steps)

library_total_time = 0:simulation_step_width:library_time_grid_points[end]             # total time incl. time solots for unknown g-values.
#library_total_g_values = [library_unknown_g_valulues..., library_g_values_interpolated...]    # zeroes until first grid point, then filled with interpolated g-values

# Working with views to avoid "OutOfStorage()" Error. Time efficient solution.
unknown_g_values_view = view(library_unknown_g_valulues, 1:length(library_time_slots_for_unknown_g_valulues))
interpolated_g_values_view = view(round.(library_g_values_interpolated, digits = 4) , 1:length(library_time))

# long timewidth g-function solution:
library_total_g_values = vcat(unknown_g_values_view, interpolated_g_values_view)

# caluclate g-function values for short timesteps with infinite linesource approach by Kelvin
g_function_short_time_step = infinite_line_source(borehole_radius, soil_diffusivity, library_total_time)

# Overlap

# find intersection between short time step solution and interpolated library values.
intersections_ils = findall(g_function_short_time_step[2:end] .== library_total_g_values[2:end])

if length(intersections_ils) > 0
    index_intersection = intersections_ils[1]
else
    start_fls = findall(library_total_g_values .> 0)
    index_intersection = start_fls[1]
end

# Final composite g-function solution.
g_function_composite = vcat(view(g_function_short_time_step, 1:index_intersection), view(library_total_g_values, (index_intersection+1):length(library_total_g_values)))

# Convert g_function_composite to a DataFrame
df = DataFrame(GValues = g_function_composite[1:round.(simulation_time/simulation_step_width)])
csv_file = "g_function.csv"

# Write the DataFrame to the CSV file
CSV.write(csv_file, df)

