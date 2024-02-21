"""
Implementation of geothermal probes.
This implementations acts as storage as it can produce and load energy.
"""

using JSON
using Interpolations
using SpecialFunctions
mutable struct GeothermalProbes <: Component
    uac::String
    controller::Controller
    sys_function::SystemFunction
    input_interfaces::InterfaceMap
    output_interfaces::InterfaceMap
    m_heat_in::Symbol
    m_heat_out::Symbol

    unloading_temperature_spread::Temperature
    loading_temperature::Temperature
    loading_temperature_spread::Temperature
    max_output_power::Float64
    max_input_power::Float64
    regeneration::Bool
    max_output_energy::Float64
    max_input_energy::Float64
    current_output_temperature::Temperature
    current_input_temperature::Temperature  
    soil_undisturbed_ground_temperature::Temperature
    soil_heat_conductivity::Float64
    soil_density::Float64
    soil_specific_heat_capacity::Float64
    borehole_thermal_resistance::Float64
    g_function::Vector{Float64}
    time_index::Int
    fluid_temperature::Temperature
    borehole_current_wall_temperature::Temperature

    energy_in_out_per_probe_meter::Vector{Float64}
    energy_in_out_difference_per_probe_meter::Vector{Float64}

    probe_depth::Float64
    number_of_probes::Int
    borehole_spacing::Float64
    probe_type::Int

    pipe_diameter_outer::Float64
    pipe_diameter_inner::Float64
    radius_pipe_inner::Float64
    radius_pipe_outer::Float64
    radius_borehole::Float64
    distance_pipe_center::Float64

    fluid_specific_heat_capacity::Float64
    fluid_density::Float64
    fluid_kinematic_viscosity::Float64
    fluid_heat_conductivity::Float64
    fluid_prandtl_number::Float64

    grout_heat_conductivity::Float64
    pipe_heat_conductivity::Float64

    borehole_diameter::Float64
    shank_spacing::Float64

    fluid_reynolds_number::Float64

    function GeothermalProbes(uac::String, config::Dict{String,Any}, sim_params::Dict{String,Any})
        m_heat_in = Symbol(default(config, "m_heat_in", "m_h_w_ht1"))
        m_heat_out = Symbol(default(config, "m_heat_out", "m_h_w_lt1"))
        register_media([m_heat_in, m_heat_out])
    
        # Current solution to get g-function values.
        # read g-function .txt file
        file = open("C:/Users/steinacker/Lokal/git_Resie/src/energy_systems/heat_sources/g_ges_vector.txt", "r") # TODO
            g_function = Vector{Float64}()
            for line in eachline(file)
                push!(g_function, parse(Float64, line))
            end
        close(file)

        return new(
            uac,                     # uac
            controller_for_strategy( # controller
                config["strategy"]["name"], config["strategy"], sim_params
            ),
            sf_storage,              # sys_function
            InterfaceMap(            # input_interfaces
                m_heat_in => nothing
            ),
            InterfaceMap(            # output_interfaces
                m_heat_out => nothing
            ),
            m_heat_in,                      # medium name of input interface
            m_heat_out,                     # medium name of output interface

            default(config, "unloading_temperature_spread", 3),   # temperature spread between forward and return flow during unloading, within one probe!
            default(config, "loading_temperature", nothing),      # nominal high temperature for loading geothermal probe storage, can also be set from other end of interface
            default(config, "loading_temperature_spread", 3),     # temperature spread between forward and return flow during loading, within one probe!
            default(config, "max_output_power", 50),  # maximum output power in W/m probe
            default(config, "max_input_power", 50),   # maximum input power in W/m probe
            default(config, "regeneration", true),    # flag if regeneration should be taken into account
            0.0,                                      # max_output_energy in every time step, calculated in control()
            0.0,                                      # max_input_energy in every time step, calculated in control()
            0.0,                                      # output temperature in current time step, calculated in control()
            0.0,                                      # input temperature in current time step, calculated in control()
            default(config, "soil_undisturbed_ground_temperature", 11.0),    # Considered as constant
            default(config, "soil_heat_conductivity", 1.5),                  # Heat conductivity of surrounding soil, homogenous and constant, in [W/(m K)]
            default(config, "soil_density", 2000.0),                         # soil density in [kg/m^3]
            default(config, "soil_specific_heat_capacity", 2400.0),          # soil specific heat capacity in [J/(kg K)]
            default(config, "borehole_thermal_resistance", 0.10),            # thermal resistance in [(m K)/W]
            g_function,                             # pre-calculated multiscale g-function. Calculated in pre-processing.
            0,                                      # index of current time step to get access on time dependent g-function values
            0.0,                                    # average fluid temperature
            default(config, "boreholewall_start_temperature", 4.0),      # boreholewall starting temperature
            
            [],                                     # vector to hold specific energy sum (in and out) per probe meter in each time step
            [],                                     # vector to hold specific energy per probe meter as difference for each step, used for g-function approach 

            default(config, "probe_depth", 150.0),  # depth (or length) of a single geothermal probe
            default(config, "number_of_probes", 1), # number of geothermal probes in the borefield
            default(config, "borehole_spacing", 1), # [m] distance between boreholes in the field, assumed to be constant. Set average spacing. 
            default(config, "probe_type", 2),       # probe type: 1: single U-pipe in one probe, 2: double U-pipe in one probe

            default(config, "pipe_diameter_outer", 0.032),  # outer pipe diameter
            default(config, "pipe_diameter_inner", 0.026),  # inner pipe diameter
            0.0,                                            # radius_pipe_inner, will be calculated in initialization
            0.0,                                            # radius_pipe_outer, will be calculated in initialization
            0.0,                                            # radius_borehole, will be calculated in initialization
            0.0,                                            # distance_pipe_center, will be calculated in initialization
                    
            default(config, "fluid_specific_heat_capacity", 3800.0), # specific heat capacity brine at 0 °C (25 % glycol 75 % water (interpolated)) 
            default(config, "fluid_density", 1045.0),                # density brine at 0 °C (25 % glycol 75 % water (interpolated))
            default(config, "fluid_kinematic_viscosity", 3.9e-6),    # viscosity brine at 0 °C (25 % glycol 75 % water (interpolated)) 
            default(config, "fluid_heat_conductivity", 0.5) ,        # heat conductivity brine at 0 °C (25 % glycol 75 % water (interpolated))
            default(config, "fluid_prandtl_number", 30.0),           # prandtl-number brine at 0 °C (25 % glycol 75 % water (interpolated)) 
            
            default(config, "grout_heat_conductivity", 2.0),         # lambda grout / filling material in W/(mK)   
            default(config, "pipe_heat_conductivity", 0.42),         # lambda of inner pipes

            default(config, "borehole_diameter", 0.15),              # borehole diameter in m.
            default(config, "shank_spacing", 0.1),                   # shank-spacing = distance between inner pipes in borehole, diagonal through borehole center. required for calculation of thermal borehole resistance.

            0.0        # Reynoldsnumber. To be calculated in function later.
            )
    end
end

function initialise!(unit::GeothermalProbes, sim_params::Dict{String,Any})
    set_storage_transfer!(
        unit.input_interfaces[unit.m_heat_in],
        default(
            unit.controller.parameter, "unload_storages " * String(unit.m_heat_in), true
        )
    )
    set_storage_transfer!(
        unit.output_interfaces[unit.m_heat_out],
        default(
            unit.controller.parameter, "load_storages " * String(unit.m_heat_out), true
        )
    )

    # calculate and initialize constant variables
    unit.energy_in_out_per_probe_meter = zeros(sim_params["number_of_time_steps"])
    unit.energy_in_out_difference_per_probe_meter = zeros(sim_params["number_of_time_steps"])

    unit.radius_pipe_inner = unit.pipe_diameter_inner / 2            # [m]
    unit.radius_pipe_outer = unit.pipe_diameter_outer / 2            # [m]
    unit.radius_borehole = unit.borehole_diameter / 2                # [m]
    unit.distance_pipe_center = unit.shank_spacing / 2               # [m]

    unit.max_output_energy = watt_to_wh(unit.max_output_power * unit.probe_depth * unit.number_of_probes)
    unit.max_input_energy = watt_to_wh(unit.max_input_power * unit.probe_depth * unit.number_of_probes)

    # calculate g-function
    # The long term g-functions will be read out from the library provided by
    # G-Function Library V1.0 by T. West, J. Cook and D. Spitler, 2021, Oklahoma State University
    # The missing short term g-function will be calculated using the infinite line source/sink theory by Kelvin
    #
    # set Parameters
    # User Input: Set key_1 and key_2 for library geometrie properties. TODO: move to input
    m = 10
    n = 10
    key1 = "$(m)_$(n)"      # Note that  m <= n!
    key2 = "1"
    libfile_path = "c:/Users/steinacker/Lokal/g-function_library_1.0/Open_configurations_5m_v1.0.json"

    # calculate g-functions
    soil_diffusivity = unit.soil_heat_conductivity / (unit.soil_density * unit.soil_specific_heat_capacity)    # [m^2/s]
    steady_state_time = unit.probe_depth^2 / (9 * soil_diffusivity)  # [s]

    # Get pre-caluclated g-function values out of g-function library json file
    # time stamps are the same for each set of grid points out of the library, as follows:
    library_time_grid_normalized = [-8.5,   -7.8,   -7.2,   -6.5,   -5.9,   -5.2,   -4.5,
                                    -3.963, -3.27,  -2.864, -2.577, -2.171, -1.884, -1.191,
                                    -0.497, -0.274, -0.051,  0.196,  0.419,  0.642,  0.873,
                                     1.112,  1.335,  1.679,  2.028,  2.275,  3.003]

    g_values_long_term_library = get_library_g_values(unit, unit.probe_depth, unit.radius_borehole, unit.borehole_spacing, key1, key2, libfile_path)

    # Change normalized time values to absolut time and round to fit into simulation step width.
    library_time_grid_absolute = round_to_simulation_step_width(ln_to_normal(library_time_grid_normalized, steady_state_time), sim_params["time_step_seconds"])

    # Time-Interpolation
    # Interpolation between grid points of library, to get g-function values for eacht time-step. Short time-step solution will be calculated later.
    simulation_end_timestamp = Int(sim_params["number_of_time_steps"] * sim_params["time_step_seconds"])
    simulation_time_grid_precalculated = collect(library_time_grid_absolute[1]:sim_params["time_step_seconds"]:simulation_end_timestamp)  # time array, where g-values are available

    itp = Interpolations.LinearInterpolation(library_time_grid_absolute, g_values_long_term_library, )
    g_values_long_term_library_interpolated = round.(itp(simulation_time_grid_precalculated), digits=4)

    # caluclate g-function values for short timesteps with infinite linesource approach by Kelvin
    simulation_time_grid_total = 0:sim_params["time_step_seconds"]:simulation_end_timestamp  # total time incl. time slots for unknown g-values.
    g_values_short_term = infinite_line_source(unit.radius_borehole, soil_diffusivity, simulation_time_grid_total)

    # fill long term g-function with missing short term values with zeros
    time_index_starting_long_term_data = length(g_values_short_term) - length(g_values_long_term_library_interpolated) + 1
    g_values_long_term = vcat(zeros(time_index_starting_long_term_data - 1), g_values_long_term_library_interpolated)

    # Find intersection between short term and long term g-function values.
    intersection_found = false
    merge_index = time_index_starting_long_term_data
    for i in time_index_starting_long_term_data:(length(g_values_long_term)-1)
        if (g_values_long_term[i] <= g_values_short_term[i] &&
            g_values_long_term[i+1] >= g_values_short_term[i+1]
        )
            intersection_found = true
            merge_index = i
            break
        end
    end

    unit.g_function = copy(g_values_short_term)
    if intersection_found
        # Merge at the intersection point
        unit.g_function[merge_index:end] = g_values_long_term[merge_index:end]
    else
        # Merge at the first non-zero point if no intersection was found
        unit.g_function[time_index_starting_long_term_data:end] = g_values_long_term[time_index_starting_long_term_data:end]
    end

    # # Convert unit.g_function to a DataFrame and write to file (add using DataFrame)
    # df = DataFrame(GValues = unit.g_function[1:sim_params["number_of_time_steps"]])
    # CSV.write("g_function.csv", df)

    # # create and save plots (add using Plots)
    # plot(1:length(g_values_long_term), g_values_short_term, label="short term g-function", linewidth=2)
    # plot!(1:length(g_values_long_term), g_values_long_term, label="long term g-function", linewidth=2)
    # plot!(1:length(g_values_long_term), unit.g_function, label="merged g-function", linewidth=2)
    # savefig("output/g_funcion.png")

    @info "Successfully calculated g-function for geothermal probe field $(unit.uac)"

end


"""
    infinite_line_source()

Calculates the g-function of an infinite line source/sink by Kelvin as short-term 
responses of a single borehole.

Inputs:
    borehole_radius::Float64                    - radius of the borehole in [m]
    soil_diffusivity::Float64                   - soil diffusivity around the borehole in [m^2/s]    
    g_time_vector::StepRange{UInt64, UInt64}    - StepRange with time index of requested short-term g-function
Output:
    g_ILS::Vector{Floats64}                     - g-function values corresponding to g_time_vector for a single borehole

"""
function infinite_line_source(borehole_radius::Float64, soil_diffusivity::Float64, g_time_vector::StepRange{UInt64,UInt64})
    g_ILS = zeros(length(g_time_vector))

    for i in eachindex(g_time_vector)
        if g_time_vector[i] == 0    # if time vector starts at 0, g-function will be set to 0.
            g_ILS[i] = 0
        else
            g_ILS[i] = round.(0.5 * expint(borehole_radius^2 / (4 * soil_diffusivity * g_time_vector[i])), digits=4)
        end
    end

    return g_ILS
end

"""
    ln_to_normal()

Calculates the absolute from normalized time values as used in library by T. West, J. Cook and D. Spitler

Inputs:
    library_time_grid_normalized::Array{Float64}    - array of normalized time from West/Cook/Spitler library, 
                                                      given as ln(absolute_time/steady_state_time) 
    steady_state_time::Float64                      - steady state time = unit.probe_depth^2/(9* soil_diffusivity) in [s]
Outputs:
    time::Array{Float64}                            - array of absolute time calculated from the normlized time as
                                                      exp(library_time_grid_normalized)*steady_state_time
"""
function ln_to_normal(library_time_grid_normalized::Array{Float64}, steady_state_time::Float64)
    n = length(library_time_grid_normalized)     # number of grid points of spitler/cook data set.
    library_time_grid_absolute = zeros(n)

    for i in 1:n
        library_time_grid_absolute[i] = round(exp(library_time_grid_normalized[i]) * steady_state_time)
    end

    return library_time_grid_absolute
end

"""
    round_to_simulation_step_width()

rounds given time steps (from West/Cook/Spitler library) to nearest multiple of the simulation time step.

Inputs:
    grid_points_time::Vector{Fl at64}             - vector of absolute time steps as given in the West/Cook/Spitler library in [s]
    simulation_time_step_width::UInt64            - desired simulation time step in [s]
Outputs:
    grid_points_time_rounded:: Vector{Float64}    - vector of time steps of grid_points_time, rounded to the nearest time that
                                                    is a multiple of the simulation_time_step_width 

"""
function round_to_simulation_step_width(grid_points_time::Vector{Float64}, simulation_time_step_width::UInt64)
    grid_points_time_rounded = zeros(length(grid_points_time))

    for i = 1:length(grid_points_time)
        grid_points_time_rounded[i] = Int(round(grid_points_time[i] / simulation_time_step_width) * simulation_time_step_width)
    end

    return grid_points_time_rounded
end

"""
    get_library_g_values()

Read long term g-function from the precalculated library provided by:
G-Function Library V1.0 by T. West, J. Cook and D. Spitler, 2021, Oklahoma State University
See: https://doi.org/10.15121/1811518

There are various sets of grid points refering to different borehole depths and borehole cofigurations.
Two sets of grid points have to be read out of one .json library-file to do interpolation, 
as they are refered to specific borehole depth.
Check, which two sets of grid points (refered to borehole depth borehole_depth_library_lower 
and borehole depth borehole_depth_library_upper) will be read out of json file later.
First, each set of grid points is refered to a default borehole radius 
(borehole_radius_library_lower, borehole_radius_library_upper). 
Later, the grid points will be corrected to the real borehole radius.

Inputs:
    unit::Component              - current unit
    borehole_depth::Float64      - borehole depth in [m]
    borehole_radius::Float64     - borehole radius in [m]
    borehole_spacing::Float64    - borehole spacing in [m], assumed to be equal for all boreholes in both directions
    key1::String                 - parameter of West/Cook/Spitler library, normaly the number of probes in each direction: "m_n" while n >= m! See dokumentation of library for details.
    key2::String                 - parameter of West/Cook/Spitler library, meaning depends of type of probefield geometry, can also be "" if none is needed. See dokumentation of library for details.
    libfile_path::String         - file path to the library. Note that different probefield geometries are stored in different files.
Outputs:
    g_values_library_corrected::Array{Float64}  - interpolated and corrected g-values read out of library for given key1 and key2

"""
function get_library_g_values(unit::Component, borehole_depth::Float64, borehole_radius::Float64, borehole_spacing::Float64, key1::String, key2::String, libfile_path::String)
    # borehole_depth is set by user. 
    # Check, which two sets of grid-points will be read out of the library-file to be interpolated later.
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

    # Read in given library file
    local library
    try
        library = JSON.parsefile(libfile_path)
    catch e
        @error "The library with precalculated g-values for the geothermal probe $(unit.uac) could not be read in. The following error occured: $e\n" *
                "Check the file located at $libfile_path."
        exit()
    end

    # Get a specific configuration from the library
    local gVals
    try
        if key2 == ""
            gVals = library[key1]["g"]
        else
            gVals = library[key1][key2]["g"]
        end
    catch e
        @error "The probe field configuration for the geothermal probe  $(unit.uac) could not be detected from the given library. The following error occured: $e\n" *
        "Check the probe field configuration given as key1 = $key1 and key2 = $key2."
        exit()
    end

    # Getting a specific G-Function for the default configuration
    borehole_spacing_library_default = 5   # Default borehole spacing for each configuration.
    g_values_library_lower = gVals["$borehole_spacing_library_default._$borehole_depth_library_lower._$borehole_radius_library_lower"]
    g_values_library_upper = gVals["$borehole_spacing_library_default._$borehole_depth_library_upper._$borehole_radius_library_upper"]

    # Linear intepolation between upper and lower default g_values on B/H-ratio. B: Borehole spacing. H: Borehole depth.
    B_H_ratio_configuration = borehole_spacing / borehole_depth
    B_H_ratio_library_default_lower = borehole_spacing_library_default / borehole_depth_library_lower
    B_H_ratio_library_default_upper = borehole_spacing_library_default / borehole_depth_library_upper

    # Substitution "interpolation_factor". Describes, if B_H_ratio_configuration is closer to B_H_ratio_library_default_lower or B_H_ratio_library_default_upper
    interpolation_factor = (B_H_ratio_configuration - B_H_ratio_library_default_lower) / (B_H_ratio_library_default_upper - B_H_ratio_library_default_lower)

    g_values_library_interpolated = g_values_library_upper - (1 - interpolation_factor) * (g_values_library_upper - g_values_library_lower)

    # g_values_library_interpolated are refered to a not existing boreholeradius (borehole_radius_interpolated), because Interpolation is based on B/H.
    # Within the interpolation based on B/H, the ratio r_b/H is interpolated, too. That's why a correction is needed to refer to the real borehole radius (which is set by user or default).
    r_b_H_ratio_library_default_lower = borehole_radius_library_lower / borehole_depth_library_lower
    r_b_H_ratio_library_default_upper = borehole_radius_library_upper / borehole_depth_library_upper

    borehole_radius_interpolated = (r_b_H_ratio_library_default_upper + (1 - interpolation_factor) *
                                                                        (r_b_H_ratio_library_default_lower - r_b_H_ratio_library_default_upper)) * borehole_depth

    g_values_library_corrected = g_values_library_interpolated .- log(borehole_radius / (borehole_radius_interpolated))

    return g_values_library_corrected
end

function control(
    unit::GeothermalProbes,
    components::Grouping,
    sim_params::Dict{String,Any}
)
    # time index, necessary for g-function approach
    unit.time_index = unit.time_index + 1 

    # get output temperature for energy output and set temperature and max_energy to output interface
    unit.current_output_temperature = unit.fluid_temperature + unit.unloading_temperature_spread/2
    set_temperature!(unit.output_interfaces[unit.m_heat_out],
                     nothing,
                     unit.current_output_temperature
                     )

    # set max_energy to output interface to provide information for connected components
    set_max_energy!(unit.output_interfaces[unit.m_heat_out], unit.max_output_energy)

    # get input temperature for energy input (regeneration) and set temperature and max_energy to input interface
    if unit.regeneration
        unit.current_input_temperature = unit.fluid_temperature - unit.loading_temperature_spread/2 # of geothermal probe field 
        set_temperature!(unit.input_interfaces[unit.m_heat_in],
                         unit.current_input_temperature,
                         nothing
                         )
        
        set_max_energy!(unit.input_interfaces[unit.m_heat_in], unit.max_input_energy)
    end

end

# function to calculate current new boreholewall temperature with g-functions
function calculate_new_boreholewall_temperature!(unit::GeothermalProbes)

    # R_B with Hellström (thermal borehole resistance)
    # calculate convective heat transfer coefficient alpha in pipie
    alpha_fluid, unit.fluid_reynolds_number = calculate_alpha_pipe(unit::GeothermalProbes)

    # calculate effective thermal borehole resistance by multipole method (Hellström 1991) depending on alpha
    sigma = (unit.grout_heat_conductivity - unit.soil_heat_conductivity) / (unit.grout_heat_conductivity + unit.soil_heat_conductivity)   # dimensionless calculation factor
    beta = 1 / (2 * pi * alpha_fluid * unit.radius_pipe_inner) + 1 / (2 * pi * unit.pipe_heat_conductivity) * log(unit.radius_pipe_outer / unit.radius_pipe_inner) # in (mK)/W

    R_1 = beta + 1 / (2 * pi * unit.grout_heat_conductivity) *
                 (log(unit.radius_borehole^2 / (2 * unit.radius_pipe_outer * unit.distance_pipe_center)) +
                 sigma * log(unit.radius_borehole^4 / (unit.radius_borehole^4 - unit.distance_pipe_center^4)) -
                 unit.radius_pipe_outer^2 / (4 * unit.distance_pipe_center^2) * (1 - sigma * 4 * unit.distance_pipe_center^4 / (unit.radius_borehole^4 - unit.distance_pipe_center^4))^2 /
                 ((1 + 2 * pi * unit.grout_heat_conductivity * beta) / (1 - 2 * pi * unit.grout_heat_conductivity * beta) +
                 unit.radius_pipe_outer^2 / (4 * unit.distance_pipe_center^2) * (1 + sigma * 16 * unit.radius_borehole^4 * unit.distance_pipe_center^4 / ((unit.radius_borehole^4 - unit.distance_pipe_center^4)^2))))

    unit.borehole_thermal_resistance = R_1 / (2 * unit.probe_type)

    # calculate new average fluid temperature with g-function approach
    if unit.time_index == 1
        unit.energy_in_out_difference_per_probe_meter[unit.time_index] = unit.energy_in_out_per_probe_meter[unit.time_index]
    else
        unit.energy_in_out_difference_per_probe_meter[unit.time_index] = unit.energy_in_out_per_probe_meter[unit.time_index] - unit.energy_in_out_per_probe_meter[unit.time_index-1]
    end

    current_temperature_difference = sum(reverse(unit.energy_in_out_difference_per_probe_meter[1:unit.time_index]) .* unit.g_function[1:unit.time_index]) / (2 * pi * unit.soil_heat_conductivity)

    unit.borehole_current_wall_temperature = unit.soil_undisturbed_ground_temperature + current_temperature_difference
    unit.fluid_temperature = unit.borehole_current_wall_temperature + unit.energy_in_out_per_probe_meter[unit.time_index] * unit.borehole_thermal_resistance

end

function calculate_alpha_pipe(unit::GeothermalProbes)
    # calculate mass flow in pipe
    power_in_out_per_pipe = wh_to_watts(abs(unit.energy_in_out_per_probe_meter[unit.time_index])) * unit.probe_depth / unit.probe_type  # W/pipe
    temperature_spread = unit.energy_in_out_per_probe_meter[unit.time_index] > 0 ? unit.loading_temperature_spread : unit.unloading_temperature_spread
    mass_flow_per_pipe = power_in_out_per_pipe / (unit.fluid_specific_heat_capacity * temperature_spread)  # kg/s

    use_dynamic_fluid_properties = false
    if use_dynamic_fluid_properties
        # calculate reynolds-number based on dynamic viscosity using dynamic temperature-dependend fluid properties, adapted from TRNSYS Type 710:
        fluid_dynamic_viscosity = 0.0000017158* unit.fluid_temperature^2 - 0.0001579079*unit.fluid_temperature+0.0048830621
        unit.fluid_heat_conductivity = 0.0010214286 * unit.fluid_temperature + 0.447
        unit.fluid_prandtl_number = fluid_dynamic_viscosity * unit.fluid_specific_heat_capacity / unit.fluid_heat_conductivity 
        fluid_reynolds_number = (4 * mass_flow_per_pipe) / (fluid_dynamic_viscosity * unit.pipe_diameter_inner * pi)
    else 
        # calculate reynolds-number based on kinematic viscosity with constant fluid properties.
        fluid_reynolds_number = (4 * mass_flow_per_pipe) / (unit.fluid_density * unit.fluid_kinematic_viscosity * unit.pipe_diameter_inner * pi)
    end
    
    if fluid_reynolds_number <= 2300  # laminar
        Nu = calculate_Nu_laminar(unit, fluid_reynolds_number)
    elseif fluid_reynolds_number > 2300 && fluid_reynolds_number <= 1e4 # transitional
        # Gielinski 1995
        factor = (fluid_reynolds_number - 2300) / (1e4 - 2300)
        Nu = (1 - factor) * calculate_Nu_laminar(unit, 2300.0) +
             factor * calculate_Nu_turbulent(unit, 10_000.0)
    else  # turbulent
        Nu = calculate_Nu_turbulent(unit, fluid_reynolds_number)
    end

    alpha = Nu * unit.fluid_heat_conductivity / unit.pipe_diameter_inner

    return alpha, fluid_reynolds_number
end
    
function calculate_Nu_laminar(unit::GeothermalProbes, fluid_reynolds_number::Float64)
    # Approach used in Ramming 2007 from Elsner, Norbert; Fischer, Siegfried; Huhn, Jörg; „Grundlagen der Technischen Thermodynamik“,  Band 2 Wärmeübertragung, Akademie Verlag, Berlin 1993. 
    k_a = 1.1 - 1 / (3.4 + 0.0667 * unit.fluid_prandtl_number)
    k_n = 0.35 + 1 / (7.825 + 2.6 * sqrt(unit.fluid_prandtl_number))

    # calculate Nu-Number
    Nu_laminar = ((k_a / (1 - k_n) * (unit.fluid_prandtl_number * unit.pipe_diameter_inner * fluid_reynolds_number / unit.probe_depth)^k_n)^3 + 4.364^3)^(1 / 3)
    return Nu_laminar
end 

function calculate_Nu_turbulent(unit::GeothermalProbes, fluid_reynolds_number::Float64)
    # Approached used from Gnielinski in: V. Gnielinski: Ein neues Berechnungsverfahren für die Wärmeübertragung im Übergangsbereich zwischen laminarer und turbulenter Rohrströmung. Forsch im Ing Wes 61:240–248, 1995. 
    zeta = (1.8 * log(fluid_reynolds_number) - 1.5)^-2
    Nu_turbulent = (zeta / 8 * fluid_reynolds_number * unit.fluid_prandtl_number) /
                   (1 + 12.7 * sqrt(zeta / 8) * (unit.fluid_prandtl_number^(2 / 3) - 1)) 
    return Nu_turbulent
end

# process function that provides energy from the geothermal probes
# according to actual delivered or received energy
function process(unit::GeothermalProbes, sim_params::Dict{String,Any})
    # get actual required energy from output interface
    outface = unit.output_interfaces[unit.m_heat_out]
    exchanges = balance_on(outface, outface.target)
    energy_demanded = balance(exchanges) +
                      energy_potential(exchanges) +
                      (outface.do_storage_transfer ? storage_potential(exchanges) : 0.0)
    energy_available = unit.max_output_energy  # is positive

    # shortcut if there is no energy demanded
    if energy_demanded >= -sim_params["epsilon"]
        set_max_energy!(unit.output_interfaces[unit.m_heat_out], 0.0)    
        return
    end

    for exchange in exchanges
        demanded_on_interface = exchange.balance +
                                exchange.energy_potential +
                                (outface.do_storage_transfer ? exchange.storage_potential : 0.0)

        if demanded_on_interface >= -sim_params["epsilon"]
            continue
        end

        if (
            exchange.temperature_min !== nothing
            && exchange.temperature_min > unit.current_output_temperature
        )
            # we can only supply energy at a temperature at or below the tank's current
            # output temperature
            continue
        end

        used_heat = abs(demanded_on_interface)

        if energy_available > used_heat
            energy_available -= used_heat
            add!(outface, used_heat, unit.current_output_temperature)
        else
            add!(outface, energy_available, unit.current_output_temperature)
            energy_available = 0.0
        end
    end

    # write output heat flux into vector
    energy_delivered = -(unit.max_output_energy - energy_available)
    unit.energy_in_out_per_probe_meter[unit.time_index] = energy_delivered  / (unit.probe_depth * unit.number_of_probes) # from total energy to specific power of one single probe.
    
end

function load(unit::GeothermalProbes, sim_params::Dict{String,Any})
    inface = unit.input_interfaces[unit.m_heat_in]
    exchanges = balance_on(inface, inface.source)
    energy_available = balance(exchanges) +
                       energy_potential(exchanges) +
                       (inface.do_storage_transfer ? storage_potential(exchanges) : 0.0)
    energy_demand = unit.max_input_energy  # is positive

    # shortcut if there is no energy to be used
    if ( energy_available <= sim_params["epsilon"] ||
         !unit.regeneration)
        set_max_energy!(unit.input_interfaces[unit.m_heat_in], 0.0)
        calculate_new_boreholewall_temperature!(unit::GeothermalProbes)
        return
    end

    for exchange in exchanges
        exchange_energy_available = exchange.balance +
                                    exchange.energy_potential +
                                    (inface.do_storage_transfer ? exchange.storage_potential : 0.0)

        if exchange_energy_available < sim_params["epsilon"]
            continue
        end

        if (
            exchange.temperature_min !== nothing
                && exchange.temperature_min > unit.current_input_temperature
            || exchange.temperature_max !== nothing
                && exchange.temperature_max < unit.current_input_temperature
        )
            # we can only take in energy if it's at a higher/equal temperature than the
            # tank's upper limit for temperatures
            continue
        end

        if energy_demand > exchange_energy_available
            energy_demand -= exchange_energy_available
            sub!(inface, exchange_energy_available, unit.current_input_temperature)
        else
            sub!(inface, energy_demand, unit.current_input_temperature)
            energy_demand = 0.0
        end
    end

    # Add loaded specific heat flux to vector
    energy_taken = unit.max_input_energy - energy_demand
    unit.energy_in_out_per_probe_meter[unit.time_index] += energy_taken / (unit.probe_depth * unit.number_of_probes)
    
    # recalculate borehole temperature for next timestep
    calculate_new_boreholewall_temperature!(unit::GeothermalProbes)
    
end

function balance_on(
    interface::SystemInterface,
    unit::GeothermalProbes
)::Vector{EnergyExchange}

caller_is_input = unit.uac == interface.target.uac

    return [EnEx(
        balance=interface.balance,
        uac=unit.uac,
        energy_potential=0.0,
        storage_potential=caller_is_input ? - unit.max_input_energy : unit.max_output_energy,  # TODO is this to be assuemd as storage_potential?
        temperature_min=interface.temperature_min,
        temperature_max=interface.temperature_max,
        pressure=nothing,
        voltage=nothing,
    )]
end

function output_values(unit::GeothermalProbes)::Vector{String}
    return [string(unit.m_heat_in)*" IN",
            string(unit.m_heat_out)*" OUT",
            "TEMPERATURE_#NodeNum",
            "borehole_temperature",
            "fluid_temperature",
            "borehole_thermal_resistance",
            "fluid_reynolds_number",
            "T_out",
            "Q_out",
            "Q_in",
            "dT_monthly",
            "dT_hourly" ]
end

function output_value(unit::GeothermalProbes, key::OutputKey)::Float64
    if key.value_key == "IN"
        return calculate_energy_flow(unit.input_interfaces[key.medium])
    elseif key.value_key == "OUT"
        return calculate_energy_flow(unit.output_interfaces[key.medium])
    elseif startswith(key.value_key, "Temperature_")
        idx = parse(Int, split(key.value_key, "_")[2])
        if !(1 <= idx <= length(unit.temperature_field)) 
            throw(ArgumentError("Index \"$idx\" of requested temperature-output of geothermal probe field exeeds the number of available temperatur datapoints.")) 
        else
            return unit.temperature_field[idx]
        end
    elseif key.value_key =="borehole_temperature"
        return unit.borehole_current_wall_temperature
    elseif key.value_key =="fluid_temperature"
        return unit.fluid_temperature
    elseif key.value_key =="borehole_thermal_resistance"
        return unit.borehole_thermal_resistance
    elseif key.value_key =="fluid_reynolds_number"
        return unit.fluid_reynolds_number
    elseif key.value_key =="dT_monthly"
        return unit.dT_monthly
    elseif key.value_key =="dT_hourly"
        return unit.dT_hourly
    end
    throw(KeyError(key.value_key))
end


export GeothermalProbes