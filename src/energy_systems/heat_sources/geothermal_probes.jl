"""
Implementation of geothermal probes.
This implementations acts as storage as it can produce and load energy.
"""

using JSON
using Interpolations
using Plots: plot, scatter, savefig
import Plots
using Roots
mutable struct GeothermalProbes <: Component
    uac::String
    controller::Controller
    sys_function::SystemFunction
    input_interfaces::InterfaceMap
    output_interfaces::InterfaceMap
    m_heat_in::Symbol
    m_heat_out::Symbol

    model_type::String

    unloading_temperature_spread::Temperature
    loading_temperature::Temperature
    loading_temperature_spread::Temperature
    max_output_power::Float64
    max_input_power::Float64
    regeneration::Bool
    max_output_energy::Float64
    current_max_output_energy::Float64
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
    output_temperatures_last_timestep::Vector{Temperature}

    energy_in_out_per_probe_meter::Vector{Float64}
    power_in_out_difference_per_probe_meter::Vector{Float64}

    probe_field_geometry::String
    number_of_probes_x::Int
    number_of_probes_y::Int
    probe_field_key_2::String
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

    probe_coordinates::Union{Nothing, Vector{Vector{Float64}}}

    function GeothermalProbes(uac::String, config::Dict{String,Any}, sim_params::Dict{String,Any})
        m_heat_in = Symbol(default(config, "m_heat_in", "m_h_w_ht1"))
        m_heat_out = Symbol(default(config, "m_heat_out", "m_h_w_lt1"))
        register_media([m_heat_in, m_heat_out])
    
        # get model type from input file
        model_type = default(config, "model_type", "simplified")
        model_type_allowed_values = ["simplified", "detailed"]
        if !(model_type in model_type_allowed_values)
            @error "Undefined model type \"$(model_type)\" of unit \"$(uac)\". Has to be one of: $(model_type_allowed_values)."
            throw(InputError)
        end

        # get probe field geometry from input file
        probe_field_geometry = default(config, "probe_field_geometry", "rectangle")
        probe_field_geometry_allowed_values = ["rectangle",
                                               "open_rectangle",
                                               "zoned_rectangle",
                                               "U_configurations",
                                               "lopsided_U_configuration",
                                               "C_configuration",
                                               "L_configuration"]
        if !(probe_field_geometry in probe_field_geometry_allowed_values)
            @error "Undefined probe field configuration \"$(probe_field_geometry)\" of unit \"$(uac)\". Has to be one of: $(probe_field_geometry_allowed_values)."
            throw(InputError)
        end    
    
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

            model_type,                                           # model type. currently "simplified" with constant thermal borehole resistance and 
                                                                  # "detailed" with calculated thermal borehole resistance in every time step are available.

            default(config, "unloading_temperature_spread", 3),   # temperature spread between forward and return flow during unloading, within one probe!
            default(config, "loading_temperature", nothing),      # nominal high temperature for loading geothermal probe storage, can also be set from other end of interface
            default(config, "loading_temperature_spread", 3),     # temperature spread between forward and return flow during loading, within one probe!
            default(config, "max_output_power", 50),              # maximum output power in W/m probe
            default(config, "max_input_power", 50),               # maximum input power in W/m probe
            default(config, "regeneration", true),                # flag if regeneration should be taken into account
            0.0,                                                  # max_output_energy in every time step, calculated in control()
            0.0,                                                  # current_max_output_energy
            0.0,                                                  # max_input_energy in every time step, calculated in control()
            0.0,                                                  # output temperature in current time step, calculated in control()
            0.0,                                                  # input temperature in current time step, calculated in control()
            default(config, "soil_undisturbed_ground_temperature", 11.0),    # Considered as constant
            default(config, "soil_heat_conductivity", 1.5),                  # Heat conductivity of surrounding soil, homogenous and constant, in [W/(m K)]
            default(config, "soil_density", 2000.0),                         # soil density in [kg/m^3]
            default(config, "soil_specific_heat_capacity", 2400.0),          # soil specific heat capacity in [J/(kg K)]
            default(config, "borehole_thermal_resistance", 0.1),             # thermal resistance in [(m K)/W]. If set, borehole_thermal_resistance is constant! (default: 0.1 (m K)/W)
            [],                                                   # pre-calculated multiscale g-function. Calculated in pre-processing.
            0,                                                    # index of current time step to get access on time dependent g-function values
            0.0,                                                  # average fluid temperature
            default(config, "boreholewall_start_temperature", 4.0),          # boreholewall starting temperature
            [],                                                   # output_temperatures_last_timestep
            
            [],                                                   # vector to hold specific energy sum (in and out) per probe meter in each time step
            [],                                                   # vector to hold specific energy per probe meter as difference for each step, used for g-function approach 

            probe_field_geometry,                                 # type of probe field geometry, can be one of: rectangle, open_rectangle, zoned_rectangle, U_configurations, lopsided_U_configuration, C_configuration, L_configuration
            default(config, "number_of_probes_x", 1),             # number of probes in x direction, corresponds so m value of g-fuction library. Note that number_of_probes_x <= number_of_probes_y!
            default(config, "number_of_probes_y", 1),             # number of probes in x direction, corresponds so m value of g-fuction library. Note that number_of_probes_x <= number_of_probes_y!
            default(config, "probe_field_key_2", "") ,            # key2 of g-fuction library. Can also be "" if non is needed. The value depends on the chosen library type.
            default(config, "probe_depth", 150.0),                # depth (or length) of a single geothermal probe. Has to be between 24 m and 384 m.
            0,                                                    # number of geothermal probes in the borefield, determined in inizialize from borehole configuration
            default(config, "borehole_spacing", 5),               # [m] distance between boreholes in the field, assumed to be constant. Set average spacing. 
            default(config, "probe_type", 2),                     # probe type: 1: single U-pipe in one probe, 2: double U-pipe in one probe

            default(config, "pipe_diameter_outer", 0.032),        # outer pipe diameter
            default(config, "pipe_diameter_inner", 0.026),        # inner pipe diameter
            0.0,                                                  # radius_pipe_inner, will be calculated in initialization
            0.0,                                                  # radius_pipe_outer, will be calculated in initialization
            0.0,                                                  # radius_borehole, will be calculated in initialization
            0.0,                                                  # distance_pipe_center, will be calculated in initialization
                    
            default(config, "fluid_specific_heat_capacity", 3800.0), # specific heat capacity brine at 0 °C (25 % glycol 75 % water (interpolated)) 
            default(config, "fluid_density", 1045.0),                # density brine at 0 °C (25 % glycol 75 % water (interpolated))
            default(config, "fluid_kinematic_viscosity", 3.9e-6),    # viscosity brine at 0 °C (25 % glycol 75 % water (interpolated)) 
            default(config, "fluid_heat_conductivity", 0.5) ,        # heat conductivity brine at 0 °C (25 % glycol 75 % water (interpolated))
            default(config, "fluid_prandtl_number", 30.0),           # prandtl-number brine at 0 °C (25 % glycol 75 % water (interpolated)) 
            
            default(config, "grout_heat_conductivity", 2.0),         # lambda grout / filling material in W/(mK)   
            default(config, "pipe_heat_conductivity", 0.42),         # lambda of inner pipes

            default(config, "borehole_diameter", 0.15),              # borehole diameter in m.
            default(config, "shank_spacing", 0.1),                   # shank-spacing = distance between inner pipes in borehole, diagonal through borehole center. required for calculation of thermal borehole resistance.

            0.0,                                                     # Reynoldsnumber. To be calculated in function later.
            nothing                                                  # probe_coordinates to save them for plotting
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
    unit.power_in_out_difference_per_probe_meter = zeros(sim_params["number_of_time_steps"])
    unit.fluid_temperature = unit.borehole_current_wall_temperature
    unit.output_temperatures_last_timestep = [unit.borehole_current_wall_temperature]

    unit.radius_pipe_inner = unit.pipe_diameter_inner / 2            # [m]
    unit.radius_pipe_outer = unit.pipe_diameter_outer / 2            # [m]
    unit.radius_borehole = unit.borehole_diameter / 2                # [m]
    unit.distance_pipe_center = unit.shank_spacing / 2               # [m]

    # calculate g-function and get number of probes in the probe field
    unit.g_function, unit.number_of_probes, unit.probe_coordinates = calculate_g_function(unit, sim_params)

    # calculate max energies
    unit.max_output_energy = watt_to_wh(unit.max_output_power * unit.probe_depth * unit.number_of_probes)
    unit.max_input_energy = watt_to_wh(unit.max_input_power * unit.probe_depth * unit.number_of_probes)

end

"""
    calculate_g_function(unit::Component, sim_params::Dict{String,Any})

Calculates the g-function values.
The long term g-functions will be read out from the library provided by G-Function Library V1.0 
by T. West, J. Cook and D. Spitler, 2021, Oklahoma State University with linear interpolation.
The missing short term g-function will be calculated using the infinite line source/sink 
theory by Kelvin.

Returns:
    `g_function::Vector{Float64}`: The g-function values in the simulation time step and duration
    `number_of_probes::Int`: The total number of probes in the geothermal probe field as defined
                             through the input parameters.
    `probe_coordinates::Vector{Vector{Float64}}`: The coordinates of all probes in the probe field in m.
                         
"""
function calculate_g_function(unit::Component, sim_params::Dict{String,Any})
    # set Parameters (user input): Set key_1 and key_2 for library geometrie properties.
    key1 = "$(unit.number_of_probes_x)_$(unit.number_of_probes_y)"      # Note that  x <= y!
    key2 = unit.probe_field_key_2
    probe_field_configurations = Dict("rectangle" => "rectangle_5m_v1.0.json",
                                      "open_rectangle"=> "Open_configurations_5m_v1.0.json",
                                      "zoned_rectangle" => "zoned_rectangle_5m_v1.0.json",
                                      "U_configurations" => "U_configurations_5m_v1.0.json",
                                      "lopsided_U_configuration" => "LopU_configurations_5m_v1.0.json",
                                      "C_configuration" => "C_configurations_5m_v1.0.json",
                                      "L_configuration" => "L_configurations_5m_v1.0.json"
                                      )

    libfile_path = "src/energy_systems/heat_sources/g-function_library_1.0/" * probe_field_configurations[unit.probe_field_geometry]
    if !isfile(libfile_path)
        @error "The library for the geothermal probe field \"$(unit.uac)\" could not be found at $libfile_path"
        throw(InputError)
    end

    # calculate g-functions
    soil_diffusivity = unit.soil_heat_conductivity / (unit.soil_density * unit.soil_specific_heat_capacity)    # [m^2/s]
    steady_state_time = unit.probe_depth^2 / (9 * soil_diffusivity)  # [s]

    # Get pre-caluclated g-function values out of g-function library json file
    # time stamps are the same for each set of grid points out of the library, as follows:
    library_time_grid_normalized = [-8.5,   -7.8,   -7.2,   -6.5,   -5.9,   -5.2,   -4.5,
                                    -3.963, -3.27,  -2.864, -2.577, -2.171, -1.884, -1.191,
                                    -0.497, -0.274, -0.051,  0.196,  0.419,  0.642,  0.873,
                                     1.112,  1.335,  1.679,  2.028,  2.275,  3.003]

    g_values_long_term_library, number_of_probes, probe_coordinates = get_library_g_values(unit, 
                                                                                           unit.probe_depth,
                                                                                           unit.radius_borehole,
                                                                                           unit.borehole_spacing,
                                                                                           key1,
                                                                                           key2,
                                                                                           libfile_path)

    # Change normalized time values to absolut time and round to fit into simulation step width.
    library_time_grid_absolute = round_to_simulation_step_width(ln_to_normal(library_time_grid_normalized, steady_state_time), 
                                                                sim_params["time_step_seconds"])

    # Interpolation between grid points of library, to get g-function values for eacht time-step.
    simulation_end_timestamp = Int(sim_params["number_of_time_steps"] * sim_params["time_step_seconds"])
    if simulation_end_timestamp < library_time_grid_absolute[2]  # make sure that the log() function can be parameterised
        simulation_end_timestamp = library_time_grid_absolute[2] + sim_params["time_step_seconds"]
    end
    simulation_time_grid_precalculated = collect(0:sim_params["time_step_seconds"]:simulation_end_timestamp)

    # linear interpolation. 
    # May change later to something else than linear interpolation between irregular grid points?
    # Though spline interpolation from Dierckx does not really work here...
    itp = interpolate((vcat(0,library_time_grid_absolute),), vcat(0,g_values_long_term_library), Gridded(Linear()))
    g_function = itp.(simulation_time_grid_precalculated)

    # change first time period of g-function to fitted f(x)=a*log(b*x) function to better represent short term effects.
    # This follows the method of an infitine line source that is widely used to calculate the short term g-functions
    # of a single probe. As for short time periods, the interferrence between the probes in the field is assumed to be
    # very small, therefore this approach is commonly assumed to be valid.
    # f(x) is fitted in a way that the function goes through the first and second node of the existing non-interpolated g-function.
    # therefore: a = y1/log(x1*exp(log(x1^y2/x2^y1)/(y1 - y2))) and b = exp(log(x1^y2/x2^y1)/(y1 - y2))
    x1 = Int(library_time_grid_absolute[1] / sim_params["time_step_seconds"]) + 1  # in time step width
    y1 = g_function[x1]
    x2 = Int(library_time_grid_absolute[2] / sim_params["time_step_seconds"]) + 1  # in time step width
    y2 = g_function[x2]
    a = y1/log(x1*exp(log(x1^y2/x2^y1)/(y1 - y2)))
    b = exp(log(x1^y2/x2^y1)/(y1 - y2))

    for x in 1:(x1-1)
        g_function[x] = max(0, a * log(b * (x+1)))
    end

    @info "Successfully calculated g-function for geothermal probe field \"$(unit.uac)\" with $(number_of_probes) probes."

    return g_function, number_of_probes, probe_coordinates
end

function plot_optional_figures(unit::GeothermalProbes, output_path::String, output_formats::Vector{Any}, sim_params::Dict{String,Any})
    # plot g-function values during simulation time
    plot(0:Int(sim_params["number_of_time_steps"]),
         unit.g_function[1:Int(sim_params["number_of_time_steps"]+1)],
         title="g-fuction values for geothermal probe field \"$(unit.uac)\"",
         xlabel="time [time step]",
         ylabel="g-function value [-]",
         legend=false,
         linewidth=6,
         gridlinewidth=1,
         size=(1800, 1200),
         titlefontsize=30,
         guidefontsize=24,
         tickfontsize=24,
         legendfontsize=24,
         grid=true,
         minorgrid=true,
         margin=15Plots.mm)
    fig_name = "probe_field_g_funcion_$(unit.uac)"
    for output_format in output_formats
        savefig(output_path * "/" *  fig_name * "." *  output_format)
    end

    # plot probe field configuration
    borehole_spacing_library_default = 5
    scatter([v[1] for v in values(unit.probe_coordinates)] * (unit.borehole_spacing / borehole_spacing_library_default),
            [v[2] for v in values(unit.probe_coordinates)] * (unit.borehole_spacing / borehole_spacing_library_default),
            title="probe field configuration for \"$(unit.uac)\"",
            xlabel = "distribution in x direction [m]",
            ylabel = "distribution in y direction [m]",
            aspect_ratio=:equal,
            legend=false,
            markersize=6,
            gridlinewidth=1,
            size=(1800, 1200),
            titlefontsize=30,
            guidefontsize=24,
            tickfontsize=24,
            legendfontsize=24,
            grid=true,
            minorgrid=true,
            margin=10Plots.mm)
    fig_name = "probe_field_geometry_$(unit.uac)"
    for output_format in output_formats
        savefig(output_path * "/" *  fig_name * "." *  output_format)
    end
    # throw not needed values away
    unit.probe_coordinates = nothing

    return true
end

"""
    ln_to_normal()

Calculates the absolute from normalized time values as used in library by T. West, J. Cook and D. Spitler

Inputs:
    library_time_grid_normalized::Array{Float64}    - array of normalized time from West/Cook/Spitler library, 
                                                      given as ln(absolute_time/steady_state_time) 
    steady_state_time::Float64                      - steady state time = unit.probe_depth^2/(9* soil_diffusivity) in [s]
Returns:
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
Returns:
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
This function checks, which two sets of grid points (refered to borehole depth borehole_depth_library_lower 
and borehole depth borehole_depth_library_upper) are read out of json file.
First, each set of grid points is refered to a default borehole radius 
(borehole_radius_library_lower, borehole_radius_library_upper). Later, the grid points are
corrected to the real borehole radius.

Inputs:
    unit::Component              - current unit
    borehole_depth::Float64      - borehole depth in [m]
    borehole_radius::Float64     - borehole radius in [m]
    borehole_spacing::Float64    - borehole spacing in [m], assumed to be equal for all boreholes in both directions
    key1::String                 - parameter of West/Cook/Spitler library, normaly the number of probes in each direction: "m_n" while n >= m! See dokumentation of library for details.
    key2::String                 - parameter of West/Cook/Spitler library, meaning depends of type of probefield geometry, can also be "" if none is needed. See dokumentation of library for details.
    libfile_path::String         - file path to the library. Note that different probefield geometries are stored in different files.
Returns:
    g_values_library_corrected::Array{Float64} - interpolated and corrected g-values read out of library for given key1 and key2
    number_of_probes::Int                      - number of probes in geothermal probe field, defined by key1 and key2
    probe_coordinates::Vector{Vector{Float64}} - The coordinates of all probes in the probe field in m.

"""
function get_library_g_values(unit::Component, borehole_depth::Float64, borehole_radius::Float64, borehole_spacing::Float64, key1::String, key2::String, libfile_path::String)
    # Check, which two sets of grid-points will be read out of the library-file to be interpolated later.
    # Define the ranges and corresponding values in an array of tuples
    borehole_ranges = [ (24, 0.075),
                        (48, 0.075),
                        (96, 0.075),
                        (192, 0.08),
                        (384, 0.0875)]  # defines the borehole depth [m] and corresponding borehole radius [m] of the library. 
                                        # They are the same for all configurations.

    depths = first.(borehole_ranges)
    upper_index = findfirst(x -> x > borehole_depth, depths)
    lower_index = upper_index - 1
    
    # If borehole_depth is less than the first depth or beyond all known ranges
    if upper_index === nothing || upper_index == 1
        @error "In geothermal probe \"$(unit.uac)\", the borehole_depth needs to be between $(depths[1]) m and $(depths[end]) m, but is $(borehole_depth) m."
        throw(InputError)
    end
    
    # Extracting the lower and upper bounds details
    borehole_depth_library_lower = borehole_ranges[lower_index][1]
    borehole_depth_library_upper = borehole_ranges[upper_index][1]
    borehole_radius_library_lower = borehole_ranges[lower_index][2]
    borehole_radius_library_upper = borehole_ranges[upper_index][2]

    # Read in given library file
    local library
    try
        library = JSON.parsefile(libfile_path)
    catch e
        @error "The library with precalculated g-values for the geothermal probe \"$(unit.uac)\" could not be read in. The following error occured: $e\n" *
                "Check the file located at $libfile_path."
        throw(InputError)
    end

    # Get a specific configuration from the library
    local gVals, probe_coordinates::Vector{Vector{Float64}}
    try
        if key2 == ""
            gVals = library[key1]["g"]
            probe_coordinates = library[key1]["bore_locations"]
        else
            gVals = library[key1][key2]["g"]
            probe_coordinates = library[key1][key2]["bore_locations"]
        end
    catch e
        @error "The probe field configuration for the geothermal probe  \"$(unit.uac)\" could not be detected from the given library. The following error occured: $e\n" *
        "Check the probe field configuration given as key1 = $key1 and key2 = $key2."
        throw(InputError)
    end

    # calculate number of probes in probe field
    number_of_probes = length(probe_coordinates)

    # Getting a specific G-Function for the default configuration
    borehole_spacing_library_default = 5   # Default borehole spacing for each configuration. Attention: Also change in plot_optional_figures!
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

    return g_values_library_corrected, number_of_probes, probe_coordinates
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

    # get requested temperature of the probe output, if they are written and present, to set the desired output 
    # temperature. If the output interface does not provide a temperature, the highest requested temperature of the last
    # timestep will be used to calculate the new max_energy. The max_energy is limited to the user-input of the max_energy
    # and has only an effect if the probe field reaches the lower limit temperature while unloading.
    exchange_temperatures = temp_min_all_non_empty(balance_on(unit.output_interfaces[unit.m_heat_out], unit.output_interfaces[unit.m_heat_out].target))
    desired_output_temperature = highest(exchange_temperatures[exchange_temperatures .< unit.current_output_temperature])

    if desired_output_temperature === nothing  # desired temperature not given in output interface, use temperature of last timestep
        desired_output_temperature = highest(unit.output_temperatures_last_timestep)
    end 
    if desired_output_temperature !== nothing
        desired_output_temperature += 0.001  # add small number to avoid problems with rough tolerances of find_zero()
    end 

    local max_energy_in_out_per_probe_meter

    # calculate the minimal temperature in the probe field that could occur in the current time step.
    unit.energy_in_out_per_probe_meter[unit.time_index] = -unit.max_output_energy /(unit.probe_depth * unit.number_of_probes)
    current_min_fluid_temperature = calculate_new_fluid_temperature(unit)
    unit.energy_in_out_per_probe_meter[unit.time_index] = 0.0 

    if desired_output_temperature !== nothing && (current_min_fluid_temperature + unit.unloading_temperature_spread/2) < desired_output_temperature
        try
            max_energy_in_out_per_probe_meter = find_zero((energy_in_out_per_probe_meter -> find_max_output_energy(energy_in_out_per_probe_meter, unit, desired_output_temperature)),
                                                        -unit.max_output_energy /(unit.probe_depth * unit.number_of_probes), 
                                                        Order0();
                                                        atol = 0.001
                                                        )
        catch # handles desired_output_temperature === nothing and non-converging results of find_zero()
            max_energy_in_out_per_probe_meter = -unit.max_output_energy /(unit.probe_depth * unit.number_of_probes)
        end
    else
        max_energy_in_out_per_probe_meter = -unit.max_output_energy /(unit.probe_depth * unit.number_of_probes)
    end                                   

    unit.current_max_output_energy = - max_energy_in_out_per_probe_meter * (unit.probe_depth * unit.number_of_probes)
    unit.current_max_output_energy = min(max(0, unit.current_max_output_energy), unit.max_output_energy)

    set_max_energy!(unit.output_interfaces[unit.m_heat_out], unit.current_max_output_energy)

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

"""
    find_max_output_energy(energy_in_out_per_probe_meter::Float64, unit::GeothermalProbes, desired_output_temperature::Float64)

Helper function that calls calculate_new_fluid_temperature() with a given energy_in_out_per_probe_meter and
returns the difference of the resulting fluid temperature after applying energy_in_out_per_probe_meter to the probe field
compared to a given desired_output_temperature. 
This function can be used to find the energy that can be taken out of the probe field while not exceeding 
the given desired_output_temperature using a find_zero() algorithm.

Note that this function resets the unit.energy_in_out_per_probe_meter[unit.time_index] to zero.

Args:
    energy_in_out_per_probe_meter::Float64:   the energy input (positve) or output (negative) of the probe field (per probe meter!)        
    unit::GeothermalProbes :                  the unit struct of the geothermal probe with all its parameters
    desired_output_temperature::Float64       the desired temperature of the fluid temperature at the probe field output
Returns:
    temperature_difference::Float64:          the temperature difference between the new fluid output temperature of the
                                              geothermal probe field after applying energy_in_out_per_probe_meter and the 
                                              desired output temperature.
"""
function find_max_output_energy(energy_in_out_per_probe_meter::Float64, unit::GeothermalProbes, desired_output_temperature::Float64)
   # set energy_in_out_per_probe_meter to given energy_in_out_per_probe_meter to as vehicle 
   # to deliver the information to calculate_new_fluid_temperature(unit)
   unit.energy_in_out_per_probe_meter[unit.time_index] = energy_in_out_per_probe_meter

   new_fluid_temperature = calculate_new_fluid_temperature(unit)

    # reset energy_in_out_per_probe_meter.
   unit.energy_in_out_per_probe_meter[unit.time_index] = 0.0 

   return new_fluid_temperature + unit.unloading_temperature_spread/2 - desired_output_temperature
end

"""
    calculate_new_fluid_temperature(unit::GeothermalProbes)::Temperature

Calculates the current new boreholewall temperature using g-functions as step response functions in 
response to a given change in energy input or output compared to the last time step.

Uses parameters stored in unit::GeothermalProbes, specifically the energy_in_out_per_probe_meter that has to be 
a global variable within the geothermal probe and that is set in load() as positive value and in process() as negative value.
Note that this function has to be called at least once in every time step to retrieve the current fluid temperature.

To calculate the fluid tempeature from the borehole wall temperature, an approach by Hellström 1991 is used to
determine the effective thermal borehole resistance using the convective heat transfer coefficient within the pipe.:
Hellström, G. Ground Heat Storage [online]. Thermal Analyses of Duct Storage Systems. Theorie, 1991.

Args:
    unit::GeothermalProbes:         the unit struct of the geothermal probe with all its parameters
Returns:
    fluid_temperature::Float64:     the new average temperature of the fluid in the geothermal probe after an energy in- or output 
                                    stored in unit.energy_in_out_per_probe_meter
"""
function calculate_new_fluid_temperature(unit::GeothermalProbes)
    if unit.model_type == "detailed"
        # R_B with Hellström (thermal borehole resistance)
        # calculate convective heat transfer coefficient alpha in pipe
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

        borehole_thermal_resistance = R_1 / (2 * unit.probe_type)
    elseif unit.model_type == "simplified"
        # constant borehole thermal resistance from user input
        borehole_thermal_resistance = unit.borehole_thermal_resistance
    end

    # calculate new average fluid temperature with g-function approach
    if unit.time_index == 1
        unit.power_in_out_difference_per_probe_meter[unit.time_index] = wh_to_watts(unit.energy_in_out_per_probe_meter[unit.time_index])
    else
        unit.power_in_out_difference_per_probe_meter[unit.time_index] = wh_to_watts(unit.energy_in_out_per_probe_meter[unit.time_index] - unit.energy_in_out_per_probe_meter[unit.time_index-1])
    end

    current_temperature_difference = sum(reverse(unit.power_in_out_difference_per_probe_meter[1:unit.time_index]) .* unit.g_function[1:unit.time_index]) / (2 * pi * unit.soil_heat_conductivity)

    borehole_current_wall_temperature = unit.soil_undisturbed_ground_temperature + current_temperature_difference

    fluid_temperature = borehole_current_wall_temperature + wh_to_watts(unit.energy_in_out_per_probe_meter[unit.time_index]) * borehole_thermal_resistance
    return fluid_temperature
end

"""
    calculate_alpha_pipe(unit::GeothermalProbes)

Calculates the convective heat transfer coefficient alpha in the pipe of the probe.
Uses parameters stored in unit::GeothermalProbes to calculate the mass flow in the pipe,
the Reynolds number, the Nusselt number for laminar or turbulent flow and finally alpha. 

Args:
    unit::GeothermalProbes:         the unit struct of the geothermal probe with all its parameters
Returns:
    alpha::Float64:                 convective heat transfer coefficient in pipe in current time step
    fluid_reynolds_number::Float64  the current Reynolds number of the fluid in the pipe
"""
function calculate_alpha_pipe(unit::GeothermalProbes)
    # calculate mass flow in pipe
    power_in_out_per_pipe = wh_to_watts(abs(unit.energy_in_out_per_probe_meter[unit.time_index])) * unit.probe_depth / unit.probe_type  # W/pipe
    temperature_spread = unit.energy_in_out_per_probe_meter[unit.time_index] > 0 ? unit.loading_temperature_spread : unit.unloading_temperature_spread
    mass_flow_per_pipe = power_in_out_per_pipe / (unit.fluid_specific_heat_capacity * temperature_spread)  # kg/s

    # calculate reynolds-number based on kinematic viscosity with constant fluid properties.
    fluid_reynolds_number = (4 * mass_flow_per_pipe) / (unit.fluid_density * unit.fluid_kinematic_viscosity * unit.pipe_diameter_inner * pi)

    # # In case that someone wants to implement temperature-dependend fluid properties in a later version, this already
    # # existing code can be used to start with. The equations are not checked and it is not known for which fluid the 
    # # parameters are valid. The reynolds-number is calculated based on dynamic viscosity using dynamic temperature-dependend 
    # # fluid properties. The method is adapted from TRNSYS Type 710:
    # fluid_dynamic_viscosity = 0.0000017158* unit.fluid_temperature^2 - 0.0001579079*unit.fluid_temperature+0.0048830621
    # unit.fluid_heat_conductivity = 0.0010214286 * unit.fluid_temperature + 0.447
    # unit.fluid_prandtl_number = fluid_dynamic_viscosity * unit.fluid_specific_heat_capacity / unit.fluid_heat_conductivity 
    # fluid_reynolds_number = (4 * mass_flow_per_pipe) / (fluid_dynamic_viscosity * unit.pipe_diameter_inner * pi)
    
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
 
"""
    calculate_Nu_laminar(unit::GeothermalProbes, fluid_reynolds_number::Float64)::Float64

Calculates the Nusselt number for laminar flow through the pipe of the geothermal probe.
Uses an approach from Ramming 2007, described in:
Elsner, Norbert; Fischer, Siegfried; Huhn, Jörg; „Grundlagen der Technischen Thermodynamik“,  Band 2 Wärmeübertragung, Akademie Verlag, Berlin 1993. 

Args:
    unit::GeothermalProbes:             the unit struct of the geothermal probe with all its parameters
    fluid_reynolds_number::Float64:     the current reynolds number in the pipe of the geothermal probe
Returns:
    Nusselt_laminar::Float64:           Nusselt number in current pipe for given Reynolds number and pipe dimensions
"""
function calculate_Nu_laminar(unit::GeothermalProbes, fluid_reynolds_number::Float64)
    k_a = 1.1 - 1 / (3.4 + 0.0667 * unit.fluid_prandtl_number)
    k_n = 0.35 + 1 / (7.825 + 2.6 * sqrt(unit.fluid_prandtl_number))

    # calculate Nu-Number
    Nu_laminar = ((k_a / (1 - k_n) * (unit.fluid_prandtl_number * unit.pipe_diameter_inner * fluid_reynolds_number / ((2 * unit.probe_type) * unit.probe_depth))^k_n)^3 + 4.364^3)^(1 / 3)
    return Nu_laminar
end 

"""
calculate_Nu_turbulent(unit::GeothermalProbes, fluid_reynolds_number::Float64)::Float64

Calculates the Nusselt number for turbulent flow through the pipe of the geothermal probe.
Uses an approach from Gnielinski described in:
V. Gnielinski: Ein neues Berechnungsverfahren für die Wärmeübertragung im Übergangsbereich zwischen laminarer und turbulenter Rohrströmung. Forsch im Ing Wes 61:240 - 248, 1995. 

Args:
    unit::GeothermalProbes:             the unit struct of the geothermal probe with all its parameters
    fluid_reynolds_number::Float64:     the current reynolds number in the pipe of the geothermal probe
Returns:
    Nusselt_turbulent::Float64:         Nusselt number in current pipe for given Reynolds number and pipe dimensions
"""
function calculate_Nu_turbulent(unit::GeothermalProbes, fluid_reynolds_number::Float64)
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
    energy_demanded = balance(exchanges) + energy_potential(exchanges)
    energy_available = unit.current_max_output_energy  # is positive

    # shortcut if there is no energy demanded
    if energy_demanded >= -sim_params["epsilon"]
        set_max_energy!(unit.output_interfaces[unit.m_heat_out], 0.0)    
        return
    end

    unit.output_temperatures_last_timestep = temp_min_all_non_empty(exchanges)

    for exchange in exchanges
        demanded_on_interface = exchange.balance + exchange.energy_potential
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
    energy_delivered = -(unit.current_max_output_energy - energy_available)
    unit.energy_in_out_per_probe_meter[unit.time_index] = energy_delivered  / (unit.probe_depth * unit.number_of_probes) # from total energy to specific power of one single probe.
    
end

function load(unit::GeothermalProbes, sim_params::Dict{String,Any})
    inface = unit.input_interfaces[unit.m_heat_in]
    exchanges = balance_on(inface, inface.source)
    energy_available = balance(exchanges) + energy_potential(exchanges)
    energy_demand = unit.max_input_energy  # is positive

    # shortcut if there is no energy to be used
    if ( energy_available <= sim_params["epsilon"] ||
         !unit.regeneration)
        set_max_energy!(unit.input_interfaces[unit.m_heat_in], 0.0)
        unit.fluid_temperature = calculate_new_fluid_temperature(unit::GeothermalProbes)
        return
    end

    for exchange in exchanges
        exchange_energy_available = exchange.balance + exchange.energy_potential
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
    unit.fluid_temperature = calculate_new_fluid_temperature(unit::GeothermalProbes)
    
end

function balance_on(
    interface::SystemInterface,
    unit::GeothermalProbes
)::Vector{EnergyExchange}

    caller_is_input = unit.uac == interface.target.uac
    purpose_uac = unit.uac == interface.target.uac ? interface.target.uac : interface.source.uac

    return [EnEx(
        balance=interface.balance,
        energy_potential=caller_is_input ? - unit.max_input_energy : unit.max_output_energy,
        purpose_uac = purpose_uac,
        temperature_min=interface.temperature_min,
        temperature_max=interface.temperature_max,
        pressure=nothing,
        voltage=nothing,
    )]
end

function output_values(unit::GeothermalProbes)::Vector{String}
    return [string(unit.m_heat_in)*" IN",
            string(unit.m_heat_out)*" OUT",
            "new_fluid_temperature",
            "current_output_temperature",
            "fluid_reynolds_number"
			]
end

function output_value(unit::GeothermalProbes, key::OutputKey)::Float64
    if key.value_key == "IN"
        return calculate_energy_flow(unit.input_interfaces[key.medium])
    elseif key.value_key == "OUT"
        return calculate_energy_flow(unit.output_interfaces[key.medium])
    elseif key.value_key =="new_fluid_temperature"
        return unit.fluid_temperature
    elseif key.value_key == "current_output_temperature"
        return unit.current_output_temperature
    elseif key.value_key =="fluid_reynolds_number"
        return unit.fluid_reynolds_number
    end
    throw(KeyError(key.value_key))
end


export GeothermalProbes