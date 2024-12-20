using GLMakie
using Roots
using Plots

"""
Implementation of geothermal heat collector.
This implementations acts as storage as is can produce and load energy.
"""

mutable struct GeothermalHeatCollector <: Component
    uac::String
    controller::Controller
    sys_function::SystemFunction
    input_interfaces::InterfaceMap
    output_interfaces::InterfaceMap
    m_heat_in::Symbol
    m_heat_out::Symbol

    ambient_temperature_profile::Union{Profile,Nothing}
    global_radiation_profile::Union{Profile,Nothing}
    infrared_sky_radiation_profile::Union{Profile,Nothing}
    unloading_temperature_spread::Temperature
    fluid_min_output_temperature::Temperature
    fluid_max_input_temperature::Temperature
    loading_temperature_spread::Temperature

    max_output_power::Union{Nothing,Float64}
    max_input_power::Union{Nothing,Float64}
    regeneration::Bool
    max_output_energy::Float64
    max_input_energy::Float64

    current_output_temperature::Temperature
    current_input_temperature::Temperature
    ambient_temperature::Temperature

    soil_specific_heat_capacity::Float64
    soil_specific_heat_capacity_frozen::Float64
    soil_density::Float64
    soil_heat_conductivity::Float64
    soil_heat_conductivity_frozen::Float64
    soil_density_vector::Vector
    soil_specific_enthalpy_of_fusion::Float64

    phase_change_upper_boundary_temperature::Temperature
    phase_change_lower_boundary_temperature::Temperature

    surface_convective_heat_transfer_coefficient::Float64
    surface_reflection_factor::Float64

    t1::Array{Float64}
    t2::Array{Float64}
    cp::Array{Float64}
    soil_weight::Array{Float64}
    fluid_node_y_idx::Int
    is_pipe_surrounding::Array{Bool}

    pipe_radius_outer::Float64
    pipe_thickness::Float64
    pipe_d_i::Float64
    pipe_d_o::Float64
    pipe_laying_depth::Float64
    pipe_length::Float64
    number_of_pipes::Float64
    pipe_spacing::Float64
    considered_soil_depth::Float64
    accuracy_mode::String
    model_type::String
    pipe_soil_thermal_resistance::Floathing

    global_radiation_power::Float64
    boltzmann_constant::Float64
    surface_emissivity::Float64

    dx::Vector{Float64}
    dy::Vector{Float64}
    dx_mesh::Vector{Float64}
    dy_mesh::Vector{Float64}
    dz::Float64
    dt::Integer

    fluid_temperature::Temperature

    collector_total_heat_energy_in_out::Float64

    pipe_heat_conductivity::Float64
    fluid_specific_heat_capacity::Float64
    fluid_prandtl_number::Float64
    fluid_density::Float64
    fluid_kinematic_viscosity::Float64
    fluid_heat_conductivity::Float64
    use_dynamic_fluid_properties::Bool
    nusselt_approach::String

    fluid_reynolds_number::Float64
    alpha_fluid_pipe::Float64
    average_temperature_adjacent_to_pipe::Float64
    undisturbed_ground_temperature::Float64

    temp_field_output::Array{Float64}
    n_internal_timesteps::Int
    sigma_lat::Float64
    t_lat::Float64
    delta_t_lat::Float64
    volume_adjacent_to_pipe::Float64

    function GeothermalHeatCollector(uac::String, config::Dict{String,Any}, sim_params::Dict{String,Any})
        m_heat_in = Symbol(default(config, "m_heat_in", "m_h_w_ht1"))
        m_heat_out = Symbol(default(config, "m_heat_out", "m_h_w_lt1"))
        register_media([m_heat_in, m_heat_out])

        ambient_temperature_profile = get_temperature_profile_from_config(config, sim_params, uac)
        global_radiation_profile = get_glob_solar_radiation_profile_from_config(config, sim_params, uac) # Wh/m^2
        infrared_sky_radiation_profile = get_infrared_sky_radiation_profile_from_config(config, sim_params, uac) # W/m^2

        # get model type from input file
        model_type = default(config, "model_type", "simplified")
        model_type_allowed_values = ["simplified", "detailed"]
        if !(model_type in model_type_allowed_values)
            @error "Undefined model type \"$(model_type)\" of unit \"$(uac)\". Has to be one of: $(model_type_allowed_values)."
            throw(InputError)
        end

        return new(uac,                                                  # uac
                   Controller(default(config, "control_parameters", nothing)),
                   sf_storage,                                           # sys_function
                   InterfaceMap(m_heat_in => nothing),                   # input_interfaces
                   InterfaceMap(m_heat_out => nothing),                  # output_interfaces
                   m_heat_in,                                            # medium name of input interface
                   m_heat_out,                                           # medium name of output interface
                   ambient_temperature_profile,                          # [°C] ambient temperature profile
                   global_radiation_profile,                             # [Wh/m^2]
                   infrared_sky_radiation_profile,                       # [Wh/m^2]
                   default(config, "unloading_temperature_spread", 3.0), # [K] temperature spread between forward and return flow during unloading            
                   default(config, "fluid_min_output_temperature", nothing),   # [°C] minimum output temperature of the fluid for unloading
                   default(config, "fluid_max_input_temperature", nothing),    # [°C] maximum input temperature of the fluid for loading
                   default(config, "loading_temperature_spread", 3.0),   # [K] temperature spread between forward and return flow during loading         
                   default(config, "max_output_power", 20),              # maximum output power in W/m^2, set by user. Depending on ground and climate localization. [VDI 4640-2.]
                   default(config, "max_input_power", 20),               # maximum input power in W/m^2, set by user. Depending on ground and climate localization. [VDI 4640-2.]
                   default(config, "regeneration", true),                # flag if regeneration should be taken into account
                   0.0,                                                  # max_output_energy [Wh] in every time step, calculated in control()
                   0.0,                                                  # max_input_energy [Wh] in every time step, calculated in control()
                   0.0,                                                  # current_output_temperature [°C] in current time step, calculated in control()
                   0.0,                                                  # current_input_temperature [°C] in current time step, calculated in control()
                   0.0,                                                  # ambient_temperature [°C] in current time step, calculated in control()
                   default(config, "soil_specific_heat_capacity", 1000), # specific heat capacity of unfrozen soil [J/(kgK)]
                   default(config, "soil_specific_heat_capacity_frozen", 900), # specific heat capacity of soil in fully frozen condition [J/(kgK)]
                   default(config, "soil_density", 2000),                # density of soil [kg/m^3]
                   default(config, "soil_heat_conductivity", 1.5),       # heat conductivity of unfrozen soil (lambda) [W/(mK)]
                   default(config, "soil_heat_conductivity_frozen", 2.0),# heat conductivity of frozen soil (lambda) [W/(mK)]
                   Array{Float64}(undef, 0),                             # soil_density_vector: vector to set soil density.
                   default(config, "soil_specific_enthalpy_of_fusion", 90000),            # specific enthalpy of fusion of soil [J/kg]
                   default(config, "phase_change_upper_boundary_temperature", -0.25),     # phase_change_upper_boundary_temperature [°C]
                   default(config, "phase_change_lower_boundary_temperature", -1),        # phase_change_lower_boundary_temperature [°C]
                   default(config, "surface_convective_heat_transfer_coefficient", 14.7), # convective heat transfer on surface [W/(m^2 K)]
                   default(config, "surface_reflection_factor", 0.25),   # reflection factor / albedo value of surface [-]
                   Array{Float64}(undef, 0, 0),                          # t1 [°C] holds the temperature of the last timestep
                   Array{Float64}(undef, 0, 0),                          # t2 [°C] holds the temperature of the current timestep
                   Array{Float64}(undef, 0, 0),                          # cp [J/(kg K)], holds the specific heat capacity for each node
                   Array{Float64}(undef, 0, 0),                          # soil_weight [kg], precalculated density * volume of the soil around each node
                   0,                                                    # fluid_node_y_idx, identifies the fluid node index in y direction
                   Array{Bool}(undef, 0, 0),                             # is_pipe_surrounding, identifies nodes surrounding the fluid-node: pipe-surrounding: true; else: false
                   default(config, "pipe_radius_outer", 0.016),          # pipe outer radius [m]
                   default(config, "pipe_thickness", 0.003),             # thickness of pipe [m]
                   0.0,                                                  # pipe_d_i: pipe inner diameter [m], calculated in initailize() 
                   0.0,                                                  # pipe_d_o: pipe outer diameter [m], calculated in initailize() 
                   default(config, "pipe_laying_depth", 1.5),            # deepth of pipe system below the ground sourface [m]
                   default(config, "pipe_length", 100),                  # pipe length of one collector pipe [m]
                   default(config, "number_of_pipes", 1),                # numbers of parallel pipes, each with a length of "pipe_length"
                   default(config, "pipe_spacing", 0.5),                 # distance between pipes of collector [m]
                   default(config, "considered_soil_depth", 10.0),       # depth of the soil considered in the simulation [m]
                   default(config, "accuracy_mode", "normal"),           # accuracy_mode. Has to be one of "very_rough", "rough", "normal", "high", "very_high"
                   model_type,                                           # model_type. currently "simplified" with constant fluid-to-soil resistance and 
                   #                                                       "detailed" with calculated fluid-to-soil resistance in every time step are available.
                   default(config, "pipe_soil_thermal_resistance", 0.1), # thermal resistance in [(m K)/W], only for model_type = simplified
                   0.0,                                                  # global_radiation_power [W/m^2]: solar global radiation on horizontal surface, to be read in from weather-profile
                   5.6697e-8,                                            # boltzmann_constant [W/(m^2 K^4)]: Stefan–Boltzmann-Constant
                   default(config, "surface_emissivity", 0.9),           # surface_emissivity [-]: emissivity on ground surface
                   Array{Float64}(undef, 0),                             # dx [m] horizontal dimension parallel to ground sourface and orthogonal to pipe
                   Array{Float64}(undef, 0),                             # dy [m] vertical dimension orthogonal to ground surface and orthgonal to pipe
                   Array{Float64}(undef, 0),                             # dx_mesh [m] this is dx between the nodes (while dx is the x-width assigned to each node)
                   Array{Float64}(undef, 0),                             # dy_mesh [m] this is dy between the nodes (while dy is the y-width assigned to each node)
                   0.0,                                                  # dz [m] horizontal dimension parallel to ground sourface and parallel to pipe (equals the length of one pipe, constant)
                   default(config, "internal_time_step", 0),             # dt [s], time step width of internal timestep, is calculated depending on ground properties and discretisation
                   0.0,                                                  # fluid_temperature [°C], is set to fluid_start_temperature at the beginning               
                   0.0,                                                  # collector_total_heat_energy_in_out, total heat flux in or out of collector
                   default(config, "pipe_heat_conductivity", 0.4),       # pipe_heat_conductivity [W/(mK)] 
                   default(config, "fluid_specific_heat_capacity", 3800),# fluid_specific_heat_capacity [J/(kg K)]
                   default(config, "fluid_prandtl_number", 30),          # prandtl number [-], preset for 30 % glycol at 0 °C 
                   default(config, "fluid_density", 1045),               # fluid density [kg/m^3], preset for 30 % glycol at 0 °C
                   default(config, "fluid_kinematic_viscosity", 3.9e-6), # fluid_kinematic_viscosity [m^2/s], preset fo 30 % glycol at 0 °C
                   default(config, "fluid_heat_conductivity", 0.5),      # fluid_heat_conductivity [W/(mK)], preset fo 30 % glycol at 0 °C
                   default(config, "use_dynamic_fluid_properties", false),# use_dynamic_fluid_properties [bool], false for constant, true for temperature-dependent fluid properties accoring to TRNSYS Type 710                   
                   default(config, "nusselt_approach", "Stephan"),       # Approach used for the caluclation of the Nußelt number, can be one of: Stephan, Ramming
                   0.0,                                                  # fluid_reynolds_number, to be calculated in function.
                   0.0,                                                  # alpha_fluid_pipe: convective heat transfer coefficient between fluid and pipe
                   default(config, "start_temperature_fluid_and_pipe", 15.5),  # average_temperature_adjacent_to_pipe [°C], used as starting temperature of fluid and soil near pipe during initialisation
                   default(config, "undisturbed_ground_temperature", 9.0),     # undisturbed ground temperature at the bottom of the simulation boundary [°C]
                   Array{Float64}(undef, 0, 0, 0),                       # temp_field_output [°C], holds temperature field of nodes for output plot
                   0,                                                    # n_internal_timesteps: number of internal time steps within one simulation time step
                   0.0,                                                  # sigma_lat; precalculated parameters for freezing function
                   0.0,                                                  # t_lat; precalculated parameters for freezing function
                   0.0,                                                  # delta_t_lat; precalculated parameters for freezing function
                   0.0)                                                  # volume_adjacent_to_pipe; precalculated parameter
    end
end

function initialise!(unit::GeothermalHeatCollector, sim_params::Dict{String,Any})
    if unit.regeneration
        set_storage_transfer!(unit.input_interfaces[unit.m_heat_in],
                              unload_storages(unit.controller, unit.m_heat_in))
    end
    set_storage_transfer!(unit.output_interfaces[unit.m_heat_out],
                          load_storages(unit.controller, unit.m_heat_out))

    # calculate diameters of pipe
    unit.pipe_d_i = 2 * unit.pipe_radius_outer - (2 * unit.pipe_thickness)
    unit.pipe_d_o = 2 * unit.pipe_radius_outer

    # calculate max_energy
    A_collector = unit.pipe_length * unit.pipe_spacing * unit.number_of_pipes
    unit.max_output_energy = watt_to_wh(unit.max_output_power * A_collector)
    if unit.regeneration
        unit.max_input_energy = watt_to_wh(unit.max_input_power * A_collector)
    else
        unit.max_input_energy = 0.0
    end

    # calculate coefficients for freezing function 
    unit.sigma_lat = (unit.phase_change_upper_boundary_temperature - unit.phase_change_lower_boundary_temperature) / 5
    unit.t_lat = (unit.phase_change_upper_boundary_temperature + unit.phase_change_lower_boundary_temperature) / 2
    unit.delta_t_lat = unit.phase_change_upper_boundary_temperature - unit.phase_change_lower_boundary_temperature

    # set fluid start temperature
    unit.fluid_temperature = unit.average_temperature_adjacent_to_pipe

    # calculate simulation mesh
    if unit.accuracy_mode == "very_rough"
        min_mesh_width = unit.pipe_d_o
        max_mesh_width = unit.pipe_d_o * 256
        expansion_factor = 2.0
    elseif unit.accuracy_mode == "rough"
        min_mesh_width = unit.pipe_d_o / 2
        max_mesh_width = unit.pipe_d_o * 128
        expansion_factor = 2.0
    elseif unit.accuracy_mode == "normal"
        min_mesh_width = unit.pipe_d_o / 4
        max_mesh_width = unit.pipe_d_o * 64
        expansion_factor = 2.0
    elseif unit.accuracy_mode == "high"
        min_mesh_width = unit.pipe_d_o / 8
        max_mesh_width = unit.pipe_d_o * 32
        expansion_factor = 2.0
    elseif unit.accuracy_mode == "very_high"
        min_mesh_width = unit.pipe_d_o / 16
        max_mesh_width = unit.pipe_d_o * 16
        expansion_factor = 2.0
    else
        @error "In geothermal collector $(unit.uac), the accuracy_mode has to be one of: very_rough, rough, " *
               "normal, high or very_high"
        throw(InputError)
    end

    # dy_mesh holds the delta between the nodes, while dy is the y-width assigned to each node
    unit.dy, y_pipe_node_num = create_mesh_y(min_mesh_width,
                                             max_mesh_width,
                                             expansion_factor,
                                             unit.pipe_laying_depth,
                                             unit.pipe_radius_outer,
                                             unit.considered_soil_depth)

    for y in 1:(length(unit.dy) - 1)
        push!(unit.dy_mesh, (unit.dy[y] + unit.dy[y + 1]) / 2)
    end

    # dx_mesh holds the delta between the nodes, while dx is the x-width assigned to each node 
    unit.dx = create_mesh_x(min_mesh_width,
                            max_mesh_width,
                            expansion_factor,
                            unit.pipe_radius_outer,
                            unit.pipe_spacing)

    for x in 1:(length(unit.dx) - 1)
        push!(unit.dx_mesh, (unit.dx[x] + unit.dx[x + 1]) / 2)
    end

    unit.dz = copy(unit.pipe_length)

    n_nodes_x = length(unit.dx)
    n_nodes_y = length(unit.dy)

    # localize fluid and adjacent nodes as they will be calculated separatly later
    unit.fluid_node_y_idx = y_pipe_node_num

    # set pipe-surrounding nodes assuming that the pipe has one central
    # and 8 surrounding nodes, of which 5 are considered in the axisymmetric grid used here.
    unit.is_pipe_surrounding = fill(false, n_nodes_y, n_nodes_x)  # pipe surrounding nodes are true, all others are false
    unit.is_pipe_surrounding[y_pipe_node_num - 1, 1] = true
    unit.is_pipe_surrounding[y_pipe_node_num - 1, 2] = true
    unit.is_pipe_surrounding[y_pipe_node_num + 0, 2] = true
    unit.is_pipe_surrounding[y_pipe_node_num + 1, 1] = true
    unit.is_pipe_surrounding[y_pipe_node_num + 1, 2] = true

    # set starting temperature distribution as follows:
    # top:         average ambient temperature, linearly interpolated to pipe surrounding
    # pipe:        6 nodes with the same starting temperature for pipe and pipe surrounding
    # lower bound: undistrubed ground temperature
    unit.t1 = fill(0.0, n_nodes_y, n_nodes_x)
    n_nodes_surface_to_pipe_surrounding = y_pipe_node_num - 3
    n_nodes_pipe_surrounding_to_ground = n_nodes_y - n_nodes_surface_to_pipe_surrounding - 5

    # set temperatures from surface to pipe surrounding with linear interpolation
    average_ambient_temperature = sum(values(unit.ambient_temperature_profile.data)) /
                                  length(values(unit.ambient_temperature_profile.data))
    step = (unit.average_temperature_adjacent_to_pipe - average_ambient_temperature) /
           sum(unit.dy_mesh[1:n_nodes_surface_to_pipe_surrounding])
    for i in 1:n_nodes_surface_to_pipe_surrounding
        unit.t1[i, :] .= average_ambient_temperature + sum(unit.dy_mesh[1:(i - 1)]) * step
    end

    # set temperatures for pipe and pipe surrounding
    unit.t1[(n_nodes_surface_to_pipe_surrounding + 1):(n_nodes_surface_to_pipe_surrounding + 5), :] .= unit.average_temperature_adjacent_to_pipe

    # set temperatures from pipe surrounding to ground with linear interpolation
    step = (unit.undisturbed_ground_temperature - unit.average_temperature_adjacent_to_pipe) /
           sum(unit.dy_mesh[(n_nodes_surface_to_pipe_surrounding + 5):(n_nodes_y - 1)])
    for i in 1:n_nodes_pipe_surrounding_to_ground
        dy_idx = i + n_nodes_surface_to_pipe_surrounding + 5
        unit.t1[dy_idx, :] = unit.t1[dy_idx - 1, :] .+ unit.dy_mesh[dy_idx - 1] * step
    end

    unit.t2 = copy(unit.t1)

    # # set starting temperature distribution (current settings for validation purpose)
    # # TODO: remove later
    # step = Int(round((y_pipe_node_num - 3) / 4))
    # unit.t1[(0 * step + 1):(1 * step), :] .= 9.5
    # unit.t1[(1 * step + 1):(2 * step), :] .= 11.5
    # unit.t1[(2 * step + 1):(3 * step), :] .= 13.5
    # unit.t1[(3 * step + 1):(y_pipe_node_num - 4), :] .= 15.5
    # unit.t1[(y_pipe_node_num - 3):(y_pipe_node_num + 3), :] .= 15.5
    # unit.t1[(y_pipe_node_num + 4):(end - 1), :] .= 15.5
    # unit.t1[end, :] .= 9

    # calculate volume around the pipe
    unit.volume_adjacent_to_pipe = ((unit.dx[1] + unit.dx[2]) *                                             # area adjacent to pipe in x-direction
                                    (unit.dy[unit.fluid_node_y_idx - 1] + unit.dy[unit.fluid_node_y_idx] +  # area adjacent to pipe in y-direction
                                     unit.dy[unit.fluid_node_y_idx + 1]) -
                                    (0.5 * pi * unit.pipe_radius_outer^2)) * unit.dz                        # 1/2 cross section of pipe

    # set soil density. currently only homogenous soil i  s considered, but more layers are possible
    unit.soil_density_vector = fill(unit.soil_density, n_nodes_y)

    # calculate soil weight of the soil around each node
    unit.soil_weight = fill(0.0, n_nodes_y, n_nodes_x)
    for h in 1:n_nodes_y
        for i in 1:n_nodes_x
            unit.soil_weight[h, i] = unit.soil_density_vector[h] * unit.dz * unit.dx[i] * unit.dy[h]
        end
    end

    # specific heat capacity for each node needed, because of apparent heat capacity method.
    unit.cp = fill(unit.soil_specific_heat_capacity, n_nodes_y, n_nodes_x)

    # internal time step according to TRNSYS Type 710 Model (Hirsch, Hüsing & Rockendorf 2017):
    if unit.dt == 0 # no dt given from input file
        unit.dt = Int(max(1,
                          floor(min(sim_params["time_step_seconds"],
                                    (min(unit.soil_specific_heat_capacity / unit.soil_heat_conductivity,
                                         unit.soil_specific_heat_capacity_frozen / unit.soil_heat_conductivity_frozen) *
                                     unit.soil_density * min(minimum(unit.dx_mesh), minimum(unit.dy_mesh))^2) / 4))))
    end
    # ensure that dt is a divisor of the simulation time step
    # find the largest divisor that is smaller that the calculated dt
    if sim_params["time_step_seconds"] % unit.dt != 0
        divisors = [i for i in 1:(unit.dt) if sim_params["time_step_seconds"] % i == 0]
        unit.dt = divisors[end]
    end
    @info "The geothermal collector $(unit.uac) will be simulated with an internal time step of $(unit.dt) s."

    unit.n_internal_timesteps = Int(floor(sim_params["time_step_seconds"] / unit.dt))

    # vector to hold the results of the temperatures for each node in each simulation time step
    unit.temp_field_output = zeros(Float64, sim_params["number_of_time_steps"], n_nodes_y, n_nodes_x)
end

""" 
    create_mesh_y(min_mesh_width, max_mesh_width, expansion_factor, pipe_laying_depth, pipe_radius_outer, total_depth_simulation_domain)

Creates a non-uniform mesh for geothermal collector in y-direction (orthogonal to ground surface and orthogonal to pipe).
| --------(ground surface)
| O       (pipe cross section)
V ........(lower bound)
y-direction

The mesh starts at the ground surface with an initial width of min_mesh_width. The width expands by a factor of 
expansion_factor for each node until reaching halfway to the pipe-laying depth. At that point, the mesh width begins
to decrease. The pipe itself is represented by one node, surrounded by one node each in the positive and negative
y-directions.
Below the pipe area, the mesh width starts again at min_mesh_width, increases up to halfway to the lower bound, and 
then decreases symmetrically.

The function determines the width of each volume element around a node point, starting with the surface layer.
The node itself is located in the centre of the respective width.

Note: The node spacing at the turning point between increasing and decreasing mesh width may deviate from the value
calculated using the expansion_factor.
"""
function create_mesh_y(min_mesh_width::Float64,
                       max_mesh_width::Float64,
                       expansion_factor::Float64,
                       pipe_laying_depth::Float64,
                       pipe_radius_outer::Float64,
                       total_depth_simulation_domain::Float64)
    function calculate_increasing_decreasing_distances(midpoint::Float64,
                                                       min_mesh_width::Float64,
                                                       max_mesh_width::Float64,
                                                       expansion_factor::Float64)
        distances = []
        current_width = min_mesh_width

        while distances == [] || sum(distances) + current_width < midpoint
            push!(distances, current_width)
            current_width = min(max_mesh_width, current_width * expansion_factor)
        end
        distance_to_midpoint = midpoint - sum(distances)
        if distance_to_midpoint > distances[end]
            push!(distances, distance_to_midpoint)
            push!(distances, distance_to_midpoint)
            mirrored_values = reverse(distances[1:(end - 2)])
        else
            push!(distances, 2 * distance_to_midpoint)
            mirrored_values = reverse(distances[1:(end - 1)])
        end
        append!(distances, mirrored_values)

        return distances
    end

    # define segments of the computing grid in y-direction 
    sy1 = 0                                         # surface
    sy2 = pipe_laying_depth - pipe_radius_outer * 3 / 2     # node above fluid node
    sy3 = pipe_laying_depth + pipe_radius_outer * 3 / 2     # node below fluid node
    sy4 = total_depth_simulation_domain             # lower simulation boundary

    # segment 1: surface to pipe 
    midpoint = (sy2 - sy1) / 2
    dy_1 = calculate_increasing_decreasing_distances(midpoint, min_mesh_width, max_mesh_width, expansion_factor)

    # detect node number of pipe 
    y_pipe_node_num = length(dy_1) + 2

    # segment 2: pipe
    dy_2 = [pipe_radius_outer, pipe_radius_outer, pipe_radius_outer]

    # segment 3: pipe to lower boundary 
    midpoint = (sy4 - sy3) / 2
    dy_3 = calculate_increasing_decreasing_distances(midpoint, min_mesh_width, max_mesh_width, expansion_factor)

    return [dy_1..., dy_2..., dy_3...], y_pipe_node_num
end

""" 
    create_mesh_x(min_mesh_width, max_mesh_width, expansion_factor, pipe_radius_outer, pipe_spacing)
  
Creates a non-uniform mesh for geothermal collector in x-direction (parallel to ground surface and orthogonal to pipe).
-----------------------  (ground surface)
O --> x-direction        (pipe cross section) 
.......................  (lower bound)

Starts with half the pipe diameter to represent the pipe nodes and increases the mesh width by the expansion_factor
for each node until half the pipe_spacing is reached.

The function determines the width of each volume element around a node point, starting with the vertical mirror axis
in which the pipe lies. The node itself is located in the centre of the respective width, except of the first node
which is located right on the mirror axis.

Note that the last dx can be smaller than calculated to meet the given boundary by pipe_spacing.
"""
function create_mesh_x(min_mesh_width::Float64,
                       max_mesh_width::Float64,
                       expansion_factor::Float64,
                       pipe_radius_outer::Float64,
                       pipe_spacing::Float64)
    # maximum width of grid in x direction
    x_bound = pipe_spacing / 2     # [m]

    # set first two dx to half the pipe radius
    dx = [pipe_radius_outer / 2, pipe_radius_outer]       # [m]

    # set next dx to the minumum width
    append!(dx, min_mesh_width)

    # add dx while x_bound is not exeeded
    while sum(dx) < x_bound
        append!(dx, min(dx[end] * expansion_factor, max_mesh_width))
    end

    # limit to x_bound
    if sum(dx) > x_bound
        dx[end] -= sum(dx) - x_bound
    end

    return dx
end

function control(unit::GeothermalHeatCollector,
                 components::Grouping,
                 sim_params::Dict{String,Any})
    update(unit.controller)

    # reset energy summarizer
    unit.collector_total_heat_energy_in_out = 0.0

    # get ambient temperature and global radiation from profile for current time step if needed
    unit.ambient_temperature = Profiles.value_at_time(unit.ambient_temperature_profile, sim_params)
    unit.global_radiation_power = Profiles.power_at_time(unit.global_radiation_profile, sim_params) # W/m^2

    unit.current_output_temperature = unit.fluid_temperature + unit.unloading_temperature_spread / 2
    unit.current_output_temperature = highest(unit.fluid_min_output_temperature, unit.current_output_temperature)
    set_temperature!(unit.output_interfaces[unit.m_heat_out],
                     nothing,
                     unit.current_output_temperature)

    set_max_energy!(unit.output_interfaces[unit.m_heat_out], unit.max_output_energy)

    # get input temperature for energy input (regeneration) and set temperature and max_energy to input interface
    if unit.regeneration
        unit.current_input_temperature = unit.fluid_temperature - unit.loading_temperature_spread / 2
        unit.current_input_temperature = lowest(unit.fluid_max_input_temperature, unit.current_input_temperature)
        set_temperature!(unit.input_interfaces[unit.m_heat_in],
                         unit.current_input_temperature,
                         nothing)

        set_max_energy!(unit.input_interfaces[unit.m_heat_in], unit.max_input_energy)
    end
end

function calculate_new_temperature_field!(unit::GeothermalHeatCollector, q_in_out::Float64, sim_params)
    function determine_soil_heat_conductivity(h1::Int, i1::Int, h2::Int, i2::Int)
        if unit.t1[h1, i1] >= unit.phase_change_upper_boundary_temperature
            soil_heat_conductivity_1 = unit.soil_heat_conductivity
        elseif unit.t1[h1, i1] <= unit.phase_change_lower_boundary_temperature
            soil_heat_conductivity_1 = unit.soil_heat_conductivity_frozen
        else # linear interpolation between frozen and unfrozen conductivity
            factor = (unit.phase_change_upper_boundary_temperature - unit.t1[h1, i1]) /
                     (unit.phase_change_upper_boundary_temperature - unit.phase_change_lower_boundary_temperature)
            soil_heat_conductivity_1 = unit.soil_heat_conductivity * (1 - factor) +
                                       unit.soil_heat_conductivity_frozen * factor
        end

        if unit.t1[h2, i2] >= unit.phase_change_upper_boundary_temperature
            soil_heat_conductivity_2 = unit.soil_heat_conductivity
        elseif unit.t1[h2, i2] <= unit.phase_change_lower_boundary_temperature
            soil_heat_conductivity_2 = unit.soil_heat_conductivity_frozen
        else # linear interpolation between frozen and unfrozen conductivity
            factor = (unit.phase_change_upper_boundary_temperature - unit.t1[h2, i2]) /
                     (unit.phase_change_upper_boundary_temperature - unit.phase_change_lower_boundary_temperature)
            soil_heat_conductivity_2 = unit.soil_heat_conductivity * (1 - factor) +
                                       unit.soil_heat_conductivity_frozen * factor
        end

        # arithmetic mean is used here for simplicity
        return (soil_heat_conductivity_1 + soil_heat_conductivity_2) / 2
    end

    if unit.model_type == "simplified"
        pipe_thermal_resistance_length_specific = unit.pipe_soil_thermal_resistance
    elseif unit.model_type == "detailed"
        # calculate heat transfer coefficient and thermal resistance
        unit.alpha_fluid_pipe, unit.fluid_reynolds_number = calculate_alpha_pipe(unit, q_in_out)

        # calculation of pipe_thermal_resistance_length_specific with the approach by TRNSYS Type 710 publication
        k = 1 / ((unit.pipe_d_o / (unit.alpha_fluid_pipe * unit.pipe_d_i) +
                  (log(unit.pipe_d_o / unit.pipe_d_i) * unit.pipe_d_o) / (2 * unit.pipe_heat_conductivity) +
                  unit.dx_mesh[2] /
                  (2 * determine_soil_heat_conductivity(unit.fluid_node_y_idx, 2, unit.fluid_node_y_idx, 2))))

        pipe_thermal_resistance_length_specific = 1 / (k * pi * unit.pipe_d_o)
    end

    unit.temp_field_output[Int(sim_params["time"] / sim_params["time_step_seconds"]) + 1, :, :] = copy(unit.t1)

    # calculate specific heat extraction for 1 pipe per length
    specific_heat_flux_pipe = wh_to_watts(q_in_out) / (unit.pipe_length * unit.number_of_pipes)   # [W/m]
    q_in_out_surrounding = specific_heat_flux_pipe * unit.pipe_length   # [W/pipe]

    # calculate effective mean sky temperature (sky radiative temperature) according to Stefan-Boltzmann law assuming 
    # the sky beeing a black body for the current time step:
    sky_temperature = (Profiles.power_at_time(unit.infrared_sky_radiation_profile, sim_params) /
                       unit.boltzmann_constant)^0.25  # [K]

    # loop over internal timesteps
    for _ in 1:(unit.n_internal_timesteps)

        # calculate fluid temperature and pipe-surrounding nodes
        for h in [unit.fluid_node_y_idx]  # vertical direction
            for i in [1]                  # horizontal direction
                unit.fluid_temperature = unit.average_temperature_adjacent_to_pipe +
                                         pipe_thermal_resistance_length_specific * specific_heat_flux_pipe # "+", because specific_heat_flux_pipe is negative during heat extraction.        

                unit.t2[h, i] = unit.fluid_temperature

                #! format: off
                # upper node
                unit.t2[h - 1, i] = unit.t1[h - 1, i] +
                                    (unit.dz * unit.dx[i] * determine_soil_heat_conductivity(h - 1, i, h - 2, i) * (unit.t1[h - 2, i] - unit.t1[h - 1, i]) / unit.dy_mesh[h - 2] +   # heat conduction in negative y-direction
                                    1 / 16 * q_in_out_surrounding) *                                                                                                       # heat source / sink
                                    unit.dt / (unit.soil_density_vector[h] * unit.cp[h - 1, i] * unit.volume_adjacent_to_pipe / 8)

                unit.t2[h - 1, i], unit.cp[h - 1, i] = freezing(unit, unit.t1[h - 1, i], unit.t2[h - 1, i], unit.cp[h - 1, i], sim_params)

                # lower node
                unit.t2[h + 1, i] = unit.t1[h + 1, i] +
                                    (unit.dz * unit.dx[i] * determine_soil_heat_conductivity(h + 1, i, h + 2, i) *  (unit.t1[h + 2, i] - unit.t1[h + 1, i]) / unit.dy_mesh[h + 1] +  # heat conduction in positive y-direction
                                    1 / 16 * q_in_out_surrounding) *                                                                                                       # heat source / sink
                                    unit.dt / (unit.soil_density_vector[h] * unit.cp[h + 1, i] * unit.volume_adjacent_to_pipe / 8)

                unit.t2[h + 1, i], unit.cp[h + 1, i] = freezing(unit, unit.t1[h + 1, i], unit.t2[h + 1, i], unit.cp[h + 1, i], sim_params)

                # right node
                unit.t2[h, i + 1] = unit.t1[h, i + 1] +
                                    (unit.dz * unit.dy[h] * determine_soil_heat_conductivity(h, i + 1, h, i + 2) * (unit.t1[h, i + 2] - unit.t1[h, i + 1]) / unit.dx_mesh[i + 1] +   # heat conduction in positive x-direction
                                    1 / 8 * q_in_out_surrounding) *                                                                                                        # heat source / sink
                                    unit.dt / (unit.soil_density_vector[h] * unit.cp[h, i + 1] * unit.volume_adjacent_to_pipe / 4)

                unit.t2[h, i + 1], unit.cp[h, i + 1] = freezing(unit, unit.t1[h, i + 1], unit.t2[h, i + 1], unit.cp[h, i + 1], sim_params)

                # upper right node
                unit.t2[h - 1, i + 1] = unit.t1[h - 1, i + 1] +
                                        (unit.dz * unit.dy[h - 1] * determine_soil_heat_conductivity(h - 1, i + 1, h - 1, i + 2) * (unit.t1[h - 1, i + 2] - unit.t1[h - 1, i + 1]) / unit.dx_mesh[i + 1] +  # heat conduction in positive x-direction
                                        unit.dz * unit.dx[i + 1] * determine_soil_heat_conductivity(h - 1, i + 1, h - 2, i + 1) * (unit.t1[h - 2, i + 1] - unit.t1[h - 1, i + 1]) / unit.dy_mesh[h - 2] +   # heat conduction in negative y-direction
                                        1 / 8 * q_in_out_surrounding) *                                                                                                                       # heat source / sink
                                        unit.dt / (unit.soil_density_vector[h] * unit.cp[h - 1, i + 1] * unit.volume_adjacent_to_pipe / 4)

                unit.t2[h - 1, i + 1], unit.cp[h - 1, i + 1] = freezing(unit, unit.t1[h - 1, i + 1], unit.t2[h - 1, i + 1], unit.cp[h - 1, i + 1], sim_params)

                # lower right node
                unit.t2[h + 1, i + 1] = unit.t1[h + 1, i + 1] +
                                        (unit.dz * unit.dx[i + 1] * determine_soil_heat_conductivity(h + 1, i + 1, h + 2, i + 1) * (unit.t1[h + 2, i + 1] - unit.t1[h + 1, i + 1]) / unit.dy_mesh[h + 1] +   # heat conduction in positive y-direction
                                        unit.dz * unit.dy[h + 1] * determine_soil_heat_conductivity(h + 1, i + 1, h + 1, i + 2) * (unit.t1[h + 1, i + 2] - unit.t1[h + 1, i + 1]) / unit.dx_mesh[i + 1] +    # heat conduction in positive x-direction
                                        1 / 8 * q_in_out_surrounding) *                                                                                                                        # heat source / sink
                                        unit.dt / (unit.soil_density_vector[h] * unit.cp[h + 1, i + 1] * unit.volume_adjacent_to_pipe / 4)

                unit.t2[h + 1, i + 1], unit.cp[h + 1, i + 1] = freezing(unit, unit.t1[h + 1, i + 1], unit.t2[h + 1, i + 1], unit.cp[h + 1, i + 1], sim_params)

                unit.average_temperature_adjacent_to_pipe = (unit.t2[h - 1, i]     +
                                                             unit.t2[h + 1, i]     +
                                                             2 * unit.t2[h, i + 1]     +
                                                             2 * unit.t2[h + 1, i + 1] +
                                                             2 * unit.t2[h - 1, i + 1]) / 8
                #! format: on 

                # average temperature of the nodes around the pipe
                unit.t2[h - 1, i] = unit.average_temperature_adjacent_to_pipe
                unit.t2[h + 1, i] = unit.average_temperature_adjacent_to_pipe
                unit.t2[h, i + 1] = unit.average_temperature_adjacent_to_pipe
                unit.t2[h + 1, i + 1] = unit.average_temperature_adjacent_to_pipe
                unit.t2[h - 1, i + 1] = unit.average_temperature_adjacent_to_pipe
            end
        end

        # loop over vertical direction
        for h in 1:(length(unit.dy_mesh))
            # loop over horizontal direction
            for i in 1:(length(unit.dx))
                # calculate temperature of adjacent nodes to fluid
                if (h == unit.fluid_node_y_idx && i == 1) || unit.is_pipe_surrounding[h, i]
                    # do nothing, because fluid temperature and surrounding nodes are already calculated before
                    continue
                elseif h == 1  # surface layer
                    #! format: off
                    # left boundary
                    if i == 1
                        unit.t2[h, i] = unit.t1[h, i] +
                                        (unit.dz * unit.dx[i] * unit.surface_convective_heat_transfer_coefficient * (unit.ambient_temperature - unit.t1[h, i]) +                            # Convection on surface
                                         unit.dz * unit.dx[i] * (1 - unit.surface_reflection_factor) * unit.global_radiation_power +                                                        # Radiation on surface
                                         unit.dz * unit.dx[i] * unit.surface_emissivity * unit.boltzmann_constant * (sky_temperature^4 - (unit.t1[h, i] + 273.15)^4) +  # Radiation out of surfac
                                         unit.dz * unit.dy[h] * determine_soil_heat_conductivity(h, i, h, i + 1) * (unit.t1[h, i + 1] - unit.t1[h, i]) / unit.dx_mesh[i] +                            # heat conduction in positve x-direction
                                         unit.dz * unit.dx[i] * determine_soil_heat_conductivity(h, i, h + 1, i) * (unit.t1[h + 1, i] - unit.t1[h, i]) / unit.dy_mesh[h]) *                           # heat conduction in positive y-direction
                                        unit.dt / (unit.cp[h, i] * unit.soil_weight[h, i])

                    # right boundary
                    elseif i == length(unit.dx)
                        unit.t2[h, i] = unit.t1[h, i] +
                                        (unit.dz * unit.dx[i] * unit.surface_convective_heat_transfer_coefficient * (unit.ambient_temperature - unit.t1[h, i]) +                            # Convection on surfac
                                         unit.dz * unit.dx[i] * (1 - unit.surface_reflection_factor) * unit.global_radiation_power +                                                        # Radiation on surface
                                         unit.dz * unit.dx[i] * unit.surface_emissivity * unit.boltzmann_constant * (sky_temperature^4 - (unit.t1[h, i] + 273.15)^4) +  # Radiation out of surfac
                                         unit.dz * unit.dy[h] * determine_soil_heat_conductivity(h, i, h, i - 1) * (unit.t1[h, i - 1] - unit.t1[h, i]) / unit.dx_mesh[i - 1] +                        # heat conduction in negative x-direction
                                         unit.dz * unit.dx[i] * determine_soil_heat_conductivity(h, i, h + 1, i) * (unit.t1[h + 1, i] - unit.t1[h, i]) / unit.dy_mesh[h]) *                           # heat conduction in positive y-direction
                                        unit.dt / (unit.cp[h, i] * unit.soil_weight[h, i])

                    else
                        unit.t2[h, i] = unit.t1[h, i] +
                                        (unit.dz * unit.dx[i] * unit.surface_convective_heat_transfer_coefficient *  (unit.ambient_temperature - unit.t1[h, i]) +                           # Convection on surface
                                         unit.dz * unit.dx[i] * (1 - unit.surface_reflection_factor) * unit.global_radiation_power +                                                        # Radiation on surface
                                         unit.dz * unit.dx[i] * unit.surface_emissivity * unit.boltzmann_constant * (sky_temperature^4 - (unit.t1[h, i] + 273.15)^4) +  # Radiation out of surface
                                         unit.dz * unit.dy[h] * determine_soil_heat_conductivity(h, i, h, i + 1) * (unit.t1[h, i + 1] - unit.t1[h, i]) / unit.dx_mesh[i] +                            # heat conduction in positve x-direction
                                         unit.dz * unit.dy[h] * determine_soil_heat_conductivity(h, i, h, i - 1) * (unit.t1[h, i - 1] - unit.t1[h, i]) / unit.dx_mesh[i - 1] +                        # heat conduction in negative x-direction
                                         unit.dz * unit.dx[i] * determine_soil_heat_conductivity(h, i, h + 1, i) * (unit.t1[h + 1, i] - unit.t1[h, i]) / unit.dy_mesh[h]) *                           # heat conduction in positive y-direction
                                        unit.dt / (unit.cp[h, i] * unit.soil_weight[h, i])

                    end
                else # other layers below surface
                    # left boundary (forward difference with Neumann condition to the left)
                    if i == 1
                        unit.t2[h, i] = unit.t1[h, i] +
                                        (unit.dz * unit.dy[h] * determine_soil_heat_conductivity(h, i, h, i + 1) * (unit.t1[h, i + 1] - unit.t1[h, i]) / unit.dx_mesh[i] +       # heat conduction in positve x-direction
                                         unit.dz * unit.dx[i] * determine_soil_heat_conductivity(h, i, h + 1, i) * (unit.t1[h + 1, i] - unit.t1[h, i]) / unit.dy_mesh[h] +       # heat conduction in positive y-direction
                                         unit.dz * unit.dx[i] * determine_soil_heat_conductivity(h, i, h - 1, i) * (unit.t1[h - 1, i] - unit.t1[h, i]) / unit.dy_mesh[h - 1]) *  # heat conduction in negative y-direction
                                        unit.dt / (unit.cp[h, i] * unit.soil_weight[h, i])

                    # right boundary (backward difference with Neumann condition to the right)
                    elseif i == length(unit.dx)
                        unit.t2[h, i] = unit.t1[h, i] +
                                        (unit.dz * unit.dy[h] * determine_soil_heat_conductivity(h, i, h, i - 1) * (unit.t1[h, i - 1] - unit.t1[h, i]) / unit.dx_mesh[i - 1] +   # heat conduction in negative x-direction
                                         unit.dz * unit.dx[i] * determine_soil_heat_conductivity(h, i, h + 1, i) * (unit.t1[h + 1, i] - unit.t1[h, i]) / unit.dy_mesh[h] +       # heat conduction in positive y-direction
                                         unit.dz * unit.dx[i] * determine_soil_heat_conductivity(h, i, h - 1, i) * (unit.t1[h - 1, i] - unit.t1[h, i]) / unit.dy_mesh[h - 1]) *  # heat conduction in negative y-direction
                                        unit.dt / (unit.cp[h, i] * unit.soil_weight[h, i])
                    
                    else
                        unit.t2[h, i] = unit.t1[h, i] +
                                        (unit.dz * unit.dy[h] * determine_soil_heat_conductivity(h, i, h, i + 1) * (unit.t1[h, i + 1] - unit.t1[h, i]) / unit.dx_mesh[i] +       # heat conduction in positve x-direction
                                         unit.dz * unit.dy[h] * determine_soil_heat_conductivity(h, i, h, i - 1) * (unit.t1[h, i - 1] - unit.t1[h, i]) / unit.dx_mesh[i - 1] +   # heat conduction in negative x-direction
                                         unit.dz * unit.dx[i] * determine_soil_heat_conductivity(h, i, h + 1, i) * (unit.t1[h + 1, i] - unit.t1[h, i]) / unit.dy_mesh[h] +       # heat conduction in positive y-direction
                                         unit.dz * unit.dx[i] * determine_soil_heat_conductivity(h, i, h - 1, i) * (unit.t1[h - 1, i] - unit.t1[h, i]) / unit.dy_mesh[h - 1]) *  # heat conduction in negative y-direction
                                        unit.dt / (unit.cp[h, i] * unit.soil_weight[h, i])

                    end
                    #! format: on
                end
                unit.t2[h, i], unit.cp[h, i] = freezing(unit, unit.t1[h, i], unit.t2[h, i], unit.cp[h, i], sim_params)
            end
        end
        unit.t1 = copy(unit.t2)
    end
end

function plot_optional_figures_begin(unit::GeothermalHeatCollector,
                                     output_path::String,
                                     output_formats::Vector{String},
                                     sim_params::Dict{String,Any})
    # plot mesh of geothermal collector
    x_positions = cumsum([0; unit.dx])
    y_positions = cumsum([0; unit.dy]) .* (-1)

    plt = plot(;
               title="Mesh of the collector at \"$(unit.accuracy_mode)\"\n" *
                     "accuracy (dx_min = $(round(minimum(unit.dx)*100;digits=1)) mm) ",
               xlabel="horizontal dimension [m]",
               ylabel="vertical dimension [m]",
               legend=false,
               linewidth=6,
               gridlinewidth=1,
               size=(900, 1200),
               titlefontsize=28,
               guidefontsize=24,
               tickfontsize=24,
               legendfontsize=24,
               margin=15Plots.mm,
               aspect_ratio=:equal)

    # Draw vertical lines
    for x in x_positions
        Plots.plot!(plt, [x, x], [y_positions[1], y_positions[end]]; color=:black, linewidth=1)
    end

    # Draw horizontal lines
    for y in y_positions
        Plots.plot!(plt, [x_positions[1], x_positions[end]], [y, y]; color=:black, linewidth=1)
    end

    fig_name = "collector_simulation_mesh_$(unit.uac)"
    for output_format in output_formats
        savefig(output_path * "/" * fig_name * "." * output_format)
    end

    return true
end

function plot_optional_figures_end(unit::GeothermalHeatCollector,
                                   sim_params::Dict{String,Any})
    # plot temperature field as 3D mesh with time-slider
    @info "Plotting time-shiftable temperature distribution of geothermal collector $(unit.uac). " *
          "Close figure to continue..."

    f = Figure()
    ax = Axis3(f[1, 1])

    ax.zlabel = "Temperature [°C]"
    ax.xlabel = "Vertical expansion (depth) [m]"
    ax.ylabel = "Horizontal expansion [m]"
    min_temp = minimum(unit.temp_field_output)
    min_temp = min_temp < 0.0 ? 1.1 * min_temp : 0.9 * min_temp
    max_temp = maximum(unit.temp_field_output)
    max_temp = max_temp < 0.0 ? 0.9 * max_temp : 1.1 * max_temp
    Makie.zlims!(ax, min_temp, max_temp)

    x_abs = [0; cumsum(unit.dx_mesh)]               # Absolute x coordinates
    y_abs = [unit.dy[1] / 2; cumsum(unit.dy_mesh)]  # Absolute y coordinates

    ## activate for equal axis ratio
    # xlims!(ax, 0, max(x_abs[end], y_abs[end]))
    # ylims!(ax, 0, max(x_abs[end], y_abs[end]))

    time = Observable(1)
    surfdata = @lift(unit.temp_field_output[$time, :, :])
    GLMakie.surface!(ax, y_abs, x_abs, surfdata)
    GLMakie.scatter!(ax, y_abs, x_abs, surfdata)
    slg = SliderGrid(f[2, 1], (; range=1:1:sim_params["number_of_time_steps"], label="Time"))

    on(slg.sliders[1].value) do v
        time[] = v
    end
    wait(display(f))
end

# function to handle freezing of the soil. Corrects the new temperature to include the enthalpy of fusion
# both for melting and freezing.
function freezing(unit::GeothermalHeatCollector,
                  t_old::Float64,
                  t_new::Float64,
                  cp::Float64,
                  sim_params::Dict{String,Any})
    function f(t_new_correct)
        function cp_freezing(t_new_correct)
            return (unit.soil_specific_enthalpy_of_fusion / (unit.sigma_lat * sqrt(2 * pi)) *
                    exp(-0.5 * (t_new_correct - unit.t_lat)^2 / unit.sigma_lat^2))
        end
        if t_new_correct > unit.phase_change_upper_boundary_temperature
            return t_old -
                   t_new_correct +
                   coeff / (cp_freezing(t_new_correct) + unit.soil_specific_heat_capacity)
        elseif t_new_correct < unit.phase_change_lower_boundary_temperature
            return t_old -
                   t_new_correct +
                   coeff / (cp_freezing(t_new_correct) + unit.soil_specific_heat_capacity_frozen)
        else
            return t_old -
                   t_new_correct +
                   coeff / (cp_freezing(t_new_correct) + unit.soil_specific_heat_capacity_frozen +
                            (unit.soil_specific_heat_capacity - unit.soil_specific_heat_capacity_frozen) *
                            (t_new_correct - (unit.t_lat - unit.delta_t_lat / 2)) / unit.delta_t_lat)
        end
    end

    if (t_old >= unit.phase_change_upper_boundary_temperature && t_new < unit.phase_change_upper_boundary_temperature ||
        t_old <= unit.phase_change_lower_boundary_temperature && t_new > unit.phase_change_lower_boundary_temperature ||
        t_old < unit.phase_change_upper_boundary_temperature && t_old > unit.phase_change_lower_boundary_temperature) &&
       abs(t_old - t_new) > sim_params["epsilon"]
        # solve t_new_corrected = t_old + coeff / cp_new(t_new_corrected) for t_new_corrected
        coeff = (t_new - t_old) * cp
        t_new_corrected = 0.0
        try
            t_new_corrected = find_zero(f, (-20, 20), Roots.A42())
        catch e
            @warn "Initial solving of freezing function of the geothermal collector $(unit.uac) with strict tolerance " *
                  "failed, retrying with rougher tolerance for the current time step..."
            try
                t_new_corrected = find_zero(f, unit.t_lat; tol=1e-3)
            catch e2
                @error "The solving algorithm for the phase change process in the geothermal collector $(unit.uac) was not successful. " *
                       "Manually decrease the \"internal_time_step\" which is currently $(unit.dt) s or change the accuracy_mode."
                throw(CalculationError)
            end
        end
        cp_corrected = coeff / (t_new_corrected - t_old)
        return t_new_corrected, cp_corrected
    elseif t_new >= unit.phase_change_upper_boundary_temperature
        return t_new, unit.soil_specific_heat_capacity
    else # only happens if t_new <= unit.phase_change_lower_boundary_temperature
        return t_new, unit.soil_specific_heat_capacity_frozen
    end
end

# function to calculate heat transfer coefficient alpha.
function calculate_alpha_pipe(unit::GeothermalHeatCollector, q_in_out::Float64)

    # calculate mass flow in pipe
    collector_power_in_out_per_pipe = wh_to_watts(abs(q_in_out)) / unit.number_of_pipes  # W/pipe
    temperature_spread = q_in_out > 0 ? unit.loading_temperature_spread : unit.unloading_temperature_spread
    collector_mass_flow_per_pipe = collector_power_in_out_per_pipe /
                                   (unit.fluid_specific_heat_capacity * temperature_spread)  # kg/s

    if unit.use_dynamic_fluid_properties
        # calculate reynolds-number based on dynamic viscosity using dynamic temperature-dependend fluid properties, 
        # adapted from TRNSYS Type 710, for 30 Vol-% ethylene glycol mix:
        fluid_dynamic_viscosity = 0.0000017158 * unit.fluid_temperature^2 -
                                  0.0001579079 * unit.fluid_temperature + 0.0048830621
        unit.fluid_heat_conductivity = 0.0010214286 * unit.fluid_temperature + 0.447
        unit.fluid_prandtl_number = fluid_dynamic_viscosity * unit.fluid_specific_heat_capacity /
                                    unit.fluid_heat_conductivity
        fluid_reynolds_number = (4 * collector_mass_flow_per_pipe) / (pi * unit.pipe_d_i * fluid_dynamic_viscosity)
    else
        # calculate reynolds-number, based on kinematic viscosity with constant fluid properties.
        fluid_reynolds_number = (4 * collector_mass_flow_per_pipe) /
                                (unit.fluid_density * unit.fluid_kinematic_viscosity * unit.pipe_d_i * pi)
    end

    if fluid_reynolds_number <= 2300  # laminar
        nusselt = calculate_Nu_laminar(unit, fluid_reynolds_number)
    elseif fluid_reynolds_number > 2300 && fluid_reynolds_number <= 1e4 # transitional
        # Gielinski 1995
        factor = (fluid_reynolds_number - 2300) / (1e4 - 2300)
        nusselt = (1 - factor) * calculate_Nu_laminar(unit, 2300.0) +
                  factor * calculate_Nu_turbulent(unit, 10_000.0)
    else  # turbulent
        nusselt = calculate_Nu_turbulent(unit, fluid_reynolds_number)
    end

    alpha = nusselt * unit.fluid_heat_conductivity / unit.pipe_d_i

    return alpha, fluid_reynolds_number
end

function calculate_Nu_laminar(unit::GeothermalHeatCollector, fluid_reynolds_number::Float64)
    if unit.nusselt_approach == "Ramming"
        # Approach used in Ramming 2007 from Elsner, Norbert; Fischer, Siegfried; Huhn, Jörg; „Grundlagen der
        # Technischen Thermodynamik“,  Band 2 Wärmeübertragung, Akademie Verlag, Berlin 1993. 
        k_a = 1.1 - 1 / (3.4 + 0.0667 * unit.fluid_prandtl_number)
        k_n = 0.35 + 1 / (7.825 + 2.6 * sqrt(unit.fluid_prandtl_number))

        # calculate Nu-Number
        nusselt_laminar = ((k_a / (1 - k_n) *
                            (unit.fluid_prandtl_number * unit.pipe_d_i * fluid_reynolds_number / unit.pipe_length)^k_n)^3 +
                           4.364^3)^(1 / 3)
    elseif unit.nusselt_approach == "Stephan"
        # Stephan
        pr_water = 13.44                # Pr Number Water 0 °C as reference
        nusselt_laminar = 3.66 +
                          (0.0677 *
                           (fluid_reynolds_number * unit.fluid_prandtl_number * unit.pipe_d_i / unit.pipe_length)^1.33) /
                          (1 +
                           0.1 * unit.fluid_prandtl_number *
                           (fluid_reynolds_number * unit.pipe_d_i / unit.pipe_length)^0.83) *
                          (unit.fluid_prandtl_number / pr_water)^0.11
    else
        @error "In geothermal collector $(unit.uac), the nusselt_approach has to be one of: Ramming, Stephan."
        throw(InputError)
    end
    return nusselt_laminar
end

function calculate_Nu_turbulent(unit::GeothermalHeatCollector, fluid_reynolds_number::Float64)
    # Approached used from Gnielinski in: V. Gnielinski: Ein neues Berechnungsverfahren für die Wärmeübertragung 
    # im Übergangsbereich zwischen laminarer und turbulenter Rohrströmung. Forsch im Ing Wes 61:240–248, 1995. 
    zeta = (1.8 * log(fluid_reynolds_number) - 1.5)^-2
    nusselt_turbulent = (zeta / 8 * fluid_reynolds_number * unit.fluid_prandtl_number) /
                        (1 + 12.7 * sqrt(zeta / 8) * (unit.fluid_prandtl_number^(2 / 3) - 1))
    return nusselt_turbulent
end

# process function that provides energy from the geothermal heat collector and calculates new temperatures 
# according to actual delivered or received energy
function process(unit::GeothermalHeatCollector, sim_params::Dict{String,Any})
    # get actual required energy from output interface
    outface = unit.output_interfaces[unit.m_heat_out]
    exchanges = balance_on(outface, outface.target)
    energy_demanded = balance(exchanges) + energy_potential(exchanges)
    energy_available = unit.max_output_energy  # is positive

    # shortcut if there is no energy demanded
    if energy_demanded >= -sim_params["epsilon"]
        set_max_energy!(unit.output_interfaces[unit.m_heat_out], 0.0)
        return
    end

    for exchange in exchanges
        demanded_on_interface = exchange.balance + exchange.energy_potential

        if demanded_on_interface >= -sim_params["epsilon"]
            continue
        end

        if (exchange.temperature_min !== nothing &&
            exchange.temperature_min > unit.current_output_temperature)
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
    unit.collector_total_heat_energy_in_out = energy_delivered
end

function load(unit::GeothermalHeatCollector, sim_params::Dict{String,Any})
    if !unit.regeneration
        # calculate new temperatures of field to account for possible ambient effects
        calculate_new_temperature_field!(unit::GeothermalHeatCollector, unit.collector_total_heat_energy_in_out,
                                         sim_params)
        return
    end

    inface = unit.input_interfaces[unit.m_heat_in]
    exchanges = balance_on(inface, inface.source)
    energy_available = balance(exchanges) + energy_potential(exchanges)
    energy_demand = unit.max_input_energy  # is positive

    if energy_available <= sim_params["epsilon"]
        # shortcut if there is no energy to be used
        set_max_energy!(unit.input_interfaces[unit.m_heat_in], 0.0)
        # calculate new temperatures of field to account for possible ambient effects
        calculate_new_temperature_field!(unit::GeothermalHeatCollector, unit.collector_total_heat_energy_in_out,
                                         sim_params)
        return
    end

    for exchange in exchanges
        exchange_energy_available = exchange.balance + exchange.energy_potential

        if exchange_energy_available < sim_params["epsilon"]
            continue
        end

        if (exchange.temperature_min !== nothing &&
            exchange.temperature_min > unit.current_input_temperature
            ||
            exchange.temperature_max !== nothing &&
            exchange.temperature_max < unit.current_input_temperature)
            # we can only take energy if it's at a higher/equal temperature than the
            # collector's current input temperature
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

    energy_taken = unit.max_input_energy - energy_demand
    unit.collector_total_heat_energy_in_out += energy_taken
    calculate_new_temperature_field!(unit::GeothermalHeatCollector, unit.collector_total_heat_energy_in_out, sim_params)
end

function balance_on(interface::SystemInterface, unit::GeothermalHeatCollector)::Vector{EnergyExchange}
    caller_is_input = unit.uac == interface.target.uac
    balance_written = interface.max_energy.max_energy[1] === nothing || interface.sum_abs_change > 0.0
    purpose_uac = unit.uac == interface.target.uac ? interface.target.uac : interface.source.uac

    return [EnEx(; balance=interface.balance,
                 energy_potential=balance_written ? 0.0 :
                                  (caller_is_input ? -unit.max_input_energy : unit.max_output_energy),
                 purpose_uac=purpose_uac,
                 temperature_min=interface.temperature_min,
                 temperature_max=interface.temperature_max,
                 pressure=nothing,
                 voltage=nothing)]
end

function output_values(unit::GeothermalHeatCollector)::Vector{String}
    output_vals = []
    if unit.regeneration
        push!(output_vals, string(unit.m_heat_in) * " IN")
    end
    if unit.model_type == "detailed"
        push!(output_vals, "fluid_reynolds_number")
        push!(output_vals, "alpha_fluid_pipe")
    end
    append!(output_vals,
            [string(unit.m_heat_out) * " OUT",
             "fluid_temperature",
             "ambient_temperature",
             "global_radiation_power"])
    # push!(output_vals, "TEMPERATURE_xNodeNum_yNodeNum")

    return output_vals
end

function output_value(unit::GeothermalHeatCollector, key::OutputKey)::Float64
    if key.value_key == "IN"
        return calculate_energy_flow(unit.input_interfaces[key.medium])
    elseif key.value_key == "OUT"
        return calculate_energy_flow(unit.output_interfaces[key.medium])
    elseif startswith(key.value_key, "Temperature_")
        splitted = split(key.value_key, "_")
        x_idx = parse(Int, splitted[2])
        y_idx = parse(Int, splitted[3])
        if !(1 <= y_idx <= length(unit.dy) + 1) || !(1 <= x_idx <= length(unit.dx) + 1)
            throw(ArgumentError("The indexes ($(x_idx), $(y_idx)) of the requested temperature-output of the " *
                                "geothermal collector $(unit.uac) exeed the number of nodes of the mesh. The maximum " *
                                "is ($(length(unit.dx)), $(length(unit.dy)))."))
        else
            return unit.t2[x_idx, y_idx]
        end
    elseif key.value_key == "fluid_temperature"
        return unit.fluid_temperature
    elseif key.value_key == "fluid_reynolds_number"
        return unit.fluid_reynolds_number
    elseif key.value_key == "ambient_temperature"
        return unit.ambient_temperature
    elseif key.value_key == "global_radiation_power"
        return unit.global_radiation_power
    elseif key.value_key == "alpha_fluid_pipe"
        return unit.alpha_fluid_pipe
    end
    throw(KeyError(key.value_key))
end

export GeothermalHeatCollector
