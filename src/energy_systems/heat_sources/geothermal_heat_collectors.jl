using GLMakie

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
    unloading_temperature_spread::Temperature
    loading_temperature::Temperature
    loading_temperature_spread::Temperature

    max_output_power::Union{Nothing,Float64}
    max_input_power::Union{Nothing,Float64}
    regeneration::Bool
    max_output_energy::Float64
    max_input_energy::Float64

    current_output_temperature::Temperature
    current_input_temperature::Temperature
    ambient_temperature::Temperature

    soil_starting_temperature::Temperature
    soil_specific_heat_capacity::Float64
    soil_density::Float64
    soil_heat_conductivity::Float64
    soil_density_vector::Vector
    soil_heat_conductivity_vector::Vector
    soil_specific_enthalpy_of_fusion::Float64

    phase_change_upper_boundary_temperature::Temperature
    phase_change_lower_boundary_temperature::Temperature

    surface_convective_heat_transfer_coefficient::Float64
    surface_reflection_factor::Float64

    t1::Array{Float64}
    t2::Array{Float64}
    phase_change_state::Array{Int}
    cp1::Array{Float64}
    cp2::Array{Float64}
    is_fluid_node::Array{Bool}
    is_pipe_surrounding::Array{Bool}

    pipe_radius_outer::Float64
    pipe_thickness::Float64
    pipe_d_i::Float64
    pipe_d_o::Float64
    pipe_laying_depth::Float64
    pipe_length::Float64
    number_of_pipes::Float64
    pipe_spacing::Float64

    global_radiation_power::Float64
    boltzmann_constant::Float64
    surface_emissivity::Float64

    dx::Vector{Float64}
    dy::Vector{Float64}
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

    fluid_reynolds_number::Float64
    average_temperature_adjacent_to_pipe::Float64

    temp_field_output::Array{Float64}
    n_internal_timesteps::Int

    function GeothermalHeatCollector(uac::String, config::Dict{String,Any}, sim_params::Dict{String,Any})
        m_heat_in = Symbol(default(config, "m_heat_in", "m_h_w_ht1"))
        m_heat_out = Symbol(default(config, "m_heat_out", "m_h_w_lt1"))
        register_media([m_heat_in, m_heat_out])

        ambient_temperature_profile = get_temperature_profile_from_config(config, sim_params, uac)
        global_radiation_profile = get_glob_solar_radiation_profile_from_config(config, sim_params, uac) # Wh/m^2

        return new(uac,                    # uac
                   Controller(default(config, "control_parameters", nothing)),
                   sf_storage,                           # sys_function
                   InterfaceMap(m_heat_in => nothing),   # input_interfaces
                   InterfaceMap(m_heat_out => nothing),  # output_interfaces
                   m_heat_in,                            # medium name of input interface
                   m_heat_out,                           # medium name of output interface
                   ambient_temperature_profile,          # [°C] ambient temperature profile
                   global_radiation_profile,             # [Wh/m^2]
                   default(config, "unloading_temperature_spread", 3),   # temperature spread between forward and return flow during unloading            
                   default(config, "loading_temperature", nothing),      # nominal high temperature for loading geothermal heat collector storage, can also be set from other end of interface TODO
                   default(config, "loading_temperature_spread", 3),     # temperature spread between forward and return flow during loading         
                   default(config, "max_output_power", 20),              # maximum output power in W/m^2, set by user. Depending on ground and climate localization. [VDI 4640-2.]
                   default(config, "max_input_power", 20),               # maximum input power in W/m^2, set by user. Depending on ground and climate localization. [VDI 4640-2.]
                   default(config, "regeneration", true),                # flag if regeneration should be taken into account
                   0.0,                        # max_output_energy [Wh] in every time step, calculated in control()
                   0.0,                        # max_input_energy [Wh] in every time step, calculated in control()
                   0.0,                        # current_output_temperature [°C] in current time step, calculated in control()
                   0.0,                        # current_input_temperature [°C] in current time step, calculated in control()
                   0.0,                        # ambient_temperature [°C] in current time step, calculated in control()
                   default(config, "soil_starting_temperature", 15.5),   # starting temperature of soil near pipe [°C]
                   default(config, "soil_specific_heat_capacity", 1000), # specific heat capacity of soil [J/(kgK)]
                   default(config, "soil_density", 2000),                # density of soil [kg/m^3]              
                   default(config, "soil_heat_conductivity", 1.5),       # heat conductivity of soil (lambda) [W/(mK)]
                   Array{Float64}(undef, 0),                  # soil_density_vector: vector to set soil density.
                   Array{Float64}(undef, 0),                  # soil_heat_conductivity_vector: vector to hold the heat conductivity of soil.
                   default(config, "soil_specific_enthalpy_of_fusion", 90000),            # specific enthalpy of fusion of soil [J/kg]
                   default(config, "phase_change_upper_boundary_temperature", -0.25),     # phase_change_upper_boundary_temperature [°C]
                   default(config, "phase_change_lower_boundary_temperature", -1),        # phase_change_lower_boundary_temperature [°C]
                   default(config, "surface_convective_heat_transfer_coefficient", 14.7), # convective heat transfer on surface [W/(m^2 K)]
                   default(config, "surface_reflection_factor", 0.25),                    # reflection factor / Albedo value of surface [-]
                   Array{Float64}(undef, 0, 0),           # t1 [°C]
                   Array{Float64}(undef, 0, 0),           # t2 [°C]
                   Array{Int}(undef, 0, 0),               # phase_change_state [-]
                   Array{Float64}(undef, 0, 0),           # cp1 [J/(kg K)], holds the specific heat capacity of previous timestep (for apparent heat capacity method)
                   Array{Float64}(undef, 0, 0),           # cp2 [J/(kg K)], holds the specific heat capacity of the current timestep (for apparent heat capacity method)
                   Array{Bool}(undef, 0, 0),              # is_fluid_node, identifies the fluid node: fuid node: true; soil-node: false
                   Array{Bool}(undef, 0, 0),              # is_pipe_surrounding, identifies nodes surrounding the fluid-node: pipe-surrounding: true; else: false
                   default(config, "pipe_radius_outer", 0.016),          # pipe outer radius [m]
                   default(config, "pipe_thickness", 0.003),             # thickness of pipe [m]
                   0.0,                                                  # pipe_d_i: pipe inner diameter [m], calculated in initailize() 
                   0.0,                                                  # pipe_d_o: pipe outer diameter [m], calculated in initailize() 
                   default(config, "pipe_laying_depth", 1.5),            # deepth of pipe system below the ground sourface [m]
                   default(config, "pipe_length", 100),                  # pipe length of one collector pipe [m]
                   default(config, "number_of_pipes", 1),                # numbers of parallel pipes, each with a length of "pipe_length"
                   default(config, "pipe_spacing", 0.5),                 # distance between pipes of collector [m]
                   0.0,                        # global_radiation_power [W/m^2]: solar global radiation on horizontal surface, to be read in from weather-profile
                   5.67e-8,                    # boltzmann_constant [W/(m^2 K^4)]: Stefan–Boltzmann-Constant
                   0.9,                        # surface_emissivity [-]: emissivity on ground surface
                   Array{Float64}(undef, 0),   # dx [m] horizontal dimension parallel to ground sourface and orthogonal to pipe
                   Array{Float64}(undef, 0),   # dy [m] vertical dimension orthogonal to ground surface and orthgonal to pipe
                   0.0,                        # dz [m] horizontal dimension parallel to ground sourface and parallel to pipe (equals the length of one pipe, constant)
                   0.0,                        # dt [s], time step width of internal timestep, is calculated depending on ground properties and discretisation
                   default(config, "fluid_start_temperature", 11.0),     # fluid_temperature [°C], is set to fluid_start_temperature at the beginning               
                   0.0,                        # collector_total_heat_energy_in_out, total heat flux in or out of collector
                   default(config, "pipe_heat_conductivity", 0.4),                 # pipe_heat_conductivity [W/(mK)] 
                   default(config, "fluid_specific_heat_capacity", 3800),          # fluid_specific_heat_capacity [J/(kg K)]
                   default(config, "fluid_prandtl_number", 30),                    # prandtl number [-], preset for 30 % glycol at 0 °C 
                   default(config, "fluid_density", 1045),                         # fluid density [kg/m^3], preset for 30 % glycol at 0 °C
                   default(config, "fluid_kinematic_viscosity", 3.9e-6),           # fluid_kinematic_viscosity [m^2/s], preset fo 30 % glycol at 0 °C
                   default(config, "fluid_heat_conductivity", 0.5),                # fluid_heat_conductivity [W/(mK)], preset fo 30 % glycol at 0 °C
                   0.0,                                                            # fluid_reynolds_number, to be calculated in function.
                   default(config, "soil_around_pipe_starting_temperature", 16.0), # average_temperature_adjacent_to_pipe [°C] TODO: Doubling with "soil_starting_temperature"?
                   Array{Float64}(undef, 0, 0, 0),                                 # temp_field_output [°C], holds temperature field of nodes for output plot
                   0)                                                              # n_internal_timesteps: number of internal time steps within one simulation time step
    end
end

function initialise!(unit::GeothermalHeatCollector, sim_params::Dict{String,Any})
    if unit.regeneration
        set_storage_transfer!(unit.input_interfaces[unit.m_heat_in],
                              unload_storages(unit.controller, unit.m_heat_in))
    end
    set_storage_transfer!(unit.output_interfaces[unit.m_heat_out],
                          load_storages(unit.controller, unit.m_heat_out))

    # Discretization and starting Temperature-Field. --> ADD PREPROCESSING! TODO
    # calculate diameters of pipe
    unit.pipe_d_i = 2 * unit.pipe_radius_outer - (2 * unit.pipe_thickness)
    unit.pipe_d_o = 2 * unit.pipe_radius_outer

    # calculate max_energy
    A_collector = unit.pipe_length * unit.pipe_spacing * unit.number_of_pipes
    unit.max_output_energy = watt_to_wh(unit.max_output_power * A_collector)
    unit.max_input_energy = watt_to_wh(unit.max_input_power * A_collector)

    # Discretization. Will be done in Pre-Processing. TODO.
    dx_R = [unit.pipe_radius_outer]
    dx_RM = [unit.pipe_radius_outer, unit.pipe_radius_outer, 2 * unit.pipe_radius_outer, 4 * unit.pipe_radius_outer]
    dx_End = unit.pipe_spacing / 2 - sum(dx_R) - sum(dx_RM)
    unit.dx = [dx_R..., dx_RM..., dx_End...] # append()

    dy_O = [unit.pipe_radius_outer / 2, unit.pipe_radius_outer, 4 * unit.pipe_radius_outer,
            8 * unit.pipe_radius_outer, 16 * unit.pipe_radius_outer]                                # dy near surface
    dy_MR = [16 * unit.pipe_radius_outer, 8 * unit.pipe_radius_outer, 4 * unit.pipe_radius_outer,
             2 * unit.pipe_radius_outer, unit.pipe_radius_outer / 2]                                # dy getting smaller untill pipe node
    dy_M = [(unit.pipe_laying_depth - (sum(dy_O) + sum(dy_MR)) - unit.pipe_radius_outer)]           # dy is wide in the middle    
    dy_R = [unit.pipe_radius_outer, unit.pipe_radius_outer]                                         # distance between nodes adjacent to fluid
    dy_RU = [unit.pipe_radius_outer / 2,
             1 * unit.pipe_radius_outer,
             1 * unit.pipe_radius_outer,
             3 / 2 * unit.pipe_radius_outer,
             2 * unit.pipe_radius_outer,
             2 * unit.pipe_radius_outer,
             4 * unit.pipe_radius_outer,
             16 * unit.pipe_radius_outer,
             32 * unit.pipe_radius_outer,
             64 * unit.pipe_radius_outer,
             64 * unit.pipe_radius_outer,
             64 * unit.pipe_radius_outer,
             128 * unit.pipe_radius_outer,
             64 * unit.pipe_radius_outer,
             32 * unit.pipe_radius_outer,
             16 * unit.pipe_radius_outer,
             8 * unit.pipe_radius_outer]  # mesh below the fluid-node until lower simulation boundary   
    unit.dy = [dy_O..., dy_M..., dy_MR..., dy_R..., dy_RU...]   #size:27     # build dy-vector

    unit.dz = copy(unit.pipe_length)

    n_nodes_x = length(unit.dx) + 1
    n_nodes_y = length(unit.dy) + 1
    n_nodes_z = length(unit.dz) + 1

    # Localize fuid and adjacent Nodes. They will be calculated separatly.
    unit.is_fluid_node = fill(false, n_nodes_y, n_nodes_x)
    unit.is_fluid_node[(length(dy_O) + length(dy_MR) + length(dy_M) + length(dy_R) + 1 - 1), 1] = true

    unit.is_pipe_surrounding = fill(false, n_nodes_y, n_nodes_x)
    unit.is_pipe_surrounding[(length(dy_O) + length(dy_MR) + length(dy_M) + length(dy_R) + 1), 1] = true
    unit.is_pipe_surrounding[(length(dy_O) + length(dy_MR) + length(dy_M) + length(dy_R) + 1 - 2), 1] = true
    unit.is_pipe_surrounding[(length(dy_O) + length(dy_MR) + length(dy_M) + length(dy_R) + 1 - 1), 2] = true

    unit.is_pipe_surrounding[(length(dy_O) + length(dy_MR) + length(dy_M) + length(dy_R) + 1 - 2), 2] = true
    unit.is_pipe_surrounding[(length(dy_O) + length(dy_MR) + length(dy_M) + length(dy_R) + 1), 2] = true

    # set starting temperature distribution (current settings for validation purpose)
    unit.t1 = fill(0.0, n_nodes_y, n_nodes_x)
    unit.t1[1:3, :] .= 9.5
    unit.t1[4:6, :] .= 11.5
    unit.t1[7:9, :] .= 13.5
    unit.t1[10:13, :] .= 15.5
    unit.t1[14:18, :] .= unit.soil_starting_temperature
    unit.t1[19:25, :] .= 15.5
    unit.t1[26:30, :] .= 15.5
    unit.t1[31, :] .= 9    # TODO: Use avg ambient temperature out of weather data set.
    unit.t2 = copy(unit.t1)

    # set soil heat conductivity. currently only homogenous soil is considered, but more layers are possible
    unit.soil_heat_conductivity_vector = fill(unit.soil_heat_conductivity, n_nodes_y)

    # set soil density. currently only homogenous soil is considered, but more layers are possible
    unit.soil_density_vector = fill(unit.soil_density, n_nodes_y)

    # specific heat capacity for each node needed, because of apparent heat capacity method.
    unit.cp1 = fill(unit.soil_specific_heat_capacity, n_nodes_y, n_nodes_x)
    unit.cp2 = copy(unit.cp1)

    # create phase_change_state vector
    unit.phase_change_state = fill(0, n_nodes_y, n_nodes_x)

    # internal time step according to TRNSYS Type 710 Model (Hirsch, Hüsing & Rockendorf 2017):
    unit.dt = Int(floor(min(sim_params["time_step_seconds"],
                            (unit.soil_density * unit.soil_specific_heat_capacity *
                             min(minimum(unit.dx), minimum(unit.dy))^2) / (4 * unit.soil_heat_conductivity))))
    # ensure that dt is a divisor of the simulation time step
    if sim_params["time_step_seconds"] % unit.dt != 0
        unit.dt = Int(sim_params["time_step_seconds"] / Int(ceil(sim_params["time_step_seconds"] / unit.dt)))
    end
    unit.n_internal_timesteps = Int(floor(sim_params["time_step_seconds"] / unit.dt))

    # vector to hold the results of the temperatures for each node in each simulation time step
    unit.temp_field_output = zeros(Float64, sim_params["number_of_time_steps"], n_nodes_y, n_nodes_x)
end

function control(unit::GeothermalHeatCollector,
                 components::Grouping,
                 sim_params::Dict{String,Any})
    update(unit.controller)

    # reset energy summarizer
    unit.collector_total_heat_energy_in_out = 0.0

    # get ambient temperature and global radiation from profile for current time step if needed
    unit.ambient_temperature = Profiles.value_at_time(unit.ambient_temperature_profile, sim_params)
    unit.global_radiation_power = wh_to_watts(Profiles.value_at_time(unit.global_radiation_profile, sim_params)) # from Wh/m^2 to W/m^2

    unit.current_output_temperature = unit.fluid_temperature + unit.unloading_temperature_spread / 2
    set_temperature!(unit.output_interfaces[unit.m_heat_out],
                     nothing,
                     unit.current_output_temperature)

    set_max_energy!(unit.output_interfaces[unit.m_heat_out], unit.max_output_energy)

    # get input temperature for energy input (regeneration) and set temperature and max_energy to input interface
    if unit.regeneration
        unit.current_input_temperature = unit.average_temperature_adjacent_to_pipe - unit.loading_temperature_spread / 2 # of geothermal heat collector 
        set_temperature!(unit.input_interfaces[unit.m_heat_in],
                         unit.current_input_temperature,
                         nothing)

        set_max_energy!(unit.input_interfaces[unit.m_heat_in], unit.max_input_energy)
    end
end

function calculate_new_temperature_field!(unit::GeothermalHeatCollector, q_in_out::Float64, sim_params)

    # calculate heat transfer coefficient and thermal resistance
    alpha_fluid, unit.fluid_reynolds_number = calculate_alpha_pipe(unit, q_in_out)

    # calculation of pipe_thermal_resistance_length_specific with the approach by TRNSYS Type 710 publication
    k = (pi / 4) / ((unit.pipe_d_o / (alpha_fluid * unit.pipe_d_i) +
                     (log(unit.pipe_d_o / unit.pipe_d_i) * unit.pipe_d_o) / (2 * unit.pipe_heat_conductivity) +
                     minimum(unit.dx) / (2 * unit.soil_heat_conductivity)))

    pipe_thermal_resistance_length_specific = 1 / (k * pi * unit.pipe_d_o)

    # calculate fluid temperature
    for h in 1:length(unit.dy)
        for i in 1:2
            if unit.is_fluid_node[h, i]
                # calculate specific heat extraction for 1 pipe per length
                specific_heat_flux_pipe = (wh_to_watts(q_in_out) / (unit.pipe_length * unit.number_of_pipes))
                unit.fluid_temperature = unit.average_temperature_adjacent_to_pipe +
                                         pipe_thermal_resistance_length_specific * specific_heat_flux_pipe # "+", because specific_heat_flux_pipe is negative during heat extraction.        

                unit.t2[h, i] = unit.fluid_temperature

                volume_adjacent_to_pipe = ((unit.dx[i] + unit.dx[i + 1] / 2) *                                        # total x-Direction
                                           (unit.dy[h - 2] / 2 + unit.dy[h - 1] + unit.dy[h] + unit.dy[h + 1] / 2) -  # total y-Direction
                                           (pi * unit.dy[h]^2 / 2)) * unit.dz                                         # 1/2 Area pipe.          

                q_in_out_surrounding = wh_to_watts(q_in_out) / (unit.number_of_pipes)

                # upper node               
                unit.t2[h - 1, i] = unit.t1[h - 1, i] +
                                    (unit.dz * unit.dx[i] / 2 * unit.soil_heat_conductivity_vector[h] *  # heat conduction in negative y-direction
                                     (unit.t1[h - 2, i] - unit.t1[h - 1, i]) / (unit.dy[h - 2]) +
                                     1 / 16 * q_in_out_surrounding) *                                    # heat source / sink
                                    unit.dt * 1 /
                                    (unit.soil_density_vector[h] * unit.cp1[h, i] * volume_adjacent_to_pipe / 5)

                # lower node
                unit.t2[h + 1, i] = unit.t1[h + 1, i] +
                                    (unit.dz * unit.dx[i] / 2 * unit.soil_heat_conductivity_vector[h] *  # heat conduction in positive y-direction
                                     (unit.t1[h + 2, i] - unit.t1[h + 1, i]) / (unit.dy[h + 1]) +
                                     1 / 16 * q_in_out_surrounding) *                                    # heat source / sink
                                    unit.dt * 1 /
                                    (unit.soil_density_vector[h] * unit.cp1[h, i] * volume_adjacent_to_pipe / 5)

                # right node
                unit.t2[h, i + 1] = unit.t1[h, i + 1] +
                                    (unit.dz * 0.5 * (unit.dy[h] + unit.dy[h + 1]) *
                                     unit.soil_heat_conductivity_vector[h] * (unit.t1[h, i + 2] - unit.t1[h, i + 1]) /
                                     unit.dx[i + 1] + 1 / 8 * q_in_out_surrounding) *                   # heat source / sink
                                    unit.dt * 1 /
                                    (unit.soil_density_vector[h] * unit.cp1[h, i] * volume_adjacent_to_pipe / 5)

                # upper right node
                unit.t2[h - 1, i + 1] = unit.t1[h - 1, i + 1] +
                                        (unit.dz * 0.5 * (unit.dy[h - 2] + unit.dy[h - 1]) *                    # heat conduction in positive x-direction
                                         unit.soil_heat_conductivity_vector[h] *
                                         (unit.t1[h - 1, i + 2] - unit.t1[h - 1, i + 1]) / unit.dx[i] +
                                         unit.dz * 0.5 * (unit.dx[i] + unit.dx[i + 1]) *                        # heat conduction in negative y-direction
                                         unit.soil_heat_conductivity_vector[h] *
                                         (unit.t1[h - 2, i + 1] - unit.t1[h - 1, i + 1]) / (unit.dy[h - 2]) +
                                         1 / 8 * q_in_out_surrounding) *                                        # heat source / sink
                                        unit.dt * 1 /
                                        (unit.soil_density_vector[h] * unit.cp1[h, i] * volume_adjacent_to_pipe / 5)

                # lower right node
                unit.t2[h + 1, i + 1] = unit.t1[h + 1, i + 1] +
                                        (unit.dz * 0.5 * (unit.dx[i] + unit.dx[i + 1]) *                        # heat conduction in positive y-direction
                                         unit.soil_heat_conductivity_vector[h] *
                                         (unit.t1[h + 2, i + 1] - unit.t1[h + 1, i + 1]) / (unit.dy[h + 1]) +
                                         unit.dz * 0.5 * (unit.dy[h] + unit.dy[h + 1]) *                        # heat conduction in positive x-direction
                                         unit.soil_heat_conductivity_vector[h] *
                                         (unit.t1[h + 1, i + 2] - unit.t1[h + 1, i + 1]) / unit.dx[i + 1] +
                                         1 / 8 * q_in_out_surrounding) *                                        # heat source / sink
                                        unit.dt * 1 /
                                        (unit.soil_density_vector[h] * unit.cp1[h, i] * volume_adjacent_to_pipe / 5)

                # adjacent nodes to fluid node have uniform temperature average_temperature_adjacent_to_pipe
                unit.average_temperature_adjacent_to_pipe = (unit.t2[h - 1, i] + unit.t2[h + 1, i] + unit.t2[h, i + 1] +
                                                             unit.t2[h + 1, i + 1] + unit.t2[h - 1, i + 1]) / 5
                unit.t2[h - 1, i] = copy(unit.average_temperature_adjacent_to_pipe)
                unit.t2[h + 1, i] = copy(unit.average_temperature_adjacent_to_pipe)
                unit.t2[h, i + 1] = copy(unit.average_temperature_adjacent_to_pipe)

                unit.t2[h + 1, i + 1] = copy(unit.average_temperature_adjacent_to_pipe)
                unit.t2[h - 1, i + 1] = copy(unit.average_temperature_adjacent_to_pipe)

                unit.t2[h - 1, i],
                unit.phase_change_state[h - 1, i],
                unit.cp2[h - 1, i] = freezing(unit,
                                              unit.t2[h - 1, i],
                                              unit.phase_change_state[h - 1, i])
                unit.t2[h + 1, i],
                unit.phase_change_state[h + 1, i],
                unit.cp2[h + 1, i] = freezing(unit,
                                              unit.t2[h + 1, i],
                                              unit.phase_change_state[h + 1, i])
                unit.t2[h, i + 1],
                unit.phase_change_state[h, i + 1],
                unit.cp2[h, i + 1] = freezing(unit,
                                              unit.t2[h, i + 1],
                                              unit.phase_change_state[h, i + 1])

                unit.t2[h, i + 1],
                unit.phase_change_state[h, i + 1],
                unit.cp2[h, i + 1] = freezing(unit,
                                              unit.t2[h + 1, i + 1],
                                              unit.phase_change_state[h + 1, i + 1])
                unit.t2[h, i + 1],
                unit.phase_change_state[h, i + 1],
                unit.cp2[h, i + 1] = freezing(unit,
                                              unit.t2[h - 1, i + 1],
                                              unit.phase_change_state[h - 1, i + 1])
            end
        end
    end

    for _ in 1:(unit.n_internal_timesteps)

        # Loop - vertical direction
        for h in 1:(length(unit.dy))

            # Loop - horizontal direction
            for i in 1:(length(unit.dx) + 1)

                # calculate temperature of adjacent nodes to fluid
                if unit.is_fluid_node[h, i]
                    # do nothing, because fuid temperature and surrounding nodes are already calculated before iteration started.

                elseif unit.is_pipe_surrounding[h, i]
                    # Skip fluid_reynolds_number-calculation. already calculated.      

                    # ground surface boundary 
                elseif h == 1
                    #! format: off
                    # left boundary
                    if i == 1
                        unit.t2[h, i] = unit.t1[h, i] +
                                        (unit.dz * 0.5 * unit.dx[i] * unit.surface_convective_heat_transfer_coefficient * (unit.ambient_temperature - unit.t1[h, i]) +
                                         unit.dz * 0.5 * unit.dx[i] * (1 - unit.surface_reflection_factor) * unit.global_radiation_power +
                                         unit.dz * 0.5 * unit.dx[i] * unit.surface_emissivity * unit.boltzmann_constant * ((unit.ambient_temperature + 273.15)^4 - (unit.t1[h, i] + 273.15)^4) +
                                         unit.dz * unit.dy[1] * unit.soil_heat_conductivity_vector[h] * (unit.t1[h, i + 1] - unit.t1[h, i]) / unit.dx[i] +
                                         unit.dz * unit.dx[i] * unit.soil_heat_conductivity_vector[h] * (unit.t1[h + 1, i] - unit.t1[h, i]) / unit.dy[h]) *
                                        unit.dt * 1 / (unit.soil_density_vector[h] * unit.cp1[h, i] * unit.dz * 0.5 * unit.dx[i] * unit.dy[h])

                        unit.t2[h, i],
                            unit.phase_change_state[h, i],
                            unit.cp2[h, i] = freezing(unit::GeothermalHeatCollector,
                                                      unit.t2[h, i],
                                                      unit.phase_change_state[h, i])

                    # right boundary
                    elseif i == (length(unit.dx) + 1)
                        unit.t2[h, i] = unit.t1[h, i] +
                                        (unit.dz * 0.5 * unit.dx[i - 1] * unit.surface_convective_heat_transfer_coefficient * (unit.ambient_temperature - unit.t1[h, i]) +                           # Convection on surface
                                         unit.dz * 0.5 * unit.dx[i - 1] * (1 - unit.surface_reflection_factor) * unit.global_radiation_power +                                                       # Radiation on surface                                 
                                         unit.dz * 0.5 * unit.dx[i - 1] * unit.surface_emissivity * unit.boltzmann_constant * ((unit.ambient_temperature + 273.15)^4 - (unit.t1[h, i] + 273.15)^4) + # Radiation out of surface
                                         unit.dz * unit.dy[h] * unit.soil_heat_conductivity_vector[h] * (unit.t1[h, i - 1] - unit.t1[h, i]) / unit.dx[i - 1] +
                                         unit.dz * unit.dx[i - 1] * unit.soil_heat_conductivity_vector[h] * (unit.t1[h + 1, i] - unit.t1[h, i]) / (unit.dy[h])) *
                                        unit.dt * 1 / (unit.soil_density_vector[h] * unit.cp1[h, i] * unit.dz * 0.5 * unit.dx[i - 1] * unit.dy[1])

                        unit.t2[h, i],
                            unit.phase_change_state[h, i],
                            unit.cp2[h, i] = freezing(unit::GeothermalHeatCollector,
                                                      unit.t2[h, i],
                                                      unit.phase_change_state[h, i])

                    else
                        unit.t2[h, i] = unit.t1[h, i] +
                                        (unit.dz * 0.5 * (unit.dx[i] + unit.dx[i - 1]) * unit.surface_convective_heat_transfer_coefficient *  (unit.ambient_temperature - unit.t1[h, i]) +                           # Convection on surface
                                         unit.dz * 0.5 * (unit.dx[i] + unit.dx[i - 1]) * (1 - unit.surface_reflection_factor) * unit.global_radiation_power +                                                        # Radiation on surface
                                         unit.dz * 0.5 * (unit.dx[i] + unit.dx[i - 1]) * unit.surface_emissivity * unit.boltzmann_constant * ((unit.ambient_temperature + 273.15)^4 - (unit.t1[h, i] + 273.15)^4) +  # Radiation out of surface
                                         unit.dz * unit.dy[h] * unit.soil_heat_conductivity_vector[h] * (unit.t1[h, i + 1] - unit.t1[h, i]) / unit.dx[i] +                                                           # heat conduction in positve x-direction
                                         unit.dz * unit.dy[h] * unit.soil_heat_conductivity_vector[h] * (unit.t1[h, i - 1] - unit.t1[h, i]) / unit.dx[i - 1] +                                                       # heat conduction in negative x-direction
                                         unit.dz * 0.5 * (unit.dx[i] + unit.dx[i - 1]) * unit.soil_heat_conductivity_vector[h] * (unit.t1[h + 1, i] - unit.t1[h, i]) / (unit.dy[1])) *                               # heat conduction in y-direction
                                        unit.dt * 1 / (unit.soil_density_vector[h] * unit.cp1[h, i] * unit.dz * 0.5 * (unit.dx[i] + unit.dx[i - 1]) * unit.dy[1])

                        unit.t2[h, i],
                            unit.phase_change_state[h, i],
                            unit.cp2[h, i] = freezing(unit::GeothermalHeatCollector,
                                                      unit.t2[h, i],
                                                      unit.phase_change_state[h, i])

                        # for validation purposes to set constant temperatures at ground surface.
                        # unit.t2[1,:] .= 9
                    end
                else
                    # left boundary
                    if i == 1 && !unit.is_pipe_surrounding[h, i]
                        unit.t2[h, i] = unit.t1[h, i] +
                                        (unit.dz * 0.5 * (unit.dy[h] + unit.dy[h - 1]) * unit.soil_heat_conductivity_vector[h] * (unit.t1[h, i + 1] - unit.t1[h, i]) / unit.dx[i] +  # heat conduction in positve x-direction
                                         unit.dz * 0.5 * unit.dx[i] * unit.soil_heat_conductivity_vector[h] * (unit.t1[h + 1, i] - unit.t1[h, i]) / (unit.dy[h]) +                   # heat conduction in positive y-direction
                                         unit.dz * 0.5 * unit.dx[i] * unit.soil_heat_conductivity_vector[h] * (unit.t1[h - 1, i] - unit.t1[h, i]) / (unit.dy[h - 1])) *              # heat conduction in negative y-direction
                                        unit.dt * 1 / (unit.soil_density_vector[h] * unit.cp1[h, i] * unit.dz * 0.5 * unit.dx[i] * 0.5 * (unit.dy[h] + unit.dy[h - 1]))              # heat conduction in negative y-direction

                        unit.t2[h, i],
                            unit.phase_change_state[h, i],
                            unit.cp2[h, i] = freezing(unit,
                                                      unit.t2[h, i],
                                                      unit.phase_change_state[h, i])

                    # right boundary 
                    elseif i == length(unit.dx) + 1
                        unit.t2[h, i] = unit.t1[h, i] +
                                        (unit.dz * 0.5 * (unit.dy[h] + unit.dy[h - 1]) * unit.soil_heat_conductivity_vector[h] * (unit.t1[h, i - 1] - unit.t1[h, i]) / unit.dx[i - 1] +  # heat conduction in negative x-direction
                                         unit.dz * 0.5 * unit.dx[i - 1] * unit.soil_heat_conductivity_vector[h] * (unit.t1[h + 1, i] - unit.t1[h, i]) / (unit.dy[h]) +                   # heat conduction in positive y-direction
                                         unit.dz * 0.5 * unit.dx[i - 1] * unit.soil_heat_conductivity_vector[h] * (unit.t1[h - 1, i] - unit.t1[h, i]) / (unit.dy[h - 1])) *              # heat conduction in negative y-direction
                                        unit.dt * 1 / (unit.soil_density_vector[h] * unit.cp1[h, i] * unit.dz * 0.5 * unit.dx[i - 1] * 0.5 * (unit.dy[h] + unit.dy[h - 1]))

                        unit.t2[h, i],
                            unit.phase_change_state[h, i],
                            unit.cp2[h, i] = freezing(unit::GeothermalHeatCollector,
                                                      unit.t2[h, i],
                                                      unit.phase_change_state[h, i])

                    else
                        unit.t2[h, i] = unit.t1[h, i] +
                                        (unit.dz * 0.5 * (unit.dy[h] + unit.dy[h - 1]) * unit.soil_heat_conductivity_vector[h] * (unit.t1[h, i + 1] - unit.t1[h, i]) / unit.dx[i] +         # heat conduction in positve x-direction
                                         unit.dz * 0.5 * (unit.dy[h] + unit.dy[h - 1]) * unit.soil_heat_conductivity_vector[h] * (unit.t1[h, i - 1] - unit.t1[h, i]) / unit.dx[i - 1] +     # heat conduction in negative x-direction
                                         unit.dz * 0.5 * (unit.dx[i] + unit.dx[i - 1]) * unit.soil_heat_conductivity_vector[h] * (unit.t1[h + 1, i] - unit.t1[h, i]) / (unit.dy[h]) +       # heat conduction in positive y-direction
                                         unit.dz * 0.5 * (unit.dx[i] + unit.dx[i - 1]) * unit.soil_heat_conductivity_vector[h] * (unit.t1[h - 1, i] - unit.t1[h, i]) / (unit.dy[h - 1])) *  # heat conduction in negative y-direction
                                        unit.dt * 1 / (unit.soil_density_vector[h] * unit.cp1[h, i] * unit.dz * 0.5 * (unit.dx[i] + unit.dx[i - 1]) * 0.5 * (unit.dy[h] + unit.dy[h - 1]))   

                        unit.t2[h, i],
                            unit.phase_change_state[h, i],
                            unit.cp2[h, i] = freezing(unit::GeothermalHeatCollector,
                                                      unit.t2[h, i],
                                                      unit.phase_change_state[h, i])
                    end
                    #! format: on
                end
            end
        end
        unit.t1 = copy(unit.t2)
        unit.cp1 = copy(unit.cp2)
    end
    unit.temp_field_output[Int(sim_params["time"] / sim_params["time_step_seconds"]) + 1, :, :] = copy(unit.t2)

    activate_plot = true
    if (Int(sim_params["time"] / sim_params["time_step_seconds"]) == sim_params["number_of_time_steps"] - 1) &&
       activate_plot
        f = Figure()
        ax = Axis3(f[1, 1])

        ax.zlabel = "Temperature [°C]"
        ax.xlabel = "Vertical expansion (depth) [m]"
        ax.ylabel = "Horizontal expansion [m]"
        min_temp = minimum(unit.temp_field_output)
        min_temp = min_temp < 0.0 ? 1.1 * min_temp : 0.9 * min_temp
        max_temp = maximum(unit.temp_field_output)
        max_temp = max_temp < 0.0 ? 0.9 * max_temp : 1.1 * max_temp
        zlims!(ax, min_temp, max_temp)

        y_abs = [0; cumsum(unit.dx)]  # Absolute x coordinates
        x_abs = [0; cumsum(unit.dy)]  # Absolute y coordinates

        time = Observable(1)
        surfdata = @lift(unit.temp_field_output[$time, :, :])
        surface!(ax, x_abs, y_abs, surfdata)
        scatter!(ax, x_abs, y_abs, surfdata)
        slg = SliderGrid(f[2, 1], (; range=1:1:sim_params["number_of_time_steps"], label="Time"))

        on(slg.sliders[1].value) do v
            time[] = v
        end
        wait(display(f))
    end
end

# function to check/set phase state and calculate specific heat capacity with apparent heat capacity method. Not validated yet.
function freezing(unit::GeothermalHeatCollector, t2_in, phase_change_state_in)
    cp = unit.soil_specific_heat_capacity

    if t2_in < unit.phase_change_upper_boundary_temperature && phase_change_state_in == 0
        t2_out = unit.phase_change_upper_boundary_temperature
        phase_change_state_out = 1
    elseif t2_in > unit.phase_change_upper_boundary_temperature && phase_change_state_in == 1
        t2_out = unit.phase_change_upper_boundary_temperature
        phase_change_state_out = 0
    elseif t2_in < unit.phase_change_lower_boundary_temperature && phase_change_state_in == 1
        t2_out = unit.phase_change_lower_boundary_temperature
        phase_change_state_out = 2
    elseif t2_in > unit.phase_change_lower_boundary_temperature && phase_change_state_in == 2
        t2_out = unit.phase_change_lower_boundary_temperature
        phase_change_state_out = 1
    else
        t2_out = t2_in
        phase_change_state_out = phase_change_state_in
    end

    if t2_in > unit.phase_change_upper_boundary_temperature
        cp = unit.soil_specific_heat_capacity
    elseif t2_in <= unit.phase_change_upper_boundary_temperature &&
           t2_in >= unit.phase_change_lower_boundary_temperature
        dt_lat = unit.phase_change_upper_boundary_temperature - unit.phase_change_lower_boundary_temperature
        sigma_lat = 1 / 5 * dt_lat
        t_lat = (unit.phase_change_upper_boundary_temperature + unit.phase_change_lower_boundary_temperature) / 2
        cp = unit.soil_specific_enthalpy_of_fusion * 1 / (sigma_lat * sqrt(2 * pi)) *
             exp(-0.5 * (t2_in - t_lat)^2 / sigma_lat^2)
    elseif t2_in < unit.phase_change_lower_boundary_temperature
        cp = unit.soil_specific_heat_capacity - 100
    end

    return (t2_out, phase_change_state_out, cp)
end

# function to calculate heat transfer coefficient alpha.
function calculate_alpha_pipe(unit::GeothermalHeatCollector, q_in_out::Float64)

    # calculate mass flow in pipe
    collector_power_in_out_per_pipe = wh_to_watts(abs(q_in_out)) / unit.number_of_pipes  # W/pipe
    temperature_spread = q_in_out > 0 ? unit.loading_temperature_spread : unit.unloading_temperature_spread
    collector_mass_flow_per_pipe = collector_power_in_out_per_pipe /
                                   (unit.fluid_specific_heat_capacity * temperature_spread)  # kg/s

    use_dynamic_fluid_properties = true
    if use_dynamic_fluid_properties
        # calculate reynolds-number based on dynamic viscosity using dynamic temperature-dependend fluid properties, adapted from TRNSYS Type 710:
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
    approach = :stephan  # can be one of :ramming or :stephan

    if approach == :ramming
        # Approach used in Ramming 2007 from Elsner, Norbert; Fischer, Siegfried; Huhn, Jörg; „Grundlagen der Technischen Thermodynamik“,  Band 2 Wärmeübertragung, Akademie Verlag, Berlin 1993. 
        k_a = 1.1 - 1 / (3.4 + 0.0667 * unit.fluid_prandtl_number)
        k_n = 0.35 + 1 / (7.825 + 2.6 * sqrt(unit.fluid_prandtl_number))

        # calculate Nu-Number
        nusselt_laminar = ((k_a / (1 - k_n) *
                            (unit.fluid_prandtl_number * unit.pipe_d_i * fluid_reynolds_number / unit.pipe_length)^k_n)^3 +
                           4.364^3)^(1 / 3)
    elseif approach == :stephan
        # Stephan
        pr_water = 13.44                # Pr Number Water 0 °C as reference
        nusselt_laminar = 3.66 +
                          (0.0677 *
                           (fluid_reynolds_number * unit.fluid_prandtl_number * unit.pipe_d_i / unit.pipe_length)^1.33) /
                          (1 +
                           0.1 * unit.fluid_prandtl_number *
                           (fluid_reynolds_number * unit.pipe_d_i / unit.pipe_length)^0.83) *
                          (unit.fluid_prandtl_number / pr_water)^0.1
    end
    return nusselt_laminar
end

function calculate_Nu_turbulent(unit::GeothermalHeatCollector, fluid_reynolds_number::Float64)
    # Approached used from Gnielinski in: V. Gnielinski: Ein neues Berechnungsverfahren für die Wärmeübertragung im Übergangsbereich zwischen laminarer und turbulenter Rohrströmung. Forsch im Ing Wes 61:240–248, 1995. 
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
            exchange.temperature_min > unit.current_input_temperature ||
            exchange.temperature_max !== nothing &&
            exchange.temperature_max < unit.current_input_temperature)
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

    energy_taken = unit.max_input_energy - energy_demand
    unit.collector_total_heat_energy_in_out += energy_taken
    calculate_new_temperature_field!(unit::GeothermalHeatCollector, unit.collector_total_heat_energy_in_out, sim_params)
end

function balance_on(interface::SystemInterface,
                    unit::GeothermalHeatCollector)::Vector{EnergyExchange}
    caller_is_input = unit.uac == interface.target.uac
    purpose_uac = unit.uac == interface.target.uac ? interface.target.uac : interface.source.uac

    return [EnEx(; balance=interface.balance,
                 energy_potential=caller_is_input ? -unit.max_input_energy : unit.max_output_energy,
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
    append!(output_vals,
            [string(unit.m_heat_out) * " OUT",
             "fluid_temperature",
             "fluid_reynolds_number",
             "ambient_temperature",
             "global_radiation_power"])
    # push!(output_vals, "TEMPERATURE_#NodeNum")

    return output_vals
end

function output_value(unit::GeothermalHeatCollector, key::OutputKey)::Float64
    if key.value_key == "IN"
        return calculate_energy_flow(unit.input_interfaces[key.medium])
    elseif key.value_key == "OUT"
        return calculate_energy_flow(unit.output_interfaces[key.medium])
    elseif startswith(key.value_key, "Temperature_")
        idx = parse(Int, split(key.value_key, "_")[2])
        if !(1 <= idx <= length(unit.dy) + 1) # unit.t2!
            throw(ArgumentError("Index \"$idx\" of requested temperature-output of geothermal heat collector exeeds the number of available temperatur datapoints."))
        else
            return unit.t2[idx, 1] # unit.t2! # TODO: Every Node
        end
    elseif key.value_key == "fluid_temperature"
        return unit.fluid_temperature
    elseif key.value_key == "fluid_reynolds_number"
        return unit.fluid_reynolds_number
    elseif key.value_key == "ambient_temperature"
        return unit.ambient_temperature
    elseif key.value_key == "global_radiation_power"
        return unit.global_radiation_power
    end
    throw(KeyError(key.value_key))
end

export GeothermalHeatCollector
