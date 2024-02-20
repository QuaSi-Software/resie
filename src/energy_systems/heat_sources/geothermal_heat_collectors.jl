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
    last_timestep_calculated::Float64

    soil_starting_temperature::Temperature
    soil_specific_heat_capacity::Float64
    soil_density::Float64
    soil_heat_conductivity::Float64
    soil_density_vector::Vector
    soil_heat_conductivity_vector::Vector
    soil_heat_fusion_energy::Float64

    phase_change_upper_boundary_temperature::Temperature
    phase_change_lower_boundary_temperature::Temperature

    surface_convective_heat_transfer_coefficient::Float64
    surface_reflection_factor::Float64

    t1::Array{Float64}
    t2::Array{Float64}
    phase_change_state::Array{Float64}
    cp1::Array{Float64}
    cp2::Array{Float64}
    fluid::Array{Float64}
    pipe_surrounding::Array{Float64}

    pipe_radius_outer::Float64
    pipe_thickness::Float64
    pipe_d_i::Float64
    pipe_d_o::Float64
    pipe_laying_depth::Float64
    pipe_length::Float64
    number_of_pipes::Float64
    pipe_spacing::Float64

    global_radiation::Float64
    boltzmann_constant::Float64
    surface_emissivity::Float64

    dx::Vector{Float64}
    dy::Vector{Float64}
    dz::Float64
    dt::Integer

    timestep_index::Integer

    fluid_temperature::Temperature
    pipe_temperature::Temperature

    collector_total_heat_energy_in_out::Float64

    pipe_heat_conductivity::Float64
    fluid_specific_heat_capacity::Float64
    fluid_prandtl_number::Float64
    fluid_density::Float64
    fluid_kinematic_viscosity::Float64
    fluid_heat_conductivity::Float64

    fluid_reynolds_number::Float64
    average_temperature_adjacent_to_pipe::Float64

    function GeothermalHeatCollector(uac::String, config::Dict{String,Any}, sim_params::Dict{String,Any})
        m_heat_in = Symbol(default(config, "m_heat_in", "m_h_w_ht1"))
        m_heat_out = Symbol(default(config, "m_heat_out", "m_h_w_lt1"))
        register_media([m_heat_in, m_heat_out])
        # read in ambient temperature profilfe (downloaded from dwd)
        ambient_temperature_profile = "ambient_temperature_profile_path" in keys(config) ?
                                      Profile(config["ambient_temperature_profile_path"], sim_params) :
                                      nothing            # TODO: add global weather file
        # read in ambient global radiation profilfe (downloaded from dwd)
        global_radiation_profile = "global_radiation_profile_path" in keys(config) ?
                                   Profile(config["global_radiation_profile_path"], sim_params) :
                                   nothing               # TODO: add global weather file                   
        
        return new(
            uac,                    # uac
            controller_for_strategy( # controller
                config["strategy"]["name"], config["strategy"], sim_params
            ),
            sf_storage,             # sys_function
            InterfaceMap(           # input_interfaces
                m_heat_in => nothing
            ),
            InterfaceMap(           # output_interfaces
                m_heat_out => nothing
            ),
            m_heat_in,                      # medium name of input interface
            m_heat_out,                     # medium name of output interface
            ambient_temperature_profile,    # ambient temperature profile
            global_radiation_profile,
            default(config, "unloading_temperature_spread", 3),   # temperature spread between forward and return flow during unloading            
            default(config, "loading_temperature", nothing),      # nominal high temperature for loading geothermal heat collector storage, can also be set from other end of interface TODO
            default(config, "loading_temperature_spread", 3),     # temperature spread between forward and return flow during loading         
            default(config, "max_output_power", 20),              # maximum output power in W/m^2, set by user. Depending on ground and climate localization. [VDI 4640-2.]
            default(config, "max_input_power", 20),               # maximum input power in W/m^2, set by user. Depending on ground and climate localization. [VDI 4640-2.]
            default(config, "regeneration", true),                # flag if regeneration should be taken into account
            0.0,                        # max_output_energy in every time step, calculated in control()
            0.0,                        # max_input_energy in every time step, calculated in control()
            0.0,                        # output temperature in current time step, calculated in control()
            0.0,                        # input temperature in current time step, calculated in control()
            0.0,                        # ambient temperature in current time step, calculated in control()
            -1.0,                       # last timestep that was calculated; used to avoid double calculation of temperature field. set to -1 for the beginning
            15.5,                       # starting temperature near pipe. Value currently set for validation purpose. 
            default(config, "soil_specific_heat_capacity", 1000),  # specific heat capacity soil, 840, in J/(KgK)
            default(config, "soil_density", 2000),  # in kg/m^3              
            default(config, "soil_heat_conductivity", 1.5), # soil heat conductivity in W/(mK), 
            zeros(59),              # vector to set soil density. Size depends on discretitaion
            zeros(59),              # lamda soil. TODO: Use Soil propertys of validation project! 1.9
            default(config, "soil_heat_fusion_energy", 90000),                  # Heat fusion Energy in J/kg. Currently out of TRNSYS Default value.
            -0.25,                  # phase_change_upper_boundary_temperature
            -1,                     # phase_change_lower_boundary_temperature
            14.7,                   # convective heat transfer on surface - TODO: Search for literature source.
            0.25,                   # Reflection factor / Albedo value of surface. 
            zeros(31,7),            # t1.                   # TODO / Attention: Size of this matrix depends on discretization, which will be done in Pre-Processing. 
            zeros(31,7),            # t2 grob:zeros(31,7)   # TODO / Attention: Size of this matrix depends on discretization, which will be done in Pre-Processing. 
            zeros(31,7),            # phase_change_state    # TODO / Attention: Size of this matrix depends on discretization, which will be done in Pre-Processing. 
            zeros(31,7),            # array to define specific heat capacity on previous timestep (for apparent heat capacity method)     # TODO / Attention: Size of this matrix depends on discretization, which will be done in Pre-Processing. 
            zeros(31,7),            # array to define specific heat capacity on next timestep (for apparent heat capacity method)         # TODO / Attention: Size of this matrix depends on discretization, which will be done in Pre-Processing. 
            zeros(31,7),            # Matrix, which identifies fuid node. fuid node: 1. Soil-node: 0.                         TODO: Need to adapt on Pre-Processing.
            zeros(31,7),            # Matrix, which identifies nodes surrounding the fluid-node. Pipe-surrounding: 1. Else: 0.  TODO: Need to adapt on Pre-Processing. 
            default(config, "pipe_radius_outer", 0.016),        # pipe outer radius
            default(config, "pipe_thickness", 0.003),           # thickness of pipe in m
            0.0,                                                # pipe inner diamater, calculated in initailize() 
            0.0,                                                # pipe outer diamater, calculated in initailize() 
            default(config, "pipe_laying_depth", 1.5),          # deepth of pipe system below the ground in m
            default(config, "pipe_length", 100),                # pipe length of one collector in m
            default(config, "number_of_pipes", 1),              # numbers of paralell pipes, each with a length of "pipe_length"        
            default(config, "pipe_spacing", 0.5),               # distance between pipes of collector in m.
            0.0,                    # global radiation on surface. to be read in by weather-profile
            5.67e-8,                # Boltzmann-Constant
            0.9,                    # Emissivity on ground surface
            zeros(6),               # size dx. depends on discretization settings.
            zeros(30),              # size dy. depends on discretization settings.
            93.5,                   # dz in m, is constant. fuid-Direction. # Depending on pipe length. 
            20.0,                   # duration time of internal time-step dt in s depending on ground properties. TODO: Calculate based on discretization and soil properties.
            0.0,                    # time step index. necessary for validation purposes.                     
            10.0,                   # set starting fluid temperature.               
            0.0,                    # pipe temperature. 
            0.0,                    # total heat flux in or out of collector. set by ReSiE Interface.
            0.5,                    # pipe heat conductivity.               
            default(config, "fluid_specific_heat_capacity", 3800),      # fluid_specific_heat_capacity in J/(kg K)
            default(config, "fluid_prandtl_number", 30),                 # prandtl number at 30 % glycol, 0 °C 
            default(config, "fluid_density", 1045),                     # fluid density at 30 % glycol, 0 °C in kg/m^3
            default(config, "fluid_kinematic_viscosity", 3.9e-6),       # fluid_kinematic_viscosity at 30 % glycol, 0 °C in m^2/s
            default(config, "fluid_heat_conductivity",0.5),             # fluid_heat_conductivity at 30 % glycol, 0 °C  in W/(mK)
            0,              # fluid_reynolds_number-Number, to be calculated in function.
            16.0            # starting temperature of pipe. set for validation purpose.
            )
    end
end

function initialise!(unit::GeothermalHeatCollector, sim_params::Dict{String,Any})
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
end

function control(
    unit::GeothermalHeatCollector,
    components::Grouping,
    sim_params::Dict{String,Any}
)

    unit.timestep_index += 1

    # reset energy summarizer
    unit.collector_total_heat_energy_in_out = 0.0

    # Discretization and starting Temperature-Field. --> ADD PREPROCESSING! TODO
    if unit.timestep_index == 1
        
        # calculate diameters of pipe
        unit.pipe_d_i = 2 * unit.pipe_radius_outer - (2 * unit.pipe_thickness)
        unit.pipe_d_o = 2 * unit.pipe_radius_outer 

        # calculate max_energy
        A_collector = unit.pipe_length * unit.pipe_spacing * (unit.number_of_pipes - 1)
        unit.max_output_energy = watt_to_wh(unit.max_output_power * A_collector)
        unit.max_input_energy = watt_to_wh(unit.max_input_power * A_collector)

        # Discretization. Will be done in Pre-Processing. TODO.
        dx_R = [unit.pipe_radius_outer]
        dx_RM = [unit.pipe_radius_outer, unit.pipe_radius_outer, 2*unit.pipe_radius_outer, 4*unit.pipe_radius_outer]
        dx_End = unit.pipe_spacing/2 - sum(dx_R) - sum(dx_RM)
        unit.dx = [dx_R..., dx_RM..., dx_End...] # append()
            
        dy_O = [unit.pipe_radius_outer/2, unit.pipe_radius_outer, 4*unit.pipe_radius_outer, 8*unit.pipe_radius_outer, 16*unit.pipe_radius_outer]     # dy near surface
        dy_MR = [16*unit.pipe_radius_outer, 8*unit.pipe_radius_outer, 4*unit.pipe_radius_outer, 2*unit.pipe_radius_outer, unit.pipe_radius_outer/2]  # dy getting smaller untill pipe node
        dy_M = [(unit.pipe_laying_depth - (sum(dy_O)+sum(dy_MR)) - unit.pipe_radius_outer)]                       # dy is wide in the middle    
        dy_R = [unit.pipe_radius_outer, unit.pipe_radius_outer]                                                   # distance between nodes adjacent to fluid
        dy_RU = [unit.pipe_radius_outer/2, 
                1*unit.pipe_radius_outer,
                1*unit.pipe_radius_outer,
                3/2*unit.pipe_radius_outer,
                2*unit.pipe_radius_outer,
                2*unit.pipe_radius_outer,
                4*unit.pipe_radius_outer,
                16*unit.pipe_radius_outer,
                32*unit.pipe_radius_outer,
                64*unit.pipe_radius_outer,
                64*unit.pipe_radius_outer,
                64*unit.pipe_radius_outer,
                128*unit.pipe_radius_outer,
                64*unit.pipe_radius_outer,
                32*unit.pipe_radius_outer,
                16*unit.pipe_radius_outer,
                8*unit.pipe_radius_outer]  # mesh below the fluid-node until lower simulation boundary   
        unit.dy = [dy_O..., dy_M..., dy_MR..., dy_R..., dy_RU...]   #size:27     # build dy-vector
        
        unit.dz = copy(unit.pipe_length)
           
        # Localize fuid and adjacent Nodes. They will be calculated separatly.
        unit.fluid[(length(dy_O)+ length(dy_MR)+ length(dy_M)+length(dy_R)+1 -1),1] = 1     
        unit.pipe_surrounding[(length(dy_O)+ length(dy_MR)+ length(dy_M)+length(dy_R)+1) ,1] = 1
        unit.pipe_surrounding[(length(dy_O)+ length(dy_MR)+ length(dy_M)+length(dy_R)+1 -2) , 1] = 1
        unit.pipe_surrounding[(length(dy_O)+ length(dy_MR)+ length(dy_M)+length(dy_R)+1 -1) , 2] = 1

        unit.pipe_surrounding[(length(dy_O)+ length(dy_MR)+ length(dy_M)+length(dy_R)+1 -2),2] = 1
        unit.pipe_surrounding[(length(dy_O)+ length(dy_MR)+ length(dy_M)+length(dy_R)+1),2] = 1
        
        # set starting temperature distribution (current settings for validation purpose)
        unit.t1[1:3,:] .= 9.5
        unit.t1[4:6,:] .= 11.5
        unit.t1[7:9,:] .= 13.5
        unit.t1[10:13,:] .= 15.5 
        unit.t1[14:18,:] .= unit.soil_starting_temperature 
        unit.t1[19:25,:] .= 15.5     
        unit.t1[26:30,:] .= 15.5
        unit.t1[length(unit.dy)+1, :] .= 9    # TODO: Use avg ambient temperature out of weather data set.
        unit.t2 = copy(unit.t1)

        # set soil heat conductivity, more layers possible. Currently, only one homogenous layer is considered.
        unit.soil_heat_conductivity_vector[:] .= unit.soil_heat_conductivity
        
        # set soil density, more layers possible.
        unit.soil_density_vector[:] .= unit.soil_density
        
        # specific heat capacity for each node needed, because of apparent heat capacity method.
        unit.cp1[:,:] .= unit.soil_specific_heat_capacity
        unit.cp2 = copy(unit.cp1)

        # internal time step according to TRNSYS Type 710 Model (Hirsch, Hüsing & Rockendorf 2017):
        # TODO: Values too high?
        # unit.dt = min(watt_to_wh(1.0)*60*60, Int((unit.soil_density * unit.soil_specific_heat_capacity * min(minimum(unit.dx), minimum(unit.dy))) / (4 * unit.soil_heat_conductivity)))

    end

    # in case there is a state machine for geothermal heat collectors
    move_state(unit, components, sim_params)

    # get ambient temperature and global radiation from profile for current time step if needed (probably only for geothermal collectors)
    # TODO: add link to global weather file
    unit.ambient_temperature = Profiles.value_at_time(unit.ambient_temperature_profile, sim_params["time"])
    unit.global_radiation = wh_to_watts(Profiles.value_at_time(unit.global_radiation_profile, sim_params["time"]))  # from Wh/m^2 to W/m^2

    unit.current_output_temperature = unit.fluid_temperature + unit.unloading_temperature_spread/2
    set_temperature!(unit.output_interfaces[unit.m_heat_out],
                     nothing,
                     unit.current_output_temperature
                     )

    set_max_energy!(unit.output_interfaces[unit.m_heat_out], unit.max_output_energy)

    # get input temperature for energy input (regeneration) and set temperature and max_energy to input interface
    if unit.regeneration
        unit.current_input_temperature = unit.average_temperature_adjacent_to_pipe - unit.loading_temperature_spread/2 # of geothermal heat collector 
        set_temperature!(unit.input_interfaces[unit.m_heat_in],
                         unit.current_input_temperature,
                         nothing
                         )

        set_max_energy!(unit.input_interfaces[unit.m_heat_in], unit.max_input_energy)
    end

end

function calculate_new_temperature_field!(unit::GeothermalHeatCollector, q_in_out::Float64)
    
    # calculate heat transfer coefficient and thermal resistance
    alpha_fluid, unit.fluid_reynolds_number = calculate_alpha_pipe(unit, q_in_out)   

    # calculation of pipe_thermal_resistance_length_specific with the approach by type 710 publication
    pipe_thermal_resistance_length_specific = (4 / pi * 
                                              (unit.pipe_d_o / (alpha_fluid * unit.pipe_d_i) + 
                                              (log(unit.pipe_d_o / unit.pipe_d_i) * unit.pipe_d_o) / (2 * unit.pipe_heat_conductivity) + 
                                              unit.dx[2] / (2 * unit.soil_heat_conductivity))) / 
                                              (2 * pi * unit.pipe_radius_outer)
    
    # calculate fluid temperature
    h = 1
    for h = 1: (length(unit.dy)) 
        for i =1:2  
            if unit.fluid[h,i] == 1         
                # calculate specific heat extraction for 1 pipe per length
                specific_heat_flux_pipe = (wh_to_watts(q_in_out) / (unit.pipe_length*unit.number_of_pipes))
                unit.fluid_temperature = unit.average_temperature_adjacent_to_pipe + pipe_thermal_resistance_length_specific * specific_heat_flux_pipe     # "+", because specific_heat_flux_pipe is negative during heat extraction.        
                    
                unit.t2[h,i] = unit.fluid_temperature
                
                Volume_adjacent_to_pipe = ((unit.dx[i]+unit.dx[i+1]/2) *                               # total x-Direction
                                          (unit.dy[h-2]/2+unit.dy[h-1]+unit.dy[h]+unit.dy[h+1]/2) -    # total y-Direction
                                          (pi*unit.dy[h]^2/2))*unit.dz                                 # 1/2 Area pipe.          
                
                q_in_out_surrounding = wh_to_watts(q_in_out) / (unit.number_of_pipes)
                
                # upper node               
                unit.t2[h-1,i] = unit.t1[h-1,i] +      
                                 (unit.dz * unit.dx[i]/2 * unit.soil_heat_conductivity_vector[h] * (unit.t1[h-2, i]-unit.t1[h-1,i])/(unit.dy[h-2]) +  # heat conduction in negative y-direction
                                 1/16* q_in_out_surrounding)  *                                                                                       # heat source / sink
                                 unit.dt * 1/(unit.soil_density_vector[h] * unit.cp1[h,i] * Volume_adjacent_to_pipe/5)                                                               
                
                # lower node
                unit.t2[h+1,i] = unit.t1[h+1,i] +    
                                (unit.dz * unit.dx[i]/2 * unit.soil_heat_conductivity_vector[h] * (unit.t1[h+2,i]-unit.t1[h+1,i])/(unit.dy[h+1]) +    # heat conduction in positive y-direction
                                1/16* q_in_out_surrounding)  *                                                                                        # heat source / sink
                                unit.dt * 1/(unit.soil_density_vector[h] * unit.cp1[h,i] * Volume_adjacent_to_pipe/5)                                                       
                    
                # right node
                unit.t2[h,i+1] = unit.t1[h,i+1] +
                                (unit.dz * 0.5 * (unit.dy[h]+unit.dy[h+1]) * unit.soil_heat_conductivity_vector[h] * (unit.t1[h,i+2] - unit.t1[h,i+1])/unit.dx[i+1] + 
                                1/8 * q_in_out_surrounding)*                                                                                          # heat source / sink
                                unit.dt * 1/(unit.soil_density_vector[h] * unit.cp1[h,i] * Volume_adjacent_to_pipe/5)  

                # upper right node
                unit.t2[h-1,i+1] = unit.t1[h-1,i+1] +
                                   (unit.dz * 0.5 * (unit.dy[h-2]+unit.dy[h-1]) * unit.soil_heat_conductivity_vector[h] * (unit.t1[h-1,i+2] - unit.t1[h-1,i+1])/unit.dx[i] +       # heat conduction in positive x-direction
                                   unit.dz * 0.5*(unit.dx[i]+unit.dx[i+1]) * unit.soil_heat_conductivity_vector[h] * (unit.t1[h-2, i+1]-unit.t1[h-1,i+1])/(unit.dy[h-2]) +         # heat conduction in negative y-direction
                                   1/8* q_in_out_surrounding)*                                                                                                                     # heat source / sink
                                   unit.dt * 1/(unit.soil_density_vector[h] * unit.cp1[h,i] * Volume_adjacent_to_pipe/5) 

                # lower right node
                unit.t2[h+1,i+1] = unit.t1[h+1,i+1] +
                                   (unit.dz * 0.5*(unit.dx[i]+unit.dx[i+1]) * unit.soil_heat_conductivity_vector[h] * (unit.t1[h+2,i+1]-unit.t1[h+1,i+1])/(unit.dy[h+1]) +         # heat conduction in positive y-direction
                                   unit.dz * 0.5 * (unit.dy[h]+unit.dy[h+1]) * unit.soil_heat_conductivity_vector[h] * (unit.t1[h+1,i+2] - unit.t1[h+1,i+1])/unit.dx[i+1] +        # heat conduction in positive x-direction
                                   1/8* q_in_out_surrounding)*                                                                                                                     # heat source / sink
                                   unit.dt * 1/(unit.soil_density_vector[h] * unit.cp1[h,i] * Volume_adjacent_to_pipe/5)  

                # adjacent nodes to fluid node have uniform temperature average_temperature_adjacent_to_pipe
                unit.average_temperature_adjacent_to_pipe = (unit.t2[h-1,i] + unit.t2[h+1,i] + unit.t2[h,i+1] + unit.t2[h+1,i+1] + unit.t2[h-1,i+1])/5
                unit.t2[h-1,i] = copy(unit.average_temperature_adjacent_to_pipe)
                unit.t2[h+1,i] = copy(unit.average_temperature_adjacent_to_pipe)
                unit.t2[h,i+1] = copy(unit.average_temperature_adjacent_to_pipe)

                unit.t2[h+1,i+1] = copy(unit.average_temperature_adjacent_to_pipe)
                unit.t2[h-1,i+1] = copy(unit.average_temperature_adjacent_to_pipe)

                unit.t2[h-1,i], unit.phase_change_state[h-1,i], unit.cp2[h-1,i] = freezing(unit::GeothermalHeatCollector, unit.t2[h-1,i],unit.phase_change_state[h-1,i]) 
                unit.t2[h+1,i], unit.phase_change_state[h+1,i], unit.cp2[h+1,i] = freezing(unit::GeothermalHeatCollector, unit.t2[h+1,i],unit.phase_change_state[h+1,i])
                unit.t2[h,i+1], unit.phase_change_state[h,i+1], unit.cp2[h,i+1] = freezing(unit::GeothermalHeatCollector, unit.t2[h,i+1],unit.phase_change_state[h,i+1])
                
                unit.t2[h,i+1], unit.phase_change_state[h,i+1], unit.cp2[h,i+1] = freezing(unit::GeothermalHeatCollector, unit.t2[h+1,i+1],unit.phase_change_state[h+1,i+1])
                unit.t2[h,i+1], unit.phase_change_state[h,i+1], unit.cp2[h,i+1] = freezing(unit::GeothermalHeatCollector, unit.t2[h-1,i+1],unit.phase_change_state[h-1,i+1])
 
            end    
        end         
    end 

    # calculate number of internal timesteps depending on internal time dt
    n_it = Int(round(watt_to_wh(1.0)/unit.dt * 3600))
    it = 0

    for it=1:n_it
        
        # Loop - vertical direction
        for h = 1: (length(unit.dy))

            # Loop - horizontal direction
            for i = 1 : (length(unit.dx)+1)
                
                # calculate temperature of adjacent nodes to fluid
                if unit.fluid[h,i] == 1             
                    # do nothing, because fuid temperature and surrounding nodes are already calculated before iteration started.
                
                elseif unit.pipe_surrounding[h,i]==1           
                    # Skip fluid_reynolds_number-calculation. already calculated.      
                
                # ground surface boundary 
                elseif h==1
                    
                    # left boundary
                    if i ==1
                        unit.t2[h,i] = unit.t1[h,i] +
                                       (unit.dz * 0.5 * unit.dx[i] * unit.surface_convective_heat_transfer_coefficient * (unit.ambient_temperature - unit.t1[h,i]) +
                                       unit.dz * 0.5 * unit.dx[i] * (1-unit.surface_reflection_factor) * unit.global_radiation +
                                       unit.dz * 0.5 * unit.dx[i] * unit.surface_emissivity * unit.boltzmann_constant * ((unit.ambient_temperature+273.15)^4 - (unit.t1[h,i]+273.15)^4) +
                                       unit.dz * unit.dy[1] * unit.soil_heat_conductivity_vector[h] * (unit.t1[h,i+1] - unit.t1[h,i])/unit.dx[i] +
                                       unit.dz * unit.dx[i] * unit.soil_heat_conductivity_vector[h] * (unit.t1[h+1,i] - unit.t1[h,i])/unit.dy[h]) * 
                                       unit.dt * 1/(unit.soil_density_vector[h] * unit.cp1[h,i] * unit.dz * 0.5 * unit.dx[i] * unit.dy[h])
                        
                        unit.t2[h,i], unit.phase_change_state[h,i], unit.cp2[h,i] = freezing(unit::GeothermalHeatCollector,unit.t2[h,i],unit.phase_change_state[h,i])
                    
                    # right boundary
                    elseif i == (length(unit.dx)+1)
                        unit.t2[h,i] = unit.t1[h,i] +
                                       (unit.dz * 0.5 * unit.dx[i-1] * unit.surface_convective_heat_transfer_coefficient * (unit.ambient_temperature - unit.t1[h,i]) +                       # Convection on surface
                                       unit.dz * 0.5 * unit.dx[i-1] * (1-unit.surface_reflection_factor)* unit.global_radiation +                                                            # Radiation on surface
                                       unit.dz * 0.5 * unit.dx[i-1] * unit.surface_emissivity * unit.boltzmann_constant * ((unit.ambient_temperature+273.15)^4 - (unit.t1[h,i]+273.15)^4) +  # Radiation out of surface
                                       unit.dz * unit.dy[h] * unit.soil_heat_conductivity_vector[h] * (unit.t1[h,i-1] - unit.t1[h,i])/unit.dx[i-1] +
                                       unit.dz * unit.dx[i-1] * unit.soil_heat_conductivity_vector[h] * (unit.t1[h+1,i]-unit.t1[h,i])/(unit.dy[h]) )* 
                                       unit.dt * 1/(unit.soil_density_vector[h]*unit.cp1[h,i]*unit.dz*0.5*unit.dx[i-1]*unit.dy[1])
                            

                        unit.t2[h,i], unit.phase_change_state[h,i], unit.cp2[h,i]= freezing(unit::GeothermalHeatCollector, unit.t2[h,i],unit.phase_change_state[h,i])
                            
                    else
                        unit.t2[h,i] = unit.t1[h,i] +
                                    (unit.dz * 0.5 * (unit.dx[i]+unit.dx[i-1]) * unit.surface_convective_heat_transfer_coefficient * (unit.ambient_temperature - unit.t1[h,i]) +                          # Convection on surface
                                    unit.dz * 0.5 * (unit.dx[i]+unit.dx[i-1]) * (1-unit.surface_reflection_factor)* unit.global_radiation +                                                               # Radiation on surface
                                    unit.dz * 0.5 * (unit.dx[i]+unit.dx[i-1]) * unit.surface_emissivity * unit.boltzmann_constant * ((unit.ambient_temperature+273.15)^4 - (unit.t1[h,i]+273.15)^4) +     # Radiation out of surface
                                    unit.dz * unit.dy[h] * unit.soil_heat_conductivity_vector[h] * (unit.t1[h,i+1] - unit.t1[h,i])/unit.dx[i] +                                                           # heat conduction in positve x-direction
                                    unit.dz * unit.dy[h] * unit.soil_heat_conductivity_vector[h] * (unit.t1[h,i-1] - unit.t1[h,i])/unit.dx[i-1] +                                                         # heat conduction in negative x-direction
                                    unit.dz * 0.5 * (unit.dx[i]+unit.dx[i-1]) * unit.soil_heat_conductivity_vector[h] *                                                                                   # heat conduction in y-direction
                                    (unit.t1[h+1,i]-unit.t1[h,i])/(unit.dy[1])) * unit.dt * 1/(unit.soil_density_vector[h] * unit.cp1[h,i] * unit.dz * 0.5*(unit.dx[i]+unit.dx[i-1])*unit.dy[1])          # heat conduction in y-direction
                            
                        unit.t2[h,i], unit.phase_change_state[h,i], unit.cp2[h,i]= freezing(unit::GeothermalHeatCollector, unit.t2[h,i],unit.phase_change_state[h,i])
                    
                        # for validation purposes to set constant temperatures at ground surface.
                        # unit.t2[1,:] .= 9
                    end 
                else
                    # left boundary
                    if i == 1 && unit.pipe_surrounding[h,i] == 0
                        unit.t2[h,i] = unit.t1[h,i] +
                                       (unit.dz * 0.5 * (unit.dy[h]+unit.dy[h-1]) * unit.soil_heat_conductivity_vector[h] * (unit.t1[h,i+1] - unit.t1[h,i])/unit.dx[i] +     # heat conduction in positve x-direction
                                       unit.dz * 0.5* unit.dx[i] * unit.soil_heat_conductivity_vector[h] * (unit.t1[h+1,i]-unit.t1[h,i])/(unit.dy[h]) +                      # heat conduction in positive y-direction
                                       unit.dz * 0.5* unit.dx[i] * unit.soil_heat_conductivity_vector[h] * (unit.t1[h-1,i]-unit.t1[h,i])/(unit.dy[h-1])) *                   # heat conduction in negative y-direction
                                       unit.dt * 1/(unit.soil_density_vector[h] * unit.cp1[h,i] * unit.dz * 0.5 * unit.dx[i] * 0.5 * (unit.dy[h]+unit.dy[h-1]))              # heat conduction in negative y-direction
                            
                        
                        unit.t2[h,i], unit.phase_change_state[h,i], unit.cp2[h,i]= freezing(unit::GeothermalHeatCollector, unit.t2[h,i],unit.phase_change_state[h,i])
                        
                    # right boundary 
                    elseif i == length(unit.dx)+1
                        unit.t2[h,i] = unit.t1[h,i] +
                                       (unit.dz * 0.5 * (unit.dy[h]+unit.dy[h-1]) * unit.soil_heat_conductivity_vector[h] * (unit.t1[h,i-1] - unit.t1[h,i])/unit.dx[i-1] +  # heat conduction in negative x-direction
                                       unit.dz * 0.5*unit.dx[i-1] * unit.soil_heat_conductivity_vector[h] * (unit.t1[h+1,i]-unit.t1[h,i])/(unit.dy[h]) +                    # heat conduction in positive y-direction
                                       unit.dz * 0.5*unit.dx[i-1] * unit.soil_heat_conductivity_vector[h] * (unit.t1[h-1,i]-unit.t1[h,i])/(unit.dy[h-1])) *                 # heat conduction in negative y-direction
                                       unit.dt * 1/(unit.soil_density_vector[h] * unit.cp1[h,i] * unit.dz * 0.5 * unit.dx[i-1]* 0.5 * (unit.dy[h]+unit.dy[h-1]))            # heat conduction in negative y-direction
                        
                        unit.t2[h,i], unit.phase_change_state[h,i], unit.cp2[h,i]= freezing(unit::GeothermalHeatCollector, unit.t2[h,i],unit.phase_change_state[h,i])
                        
                    else
                        unit.t2[h,i] = unit.t1[h,i] +
                                       (unit.dz * 0.5 * (unit.dy[h]+unit.dy[h-1]) * unit.soil_heat_conductivity_vector[h] * (unit.t1[h,i+1] - unit.t1[h,i])/unit.dx[i] +       # heat conduction in positve x-direction
                                       unit.dz * 0.5 * (unit.dy[h]+unit.dy[h-1]) * unit.soil_heat_conductivity_vector[h] * (unit.t1[h,i-1] - unit.t1[h,i])/unit.dx[i-1] +      # heat conduction in negative x-direction
                                       unit.dz * 0.5 * (unit.dx[i]+unit.dx[i-1]) * unit.soil_heat_conductivity_vector[h] * (unit.t1[h+1,i]-unit.t1[h,i])/(unit.dy[h])+         # heat conduction in positive y-direction
                                       unit.dz * 0.5 * (unit.dx[i]+unit.dx[i-1]) * unit.soil_heat_conductivity_vector[h] * (unit.t1[h-1,i]-unit.t1[h,i])/(unit.dy[h-1])) *     # heat conduction in negative y-direction
                                       unit.dt * 1/(unit.soil_density_vector[h] * unit.cp1[h,i] * unit.dz * 0.5*(unit.dx[i]+unit.dx[i-1])* 0.5 * (unit.dy[h]+unit.dy[h-1]))    # heat conduction in negative y-direction
                            
                        unit.t2[h,i], unit.phase_change_state[h,i], unit.cp2[h,i]= freezing(unit::GeothermalHeatCollector, unit.t2[h,i],unit.phase_change_state[h,i]) 
                    end 
                end
            end 
        end
        unit.t1 = copy(unit.t2)
        unit.cp1 = copy(unit.cp2)
    end
end 
       

# function to check/set phase state and calculate specific heat capacity with apparent heat capacity method. Not validated yet.
function freezing(unit::GeothermalHeatCollector,t2_in, phase_change_state_in)
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
    elseif t2_in<= unit.phase_change_upper_boundary_temperature && t2_in >= unit.phase_change_lower_boundary_temperature
        dt_lat = unit.phase_change_upper_boundary_temperature - unit.phase_change_lower_boundary_temperature
        sigma_lat = 1/5 * dt_lat
        t_lat = (unit.phase_change_upper_boundary_temperature + unit.phase_change_lower_boundary_temperature)/2
        cp = unit.soil_heat_fusion_energy * 1/(sigma_lat* sqrt(2*pi)) * exp(-0.5 * (t2_in-t_lat)^2/sigma_lat^2)
    elseif t2_in<unit.phase_change_lower_boundary_temperature
        cp = unit.soil_specific_heat_capacity - 100
    end

    return (t2_out, phase_change_state_out, cp)
end 

# function to calculate heat transfer coefficient alpha.
function calculate_alpha_pipe(unit::GeothermalHeatCollector, q_in_out::Float64)

    # calculate mass flow in pipe
    collector_power_in_out_per_pipe = wh_to_watts(abs(q_in_out)) / unit.number_of_pipes  # W/pipe
    temperature_spread = q_in_out > 0 ? unit.loading_temperature_spread : unit.unloading_temperature_spread
    collector_mass_flow_per_pipe = collector_power_in_out_per_pipe / (unit.fluid_specific_heat_capacity * temperature_spread)  # kg/s

    use_dynamic_fluid_properties = false
    if use_dynamic_fluid_properties
        # calculate reynolds-number based on dynamic viscosity using dynamic temperature-dependend fluid properties, adapted from TRNSYS Type 710:
        fluid_dynamic_viscosity = 0.0000017158* unit.fluid_temperature^2 - 0.0001579079*unit.fluid_temperature+0.0048830621
        unit.fluid_heat_conductivity = 0.0010214286 * unit.fluid_temperature + 0.447
        unit.fluid_prandtl_number = fluid_dynamic_viscosity * unit.fluid_specific_heat_capacity / unit.fluid_heat_conductivity 
        fluid_reynolds_number = (unit.collector_total_heat_flux_in_out)/(unit.fluid_specific_heat_capacity * unit.unloading_temperature_spread * pi / 4 * d_i_pipe * fluid_dynamic_viscosity)
    else
        # calculate reynolds-number, based on kinematic viscosity with constant fluid properties.
        fluid_reynolds_number = (4 * collector_mass_flow_per_pipe) / (unit.fluid_density * unit.fluid_kinematic_viscosity * unit.pipe_d_i * pi)
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

    alpha = Nu * unit.fluid_heat_conductivity / unit.pipe_d_i

    return alpha, fluid_reynolds_number
end

function calculate_Nu_laminar(unit::GeothermalHeatCollector, fluid_reynolds_number::Float64)
    # Approach used in Ramming 2007 from Elsner, Norbert; Fischer, Siegfried; Huhn, Jörg; „Grundlagen der Technischen Thermodynamik“,  Band 2 Wärmeübertragung, Akademie Verlag, Berlin 1993. 
    k_a = 1.1 - 1 / (3.4 + 0.0667 * unit.fluid_prandtl_number)
    k_n = 0.35 + 1 / (7.825 + 2.6 * sqrt(unit.fluid_prandtl_number))

    # calculate Nu-Number
    Nu_laminar = ((k_a / (1 - k_n) * (unit.fluid_prandtl_number * unit.pipe_d_i * fluid_reynolds_number / unit.pipe_length)^k_n)^3 + 4.364^3)^(1 / 3)
    return Nu_laminar
end 

function calculate_Nu_turbulent(unit::GeothermalHeatCollector, fluid_reynolds_number::Float64)
    # Approached used from Gnielinski in: V. Gnielinski: Ein neues Berechnungsverfahren für die Wärmeübertragung im Übergangsbereich zwischen laminarer und turbulenter Rohrströmung. Forsch im Ing Wes 61:240–248, 1995. 
    zeta = (1.8 * log(fluid_reynolds_number) - 1.5)^-2
    Nu_turbulent = (zeta / 8 * fluid_reynolds_number * unit.fluid_prandtl_number) /
                   (1 + 12.7 * sqrt(zeta / 8) * (unit.fluid_prandtl_number^(2 / 3) - 1)) 
    return Nu_turbulent
end

# process function that provides energy from the geothermal heat collector and calculates new temperatures 
# according to actual delivered or received energy
function process(unit::GeothermalHeatCollector, sim_params::Dict{String,Any})
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
    unit.collector_total_heat_energy_in_out = energy_delivered
end

function load(unit::GeothermalHeatCollector, sim_params::Dict{String,Any})   
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
        calculate_new_temperature_field!(unit::GeothermalHeatCollector, unit.collector_total_heat_energy_in_out)  # call calculate_new_heat_collector_temperatures!() to calculate new temperatures of field to account for possible ambient effects
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

    energy_taken = unit.max_input_energy - energy_demand
    unit.collector_total_heat_energy_in_out += energy_taken
    calculate_new_temperature_field!(unit::GeothermalHeatCollector, unit.collector_total_heat_energy_in_out) 
end

function balance_on(
    interface::SystemInterface,
    unit::GeothermalHeatCollector
)::Vector{EnergyExchange}

caller_is_input = unit.uac == interface.target.uac

    return [EnEx(
        balance=interface.balance,
        uac=unit.uac,
        energy_potential=0.0,
        storage_potential=caller_is_input ? - unit.max_input_energy : unit.max_output_energy,   # TODO is this to be assuemd as storage_potential?
        temperature_min=interface.temperature_min,
        temperature_max=interface.temperature_max,
        pressure=nothing,
        voltage=nothing,
    )]
end

function output_values(unit::GeothermalHeatCollector)::Vector{String}
    return [string(unit.m_heat_in)*" IN",
            string(unit.m_heat_out)*" OUT",
            "TEMPERATURE_#NodeNum",
            "fluid_temperature"]
end

function output_value(unit::GeothermalHeatCollector, key::OutputKey)::Float64
    if key.value_key == "IN"
        return calculate_energy_flow(unit.input_interfaces[key.medium])
    elseif key.value_key == "OUT"
        return calculate_energy_flow(unit.output_interfaces[key.medium])
    elseif startswith(key.value_key, "Temperature_")
        idx = parse(Int, split(key.value_key, "_")[2])
        if !(1 <= idx <= length(unit.dy)+1) # unit.t2!
            throw(ArgumentError("Index \"$idx\" of requested temperature-output of geothermal heat collector exeeds the number of available temperatur datapoints.")) 
        else
            return unit.t2[idx,1] # unit.t2! # TODO: Every Node
        end
    elseif key.value_key =="fluid_temperature"
        return unit.fluid_temperature
    elseif key.value_key =="fluid_reynolds_number"
        return unit.fluid_reynolds_number
    elseif key.value_key =="ambient_temperature"
        return unit.ambient_temperature
    end
    throw(KeyError(key.value_key))
end

export GeothermalHeatCollector

