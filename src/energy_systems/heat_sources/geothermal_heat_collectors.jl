"""
Implementation of geothermal heat collector.
This implementations acts as storage as is can produce and load energy.
"""
# only for validation purposes. 
file = open("C:/Users/vollmer/Documents/GitHub/resie/src/energy_systems/heat_sources/q_in_out_nov_h.txt", "r")
    q_nov_wgg = Vector{Float64}()
    for line in eachline(file)
        push!(q_nov_wgg, parse(Float64, line))
    end
close(file)

mutable struct GeothermalHeatCollector <: ControlledComponent
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
    max_output_power::Float64
    max_input_power::Float64
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
    t1::Array
    t2::Array
    phase_change_state::Array
    cp1::Array
    cp2::Array
    fluid::Array
    pipe_surrounding::Array
    pipe_radius_outer::Float64
    pipe_thickness::Float64
    pipe_laying_depth::Float64
    pipe_length::Float64
    pipe_spacing::Float64
    global_radiation::Float64
    boltzmann_constant::Float64
    surface_emissivity::Float64
    dx::Vector
    dy::Vector
    dz::Float64
    dt::Integer
    timestep_index::Integer
    fluid_temperature::Temperature
    pipe_temperature::Temperature
    collector_total_heat_flux_in_out::Float64
    pipe_heat_conductivity::Float64
    fluid_specific_heat_capacity::Float64
    fluid_prantl_number::Float64
    fluid_density::Float64
    fluid_kinematic_viscosity::Float64
    fluid_heat_conductivity::Float64
    specific_heat_flux_in_out::Float64
    wh_to_w_timestep::Float64
    fluid_reynolds_number::Float64
    average_temperature_adjacent_to_pipe::Float64
    q_in_out_nov::Vector

    function GeothermalHeatCollector(uac::String, config::Dict{String,Any})
        m_heat_in = Symbol(default(config, "m_heat_in", "m_h_w_ht1"))
        m_heat_out = Symbol(default(config, "m_heat_out", "m_h_w_lt1"))
        register_media([m_heat_in, m_heat_out])
        # read in ambient temperature profilfe (downloaded from dwd)
        ambient_temperature_profile = "ambient_temperature_profile_path" in keys(config) ?
                                      Profile(config["ambient_temperature_profile_path"]) :
                                      nothing
        # read in ambient global radiation profilfe (downloaded from dwd)
        global_radiation_profile = "global_radiation_profile_path" in keys(config) ?
                                      Profile(config["global_radiation_profile_path"]) :
                                      nothing                              

        return new(
            uac,                    # uac
            controller_for_strategy( # controller
                config["strategy"]["name"], config["strategy"]
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
            default(config, "loading_temperature", nothing),        # nominal high temperature for loading geothermal heat collector storage, can also be set from other end of interface
            default(config, "loading_temperature_spread", 3),     # temperature spread between forward and return flow during loading         
            config["max_output_power"],     # maximum output power set by user, may change later to be calculated from other inputs like specific heat transfer rate
            config["max_input_power"],      # maximum input power set by user, may change later to be calculated from other inputs like specific heat transfer rate
            default(config, "regeneration", true),   # flag if regeneration should be taken into account
            0.0,                        # max_output_energy in every time step, calculated in control()
            0.0,                        # max_input_energy in every time step, calculated in control()
            0.0,                        # output temperature in current time step, calculated in control()
            0.0,                        # input temperature in current time step, calculated in control()
            0.0,                        # ambient temperature in current time step, calculated in control()
            -1.0,                       # last timestep that was calculated; used to avoid double calculation of temperature field. set to -1 for the beginning
            15.5,                       # starting temperature near pipe. Value currently set for validation purpose. 
            default(config, "soil_specific_heat_capacity", 1000),  # specific heat capacity soil, 840
            default(config, "soil_density", 2000),                
            default(config, "soil_heat_conductivity", 1.5), # soil heat conductivity, 
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
            default(config, "pipe_laying_depth", 1.5),          # deepth of pipe system below the ground in m
            default(config, "pipe_length", 100),                # pipe length of one collector in m
            default(config, "pipe_spacing", 0.5),               # distance between pipes of collector in m.
            0,                      # global radiation on surface. to be read in by weather-profile
            5.67e-8,                # Boltzmann-Constant
            0.9,                    # Emissivity on ground surface
            zeros(6),               # size dx. depends on discretization settings.
            zeros(30),              # size dy. depends on discretization settings.
            93.5,                   # dz in m, is constant. fuid-Direction. # Depending on pipe length. 
            20,                     # duration time of internal time-step dt in s depending on ground properties. TODO: Calculate based on discretization and soil properties.
            0,                      # time step index. necessary for validation purposes.                     
            10,                     # set starting fluid temperature.               
            0,                      # pipe temperature. 
            0,                      # total heat flux in or out of collector. set by ReSiE Interface.
            0.5,                    # pipe heat conductivity.               
            default(config, "fluid_specific_heat_capacity", 3800),      # fluid_specific_heat_capacity in J/(kg K)
            default(config, "fluid_prantl_number", 30),                 # prandtl number at 30 % glycol, 0 °C 
            default(config, "fluid_density", 1045),                     # fluid density at 30 % glycol, 0 °C
            default(config, "fluid_kinematic_viscosity", 3.9e-6),       # fluid_kinematic_viscosity at 30 % glycol, 0 °C
            default(config, "fluid_heat_conductivity",0.5),             # fluid_heat_conductivity at 30 % glycol, 0 °C
            default(config, "specific_heat_flux_in_out",20),            # max. specific heat flux extractable out of soil in W/m^2. Depending on ground and climate localization. [VDI 4640-2.]
            1,              # time step of simulation in h to calculate power (W) out of energy (Wh) (wh_to_w_timestep), depending on simulation time step.
            0,              # fluid_reynolds_number-Number, to be calculated in function.
            16.0,           # starting temperature of pipe. set for validation purpose.
            q_nov_wgg       # provisorisch!

            )
    end
end

function control(
    unit::GeothermalHeatCollector,
    components::Grouping,
    parameters::Dict{String,Any}
)

    unit.timestep_index += 1

    # Discretization and starting Temperature-Field. --> ADD PREPROCESSING!
    if unit.timestep_index == 1
        
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

    end

    # in case there is a state machine for geothermal heat collectors
    move_state(unit, components, parameters)

    # get ambient temperature and global radiation from profile for current time step if needed (probably only for geothermal collectors)
    unit.ambient_temperature = Profiles.value_at_time(unit.ambient_temperature_profile, parameters["time"])
    unit.global_radiation = (Profiles.value_at_time(unit.global_radiation_profile, parameters["time"])) / unit.wh_to_w_timestep * 100/36# /0.25 --> from Wh/m^2 to W/m^2

    unit.current_output_temperature = unit.fluid_temperature + unit.unloading_temperature_spread/2
    unit.output_interfaces[unit.m_heat_out].temperature = highest_temperature(
                                                                            unit.current_output_temperature, 
                                                                            unit.output_interfaces[unit.m_heat_out].temperature
                                                                            )

    # get input temperature for energy input (regeneration) and set temperature to input interface
    if unit.regeneration
        unit.current_input_temperature = current_input_temperature(unit) # of geothermal heat collector 
        unit.input_interfaces[unit.m_heat_in].temperature = highest_temperature(
                                                                                unit.current_input_temperature,
                                                                                unit.input_interfaces[unit.m_heat_in].temperature
                                                                                )
    end                                                                        

    # calculate maximum input and output energy that is possible in current time step
        # sets max_energy to zero if requested/available temperature does not fit to temperature of geothermal heat collector.
        # This works as the control step of transformers is always calculated earlier than the one of storages. If temperatures
        # are written to the connected interface by a transformer, this is already done at this point.
    if unit.output_interfaces[unit.m_heat_out].temperature > unit.current_output_temperature
        unit.max_output_energy = 0.0  # no energy can be provided if requested temperature is higher than max. temperature of heat collector
    else
        unit.max_output_energy = watt_to_wh(get_max_output_power(unit))  
    end

    if unit.regeneration
        if unit.input_interfaces[unit.m_heat_in].temperature < unit.current_input_temperature
            unit.max_input_energy = 0.0 # no energy can be taken if available temperature is less than minimum possible temperature to load the heat collector
        else
            unit.max_input_energy = watt_to_wh(get_max_input_power(unit)) 
        end
    end

    # set max_energy to interfaces to provide information for connected components
    set_max_energy!(unit.output_interfaces[unit.m_heat_out], unit.max_output_energy)
    if unit.regeneration
        set_max_energy!(unit.input_interfaces[unit.m_heat_in], unit.max_input_energy)
    end

end

# function that calculates current (highest possible) output temperature of heat collector that can be provided. (lower temp. is always possible!)
function current_output_temperature(unit::GeothermalHeatCollector)::Temperature
    highest_outer_borehole_temp =  unit.average_temperature_adjacent_to_pipe - unit.unloading_temperature_spread/2
    # could be the highest outer borewhole temperature of the heat collector in the last timestep
    return highest_outer_borehole_temp
end

# function that calculates current (minimum possible) input temperature of heat collector that should be provided for regeneration 
# (higher temperature is always possible, here the minimum should be calculated!)
function current_input_temperature(unit::GeothermalHeatCollector)::Temperature
    input_temperature =  unit.average_temperature_adjacent_to_pipe - unit.loading_temperature_spread/2
    # could be the lowest outer borewhole temperature plot the given temperature spread or a user-specified loading_temperature if given
    return input_temperature
end

# function to calculate the maximum possible output energy in current time step
function get_max_output_power(unit::GeothermalHeatCollector)::Float64
    
    # TODO. Calculate without n_pipe.
    n_pipe = 47
    A_collector = unit.pipe_length * unit.pipe_spacing*(n_pipe  -1)
    unit.max_output_power = A_collector * unit.specific_heat_flux_in_out
    
    return unit.max_output_power  # could be calculated by maximum specific unloading energy or by maximal mass flow or...
    
end

# function to calculate the maximum possible input energy in current time step
function get_max_input_power(unit::GeothermalHeatCollector)::Float64
    
    # TODO. Calculate without n_pipe.
    n_pipe = 47
    A_collector = unit.pipe_length * unit.pipe_spacing*(n_pipe  -1)
    unit.max_input_power = A_collector * unit.specific_heat_flux_in_out    # for validation only
    return unit.max_input_power  # could be calculated by maximum specific loading energy or by maximal mass flow or...
end

function calculate_new_temperature_field(unit::GeothermalHeatCollector)
    
    # calculate heat transfer coefficient and thermal resistance
    alpha_fluid = calculate_alpha_pipe(unit::GeothermalHeatCollector)   

    # calculation of pipe_thermal_resistance_length_specific with the approach by type 710 publication
    d_i_pipe = 2*unit.pipe_radius_outer - (2*unit.pipe_thickness)
    d_o_pipe = 2*unit.pipe_radius_outer
    pipe_thermal_resistance_length_specific = (4/pi * (d_o_pipe/(alpha_fluid*d_i_pipe)+(log(d_o_pipe/d_i_pipe)*d_o_pipe)/(2*unit.pipe_heat_conductivity)+unit.dx[2]/(2*unit.soil_heat_conductivity))) * 1/(2*pi*unit.pipe_radius_outer)
    # R_rohr_old = 1 / (2* pi) * (2 / (alpha_fluid * (2 * (unit.pipe_radius_outer - unit.pipe_thickness))) + 1 / unit.pipe_heat_conductivity * log(unit.pipe_radius_outer/(unit.pipe_radius_outer - unit.pipe_thickness)))    
        
    # calculate fluid temperature
    h = 1
    for h = 1: (length(unit.dy)) 
        
        for i =1:2  
            if unit.fluid[h,i] == 1 
                n_rohr = 47 # number of pipes TODO
                
                # calculate specific heat extraction for 1 pipe per length
                specific_heat_flux_pipe = (-1* unit.q_in_out_nov[unit.timestep_index]*1000/ (unit.pipe_length*n_rohr))
                unit.fluid_temperature = unit.average_temperature_adjacent_to_pipe + pipe_thermal_resistance_length_specific * specific_heat_flux_pipe     # "+", because specific_heat_flux_pipe is negative during heat extraction.        
                    
                unit.t2[h,i] = unit.fluid_temperature
                
                Volume_adjacent_to_pipe = ((unit.dx[i]+unit.dx[i+1]/2) *                               # total x-Direction
                                          (unit.dy[h-2]/2+unit.dy[h-1]+unit.dy[h]+unit.dy[h+1]/2) -    # total y-Direction
                                          (pi*unit.dy[h]^2/2))*unit.dz                                 # 1/2 Area pipe.          
                
                # TODO: Connection to ReSiE value from process function. Currently, an external input file is read in. 
                q_in_out_surrounding =(-1) * unit.q_in_out_nov[unit.timestep_index]*1000 / (n_rohr)
                
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
    n_it = Int(unit.wh_to_w_timestep/unit.dt   * 3600) 
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
function freezing(unit::GeothermalHeatCollector,T2_in, phase_change_state_in)
    cp = unit.soil_specific_heat_capacity

    if T2_in < unit.phase_change_upper_boundary_temperature && phase_change_state_in == 0
        T2_out = unit.phase_change_upper_boundary_temperature
        phase_change_state_out = 1
    
    elseif T2_in > unit.phase_change_upper_boundary_temperature && phase_change_state_in==1
        T2_out = unit.phase_change_upper_boundary_temperature
        phase_change_state_out = 0
    
    elseif T2_in < unit.phase_change_lower_boundary_temperature && phase_change_state_in==1
        T2_out = unit.phase_change_lower_boundary_temperature
        phase_change_state_out = 2
    
    elseif T2_in > unit.phase_change_lower_boundary_temperature && phase_change_state_in==2
        T2_out = unit.phase_change_lower_boundary_temperature
        phase_change_state_out = 1
    
    else 
        T2_out = T2_in
        phase_change_state_out = phase_change_state_in
          
    end
    
    if T2_in > unit.phase_change_upper_boundary_temperature
        cp = unit.soil_specific_heat_capacity
    elseif T2_in<= unit.phase_change_upper_boundary_temperature && T2_in >= unit.phase_change_lower_boundary_temperature
        dt_lat = unit.phase_change_upper_boundary_temperature - unit.phase_change_lower_boundary_temperature
        sigma_lat = 1/5 * dt_lat
        t_lat = (unit.phase_change_upper_boundary_temperature + unit.phase_change_lower_boundary_temperature)/2
        cp = unit.soil_heat_fusion_energy * 1/(sigma_lat* sqrt(2*pi)) * exp(-0.5 * (T2_in-t_lat)^2/sigma_lat^2)
    
    elseif T2_in<unit.phase_change_lower_boundary_temperature
        cp = unit.soil_specific_heat_capacity - 100
    
    
    end

    return (T2_out, phase_change_state_out, cp)

end 

# function to calculate heat transfer coefficient alpha.
function calculate_alpha_pipe(unit::GeothermalHeatCollector)
    # calculate reynolds-number to choose right Nusselt-calculation method
    n_rohr = 47 # TODO
    d_i_pipe = 2*unit.pipe_radius_outer - (2*unit.pipe_thickness)

    # calculate temeprature-dependant Prantl-Number of the heat carrier fluid.
    fluid_dynamic_viscosity = 0.0000017158* unit.fluid_temperature^2 - 0.0001579079*unit.fluid_temperature+0.0048830621
    unit.fluid_heat_conductivity = 0.0010214286 * unit.fluid_temperature + 0.447
    unit.fluid_prantl_number = fluid_dynamic_viscosity * unit.fluid_specific_heat_capacity / unit.fluid_heat_conductivity 

    # only for validation purpose! 
    unit.collector_total_heat_flux_in_out = unit.q_in_out_nov[unit.timestep_index]*1000 / n_rohr

    # calculate reynolds-number.
    unit.fluid_reynolds_number = (unit.collector_total_heat_flux_in_out) /
                                 (unit.fluid_specific_heat_capacity * unit.unloading_temperature_spread * 
                                 pi / 4 * d_i_pipe * fluid_dynamic_viscosity)

    fluid_reynolds_number = copy(unit.fluid_reynolds_number)
    # check for laminar flow
    if fluid_reynolds_number <= 2300
        Nu = calculate_Nu_laminar(unit::GeothermalHeatCollector, d_i_pipe, fluid_reynolds_number)
        # check for transition flow

    elseif fluid_reynolds_number > 2300
        # Gielinski
        factor = (fluid_reynolds_number - 2300)/(1e4-2300)
        Nu = (1- factor)*
        calculate_Nu_laminar(unit::GeothermalHeatCollector, d_i_pipe, 2300) +
        factor * calculate_Nu_turbulent(unit::GeothermalHeatCollector, d_i_pipe, 1e4)

    end

    alpha = Nu * unit.fluid_heat_conductivity / d_i_pipe

    return alpha
end

function calculate_Nu_laminar(unit::GeothermalHeatCollector, d_i_pipe, fluid_reynolds_number)
    # # Stephan
    # Pr_water = 13.44                # Pr Number Water 0 °C
    # Nu = 3.66 + (0.0677 * (fluid_reynolds_number * unit.fluid_prantl_number * d_i_pipe/unit.pipe_length)^1.33) /
    #     (1+0.1* unit.fluid_prantl_number * (fluid_reynolds_number * d_i_pipe/unit.pipe_length)^0.83) *
    #     (unit.fluid_prantl_number/Pr_water)^ 0.1

    d_i_pipe = 2*unit.pipe_radius_outer - (2*unit.pipe_thickness)
    
    # Approach used in Ramming 2007 from Elsner, Norbert; Fischer, Siegfried; Huhn, Jörg; „Grundlagen der Technischen Thermodynamik“,  Band 2 Wärmeübertragung, Akademie Verlag, Berlin 1993. 
    k_a = 0.0         # initializing k_a
    k_n = 0.0        # initializing k_n
    # substitutions: 
    k_a = 1.1 - 1 / (3.4+0.0667*unit.fluid_prantl_number)
    k_n = 0.35 + 1 / (7.825 + 2.6 * sqrt(unit.fluid_prantl_number))

    # calculate Nu-Number
    Nu = ((k_a/(1-k_n)*(unit.fluid_prantl_number*d_i_pipe*fluid_reynolds_number/unit.pipe_length)^k_n)^3+4.364^3)^(1/3)

    return Nu
end 

function calculate_Nu_turbulent(unit::GeothermalHeatCollector, d_i_pipe, fluid_reynolds_number)
    zeta = (1.8*log(fluid_reynolds_number)-1.5)^-2
    # Gielinski
    Nu = (zeta/8 * fluid_reynolds_number * unit.fluid_prantl_number)/
        (1 + 12.7 * sqrt(zeta/8) * (unit.fluid_prantl_number^(2/3)-1))
    return Nu
end

# process function that provides energy from the geothermal heat collector and calculates new temperatures 
# according to actual delivered or received energy
function process(unit::GeothermalHeatCollector, parameters::Dict{String,Any})
    # get actual required energy from output interface
    outface = unit.output_interfaces[unit.m_heat_out]  # output interface
    exchange = balance_on(outface, outface.target)     # gather information of output interface
    demand_temp = exchange.temperature                 # get temperature requested by demand (equals max(unit.temperature_field) 
            	                                       # from control-step if demand is not requesting a temperature
    # check if temperature can be met
    if demand_temp !== nothing && demand_temp > unit.current_output_temperature
        # no calculate_new_heat_collector_temperatures() as this will be done in load() step to avoid double calling!
        return  # no energy delivery possible as requested temperature can not be provided
    end

    # calculate energy demand with respect to the defined control stratey
    if unit.controller.parameter["name"] == "default"
        energy_demand = exchange.balance
    elseif unit.controller.parameter["name"] == "extended_storage_control"
        if unit.controller.parameter["load_any_storage"]
            energy_demand = exchange.balance + exchange.storage_potential
        else
            energy_demand = exchange.balance
        end
    else
        energy_demand = exchange.balance
    end

    if energy_demand >= 0.0
        # no calculate_new_heat_collector_temperatures() as this will be done in load() step to avoid double calling
        return # process is only concerned with moving energy to the target
    end

    energy_demand = min(unit.max_output_energy, energy_demand)
    unit.collector_total_heat_flux_in_out = energy_demand/(unit.wh_to_w_timestep)    
    # calcute energy that acutally can be delivered and set it to the output interface 
    # no other limits are present as max_energy for geothermal heat collector was written in control-step
    add!(outface, abs(energy_demand), unit.current_output_temperature)
end

function load(unit::GeothermalHeatCollector, parameters::Dict{String,Any})   
     # get actual delivered energy from input interface
     inface = unit.input_interfaces[unit.m_heat_in]  # input interface
     exchange = balance_on(inface, inface.source)    # gather information of input interface
     supply_temp = exchange.temperature              # get temperature delivered by source 
     energy_available = exchange.balance             # get energy that is provided

    if  (!unit.regeneration ||                                                      # no energy available if regeneration is turned off
        energy_available <= 0.0 ||                                                 # no energy available for loading as load is only concerned when receiving energy from the target
        (supply_temp !== nothing && supply_temp < unit.current_input_temperature)  # we can only take in energy if it's at a higher temperature than the heat collector lowest temperature
    )
        # recalculate heat collector temperatures for next timestep (function checks if is has already been calculated in the current timestep)
        calculate_new_temperature_field(unit::GeothermalHeatCollector)  # call calculate_new_heat_collector_temperatures() to calculate new temperatures of field to account for possible ambient effects
        return
    end
    
    energy_available = min(unit.max_input_energy, energy_available)
    unit.collector_total_heat_flux_in_out += energy_available/unit.wh_to_w_timestep 

    # calcute energy that acutally has beed delivered for regeneration and set it to interface 
    # no other limits are present as max_energy for geothermal heat collector was written in control-step!
    sub!(inface, energy_available, unit.current_input_temperature)
    calculate_new_temperature_field(unit::GeothermalHeatCollector) 
end

function balance_on(
    interface::SystemInterface,
    unit::GeothermalHeatCollector
)::NamedTuple{}
    # check if interface is input or output on unit
    input_sign = unit.uac == interface.target.uac ? -1 : +1
    # check if a balance was already written --> if yes, storage potential will be set to zero as storage was already processed/loaded
    balance_written = interface.max_energy === nothing || interface.sum_abs_change > 0.0

    return (
            balance = interface.balance,
            storage_potential = balance_written ? 0.0 : input_sign * interface.max_energy,  # geothermal heat collector are handled as storages currently!
            energy_potential = 0.0,
            temperature = interface.temperature
            )
end

function output_values(unit::GeothermalHeatCollector)::Vector{String}
    return ["IN", "OUT", "TEMPERATURE_#NodeNum","fluid_temperature"]
end

function output_value(unit::GeothermalHeatCollector, key::OutputKey)::Float64
    if key.value_key == "IN"
        return calculate_energy_flow(unit.input_interfaces[key.medium])
    elseif key.value_key == "OUT"
        return calculate_energy_flow(unit.output_interfaces[key.medium])
    elseif startswith(key.value_key, "Temperature_")
        idx = parse(Int, split(key.value_key, "_")[2])
        if !(1 <= idx <= length(unit.dy)+1) # unit.T2!
            throw(ArgumentError("Index \"$idx\" of requested temperature-output of geothermal heat collector exeeds the number of available temperatur datapoints.")) 
        else
            return unit.t2[idx,1] # unit.T2! # TODO: Every Node
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

