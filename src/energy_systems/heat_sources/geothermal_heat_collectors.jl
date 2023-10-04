"""
Implementation of geothermal heat collector.
This implementations acts as storage as is can produce and load energy.
"""

# read in heat flux profile from measurement data (only for validation purposes)
file = open("C:/Users/vollmer/Documents/GitHub/resie_FA_AdVo/src/energy_systems/heat_sources/q_in_out_nov_h.txt", "r")
heat_flux_profile = Vector{Float64}()
for line in eachline(file)
    push!(heat_flux_profile, parse(Float64, line))
end
close(file)

# read in discretitaion vector (only for validation purposes)
file = open("C:/Users/vollmer/Documents/GitHub/resie_FA_AdVo/src/energy_systems/heat_sources/dy_trnsys.txt", "r")
dy_type710 = Vector{Float64}()

for line in eachline(file)
    push!(dy_type710, parse(Float64, line))
end
close(file)

# read in inlet temperature profile from measurement data (only for validation purposes)
file = open("C:/Users/vollmer/Documents/GitHub/resie_FA_AdVo/src/energy_systems/heat_sources/t_in_validation.txt", "r")
file_t_in_validation = Vector{Float64}()
for line in eachline(file)
    push!(file_t_in_validation, parse(Float64, line))
end
close(file)

# read in massflow profile from measurement data (only for validation purposes)
file = open("C:/Users/vollmer/Documents/GitHub/resie_FA_AdVo/src/energy_systems/heat_sources/m_in_validation.txt", "r")
file_m_in_validation = Vector{Float64}()
for line in eachline(file)
    push!(file_m_in_validation, parse(Float64, line))
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
    T1::Array
    T2::Array
    phase_change_state::Array
    CP1::Array
    CP2::Array
    Fluid::Array
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
    Re::Float64
    q_in_out_nov::Vector
    discretization_mode::Int
    validation_mode::Int
    t_in_validation::Vector
    t_out_validation::Vector
    m_in_validation::Vector
    p_out_validation::Vector
    t_avg_rand::Float64

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
            config["max_output_power"],  # maximum output power set by user, may change later to be calculated from other inputs like specific heat transfer rate
            config["max_input_power"],   # maximum input power set by user, may change later to be calculated from other inputs like specific heat transfer rate
            default(config, "regeneration", true),   # flag if regeneration should be taken into account
            0.0,                         # max_output_energy in every time step, calculated in control()
            0.0,                         # max_input_energy in every time step, calculated in control()
            0.0,                         # output temperature in current time step, calculated in control()
            0.0,                         # input temperature in current time step, calculated in control()
            0.0,                         # ambient temperature in current time step, calculated in control()
            -1.0,                       # last timestep that was calculated; used to avoid double calculation of temperature field. set to -1 for the beginning
            15.5,                       # starting temperature near pipe. Value currently set for validation purpose. 
            default(config, "soil_specific_heat_capacity", 1000), 
            default(config, "soil_density", 2000),                 # specific heat capacity soil, 840
            default(config, "soil_heat_conductivity", 1.5), # soil heat conductivity, 
            zeros(59),              # vector to set soil density. Size depends on discretitaion
            zeros(59),              # lamda soil. TO DO: Use Soil propertys of validation project! 1.9
            default(config, "soil_heat_fusion_energy", 90000),                  # Heat fusion Energy in J/kg. Currently out of TRNSYS Default value.
            -0.25,                  # phase_change_upper_boundary_temperature
            -1,                     # phase_change_lower_boundary_temperature
            14.7,                   # convective heat transfer on surface - TO DO: Search for literature source.
            0.25,                   # Reflection factor / Albedo value of surface. 
            zeros(31,7),            # T1
            zeros(31,7),            # T2 grob:zeros(31,7) fein: zeros(59,11)
            zeros(31,7),            # phase_change_state
            zeros(31,7),            # array to define specific heat capacity on previous timestep (for apparent heat capacity method)
            zeros(31,7),            # array to define specific heat capacity on next timestep (for apparent heat capacity method)
            zeros(31,7),            # Fluid 
            zeros(31,7),            # pipe_surrounding
            default(config, "pipe_radius_outer", 0.016),                   # pipe outer radius
            default(config, "pipe_thickness", 0.003),                 # thickness of pipe in m
            default(config, "pipe_laying_depth", 1.5),                      # deepth of pipe system below the ground in m
            default(config, "pipe_length", 100),                   # pipe length of one collector in m
            default(config, "pipe_spacing", 0.5),                   # distance between pipes of collector in m.
            0,                      # global radiation on surface. to be read in by weather-profile
            5.67e-8,                # Boltzmann-Constant
            0.9,                    # Emissivity on ground surface
            zeros(6),               # size dx. depends on discretization settings.
            zeros(30),              # size dy. depends on discretization settings.
            93.5,                   # dz in m, is constant. Fluid-Direction.
            20,                     # duration time of internal time-step dt in s depending on ground properties 
            0,                      # time step index. necessary for validation purposes.                     
            10,                     # set starting fluid temperature.               
            0,                      # pipe temperature. 
            0,                      # total heat flux in or out of collector. set by ReSiE Interface.
            0.5,                    # pipe heat conductivity.               
            default(config, "fluid_specific_heat_capacity", 3800),  # fluid_specific_heat_capacity in J/(kg K)
            default(config, "fluid_prantl_number", 30),             # prandtl number at 30 % glycol, 0 °C 
            default(config, "fluid_density", 1045),                 # fluid density at 30 % glycol, 0 °C
            default(config, "fluid_kinematic_viscosity", 3.9e-6),     # fluid_kinematic_viscosity at 30 % glycol, 0 °C
            default(config, "fluid_heat_conductivity",0.5),      # fluid_heat_conductivity at 30 % glycol, 0 °C
            default(config, "specific_heat_flux_in_out",20),          # max. specific heat flux extractable out of soil in W/m^2. Depending on ground and climate localization. [VDI 4640-2.]
            1,           # time step of simulation in h to calculate power (W) out of energy (Wh) (wh_to_w_timestep), depending on simulation time step.
            0,           # Re-Number, to be calculated 
            heat_flux_profile,    # for validation purpose
            0,      # discretization-mode. 0: Default. 1: read in dy-vector from Type 710 for validation purpose
            0,      # 1: Validation-Mode (Input: t_in, m_in, output: p_out, t_out), 0: Resi-Mode (input: p_out, output: t_avg)
            file_t_in_validation,   # inlet temperature vector. only for validation purpose
            zeros(721),     # outlet temperature vector. only for validation purpose
            file_m_in_validation, # massflow vector. only for validation purpose
            zeros(721), # heat extraction vector. only for validation purpose
            16.0    # starting temperature of pipe. set for validation purpose.
            )
    end
end

function control(
    unit::GeothermalHeatCollector,
    components::Grouping,
    parameters::Dict{String,Any}
)

    unit.timestep_index += 1

    # Discretization and starting Temperature-Field.
    if unit.timestep_index == 1
        if unit.discretization_mode ==0
        # Discretization
        dx_R = [unit.pipe_radius_outer]
        dx_RM = [unit.pipe_radius_outer, unit.pipe_radius_outer, 2*unit.pipe_radius_outer, 4*unit.pipe_radius_outer]
        dx_End = unit.pipe_spacing/2 - sum(dx_R) - sum(dx_RM)
        unit.dx = [dx_R..., dx_RM..., dx_End...] # append()
            
        dy_O = [unit.pipe_radius_outer/2, unit.pipe_radius_outer, 4*unit.pipe_radius_outer, 8*unit.pipe_radius_outer, 16*unit.pipe_radius_outer]   # dy near surface
        dy_MR = [16*unit.pipe_radius_outer, 8*unit.pipe_radius_outer, 4*unit.pipe_radius_outer, 2*unit.pipe_radius_outer, unit.pipe_radius_outer/2]  # dy getting smaller untill pipe node
        dy_M = [(unit.pipe_laying_depth - (sum(dy_O)+sum(dy_MR)) - unit.pipe_radius_outer)]                       # dy is wide in the middle    
        dy_R = [unit.pipe_radius_outer, unit.pipe_radius_outer]                                                   # distance between nodes adjacent to fluid
        dy_RU = [unit.pipe_radius_outer/2, 1*unit.pipe_radius_outer, 1*unit.pipe_radius_outer, 3/2*unit.pipe_radius_outer, 2*unit.pipe_radius_outer,2*unit.pipe_radius_outer, 4*unit.pipe_radius_outer,16*unit.pipe_radius_outer,32*unit.pipe_radius_outer, 64*unit.pipe_radius_outer, 64*unit.pipe_radius_outer, 64*unit.pipe_radius_outer, 128*unit.pipe_radius_outer, 64*unit.pipe_radius_outer, 32*unit.pipe_radius_outer, 16*unit.pipe_radius_outer, 8*unit.pipe_radius_outer]  # mesh below the fluid-node until lower simulation boundary   
        unit.dy = [dy_O..., dy_M..., dy_MR..., dy_R..., dy_RU...]   #s ize:27                         # build dy-vector
        
        unit.dz = copy(unit.pipe_length)
           
        # Localize Fluid and adjacent Nodes
        unit.Fluid[(length(dy_O)+ length(dy_MR)+ length(dy_M)+length(dy_R)+1 -1),1] = 1     
        unit.pipe_surrounding[(length(dy_O)+ length(dy_MR)+ length(dy_M)+length(dy_R)+1) ,1] = 1
        unit.pipe_surrounding[(length(dy_O)+ length(dy_MR)+ length(dy_M)+length(dy_R)+1 -2) , 1] = 1
        unit.pipe_surrounding[(length(dy_O)+ length(dy_MR)+ length(dy_M)+length(dy_R)+1 -1) , 2] = 1
        
        # set starting temperature distribution (current settings for validation purpose)
        unit.T1[1:3,:] .= 9.5
        unit.T1[4:6,:] .= 11.5
        unit.T1[7:9,:] .= 13.5
        unit.T1[10:13,:] .= 15.5 
        unit.T1[14:30,:] .= unit.soil_starting_temperature   
        unit.T1[14:30,:] .= unit.soil_starting_temperature     
        unit.T1[length(unit.dy)+1, :] .= 9    # TO DO: Use avg temperature out of weather data set.
        unit.T2 = copy(unit.T1)


        # set soil heat conductivity, more layers possible.
        unit.soil_heat_conductivity_vector[:] .= unit.soil_heat_conductivity
        
        # set soil density, more layers possible.
        unit.soil_density_vector[:] .= unit.soil_density
        
        # only for validation purposes
        elseif unit.discretization_mode ==1 # read dy-vector from type 710 for validation purposes
            dx_R = [0.0052, 0.0052, 0.0052, 0.01, unit.pipe_radius_outer]
            dx_RM = [unit.pipe_radius_outer, unit.pipe_radius_outer, 2*unit.pipe_radius_outer, 4*unit.pipe_radius_outer]
            dx_End = unit.pipe_spacing/2 - sum(dx_R) - sum(dx_RM)
            unit.dx = [dx_R..., dx_RM..., dx_End...] # append()
            
            # set starting temperature distribution (current settings for validation purpose)
            unit.T1[1:7,:] .= 9.5
            unit.T1[8:13,:] .= 11.5
            unit.T1[14:20,:] .= 13.5
            unit.T1[21:38,:] .= unit.soil_starting_temperature 
            unit.T1[39:43,:] .= unit.soil_starting_temperature 
            unit.T1[44:48,:] .= 13.5 
            unit.T1[49:53,:] .= 11.5
            unit.T1[54:58,:] .= 9.5
            unit.T1[length(unit.dy)+1, :] .= 9    
            unit.T2 = copy(unit.T1)

            unit.Fluid[39,1] = 1     
            unit.pipe_surrounding[38,1] = 1
            unit.pipe_surrounding[40, 1] = 1
            unit.pipe_surrounding[39, 2] = 1

            unit.soil_heat_conductivity_vector[:] .= unit.soil_heat_conductivity
            unit.soil_density_vector[:] .= unit.soil_density
        end
        
        unit.CP1[:,:] .= unit.soil_specific_heat_capacity
        unit.CP2 = copy(unit.CP1)
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

    # for Validation purposes: Input t_in, m_in; Output: t_out, p_out.
    if unit.validation_mode == 1
        R_rohr = 0  # initatisation of variable.

        alpha_fluid = calculate_alpha_pipe(unit::GeothermalHeatCollector)   # To DO: Calculate out of mass-flow 
        # Ramming längenbezogener Wärmedurchgangswiderstand
        R_rohr = 1 / (2* pi) * (2 / (alpha_fluid * (2 * (unit.pipe_radius_outer - unit.pipe_thickness))) + 1 / unit.pipe_heat_conductivity * log(unit.pipe_radius_outer/(unit.pipe_radius_outer - unit.pipe_thickness)))
        
        # calculate fluid output temperature segment-wise 
        n_segments = 94
        t_in_neu = unit.t_in_validation[unit.timestep_index]
        t_out = unit.t_in_validation[unit.timestep_index]
        factor = unit.pipe_length/(n_segments * R_rohr*unit.m_in_validation[unit.timestep_index]*unit.fluid_specific_heat_capacity)
        
            for s=1:n_segments

                if s == n_segments
                    unit.t_out_validation[unit.timestep_index] = (factor * (unit.t_avg_rand-t_in_neu/2)+t_in_neu) * 1/(1+factor/2)
                    #unit.t_out_validation[unit.timestep_index] = unit.pipe_length/(n_segments*R_rohr*unit.m_in_validation[unit.timestep_index]*unit.fluid_specific_heat_capacity)*(unit.t_avg_rand-t_in_neu) + t_in_neu
                end
                t_out = (factor * (unit.t_avg_rand-t_in_neu/2)+t_in_neu) * 1/(1+factor/2)
                #t_out = unit.pipe_length/(n_segments*R_rohr*unit.m_in_validation[unit.timestep_index]*unit.fluid_specific_heat_capacity)*(unit.t_avg_rand-t_in_neu) + t_in_neu
                t_in_neu = copy(t_out)

            end
        
        n_rohr = 47 # number of pipes  
        # p_out (only heat extraction)
        if unit.t_in_validation[unit.timestep_index] >= unit.t_out_validation[unit.timestep_index]
            unit.p_out_validation[unit.timestep_index] = 0
            
        elseif unit.m_in_validation[unit.timestep_index] <= 0
            unit.p_out_validation[unit.timestep_index] = 0
        else
            unit.p_out_validation[unit.timestep_index] = unit.m_in_validation[unit.timestep_index]*unit.fluid_specific_heat_capacity*(unit.t_out_validation[unit.timestep_index]-unit.t_in_validation[unit.timestep_index])*n_rohr

        end

        for h = 1: (length(unit.dy)) 
            for i =1:2  
                if unit.Fluid[h,i] == 1 
                    n_rohr = 47 
                    unit.fluid_temperature = (unit.t_in_validation[unit.timestep_index]+unit.t_out_validation[unit.timestep_index])/2         # "+", because specific_heat_flux_pipe is negative during heat extraction.        
                    unit.T2[h,i] = unit.fluid_temperature    
            
                end    
            end
        end
    end



end

# function that calculates current (highest possible) output temperature of heat collector that can be provided. (lower temp. is always possible!)
function current_output_temperature(unit::GeothermalHeatCollector)::Temperature
    highest_outer_borehole_temp = unit.pipe_temperature - unit.unloading_temperature_spread/2
    # could be the highest outer borewhole temperature of the heat collector in the last timestep
    return highest_outer_borehole_temp
end

# function that calculates current (minimum possible) input temperature of heat collector that should be provided for regeneration 
# (higher temperature is always possible, here the minimum should be calculated!)
function current_input_temperature(unit::GeothermalHeatCollector)::Temperature
    input_temperature = unit.pipe_temperature + unit.loading_temperature_spread/2
    # could be the lowest outer borewhole temperature plot the given temperature spread or a user-specified loading_temperature if given
    return input_temperature
end

# function to calculate the maximum possible output energy in current time step
function get_max_output_power(unit::GeothermalHeatCollector)::Float64
    
    n_pipe = 47
    A_collector = unit.pipe_length * unit.pipe_spacing*(n_pipe  -1)
    unit.max_output_power = A_collector * unit.specific_heat_flux_in_out
    
    return unit.max_output_power  # could be calculated by maximum specific unloading energy or by maximal mass flow or...
end

# function to calculate the maximum possible input energy in current time step
function get_max_input_power(unit::GeothermalHeatCollector)::Float64
    n_pipe = 47
    A_collector = unit.pipe_length * unit.pipe_spacing*(n_pipe  -1)
    unit.max_input_power = A_collector * unit.specific_heat_flux_in_out    # for validation only
    return unit.max_input_power  # could be calculated by maximum specific loading energy or by maximal mass flow or...
end

function calculate_new_temperature_field(unit::GeothermalHeatCollector)
    if unit.validation_mode == 0
        # random output to avoid misunderstandings while reading output file.
        unit.p_out_validation[unit.timestep_index]=99 
    end
    
    # calculate heat transfer coefficient and thermal resistance
    alpha_fluid = calculate_alpha_pipe(unit::GeothermalHeatCollector)   
    R_rohr = 1 / (2* pi) * (2 / (alpha_fluid * (2 * (unit.pipe_radius_outer - unit.pipe_thickness))) + 1 / unit.pipe_heat_conductivity * log(unit.pipe_radius_outer/(unit.pipe_radius_outer - unit.pipe_thickness)))    
    
    # calculate fluid temperature
    if  unit.validation_mode == 0
        for h = 1: (length(unit.dy)) 
         
            for i =1:2  
                if unit.Fluid[h,i] == 1 
                    n_rohr = 47 # number of pipes
                    
                    # calculate specific heat extraction for 1 pipe per length
                    specific_heat_flux_pipe = (-1* unit.q_in_out_nov[unit.timestep_index]*1000/ (unit.pipe_length*n_rohr))
                    unit.t_avg_rand = (unit.T1[h-1,i] + unit.T1[h+1,i] + unit.T1[h,i+1]) / 3 #deactivated during validation
                    unit.fluid_temperature = unit.t_avg_rand + R_rohr * specific_heat_flux_pipe          # "+", because specific_heat_flux_pipe is negative during heat extraction.        
                        
                    unit.T2[h,i] = unit.fluid_temperature    
                
                        
                end    
            end
        end 
    end
 
    # calculate number of internal timesteps depending on internal time dt
    n_it = unit.wh_to_w_timestep/unit.dt   * 3600      
for it=1:n_it
n_rohr = 47

# Loop - vertical direction
for h = 1: (length(unit.dy))

# Loop - horizontal direction
for i = 1 : (length(unit.dx)+1)
    
    # calculate temperature of adjacent nodes to fluid
    if unit.Fluid[h,i] == 1                  

        # upper node
        unit.T2[h-1,i] = unit.T1[h-1,i] +      
        (unit.dz * 0.5 * (unit.dy[h-1]+unit.dy[h-2]) * unit.soil_heat_conductivity_vector[h] * (unit.T1[h-1,i+2] - unit.T1[h-1,i+1])/unit.dx[i] +    # WL in positve x-Richtung
        unit.dz * unit.dx[i] * unit.soil_heat_conductivity_vector[h] * (unit.T1[h-2, i]-unit.T1[h-1,i])/(unit.dy[h-2]) +                       # WL in negatove y-Richtung
        1/(R_rohr*unit.dz)* unit.dz * unit.dx[i] * (unit.fluid_temperature - unit.T1[h-1,i]))  * 
        unit.dt * 1/(unit.soil_density_vector[h] * unit.CP1[h,i] * unit.dz * 0.5 * unit.dx[i] * 0.5 * (unit.dy[h]+unit.dy[h-1]))                                                               # Quellterm
        
        # lower node
        unit.T2[h+1,i] = unit.T1[h+1,i] +    
        (unit.dz * 0.5 * (unit.dy[h]+unit.dy[h+1]) * unit.soil_heat_conductivity_vector[h] * (unit.T1[h+1, i+2] - unit.T1[h+1,i+1])/unit.dx[i] +     # WL in positve x-Richtung
        unit.dz * unit.dx[i] * unit.soil_heat_conductivity_vector[h] * (unit.T1[h+2,i]-unit.T1[h+1,i])/(unit.dy[h]) +                         # WL in positive y-Richtung
        1/(R_rohr*unit.dz)* unit.dz * unit.dx[i] * (unit.fluid_temperature - unit.T1[h+1,i]))*
        unit.dt * 1/(unit.soil_density_vector[h] * unit.CP1[h,i] * unit.dz * 0.5 * unit.dx[i] * 0.5 * (unit.dy[h]+unit.dy[h-1]))                                                        # Quellterm
            
        # right node
        dq_x_p = unit.dz * 0.5 * (unit.dy[h]+unit.dy[h+1]) * unit.soil_heat_conductivity_vector[h] * (unit.T1[h,i+2] - unit.T1[h,i+1])/unit.dx[i]
        dq_y_p = unit.dz * unit.dx[i] * unit.soil_heat_conductivity_vector[h] * (unit.T1[h+1,i+1]-unit.T1[h,i+1])/(unit.dy[h])
        dq_y_n = unit.dz * unit.dx[i] * unit.soil_heat_conductivity_vector[h] * (unit.T1[h-1,i+1]-unit.T1[h,i+1])/(unit.dy[h-1])

        unit.T2[h,i+1] = unit.T1[h,i+1] +    
        (dq_x_p +   # WL in positve x-Richtung
        dq_y_n +    # WL in negative y-Richtung
        dq_y_p +    # WL in positive y-Richtung
        1/(R_rohr*unit.dz)* unit.dz * unit.dy[i] * (unit.fluid_temperature - unit.T1[h,i+1]))* 
        unit.dt * 1/(unit.soil_density_vector[h] * unit.CP1[h,i] * unit.dz * unit.dx[i] * 0.5 * (unit.dy[h]+unit.dy[h-1]))                                                   # Quellterm

    # adjacent nodes to fluid node have uniform temperature
    unit.t_avg_rand = (unit.T2[h-1,i] + unit.T2[h+1,i] + unit.T2[h,i+1])/3
    
    unit.T2[h-1,i] = copy(unit.t_avg_rand)
    unit.T2[h+1,i] = copy(unit.t_avg_rand)
    unit.T2[h,i+1] = copy(unit.t_avg_rand)
    unit.pipe_temperature = copy(unit.t_avg_rand)

    unit.T2[h-1,i], unit.phase_change_state[h-1,i], unit.CP2[h-1,i] = freezing(unit::GeothermalHeatCollector, unit.T2[h-1,i],unit.phase_change_state[h-1,i]) 
    unit.T2[h+1,i], unit.phase_change_state[h+1,i], unit.CP2[h+1,i] = freezing(unit::GeothermalHeatCollector, unit.T2[h+1,i],unit.phase_change_state[h+1,i])
    unit.T2[h,i+1], unit.phase_change_state[h,i+1], unit.CP2[h,i+1] = freezing(unit::GeothermalHeatCollector, unit.T2[h,i+1],unit.phase_change_state[h,i+1])

    
    elseif unit.pipe_surrounding[h,i]==1           
        unit.T2[h,i] = unit.T2[h,i]          # Skip re-calculation. already calculated.      
    
    # ground surface boundary 
    elseif h==1
        
        # left boundary
        if i ==1
            unit.T2[h,i] = unit.T1[h,i] +
            (unit.dz * 0.5 * unit.dx[i] * unit.surface_convective_heat_transfer_coefficient * (unit.ambient_temperature - unit.T1[h,i]) +
            unit.dz * 0.5 * unit.dx[i] * (1-unit.surface_reflection_factor) * unit.global_radiation +
            unit.dz * 0.5 * unit.dx[i] * unit.surface_emissivity * unit.boltzmann_constant * ((unit.ambient_temperature+273.15)^4 - (unit.T1[h,i]+273.15)^4) +
            unit.dz * unit.dy[1] * unit.soil_heat_conductivity_vector[h] * (unit.T1[h,i+1] - unit.T1[h,i])/unit.dx[i] +
            unit.dz * unit.dx[i] * unit.soil_heat_conductivity_vector[h] * (unit.T1[h+1,i] - unit.T1[h,i])/unit.dy[h]) * 
            unit.dt * 1/(unit.soil_density_vector[h] * unit.CP1[h,i] * unit.dz * 0.5 * unit.dx[i] * unit.dy[h])
            
         
        unit.T2[h,i], unit.phase_change_state[h,i], unit.CP2[h,i] = freezing(unit::GeothermalHeatCollector,unit.T2[h,i],unit.phase_change_state[h,i])
        
        # right boundary
        elseif i == (length(unit.dx)+1)
        unit.T2[h,i] = unit.T1[h,i] +
             (unit.dz * 0.5 * unit.dx[i-1] * unit.surface_convective_heat_transfer_coefficient * (unit.ambient_temperature - unit.T1[h,i]) +                                  # Convection on surface
             unit.dz * 0.5 * unit.dx[i-1] * (1-unit.surface_reflection_factor)* unit.global_radiation +                                            # Radiation on surface
             unit.dz * 0.5 * unit.dx[i-1] * unit.surface_emissivity * unit.boltzmann_constant * ((unit.ambient_temperature+273.15)^4 - (unit.T1[h,i]+273.15)^4)+   # Radiation out of surface
             unit.dz * unit.dy[h] * unit.soil_heat_conductivity_vector[h] * (unit.T1[h,i-1] - unit.T1[h,i])/unit.dx[i-1]+ 
             unit.dz * unit.dx[i-1] * unit.soil_heat_conductivity_vector[h] * (unit.T1[h+1,i]-unit.T1[h,i])/(unit.dy[h]) )* 
             unit.dt * 1/(unit.soil_density_vector[h]*unit.CP1[h,i]*unit.dz*0.5*unit.dx[i-1]*unit.dy[1])
            

        unit.T2[h,i], unit.phase_change_state[h,i], unit.CP2[h,i]= freezing(unit::GeothermalHeatCollector, unit.T2[h,i],unit.phase_change_state[h,i])
            
        else
        unit.T2[h,i] = unit.T1[h,i] +
            (unit.dz * 0.5 * (unit.dx[i]+unit.dx[i-1]) * unit.surface_convective_heat_transfer_coefficient * (unit.ambient_temperature - unit.T1[h,i])+                                  # Convection on surface
            unit.dz * 0.5 * (unit.dx[i]+unit.dx[i-1]) * (1-unit.surface_reflection_factor)* unit.global_radiation+                                            # Radiation on surface
            unit.dz * 0.5 * (unit.dx[i]+unit.dx[i-1]) * unit.surface_emissivity * unit.boltzmann_constant * ((unit.ambient_temperature+273.15)^4 - (unit.T1[h,i]+273.15)^4)+   # Radiation out of surface
            unit.dz * unit.dy[h] * unit.soil_heat_conductivity_vector[h] * (unit.T1[h,i+1] - unit.T1[h,i])/unit.dx[i]+                                           # WL in positve x-Richtung
            unit.dz * unit.dy[h] * unit.soil_heat_conductivity_vector[h] * (unit.T1[h,i-1] - unit.T1[h,i])/unit.dx[i-1]+                                         # WL in negative x-Richtung
            unit.dz * 0.5 * (unit.dx[i]+unit.dx[i-1]) * 
            unit.soil_heat_conductivity_vector[h] * (unit.T1[h+1,i]-unit.T1[h,i])/(unit.dy[1])) * unit.dt * 1/(unit.soil_density_vector[h] * unit.CP1[h,i] * unit.dz * 0.5*(unit.dx[i]+unit.dx[i-1])*unit.dy[1])                          # WL in y-Richtung
            

        unit.T2[h,i], unit.phase_change_state[h,i], unit.CP2[h,i]= freezing(unit::GeothermalHeatCollector, unit.T2[h,i],unit.phase_change_state[h,i])
            
        end 
    else
        # left boundary
        if i == 1 && unit.pipe_surrounding[h,i] == 0
        unit.T2[h,i] = unit.T1[h,i] +
            (unit.dz * 0.5 * (unit.dy[h]+unit.dy[h-1]) * unit.soil_heat_conductivity_vector[h] * (unit.T1[h,i+1] - unit.T1[h,i])/unit.dx[i]+      # WL in positve x-Richtung
            unit.dz * 0.5* unit.dx[i] * unit.soil_heat_conductivity_vector[h] * (unit.T1[h+1,i]-unit.T1[h,i])/(unit.dy[h])+                  # WL in positive y-Richtung
            unit.dz * 0.5* unit.dx[i] * unit.soil_heat_conductivity_vector[h] * (unit.T1[h-1,i]-unit.T1[h,i])/(unit.dy[h-1])) * 
            unit.dt * 1/(unit.soil_density_vector[h] * unit.CP1[h,i] * unit.dz * 0.5 * unit.dx[i] * 0.5 * (unit.dy[h]+unit.dy[h-1]))                # WL in negatove y-Richtung
            
        
        unit.T2[h,i], unit.phase_change_state[h,i], unit.CP2[h,i]= freezing(unit::GeothermalHeatCollector, unit.T2[h,i],unit.phase_change_state[h,i])
            

        # right boundary 
        elseif i == length(unit.dx)+1
        unit.T2[h,i] = unit.T1[h,i] +
            (unit.dz * 0.5 * (unit.dy[h]+unit.dy[h-1]) * unit.soil_heat_conductivity_vector[h] * (unit.T1[h,i-1] - unit.T1[h,i])/unit.dx[i-1]+   # WL in negative x-Richtung
            unit.dz * 0.5*unit.dx[i-1] * unit.soil_heat_conductivity_vector[h] * (unit.T1[h+1,i]-unit.T1[h,i])/(unit.dy[h])+                 # WL in positive y-Richtung
            unit.dz * 0.5*unit.dx[i-1] * unit.soil_heat_conductivity_vector[h] * (unit.T1[h-1,i]-unit.T1[h,i])/(unit.dy[h-1])) * 
            unit.dt * 1/(unit.soil_density_vector[h] * unit.CP1[h,i] * unit.dz * 0.5 * unit.dx[i-1]* 0.5 * (unit.dy[h]+unit.dy[h-1]))               # WL in negatove y-Richtung
        
        unit.T2[h,i], unit.phase_change_state[h,i], unit.CP2[h,i]= freezing(unit::GeothermalHeatCollector, unit.T2[h,i],unit.phase_change_state[h,i])
            
        else
        unit.T2[h,i] = unit.T1[h,i] +
            (unit.dz * 0.5 * (unit.dy[h]+unit.dy[h-1]) * unit.soil_heat_conductivity_vector[h] * (unit.T1[h,i+1] - unit.T1[h,i])/unit.dx[i]+        # WL in positve x-Richtung
            unit.dz * 0.5 * (unit.dy[h]+unit.dy[h-1]) * unit.soil_heat_conductivity_vector[h] * (unit.T1[h,i-1] - unit.T1[h,i])/unit.dx[i-1]+       # WL in negative x-Richtung
            unit.dz * 0.5 * (unit.dx[i]+unit.dx[i-1]) * unit.soil_heat_conductivity_vector[h] * (unit.T1[h+1,i]-unit.T1[h,i])/(unit.dy[h])+        # WL in positive y-Richtung
            unit.dz * 0.5 * (unit.dx[i]+unit.dx[i-1]) * unit.soil_heat_conductivity_vector[h] * (unit.T1[h-1,i]-unit.T1[h,i])/(unit.dy[h-1])) *  
            unit.dt * 1/(unit.soil_density_vector[h] * unit.CP1[h,i] * unit.dz * 0.5*(unit.dx[i]+unit.dx[i-1])* 0.5 * (unit.dy[h]+unit.dy[h-1]))       # WL in negatove y-Richtung
            
        unit.T2[h,i], unit.phase_change_state[h,i], unit.CP2[h,i]= freezing(unit::GeothermalHeatCollector, unit.T2[h,i],unit.phase_change_state[h,i])
            
        end 
    end


end 
end
unit.T1 = copy(unit.T2)
unit.CP1 = copy(unit.CP2)
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
n_rohr = 47
d_i_pipe = 2*unit.pipe_radius_outer - (2*unit.pipe_thickness)

if unit.validation_mode == 1
    unit.Re = (4 * unit.m_in_validation[unit.timestep_index])/
    (pi  * d_i_pipe * unit.fluid_kinematic_viscosity * unit.fluid_density)
else 

# only for validation purpose! 
unit.collector_total_heat_flux_in_out = unit.q_in_out_nov[unit.timestep_index]*1000

# calculate reynolds number
unit.Re = (unit.collector_total_heat_flux_in_out/n_rohr)/
        (unit.fluid_specific_heat_capacity * unit.unloading_temperature_spread * pi / 4 * d_i_pipe * unit.fluid_kinematic_viscosity * unit.fluid_density)
end
    

Re = copy(unit.Re)


    # check for laminar flow
    if Re <= 2300
        Nu = calculate_Nu_laminar(unit::GeothermalHeatCollector, d_i_pipe, Re)
        # check for transition flow

    elseif Re > 2300
        # Gielinski
        factor = (Re - 2300)/(1e4-2300)
        Nu = (1- factor)*
        calculate_Nu_laminar(unit::GeothermalHeatCollector, d_i_pipe, 2300) +
        factor * calculate_Nu_turbulent(unit::GeothermalHeatCollector, d_i_pipe, 1e4)

    end

    alpha = Nu * unit.fluid_heat_conductivity / d_i_pipe

    return alpha

end


function calculate_Nu_laminar(unit::GeothermalHeatCollector, d_i_pipe, Re)
    # Stephan
    Pr_water = 13.44                # Pr Number Water 0 °C
    Nu = 3.66 + (0.0677 * (Re * unit.fluid_prantl_number * d_i_pipe/unit.pipe_length)^1.33) /
        (1+0.1* unit.fluid_prantl_number * (Re * d_i_pipe/unit.pipe_length)^0.83) *
        (unit.fluid_prantl_number/Pr_water)^ 0.1

    return Nu
end 


function calculate_Nu_turbulent(unit::GeothermalHeatCollector, d_i_pipe, Re)
    zeta = (1.8*log(Re)-1.5)^-2

    # Gielinski
    Nu = (zeta/8 * Re * unit.fluid_prantl_number)/
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

    if !unit.regeneration
        # recalculate heat collector temperatures for next timestep (function checks if is has already been calculated in the current timestep)
        calculate_new_temperature_field(unit::GeothermalHeatCollector)  # call calculate_new_heat_collector_temperatures() to calculate new temperatures of field to account for possible ambient effects
        return
    end

    # no energy available for loading as load is only concerned when receiving energy from the target
    if energy_available <= 0.0
        # recalculate heat collector temperatures for next timestep (function checks if is has already been calculated in the current timestep)
        calculate_new_temperature_field(unit::GeothermalHeatCollector)  # call calculate_new_heat_collector_temperatures() to calculate new temperatures of field to account for possible ambient effects
        return
    end

    # we can only take in energy if it's at a higher temperature than the heat collector lowest temperature
    if supply_temp !== nothing && supply_temp < unit.current_input_temperature
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
    return ["IN", "OUT", "TEMPERATURE_#NodeNum","fluid_temperature","t_out","p_out"]
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
            return unit.T2[idx,1] # unit.T2! # TO DO: Every Node
        end
    elseif key.value_key =="fluid_temperature"
        return unit.fluid_temperature
    elseif key.value_key =="t_out"
        return unit.t_out_validation[unit.timestep_index]
    elseif key.value_key =="p_out"
        return unit.p_out_validation[unit.timestep_index]
    elseif key.value_key =="Re"
        return unit.Re

    elseif key.value_key =="ambient_temperature"
    return unit.ambient_temperature

    end
    throw(KeyError(key.value_key))
end

export GeothermalHeatCollector

