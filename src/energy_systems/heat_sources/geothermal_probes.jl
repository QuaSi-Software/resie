"""
Implementation of geothermal probes.
This implementations acts as storage as it can produce and load energy.
"""

# Current solution to get g-function values
# read g-function .txt file
file = open("C:/Users/vollmer/Documents/GitHub/resie_FA_AdVo/src/energy_systems/heat_sources/g_ges_25y_3600s.txt", "r")

    g_function_values = Vector{Float64}()

    for line in eachline(file)
        push!(g_function_values, parse(Float64, line))
    end
close(file)


# Read inlet temperatures .txt data-file (only for validation purposes)
file = open("C:/Users/vollmer/Documents/GitHub/resie_FA_AdVo/src/energy_systems/heat_sources/T_in_GEW_1h.txt", "r")

    probes_temperature_inlet = Vector{Float64}()

    for line in eachline(file)
        push!(probes_temperature_inlet, parse(Float64, line))
    end
close(file)


# Read massflow .txt data-file (only for validation purposes)
file = open("C:/Users/vollmer/Documents/GitHub/resie_FA_AdVo/src/energy_systems/heat_sources/m_in_GEW_1h.txt", "r")

    file_m_in = Vector{Float64}()

    for line in eachline(file)
        # Konvertieren der Zeichenkette in einen Float64-Wert und Hinzufügen zum Vektor
        push!(file_m_in, parse(Float64, line))
    end
close(file)


mutable struct GeothermalProbes <: ControlledComponent
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
    undisturbed_ground_temperature::Temperature
    soil_heat_conductivity::Float64
    thermal_resistance::Float64
    g_function::Vector
    time_index::Int
    average_fluid_temperature::Temperature
    current_borehole_wall_temperature::Temperature

    specific_heat_flux_in_out_absolut::Vector
    specific_heat_flux_in_out_step::Vector

    probe_depth::Float64
    number_of_probes::Float64

    diameter_pipe_outer::Float64
    diameter_pipe_inner::Float64

    fluid_heat_capacity::Float64
    fluid_density::Float64
    fluid_viscosity::Float64
    fluid_lambda::Float64
    fluid_prandtl::Float64

    grout_lambda::Float64
    pipe_lambda::Float64

    diameter_borehole::Float64
    shank_spacing::Float64

    Re::Float64

    diameter_pipe_eq_inner::Float64
    diameter_pipe_eq_outer::Float64

    t_in_validation::Vector
    m_in_validation::Vector
    t_out_validation::Vector
    q_in_out_validation::Vector
    q_out_abs_validation::Vector
    q_in_abs_validation::Vector   
    
    set_validation_mode::Int

    function GeothermalProbes(uac::String, config::Dict{String,Any})
        m_heat_in = Symbol(default(config, "m_heat_in", "m_h_w_ht1"))
        m_heat_out = Symbol(default(config, "m_heat_out", "m_h_w_lt1"))
        register_media([m_heat_in, m_heat_out])
    

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
            default(config, "unloading_temperature_spread", 3),   # temperature spread between forward and return flow during unloading            
            default(config, "loading_temperature", nothing),        # nominal high temperature for loading geothermal probe storage, can also be set from other end of interface
            default(config, "loading_temperature_spread", 3),     # temperature spread between forward and return flow during loading         
            config["max_output_power"],  # maximum output power set by user, may change later to be calculated from other inputs like specific heat transfer rate
            config["max_input_power"],   # maximum input power set by user, may change later to be calculated from other inputs like specific heat transfer rate
            default(config, "regeneration", true),   # flag if regeneration should be taken into account
            0.0,                         # max_output_energy in every time step, calculated in control()
            0.0,                         # max_input_energy in every time step, calculated in control()
            0.0,                         # output temperature in current time step, calculated in control()
            0.0,                         # input temperature in current time step, calculated in control()
            default(config, "undisturbed_ground_temperature", 11.0),    # Considered as constant
            default(config, "soil_heat_conductivity", 1.5),        # Heat conductivity of surrounding soil, homogenous and constant
            default(config, "thermal_resistance", 0.10), # thermal resistance in (m K)/W
            g_function_values,    # pre-calculated multiscale g-function
            0,              # index of current time step to get access on time dependent g-function values
            0.0,            # average fluid temperature
            4,              # set boreholewall-starting-temperature
            zeros(219000), # vector to write specific heat flux in eacht time step
            zeros(219000), # vector to write specific heat flux differeces in eacht time step for g-function approach
            default(config, "probe_depth", 150),    # depth (or length) of a single geothermal probe
            36,                         # number of geothermal probes in the borefield

            default(config, "diameter_pipe_outer", 0.032),  # outer pipe diameter
            default(config, "diameter_pipe_inner", 0.026),  # inner pipe diameter
           
            default(config, "fluid_heat_capacity", 3800),   # specific heat capacity brine at 0 °C (25 % glycol 75 % water (interpolated)) 
            default(config, "fluid_density", 1045),  # density brine at 0 °C (25 % glycol 75 % water (interpolated))
            default(config, "fluid_viscosity", 3.9e-6), # viscosity brine at 0 °C (25 % glycol 75 % water (interpolated)) 
            default(config, "fluid_lambda", 0.5) ,  # heat conductivity brine at 0 °C (25 % glycol 75 % water (interpolated))
            default(config, "fluid_prandtl", 30),   # prandtl-number brine at 0 °C (25 % glycol 75 % water (interpolated)) 
            
            default(config, "grout_lambda", 2),      # lambda grout / filling material in W/(mK)   
            default(config, "pipe_lambda", 0.42),   # lambda of inner pipes
            default(config, "diameter_borehole", 0.15),    # borehole diameter in m.
            0.1,     # shank-spacing = distance between inner pipes in borehole.

            0,      # Reynoldsnumber. To be calculated in Function later.
            0,      # eq diameter inner
            0,      # eq diameter outer

            probes_temperature_inlet,  # read in as input parameter for validation
            file_m_in,   # read in as input parameter for validation
            zeros(8760),
            zeros(8760),
            zeros(8760),
            zeros(8760),

            1       # set validation mode. if 1: input inlet Temperature and Massflow Profile. if 0: Input Demand and Regeneration profiles.
            )
    end
end

function control(
    unit::GeothermalProbes,
    components::Grouping,
    parameters::Dict{String,Any}
)
    # time index, necessary for g-function approach
    unit.time_index = unit.time_index + 1 

    # get input temperature for energy input (regeneration) and set temperature to input interface
    if unit.regeneration
        unit.current_input_temperature = current_input_temperature(unit) # of geothermal probe field 
        unit.input_interfaces[unit.m_heat_in].temperature = highest_temperature(
                                                                                unit.current_input_temperature,
                                                                                unit.input_interfaces[unit.m_heat_in].temperature
                                                                                )
    end                                                                        

    # calculate maximum input and output energy that is possible in current time step
        # sets max_energy to zero if requested/available temperature does not fit to temperature of geothermal probe field.
        # This works as the control step of transformers is always calculated earlier than the one of storages. If temperatures
        # are written to the connected interface by a transformer, this is already done at this point.
    if unit.output_interfaces[unit.m_heat_out].temperature > unit.current_output_temperature
        unit.max_output_energy = 0.0  # no energy can be provided if requested temperature is higher than max. temperature of probe field
    else
        unit.max_output_energy = watt_to_wh(get_max_output_power(unit))  
    end

    if unit.regeneration
        if unit.input_interfaces[unit.m_heat_in].temperature < unit.current_input_temperature
            unit.max_input_energy = 0.0 # no energy can be taken if available temperature is less than minimum possible temperature to load the probe field
        else
            unit.max_input_energy = watt_to_wh(get_max_input_power(unit)) 
        end
    end

    # set max_energy to interfaces to provide information for connected components
    set_max_energy!(unit.output_interfaces[unit.m_heat_out], unit.max_output_energy)
    if unit.regeneration
        set_max_energy!(unit.input_interfaces[unit.m_heat_in], unit.max_input_energy)
    end

    # Only for Validation purposes
    if unit.set_validation_mode == 1
        
        # calculate thermal resistance with input massflow
            radius_pipe_inner = unit.diameter_pipe_inner  / 2
            radius_pipe_outer = unit.diameter_pipe_outer / 2
            radius_borehole = unit.diameter_borehole / 2
            
            
            # set factor variable for eq. radius calculation (depending on user settings)
            # R_B with Hellström - factor = 1; R_B with radius_eq = 4*r_pipe - factor = 4; R_B with radius_eq = 2*r_pipe  - factor = 2
            factor = 1
        
        # calculate equivalent radius (depending on factor-calculation) for thermal borehole resistance calculation
        radius_pipe_eq_inner = factor*radius_pipe_inner
        radius_pipe_eq_outer = factor*radius_pipe_outer
        unit.diameter_pipe_eq_inner = 2*radius_pipe_eq_inner
        unit.diameter_pipe_eq_outer = 2*radius_pipe_eq_outer
        alpha_fluid = calculate_alpha_pipe(unit::GeothermalProbes)   
        
        # Calculate thermal borehole resistance (Hellström)
        sigma = (unit.grout_lambda-unit.soil_heat_conductivity)/(unit.grout_lambda+unit.soil_heat_conductivity)   # dimensionless calculation factor
        distance_pipe_center = unit.shank_spacing / 2
        beta = 1/(2*pi*alpha_fluid*radius_pipe_inner) + 1/(2*pi * unit.pipe_lambda) * log(radius_pipe_outer/radius_pipe_inner) # in (mK)/W

        R_1 = beta + 1/(2*pi*unit.grout_lambda)*
        (log(radius_borehole^2/(2*radius_pipe_outer*distance_pipe_center))+
        sigma*log(radius_borehole^4/(radius_borehole^4-distance_pipe_center^4))-
        radius_pipe_outer^2/(4*distance_pipe_center^2)*(1-sigma*4*distance_pipe_center^4/(radius_borehole^4-distance_pipe_center^4))^2 /
        ((1+2*pi*unit.grout_lambda*beta)/(1-2*pi*unit.grout_lambda*beta)+
        radius_pipe_outer^2/(4*distance_pipe_center^2)*(1+sigma*16*radius_borehole^4*distance_pipe_center^4/((radius_borehole^4-distance_pipe_center^4)^2))
        )
        )

        unit.thermal_resistance = R_1/4

        # calculate output temperature depending on massflow and inlet temperature
        # divide probe into segments.
        n_segments = 40     # 10 Segments (l_segment = 15 m ) down and up the probe.
        t_in_neu = copy(unit.t_in_validation[unit.time_index])
        t_out = copy(unit.t_in_validation[unit.time_index])
        factor = unit.probe_depth/(n_segments*unit.thermal_resistance*unit.m_in_validation[unit.time_index]*unit.fluid_heat_capacity)

        if unit.m_in_validation[unit.time_index] == 0
            unit.t_out_validation[unit.time_index]=unit.t_out_validation[unit.time_index-1]
        else
            for i=1:n_segments
                if i == n_segments
                    #unit.t_out_validation[unit.time_index] = unit.probe_depth/(n_segments*unit.thermal_resistance*unit.m_in_validation[unit.time_index]*unit.fluid_heat_capacity)*(unit.current_borehole_wall_temperature-t_in_neu) + t_in_neu
                    unit.t_out_validation[unit.time_index] = (factor * (unit.current_borehole_wall_temperature-t_in_neu/2)+t_in_neu) * 1/(1+factor/2)
                else
                
                # t_out = unit.probe_depth/(n_segments*unit.thermal_resistance*unit.m_in_validation[unit.time_index]*unit.fluid_heat_capacity)*(unit.current_borehole_wall_temperature-t_in_neu) + t_in_neu
                t_out = (factor * (unit.current_borehole_wall_temperature-t_in_neu/2)+t_in_neu) * 1/(1+factor/2)
                t_in_neu = copy(t_out)
                end

            end
        end 


        # calculate q_in_out_validation
        q = 0
        q = unit.m_in_validation[unit.time_index] * unit.fluid_heat_capacity * (unit.t_in_validation[unit.time_index]-unit.t_out_validation[unit.time_index])*2
        unit.q_in_out_validation[unit.time_index] = q/(unit.probe_depth)
        if q > 0
            unit.q_in_abs_validation[unit.time_index] = q * unit.number_of_probes
            unit.q_out_abs_validation[unit.time_index] = 0
        elseif q < 0
            unit.q_in_abs_validation[unit.time_index] = 0
            unit.q_out_abs_validation[unit.time_index] = q * unit.number_of_probes
        elseif unit.m_in_validation[unit.time_index] == 0
            unit.q_in_abs_validation[unit.time_index] = 0
            unit.q_out_abs_validation[unit.time_index] = 0
        else 
            unit.q_in_abs_validation[unit.time_index] = 0
            unit.q_out_abs_validation[unit.time_index] = 0
        end


        # calculate new boreholewall temperature for next timestep
        sum = 0
        if unit.time_index == 1
            unit.specific_heat_flux_in_out_step[unit.time_index] = unit.q_in_out_validation[unit.time_index]
        else 
            unit.specific_heat_flux_in_out_step[unit.time_index] = unit.q_in_out_validation[unit.time_index] - unit.q_in_out_validation[unit.time_index-1]
        
        end
        
        for i in 0:(unit.time_index - 1)
            sum += unit.specific_heat_flux_in_out_step[unit.time_index-i] * unit.g_function[i+1] / (2*pi*unit.soil_heat_conductivity)
        end
    
    
        unit.current_borehole_wall_temperature = unit.undisturbed_ground_temperature + sum

    end





end

# function that calculates current (highest possible) output temperature of probe field that can be provided. (lower temp. is always possible!)
function current_output_temperature(unit::GeothermalProbes)::Temperature
        # new:
    # max average fluidtemperature == current_borehole_wall_temperature - 0.5 * temperature spread.
    highest_outer_borehole_temp = unit.average_fluid_temperature

    return highest_outer_borehole_temp
end

# function that calculates current (minimum possible) input temperature of probe field that should be provided for regeneration 
# (higher temperature is always possible, here the minimum should be calculated!)
# currenlty only a dummy implementation!!
function current_input_temperature(unit::GeothermalProbes)::Temperature
    # new: Min. current boreholewall temperature necessary to regenerate.
    input_temperature = unit.average_fluid_temperature
    return input_temperature
end

# function to calculate the maximum possible output energy in current time step
# currenlty only a dummy implementation with simple user-input!!
function get_max_output_power(unit::GeothermalProbes)::Float64
    # new: calculate max. outputpower based on max. temperature difference between borehallwall and min. average fluid temperature (VDI 4640-2 : Input should be >= 0°C most of the time)
    
    max_output_power = unit.specific_heat_flux_in_out_absolut[unit.time_index]*unit.probe_depth*unit.number_of_probes
    return max_output_power  
end

# function to calculate the maximum possible input energy in current time step
function get_max_input_power(unit::GeothermalProbes)::Float64
    # Max. input temperature is 15 K above undistrubed ground temperature 
    max_input_power = unit.specific_heat_flux_in_out_absolut[unit.time_index]*unit.probe_depth*unit.number_of_probes


    return max(0,max_input_power)  
end

# function to calculate current new boreholewall temperature with g-functions
function calculate_new_boreholewall_temperature(unit::GeothermalProbes)::Temperature
 # If Validation-Mode is set on 1, no return of this function to save calculation time.
    if unit.set_validation_mode == 1
        
        return 99
    
    else

        radius_pipe_inner = unit.diameter_pipe_inner  / 2
        radius_pipe_outer = unit.diameter_pipe_outer / 2
        radius_borehole = unit.diameter_borehole / 2
        
        # R_B with Hellström - factor = 1; R_B with radius_eq = 4*r_pipe - factor = 4; R_B with radius_eq = 2*r_pipe  - factor = 2
        factor = 1

        
        radius_pipe_eq_inner = factor*radius_pipe_inner
        radius_pipe_eq_outer = factor*radius_pipe_outer
        unit.diameter_pipe_eq_inner = 2*radius_pipe_eq_inner
        unit.diameter_pipe_eq_outer = 2*radius_pipe_eq_outer
        
        # define variable for later calculations
        sum = 0

        # calculate convective heat transfer coefficient alpha in pipie
        # Attention: currently only a copy of geothermal heat collector modell!!!
        alpha_fluid = calculate_alpha_pipe(unit::GeothermalProbes)   # To DO: Calculate out of mass-flow 
            
        
        # calculate effective thermal borehole resistance by multipole method (Hellström 1991) depending on alpha
        sigma = (unit.grout_lambda-unit.soil_heat_conductivity)/(unit.grout_lambda+unit.soil_heat_conductivity)   # dimensionless calculation factor

        distance_pipe_center = unit.shank_spacing / 2

        beta = 1/(2*pi*alpha_fluid*radius_pipe_inner) + 1/(2*pi * unit.pipe_lambda) * log(radius_pipe_outer/radius_pipe_inner) # in (mK)/W

        R_1 = beta + 1/(2*pi*unit.grout_lambda)*
        (log(radius_borehole^2/(2*radius_pipe_outer*distance_pipe_center))+
        sigma*log(radius_borehole^4/(radius_borehole^4-distance_pipe_center^4))-
        radius_pipe_outer^2/(4*distance_pipe_center^2)*(1-sigma*4*distance_pipe_center^4/(radius_borehole^4-distance_pipe_center^4))^2 /
        ((1+2*pi*unit.grout_lambda*beta)/(1-2*pi*unit.grout_lambda*beta)+
        radius_pipe_outer^2/(4*distance_pipe_center^2)*(1+sigma*16*radius_borehole^4*distance_pipe_center^4/((radius_borehole^4-distance_pipe_center^4)^2))
        )
        )

        unit.thermal_resistance = R_1/4
        
        # g-function approach
        if unit.time_index == 1
            unit.specific_heat_flux_in_out_step[unit.time_index] = unit.specific_heat_flux_in_out_absolut[unit.time_index]
        else 
            unit.specific_heat_flux_in_out_step[unit.time_index] = unit.specific_heat_flux_in_out_absolut[unit.time_index] - unit.specific_heat_flux_in_out_absolut[unit.time_index - 1]
        
        end
        
        for i in 0:(unit.time_index - 1)
            sum += unit.specific_heat_flux_in_out_step[unit.time_index-i] * unit.g_function[i+1] / (2*pi*unit.soil_heat_conductivity)
        end

        unit.current_borehole_wall_temperature = unit.undisturbed_ground_temperature + sum
        unit.average_fluid_temperature = unit.current_borehole_wall_temperature + unit.specific_heat_flux_in_out_absolut[unit.time_index] * unit.thermal_resistance    
    end 
end 

function calculate_alpha_pipe(unit::GeothermalProbes)
    # calculate reynolds-number to choose right Nusselt-calculation method
    

    if unit.set_validation_mode == 1
        unit.Re = (4 * unit.m_in_validation[unit.time_index])/
        (pi  * unit.diameter_pipe_eq_inner * unit.fluid_viscosity * unit.fluid_density)

    else 
        unit.Re = (4 * abs(unit.specific_heat_flux_in_out_absolut[unit.time_index])/2 * unit.probe_depth)/
        (unit.fluid_heat_capacity * unit.unloading_temperature_spread * pi  * unit.diameter_pipe_eq_inner * unit.fluid_viscosity * unit.fluid_density)
    end   

            
    # check for laminar flow
    if unit.Re <= 2300
        Nu = calculate_Nu_laminar(unit::GeothermalProbes,unit.Re)
        
    # check for transition flow
    elseif unit.Re > 2300
        # Gielinski
        factor = (unit.Re - 2300)/(1e4-2300)
        Nu = (1- factor)*
        calculate_Nu_laminar(unit::GeothermalProbes, 2300) +
        factor * calculate_Nu_turbulent(unit::GeothermalProbes, 1e4)
        end
    
    alpha = Nu * unit.fluid_lambda / unit.diameter_pipe_eq_inner
    
    return alpha
    
    end
    
    
function calculate_Nu_laminar(unit::GeothermalProbes,Re)
    # Stephan
    Pr_water = 13.44 # Pr Number Water 0 °C
    
    l_rohr = 2*unit.probe_depth # length down and up -> 2x probe depth.
    Nu = 3.66 + (0.0677 * (Re * unit.fluid_prandtl * unit.diameter_pipe_eq_inner/l_rohr)^1.33) /
        (1+0.1* unit.fluid_prandtl * (Re * unit.diameter_pipe_eq_inner/l_rohr)^0.83) *
        (unit.fluid_prandtl/Pr_water)^ 0.1
    
    return Nu
    end 
    
function calculate_Nu_turbulent(unit::GeothermalProbes,Re)
    zeta = (1.8*log(Re)-1.5)^-2
    
    # Gielinski
    Nu = (zeta/8 * Re * unit.fluid_prandtl)/
        (1 + 12.7 * sqrt(zeta/8) * (unit.fluid_prandtl^(2/3)-1))
    
    return Nu
end

# process function that provides energy from the geothermal probes
# according to actual delivered or received energy
function process(unit::GeothermalProbes, parameters::Dict{String,Any})
    # get actual required energy from output interface
    outface = unit.output_interfaces[unit.m_heat_out]  # output interface
    exchange = balance_on(outface, outface.target)     # gather information of output interface
    demand_temp = exchange.temperature                 # get temperature requested by demand (equals max(unit.temperature_field) 
            	                                       # from control-step if demand is not requesting a temperature)

    # check if temperature can be met
    if demand_temp !== nothing && demand_temp > unit.current_output_temperature
        # no calculate_new_probe_field() as this will be done in load() step to avoid double calling!
        return  # no energy delivery possible as requested temperature can not be provided!
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
        # no calculate_new_probe_field() as this will be done in load() step to avoid double calling!
        return # process is only concerned with moving energy to the target
    end
    # write output heat flux into vector
    unit.specific_heat_flux_in_out_absolut[unit.time_index] = energy_demand  / (unit.probe_depth * unit.number_of_probes * 1) # from total energy to specific power of one single probe.
    add!(outface, abs(energy_demand), unit.current_output_temperature)
    
end

function load(unit::GeothermalProbes, parameters::Dict{String,Any})
  
    inface = unit.input_interfaces[unit.m_heat_in]  # input interface
    exchange = balance_on(inface, inface.source)    # gather information of input interface
    supply_temp = exchange.temperature              # get temperature delivered by source 
    energy_available = exchange.balance             # get energy that is provided  
   
    if !unit.regeneration
        # recalculate borehole temperature for next timestep
        calculate_new_boreholewall_temperature(unit::GeothermalProbes)
        return
    end

    # no energy available for loading as load is only concerned when receiving energy from the target
    if energy_available <= 0.0
        # recalculate borehole temperature for next timestep    
        calculate_new_boreholewall_temperature(unit::GeothermalProbes)
        return
    end

    # we can only take in energy if it's at a higher temperature than the probe fields lowest temperature
    if supply_temp !== nothing && supply_temp < unit.current_input_temperature
        calculate_new_boreholewall_temperature(unit::GeothermalProbes)
        return
    end
    # Add loaded specific heat flux to vector
    unit.specific_heat_flux_in_out_absolut[unit.time_index] += energy_available / (unit.probe_depth * unit.number_of_probes * 1) # TO DO: Add Global Time Step!

    # calcute energy that acutally has beed delivered for regeneration and set it to interface 
    # no other limits are present as max_energy for geothermal probes was written in control-step!
    sub!(inface, energy_available, unit.current_input_temperature)
    
    # recalculate borehole temperature for next timestep
    calculate_new_boreholewall_temperature(unit::GeothermalProbes)
    
   
end

function balance_on(
    interface::SystemInterface,
    unit::GeothermalProbes
)::NamedTuple{}
    # check if interface is input or output on unit
    input_sign = unit.uac == interface.target.uac ? -1 : +1
    # check if a balance was already written --> if yes, storage potential will be set to zero as storage was already processed/loaded
    balance_written = interface.max_energy === nothing || interface.sum_abs_change > 0.0

    return (
            balance = interface.balance,
            storage_potential = balance_written ? 0.0 : input_sign * interface.max_energy,  # geothermal probes are handled as storages currently!
            energy_potential = 0.0,
            temperature = interface.temperature
            )
end

function output_values(unit::GeothermalProbes)::Vector{String}
    return ["IN", "OUT", "TEMPERATURE_#NodeNum", "borehole_temperature","average_fluid_temperature","thermal_resistance","Re", "T_out", "Q_out", "Q_in"]
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
        return unit.current_borehole_wall_temperature

    elseif key.value_key =="average_fluid_temperature"
        return unit.average_fluid_temperature

    elseif key.value_key =="thermal_resistance"
        return unit.thermal_resistance
    
    elseif key.value_key =="Re"
        return unit.Re

    elseif key.value_key =="T_out"
        return unit.t_out_validation[unit.time_index]
    elseif key.value_key =="Q_out"
        return unit.q_out_abs_validation[unit.time_index]
    elseif key.value_key =="Q_in"
        return unit.q_in_abs_validation[unit.time_index]
    end
    throw(KeyError(key.value_key))
end

  
     

export GeothermalProbes