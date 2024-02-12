"""
Implementation of geothermal probes.
This implementations acts as storage as it can produce and load energy.
"""


# Current solution to get g-function values.
# read g-function .txt file
file = open("C:/Users/steinacker/Lokal/git_Resie/src/energy_systems/heat_sources/g_ges_vector.txt", "r") # TODO
    g_function = Vector{Float64}()
    for line in eachline(file)
        push!(g_function, parse(Float64, line))
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
    soil_undisturbed_ground_temperature::Temperature
    soil_heat_conductivity::Float64
    borehole_thermal_resistance::Float64
    g_function::Vector
    time_index::Int
    fluid_temperature::Temperature
    borehole_current_wall_temperature::Temperature

    specific_heat_flux_in_out_absolut::Vector
    specific_heat_flux_in_out_step::Vector

    probe_depth::Float64
    probe_number::Float64

    pipe_diameter_outer::Float64
    pipe_diameter_inner::Float64

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
            default(config, "loading_temperature", nothing),      # nominal high temperature for loading geothermal probe storage, can also be set from other end of interface
            default(config, "loading_temperature_spread", 3),     # temperature spread between forward and return flow during loading         
            config["max_output_power"],  # maximum output power set by user, may change later to be calculated from other inputs like specific heat transfer rate
            config["max_input_power"],   # maximum input power set by user, may change later to be calculated from other inputs like specific heat transfer rate
            default(config, "regeneration", true),   # flag if regeneration should be taken into account
            0.0,                         # max_output_energy in every time step, calculated in control()
            0.0,                         # max_input_energy in every time step, calculated in control()
            0.0,                         # output temperature in current time step, calculated in control()
            0.0,                         # input temperature in current time step, calculated in control()
            default(config, "soil_undisturbed_ground_temperature", 11.0),    # Considered as constant
            default(config, "soil_heat_conductivity", 1.5),                  # Heat conductivity of surrounding soil, homogenous and constant
            default(config, "borehole_thermal_resistance", 0.10),            # thermal resistance in (m K)/W
            g_function,                             # pre-calculated multiscale g-function. Calculated in pre-processing.
            0,                                      # index of current time step to get access on time dependent g-function values
            0.0,                                    # average fluid temperature
            4,                                      # set boreholewall-starting-temperature
            zeros(219000),                          # vector to write specific heat flux in eacht time step
            zeros(219000),                          # vector to write specific heat flux differeces in eacht time step for g-function approach
            default(config, "probe_depth", 150),    # depth (or length) of a single geothermal probe
            36,                                     # number of geothermal probes in the borefield

            default(config, "pipe_diameter_outer", 0.032),  # outer pipe diameter
            default(config, "pipe_diameter_inner", 0.026),  # inner pipe diameter
           
            default(config, "fluid_specific_heat_capacity", 3800),  # specific heat capacity brine at 0 °C (25 % glycol 75 % water (interpolated)) 
            default(config, "fluid_density", 1045),                 # density brine at 0 °C (25 % glycol 75 % water (interpolated))
            default(config, "fluid_kinematic_viscosity", 3.9e-6),   # viscosity brine at 0 °C (25 % glycol 75 % water (interpolated)) 
            default(config, "fluid_heat_conductivity", 0.5) ,       # heat conductivity brine at 0 °C (25 % glycol 75 % water (interpolated))
            default(config, "fluid_prandtl_number", 30),            # prandtl-number brine at 0 °C (25 % glycol 75 % water (interpolated)) 
            
            default(config, "grout_heat_conductivity", 2),          # lambda grout / filling material in W/(mK)   
            default(config, "pipe_heat_conductivity", 0.42),        # lambda of inner pipes
            default(config, "borehole_diameter", 0.15),             # borehole diameter in m.
            0.1,     # shank-spacing = distance between inner pipes in borehole. Needed for calculation of thermal borehole resistance.
            0        # Reynoldsnumber. To be calculated in Function later.
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
end

# function that calculates current (highest possible) output temperature of probe field that can be provided. (lower temp. is always possible!)
function current_output_temperature(unit::GeothermalProbes)::Temperature
    highest_outer_borehole_temp = unit.fluid_temperature + unit.unloading_temperature_spread/2
    return highest_outer_borehole_temp
end

# function that calculates current (minimum possible) input temperature of probe field that should be provided for regeneration 
# (higher temperature is always possible, here the minimum should be calculated!)
function current_input_temperature(unit::GeothermalProbes)::Temperature
    # new: Min. current boreholewall temperature necessary to regenerate.
    input_temperature = unit.fluid_temperature - unit.loading_temperature_spread/2
    return input_temperature
end

# function to calculate the maximum possible output energy in current time step
# currenlty only a dummy implementation with simple user-input!!
function get_max_output_power(unit::GeothermalProbes)::Float64
    # new: calculate max. outputpower based on max. temperature difference between borehallwall and min. average fluid temperature 
    # (VDI 4640-2 : Input should be >= 0°C most of the time)
    max_output_power = unit.specific_heat_flux_in_out_absolut[unit.time_index]*unit.probe_depth*unit.probe_number
    return max_output_power  
end

# function to calculate the maximum possible input energy in current time step
function get_max_input_power(unit::GeothermalProbes)::Float64
    # Max. input temperature is 15 K above undistrubed ground temperature 
    max_input_power = unit.specific_heat_flux_in_out_absolut[unit.time_index]*unit.probe_depth*unit.probe_number
    return max(0,max_input_power)  
end

# function to calculate current new boreholewall temperature with g-functions
function calculate_new_boreholewall_temperature(unit::GeothermalProbes)::Temperature
    
    # calculate radius values
    radius_pipe_inner = unit.pipe_diameter_inner  / 2
    radius_pipe_outer = unit.pipe_diameter_outer / 2
    radius_borehole = unit.borehole_diameter / 2   
    
    # R_B with Hellström 
    radius_pipe_eq_inner = radius_pipe_inner
    unit.pipe_diameter_inner = 2*radius_pipe_eq_inner
   
    # define variable for later calculations
    sum = 0

    # calculate convective heat transfer coefficient alpha in pipie
    alpha_fluid = calculate_alpha_pipe(unit::GeothermalProbes)    
    
    # calculate effective thermal borehole resistance by multipole method (Hellström 1991) depending on alpha
    sigma = (unit.grout_heat_conductivity-unit.soil_heat_conductivity)/(unit.grout_heat_conductivity+unit.soil_heat_conductivity)   # dimensionless calculation factor
    distance_pipe_center = unit.shank_spacing / 2
    beta = 1/(2*pi*alpha_fluid*radius_pipe_inner) + 1/(2*pi * unit.pipe_heat_conductivity) * log(radius_pipe_outer/radius_pipe_inner) # in (mK)/W
    
    R_1 = beta + 1/(2*pi*unit.grout_heat_conductivity)*
          (log(radius_borehole^2/(2*radius_pipe_outer*distance_pipe_center))+
          sigma*log(radius_borehole^4/(radius_borehole^4-distance_pipe_center^4))-
          radius_pipe_outer^2/(4*distance_pipe_center^2)*(1-sigma*4*distance_pipe_center^4/(radius_borehole^4-distance_pipe_center^4))^2 /
          ((1+2*pi*unit.grout_heat_conductivity*beta)/(1-2*pi*unit.grout_heat_conductivity*beta)+
          radius_pipe_outer^2/(4*distance_pipe_center^2)*(1+sigma*16*radius_borehole^4*distance_pipe_center^4/((radius_borehole^4-distance_pipe_center^4)^2))))

    unit.borehole_thermal_resistance = R_1/4
    
    # g-function approach
    if unit.time_index == 1
        unit.specific_heat_flux_in_out_step[unit.time_index] = unit.specific_heat_flux_in_out_absolut[unit.time_index]
    else 
        unit.specific_heat_flux_in_out_step[unit.time_index] = unit.specific_heat_flux_in_out_absolut[unit.time_index] - unit.specific_heat_flux_in_out_absolut[unit.time_index - 1]
    end
        
    for i in 0:(unit.time_index - 1)    
        sum += unit.specific_heat_flux_in_out_step[unit.time_index-i] * unit.g_function[i+1] / (2*pi*unit.soil_heat_conductivity)
    end
    unit.borehole_current_wall_temperature = unit.soil_undisturbed_ground_temperature + sum
    unit.fluid_temperature = unit.borehole_current_wall_temperature + unit.specific_heat_flux_in_out_absolut[unit.time_index] * unit.borehole_thermal_resistance    

end 

function calculate_alpha_pipe(unit::GeothermalProbes)
    # calculate temeprature-dependant Prantl-Number of the heat carrier fluid.
    fluid_dynamic_viscosity = 0.0000017158* unit.fluid_temperature^2 - 0.0001579079*unit.fluid_temperature+0.0048830621
    unit.fluid_heat_conductivity = 0.0010214286 * unit.fluid_temperature + 0.447
    unit.fluid_prandtl_number = fluid_dynamic_viscosity * unit.fluid_specific_heat_capacity / unit.fluid_heat_conductivity 

    probe_total_heat_flux_in_out = abs(unit.specific_heat_flux_in_out_absolut[unit.time_index])/2 * unit.probe_depth
    # calculate reynolds-number (based on dynamic viscosity)
    unit.fluid_reynolds_number = (probe_total_heat_flux_in_out)/
                                 (unit.fluid_specific_heat_capacity * unit.unloading_temperature_spread * pi / 4 * unit.pipe_diameter_inner * fluid_dynamic_viscosity)
    
    # # old: Reynoldsnumber calculation based on kinematic viscosity.
    # # calculate reynolds-number to choose right Nusselt-calculation method
    #     re_old = (4 * abs(unit.specific_heat_flux_in_out_absolut[unit.time_index])/2 * unit.probe_depth)/
    #     (unit.fluid_specific_heat_capacity * unit.unloading_temperature_spread * pi  * unit.pipe_diameter_inner * unit.fluid_kinematic_viscosity * unit.fluid_density)
     
    # check for laminar flow
    if unit.fluid_reynolds_number <= 2300
        Nu = calculate_Nu_laminar(unit::GeothermalProbes,unit.fluid_reynolds_number)
        
    # check for transition flow
    elseif unit.fluid_reynolds_number > 2300
        # Gielinski
        factor = (unit.fluid_reynolds_number - 2300)/(1e4-2300)
        Nu = (1- factor)*
             calculate_Nu_laminar(unit::GeothermalProbes, 2300) +
             factor * calculate_Nu_turbulent(unit::GeothermalProbes, 1e4)
    end
    
    alpha = Nu * unit.fluid_heat_conductivity / unit.pipe_diameter_inner
    
    return alpha
end
    
function calculate_Nu_laminar(unit::GeothermalProbes,fluid_reynolds_number)
    
    # Approach used in Ramming 2007 from Elsner, Norbert; Fischer, Siegfried; Huhn, Jörg; „Grundlagen der Technischen Thermodynamik“,  Band 2 Wärmeübertragung, Akademie Verlag, Berlin 1993. 
    # TODO: Need to be changed in documentation.
    k_a = 0.0         # initializing k_a
    k_n = 0.0        # initializing k_n
    # substitutions: 
    k_a = 1.1 - 1 / (3.4+0.0667*unit.fluid_prandtl_number)
    k_n = 0.35 + 1 / (7.825 + 2.6 * sqrt(unit.fluid_prandtl_number))
    
    # calculate Nu-Number
    Nu = ((k_a/(1-k_n)*(unit.fluid_prandtl_number*unit.pipe_diameter_inner*unit.fluid_reynolds_number/(unit.probe_depth*2))^k_n)^3+4.364^3)^(1/3)
    
    return Nu
end 
    
function calculate_Nu_turbulent(unit::GeothermalProbes,fluid_reynolds_number)
    zeta = (1.8*log(fluid_reynolds_number)-1.5)^-2
    
    # Gielinski
    Nu = (zeta/8 * fluid_reynolds_number * unit.fluid_prandtl_number)/
        (1 + 12.7 * sqrt(zeta/8) * (unit.fluid_prandtl_number^(2/3)-1))
    
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
    unit.specific_heat_flux_in_out_absolut[unit.time_index] = energy_demand  / (unit.probe_depth * unit.probe_number * 1) # from total energy to specific power of one single probe.
    add!(outface, abs(energy_demand), unit.current_output_temperature)
    
end

function load(unit::GeothermalProbes, parameters::Dict{String,Any})
  
    inface = unit.input_interfaces[unit.m_heat_in]  # input interface
    exchange = balance_on(inface, inface.source)    # gather information of input interface
    supply_temp = exchange.temperature              # get temperature delivered by source 
    energy_available = exchange.balance             # get energy that is provided  
   
    if (!unit.regeneration ||                                                        # no energy available if regeneration is turned off
        energy_available <= 0.0 ||                                                   # no energy available for loading as load is only concerned when receiving energy from the target
        (supply_temp !== nothing && supply_temp < unit.current_input_temperature)    # we can only take in energy if it's at a higher temperature than the probe fields lowest temperature
    )
        # recalculate borehole temperature for next timestep
        calculate_new_boreholewall_temperature(unit::GeothermalProbes)
        return
    end

    # Add loaded specific heat flux to vector
    unit.specific_heat_flux_in_out_absolut[unit.time_index] += energy_available / (unit.probe_depth * unit.probe_number * 1) # TODO: Add Global Time Step!

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
    return ["IN", "OUT", "TEMPERATURE_#NodeNum", "borehole_temperature","fluid_temperature","borehole_thermal_resistance","fluid_reynolds_number", "T_out", "Q_out", "Q_in", "dT_monthly","dT_hourly" ]
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