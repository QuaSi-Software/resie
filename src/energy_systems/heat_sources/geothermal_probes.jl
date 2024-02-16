"""
Implementation of geothermal probes.
This implementations acts as storage as it can produce and load energy.
"""

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
    borehole_thermal_resistance::Float64
    g_function::Vector{Float64}
    time_index::Int
    fluid_temperature::Temperature
    borehole_current_wall_temperature::Temperature

    energy_in_out_per_probe_meter::Vector{Float64}
    energy_in_out_difference_per_probe_meter::Vector{Float64}

    probe_depth::Float64
    number_of_probes::Int
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
            default(config, "soil_heat_conductivity", 1.5),                  # Heat conductivity of surrounding soil, homogenous and constant
            default(config, "borehole_thermal_resistance", 0.10),            # thermal resistance in (m K)/W
            g_function,                             # pre-calculated multiscale g-function. Calculated in pre-processing.
            0,                                      # index of current time step to get access on time dependent g-function values
            0.0,                                    # average fluid temperature
            4.0,                                    # set boreholewall-starting-temperature TODO
            
            [],                                     # vector to hold specific energy sum (in and out) per probe meter in each time step
            [],                                     # vector to hold specific energy per probe meter as difference for each step, used for g-function approach 

            default(config, "probe_depth", 150.0),  # depth (or length) of a single geothermal probe
            36,                                     # number of geothermal probes in the borefield TODO
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
            default(config, "shank_spacing", 0.1),                   # shank-spacing = distance between inner pipes in borehole, diagonal through borehole center. Needed for calculation of thermal borehole resistance.

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
end

function control(
    unit::GeothermalProbes,
    components::Grouping,
    sim_params::Dict{String,Any}
)
    # time index, necessary for g-function approach
    unit.time_index = unit.time_index + 1 

    # Initialization
    if unit.time_index == 1
        # TODO add preprocessing and move to init()

        unit.energy_in_out_per_probe_meter = zeros(219000)  # TODO change to global number of time steps 
        unit.energy_in_out_difference_per_probe_meter = zeros(219000) # TODO change to global number of time steps 

        unit.radius_pipe_inner = unit.pipe_diameter_inner / 2
        unit.radius_pipe_outer = unit.pipe_diameter_outer / 2
        unit.radius_borehole = unit.borehole_diameter / 2
        unit.distance_pipe_center = unit.shank_spacing / 2
    
        unit.max_output_energy = watt_to_wh(unit.max_output_power * unit.probe_depth * unit.number_of_probes)
        unit.max_input_energy = watt_to_wh(unit.max_input_power * unit.probe_depth * unit.number_of_probes)  
    end

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
        return calculate_energy_flow(unit.input_interfaces[key.medium], unit)
    elseif key.value_key == "OUT"
        return calculate_energy_flow(unit.output_interfaces[key.medium], unit)
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