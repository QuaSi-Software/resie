"""
Implementation of geothermal probes.

This implementations acts as storage as is can produce and load energy.
## ATTENTION: Geothermal probes are currently work in progress and not completed!!
## ATTENTION: The physical implementation is party nonsense!!
"""
mutable struct GeothermalProbes <: ControlledSystem
    uac::String
    controller::Controller
    sys_function::SystemFunction

    input_interfaces::InterfaceMap
    output_interfaces::InterfaceMap

    m_heat_in::Symbol
    m_heat_out::Symbol

    ambient_temperature_profile::Union{Profile,Nothing}

    unloading_temperature_spread::Temperature
    loading_temperature::Temperature
    loading_temperature_spread::Temperature

    max_output_power::Float64
    max_input_power::Float64

    regeneration::Bool

    temperature_field::Vector
    max_output_energy::Float64
    max_input_energy::Float64
    current_output_temperature::Temperature
    current_input_temperature::Temperature
    ambient_temperature::Temperature
    last_timestep_calculated::Float64

    function GeothermalProbes(uac::String, config::Dict{String,Any})
        m_heat_in = Symbol(default(config, "m_heat_in", "m_h_w_ht1"))
        m_heat_out = Symbol(default(config, "m_heat_out", "m_h_w_lt1"))
        register_media([m_heat_in, m_heat_out])

        ambient_temperature_profile = "ambient_temperature_profile_path" in keys(config) ?
                                      Profile(config["ambient_temperature_profile_path"]) :
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
            default(config, "unloading_temperature_spread", 5.0),   # temperature spread between forward and return flow during unloading            
            default(config, "loading_temperature", nothing),        # nominal high temperature for loading geothermal probe storage, can also be set from other end of interface
            default(config, "loading_temperature_spread", 5.0),     # temperature spread between forward and return flow during loading         
            config["max_output_power"],  # maximum output power set by user, may change later to be calculated from other inputs like specific heat transfer rate
            config["max_input_power"],   # maximum input power set by user, may change later to be calculated from other inputs like specific heat transfer rate
            default(config, "regeneration", true),   # flag if regeneration should be taken into account
            fill(15, 100),               # temperature_field of probe field with initially temperature of equally 8.5°, acts currently only as dummy!
            0.0,                         # max_output_energy in every time step, calculated in control()
            0.0,                         # max_input_energy in every time step, calculated in control()
            0.0,                         # output temperature in current time step, calculated in control()
            0.0,                         # input temperature in current time step, calculated in control()
            0.0,                         # ambient temperature in current time step, calculated in control()
            -1.0                         # last timestep that was calculated; used to avoid double calculation of temperature field. set to -1 for the beginning
        )
    end
end

function control(
    unit::GeothermalProbes,
    systems::Grouping,
    parameters::Dict{String,Any}
)
    # in case there is a state machine for geothermal probes
    move_state(unit, systems, parameters)

    # get ambient temperature from profile for current time step if needed (probably only for geothermal collectors)
    unit.ambient_temperature = Profiles.value_at_time(unit.ambient_temperature_profile, parameters["time"])

    ## Assumption:
        # Temperatur output of geothermal probes is independend of the maximum thermal energy output in the current timestep.
        # max. temperature could be assumed to be the highest outer borehole temperature (constant within one timestep)
        # this assumption is probably only valid if the timestep is not too long!

    # calculate maximum possible output temperatures for energy output and set temperature to output interface 
    # (make sure, that a possibly higher temeprature that is already written in the interface is not overwritten)
    unit.current_output_temperature = current_output_temperature(unit) # of geothermal probe field
    unit.output_interfaces[unit.m_heat_out].temperature = highest_temperature(
                                                                            unit.current_output_temperature, 
                                                                            unit.output_interfaces[unit.m_heat_out].temperature
                                                                            )
   
    # get input temperature for energy input (regeneration) and set temperature to input interface
    # (currently from user input or from lowest temperature in probe field...)
    # (make sure, that a possibly higher temperature that is already written in the interface is not overwritten)
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
        # ToDo: what if GTP is connected to Bus?!
    if unit.output_interfaces[unit.m_heat_out].temperature > unit.current_output_temperature
        unit.max_output_energy = 0.0  # no energy can be provided if requested temperature is higher than max. temperature of probe field
    else
        unit.max_output_energy = get_max_output_power(unit)  # Attention: needs watt_to_wh function!! ToDo
    end

    if unit.regeneration
        if unit.input_interfaces[unit.m_heat_in].temperature < unit.current_input_temperature
            unit.max_input_energy = 0.0 # no energy can be taken if available temperature is less than minimum possible temperature to load the probe field
        else
            unit.max_input_energy = get_max_input_power(unit) # Attention: needs watt_to_wh function!! ToDo
        end
    end

    # set max_energy to interfaces to provide information for connected systems
    set_max_energy!(unit.output_interfaces[unit.m_heat_out], unit.max_output_energy)
    if unit.regeneration
        set_max_energy!(unit.input_interfaces[unit.m_heat_in], unit.max_input_energy)
    end

end

# function that calculates current (highest possible) output temperature of probe field that can be provided. (lower temp. is always possible!)
# currenlty only a dummy implementation!!
function current_output_temperature(unit::GeothermalProbes)::Temperature
    highest_outer_borehole_temp = maximum(unit.temperature_field)
    # could be the highest outer borewhole temperature of the probe field in the last timestep
    return highest_outer_borehole_temp
end

# function that calculates current (minimum possible) input temperature of probe field that should be provided for regeneration 
# (higher temperature is always possible, here the minimum should be calculated!)
# currenlty only a dummy implementation!!
function current_input_temperature(unit::GeothermalProbes)::Temperature
    input_temperature = highest_temperature(minimum(unit.temperature_field)+unit.loading_temperature_spread, unit.loading_temperature ) 
    # could be the lowest outer borewhole temperature plot the given temperature spread or a user-specified loading_temperature if given
    return input_temperature
end

# function to calculate the maximum possible output energy in current time step
# currenlty only a dummy implementation with simple user-input!!
function get_max_output_power(unit::GeothermalProbes)::Float64
    return unit.max_output_power  # could be calculated by maximum specific unloading energy or by maximal mass flow or...
end

# function to calculate the maximum possible input energy in current time step
# currenlty only a dummy implementation with simple user-input!!
function get_max_input_power(unit::GeothermalProbes)::Float64
    return unit.max_input_power  # could be calculated by maximum specific loading energy or by maximal mass flow or...
end

# function that calculates the new state of the geothermal probe field for the next timestep
# currenlty only a dummy implementation that calculates physical nonsense!!
function calculate_new_probe_field_temperatures!(unit::GeothermalProbes, parameters::Dict{String,Any}, energy::Float64)
    # Note: energy is positive for loading (regeneration) and negative for unloading the probe field as this point
    current_timestep = parameters["time"]
    if current_timestep == unit.last_timestep_calculated # check if current timestep has already been simulated
        if !iszero(energy) # this actually should never be the case if the GTP are not loaded AND unloaded within one timestep
            # not using ambient temperature here as this has already been calculated in the current timestep, but energy has changed, some
            # a recalculation is necessary.
            unit.temperature_field = unit.temperature_field .+ energy/200 
        end
    else
        # timestep has not been calculated yet (normal case)
        unit.temperature_field = unit.temperature_field .+ energy/1000 - (unit.temperature_field .- unit.ambient_temperature)./200
    end

    # update information if the current timestep has already been calculated
    unit.last_timestep_calculated = current_timestep

end

# produce function that provides energy from the geothermal probes and calculates new temperatures 
# according to actual delivered or received energy
function produce(unit::GeothermalProbes, parameters::Dict{String,Any}, watt_to_wh::Function)
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
        return # produce is only concerned with moving energy to the target
    end

    # calcute energy that acutally can be delivered and set it to the output interface 
    # no other limits are present as max_energy for geothermal probes was written in control-step!
    add!(outface, abs(energy_demand), unit.current_output_temperature)
    # recalculate probe field temperatures for next timestep
    calculate_new_probe_field_temperatures!(unit, parameters, energy_demand) # energy_demand is negative here and has to be negative in calculate_new_probe_field_temperatures!() for unloading

end

function load(unit::GeothermalProbes, parameters::Dict{String,Any}, watt_to_wh::Function)
    # we can assume that energy will be either be taken from or fed into the geothermal probe field within one time step,
    # but not both within one time step - right? Therefore we can call calculate_new_probe_field_temperatures!() from produce
    # and from load, but it will never be calculated twice as only one call will be performed with energy>0. If both load
    # and produce happens at the same time step with an energy flow>0, the temperatures of the GTP field will be calculated twice 
    # in each timestep with the current implementation. This should work, but is not computationally effiecient.
    if !unit.regeneration
        # recalculate probe field temperatures for next timestep (function checks if is has already been calculated in the current timestep)
        calculate_new_probe_field_temperatures!(unit, parameters, 0.0)  # call calculate_new_probe_field() to calculate new temperatures of field to account for possible ambient effects
        return
    end

    # get actual delivered energy from input interface
    inface = unit.input_interfaces[unit.m_heat_in]  # input interface
    exchange = balance_on(inface, inface.source)    # gather information of input interface
    supply_temp = exchange.temperature              # get temperature delivered by source 
    energy_available = exchange.balance             # get energy that is provided

    # no energy available for loading as load is only concerned when receiving energy from the target
    if energy_available <= 0.0
        # recalculate probe field temperatures for next timestep (function checks if is has already been calculated in the current timestep)
        calculate_new_probe_field_temperatures!(unit, parameters, 0.0)  # call calculate_new_probe_field() to calculate new temperatures of field to account for possible ambient effects
        return
    end

    # we can only take in energy if it's at a higher temperature than the probe fields lowest temperature
    if supply_temp !== nothing && supply_temp < unit.current_input_temperature
        # recalculate probe field temperatures for next timestep (function checks if is has already been calculated in the current timestep)
        calculate_new_probe_field_temperatures!(unit, parameters, 0.0)  # call calculate_new_probe_field() to calculate new temperatures of field to account for possible ambient effects
        return
    end
    
    # calcute energy that acutally has beed delivered for regeneration and set it to interface 
    # no other limits are present as max_energy for geothermal probes was written in control-step!
    sub!(inface, energy_available, unit.current_input_temperature)
    # recalculate probe field temperatures for next timestep
    calculate_new_probe_field_temperatures!(unit, parameters, energy_available)  # energy_availability is positive here
        
end

function balance_on(
    interface::SystemInterface,
    unit::GeothermalProbes
)::NamedTuple{}
    # check if interface is input or output on unit
    input_sign = unit.uac == interface.target.uac ? -1 : +1
    # check if a balance was already written --> if yes, storage potential will be set to zero as storage was already produced/loaded
    balance_written = interface.max_energy === nothing || interface.sum_abs_change > 0.0

    return (
            balance = interface.balance,
            storage_potential = balance_written ? 0.0 : input_sign * interface.max_energy,  # geothermal probes are handled as storages currently!
            energy_potential = 0.0,
            temperature = interface.temperature
            )
end

function output_values(unit::GeothermalProbes)::Vector{String}
    return ["IN", "OUT", "TEMPERATURE_#NodeNum"]
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
    end
    throw(KeyError(key.value_key))
end

export GeothermalProbes