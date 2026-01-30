"""
Implementation of a battery component holding electric charge.

For the moment the implementation remains simple with only one state (its charge) and one
parameter (its capacity).
"""
Base.@kwdef mutable struct Battery <: Component
    uac::String
    controller::Controller
    sys_function::SystemFunction

    input_interfaces::InterfaceMap
    output_interfaces::InterfaceMap

    medium::Symbol

    capacity::Float64
    load::Float64
    load_end_of_last_timestep::Float64
    losses::Float64

    max_charge::Float64
    max_discharge::Float64

    process_done::Bool
    load_done::Bool

    function Battery(uac::String, config::Dict{String,Any}, sim_params::Dict{String,Any})
        medium = Symbol(default(config, "medium", "m_e_ac_230v"))
        register_media([medium])

        return new(uac, # uac
                   Controller(default(config, "control_parameters", nothing)),
                   sf_storage, # sys_function
                   InterfaceMap(medium => nothing), # input_interfaces
                   InterfaceMap(medium => nothing), # output_interfaces
                   medium,
                   config["capacity"], # capacity
                   config["capacity"] * config["load"],     # load
                   0.0, # load_end_of_last_timestep
                   0.0, # losses
                   0.0, # max_charge
                   0.0, # max_discharge
                   false,  # process_done, bool indicating if the process step has already been performed in the current time step
                   false)  # load_done, bool indicating if the load step has already been performed in the current time step
    end
end

function initialise!(unit::Battery, sim_params::Dict{String,Any})
    set_storage_transfer!(unit.input_interfaces[unit.medium],
                          unload_storages(unit.controller, unit.medium))
    set_storage_transfer!(unit.output_interfaces[unit.medium],
                          load_storages(unit.controller, unit.medium))

    unit.load_end_of_last_timestep = copy(unit.load)
end

function reset(unit::Battery)
    invoke(reset, Tuple{Component}, unit)

    unit.max_charge = 0.0
    unit.max_discharge = 0.0
end

function control(unit::Battery,
                 components::Grouping,
                 sim_params::Dict{String,Any})
    update(unit.controller)

    if discharge_is_allowed(unit.controller, sim_params)
        unit.max_discharge = unit.load
    else
        unit.max_discharge = 0.0
    end
    set_max_energy!(unit.output_interfaces[unit.medium], unit.max_discharge)

    if charge_is_allowed(unit.controller, sim_params)
        unit.max_charge = unit.capacity - unit.load
    else
        unit.max_charge = 0.0
    end
    set_max_energy!(unit.input_interfaces[unit.medium], unit.max_charge)
end

function process(unit::Battery, sim_params::Dict{String,Any})
    if unit.max_discharge < sim_params["epsilon"]
        handle_component_update!(unit, "process", sim_params)
        set_max_energy!(unit.output_interfaces[unit.medium], 0.0)
        return
    end

    outface = unit.output_interfaces[unit.medium]
    exchanges = balance_on(outface, outface.target)
    energy_demand = balance(exchanges) + energy_potential(exchanges)

    if energy_demand >= 0.0
        handle_component_update!(unit, "process", sim_params)
        set_max_energy!(unit.output_interfaces[unit.medium], 0.0)
        return # process is only concerned with moving energy to the target
    end

    if unit.load > abs(energy_demand)
        unit.load += energy_demand
        add!(outface, abs(energy_demand))
    else
        add!(outface, unit.load)
        unit.load = 0.0
    end

    handle_component_update!(unit, "process", sim_params)
end

function handle_component_update!(unit::Battery, step::String, sim_params::Dict{String,Any})
    if step == "process"
        unit.process_done = true
    elseif step == "load"
        unit.load_done = true
    end
    if unit.process_done && unit.load_done
        # update component
        unit.load_end_of_last_timestep = copy(unit.load)
        # reset 
        unit.process_done = false
        unit.load_done = false
    end
end

function load(unit::Battery, sim_params::Dict{String,Any})
    if unit.max_charge < sim_params["epsilon"]
        handle_component_update!(unit, "load", sim_params)
        set_max_energy!(unit.output_interfaces[unit.medium], 0.0)
        return
    end

    inface = unit.input_interfaces[unit.medium]
    exchanges = balance_on(inface, inface.source)
    energy_available = balance(exchanges) + energy_potential(exchanges)

    if energy_available <= 0.0
        handle_component_update!(unit, "load", sim_params)
        set_max_energy!(unit.input_interfaces[unit.medium], 0.0)
        return # load is only concerned with receiving energy from the source
    end

    diff = unit.capacity - unit.load
    if diff > energy_available
        unit.load += energy_available
        sub!(inface, energy_available)
    else
        unit.load = unit.capacity
        sub!(inface, diff)
    end

    handle_component_update!(unit, "load", sim_params)
end

function output_values(unit::Battery)::Vector{String}
    return [string(unit.medium) * ":IN",
            string(unit.medium) * ":OUT",
            "Load",
            "Load%",
            "Capacity",
            "LossesGains"]
end

function output_value(unit::Battery, key::OutputKey)::Float64
    if key.value_key == "IN"
        return calculate_energy_flow(unit.input_interfaces[key.medium])
    elseif key.value_key == "OUT"
        return calculate_energy_flow(unit.output_interfaces[key.medium])
    elseif key.value_key == "Load"
        return unit.load
    elseif key.value_key == "Load%"
        return 100 * unit.load / unit.capacity
    elseif key.value_key == "Capacity"
        return unit.capacity
    elseif key.value_key == "LossesGains"
        return -unit.losses
    end
    throw(KeyError(key.value_key))
end

export Battery
