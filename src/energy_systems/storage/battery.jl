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

    model_type::String

    capacity::Float64
    load::Float64
    
    charge_efficiency::Float64
    discharge_efficiency::Float64
    self_discharge_rate::Float64

    load_end_of_last_timestep::Float64
    losses::Float64
    max_charge::Float64
    max_discharge::Float64

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
                   config["load"], # load
                   config["charge_efficiency"],
                   config["discharge_efficiency"],
                   config["self_discharge_rate"],
                   0.0, # load_end_of_last_timestep
                   0.0, # losses
                   0.0, # max_charge
                   0.0) # max_discharge
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

    unit.losses = unit.load * unit.self_discharge_rate
    unit.load -= unit.losses

    if discharge_is_allowed(unit.controller, sim_params)
        unit.max_discharge = unit.load * unit.discharge_efficiency
    else
        unit.max_discharge = 0.0
    end
    set_max_energy!(unit.output_interfaces[unit.medium], unit.max_discharge)

    if charge_is_allowed(unit.controller, sim_params)
        unit.max_charge = (unit.capacity - unit.load) / unit.charge_efficiency
    else
        unit.max_charge = 0.0
    end
    set_max_energy!(unit.input_interfaces[unit.medium], unit.max_charge)
end

function process(unit::Battery, sim_params::Dict{String,Any})
    if unit.max_discharge < sim_params["epsilon"]
        set_max_energy!(unit.output_interfaces[unit.medium], 0.0)
        return
    end

    outface = unit.output_interfaces[unit.medium]
    exchanges = balance_on(outface, outface.target)
    energy_demand = balance(exchanges) + energy_potential(exchanges)

    if energy_demand >= 0.0
        set_max_energy!(unit.output_interfaces[unit.medium], 0.0)
        return # process is only concerned with moving energy to the target
    end

    if unit.load > abs(energy_demand) / unit.discharge_efficiency
        unit.losses += abs(energy_demand) * (1.0 / unit.discharge_efficiency - 1)
        unit.load += energy_demand / unit.discharge_efficiency
        add!(outface, abs(energy_demand))
    else
        unit.losses += unit.load * (1.0 - unit.discharge_efficiency)
        add!(outface, unit.load * unit.discharge_efficiency)
        unit.load = 0.0
    end
end

function load(unit::Battery, sim_params::Dict{String,Any})
    if unit.max_charge < sim_params["epsilon"]
        set_max_energy!(unit.input_interfaces[unit.medium], 0.0)
        unit.load_end_of_last_timestep = copy(unit.load)
        return
    end

    inface = unit.input_interfaces[unit.medium]
    exchanges = balance_on(inface, inface.source)
    energy_available = balance(exchanges) + energy_potential(exchanges)

    if energy_available <= 0.0
        set_max_energy!(unit.input_interfaces[unit.medium], 0.0)
        unit.load_end_of_last_timestep = copy(unit.load)
        return # load is only concerned with receiving energy from the source
    end

    diff = unit.capacity - unit.load
    if diff > energy_available * unit.charge_efficiency
        unit.load += energy_available * unit.charge_efficiency
        sub!(inface, energy_available * unit.charge_efficiency)
    else
        unit.load = unit.capacity
        sub!(inface, diff * unit.charge_efficiency)
    end

    unit.load_end_of_last_timestep = copy(unit.load)
end

function output_values(unit::Battery)::Vector{String}
    return [string(unit.medium) * " IN",
            string(unit.medium) * " OUT",
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

