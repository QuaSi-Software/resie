"""
Implementation of a battery component holding electric charge.

For the moment the implementation remains simple with only one state (its charge) and one
parameters (its capacity). However the default operation strategy is more complex and
toggles the processing of the battery dependant on available PV power and its own charge.
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
    losses::Float64

    function Battery(uac::String, config::Dict{String,Any}, sim_params::Dict{String,Any})
        medium = Symbol(default(config, "medium", "m_e_ac_230v"))
        register_media([medium])

        return new(
            uac, # uac
            controller_for_strategy( # controller
                config["strategy"]["name"], config["strategy"], sim_params
            ),
            sf_storage, # sys_function
            InterfaceMap( # input_interfaces
                medium => nothing
            ),
            InterfaceMap( # output_interfaces
                medium => nothing
            ),
            medium,
            config["capacity"], # capacity
            config["load"], # load
            0.0  # losses
        )
    end
end

function initialise!(unit::Battery, sim_params::Dict{String,Any})
    set_storage_transfer!(
        unit.input_interfaces[unit.medium],
        default(
            unit.controller.parameter, "unload_storages " * String(unit.medium), true
        )
    )
    set_storage_transfer!(
        unit.output_interfaces[unit.medium],
        default(
            unit.controller.parameter, "load_storages " * String(unit.medium), true
        )
    )
end

function control(
    unit::Battery,
    components::Grouping,
    sim_params::Dict{String,Any}
)
    move_state(unit, components, sim_params)

    set_max_energy!(unit.input_interfaces[unit.medium], unit.capacity - unit.load)
    set_max_energy!(unit.output_interfaces[unit.medium], unit.load)

end


function balance_on(
    interface::SystemInterface,
    unit::Battery
)::Vector{EnergyExchange}
    caller_is_input = unit.uac == interface.target.uac

    return [EnEx(
        balance=interface.balance,
        uac=unit.uac,
        energy_potential=0.0,
        storage_potential=caller_is_input ? -(unit.capacity - unit.load) : unit.load,
        temperature_min=interface.temperature_min,
        temperature_max=interface.temperature_max,
        pressure=nothing,
        voltage=nothing,
    )]
end

function process(unit::Battery, sim_params::Dict{String,Any})
    if unit.controller.state_machine.state != 2
        set_max_energy!(unit.output_interfaces[unit.medium], 0.0)    
        return
    end

    outface = unit.output_interfaces[unit.medium]
    exchanges = balance_on(outface, outface.target)

    energy_demand = balance(exchanges) +
        energy_potential(exchanges) +
        (outface.do_storage_transfer ? storage_potential(exchanges) : 0.0)

    if energy_demand >= 0.0
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
end

function load(unit::Battery, sim_params::Dict{String,Any})
    if unit.controller.state_machine.state != 1
        set_max_energy!(unit.input_interfaces[unit.medium], 0.0)
        return
    end

    inface = unit.input_interfaces[unit.medium]
    exchanges = balance_on(inface, inface.source)
    energy_available = balance(exchanges) +
        energy_potential(extanges) +
        (inface.do_storage_transfer ? storage_potential(exchanges) : 0.0)

    if energy_available <= 0.0
        set_max_energy!(unit.input_interfaces[unit.medium], 0.0)
        return # load is only concerned with receiving energy from the target
    end

    diff = unit.capacity - unit.load
    if diff > energy_available
        unit.load += energy_available
        sub!(inface, energy_available)
    else
        unit.load = unit.capacity
        sub!(inface, diff)
    end
end

function output_values(unit::Battery)::Vector{String}
    return [string(unit.medium)*" IN",
            string(unit.medium)*" OUT",
            "Load",
            "Load%",
            "Capacity",
            "Losses"]
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
    elseif key.value_key == "Losses"
        return unit.losses
    end
    throw(KeyError(key.value_key))
end

export Battery