"""
Implementation of a battery component holding electric charge.

For the moment the implementation remains simple with only one state (its charge) and one
parameters (its capacity). However the default operation strategy is more complex and
toggles the processing of the battery dependant on available PV power and its own charge.
"""
Base.@kwdef mutable struct Battery <: ControlledComponent
    uac::String
    controller::Controller
    sys_function::SystemFunction

    input_interfaces::InterfaceMap
    output_interfaces::InterfaceMap

    medium::Symbol

    capacity::Float64
    load::Float64

    function Battery(uac::String, config::Dict{String,Any})
        medium = Symbol(default(config, "medium", "m_e_ac_230v"))
        register_media([medium])

        return new(
            uac, # uac
            controller_for_strategy( # controller
                config["strategy"]["name"], config["strategy"]
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
            config["load"] # load
        )
    end
end

function balance_on(
    interface::SystemInterface,
    unit::Battery
)::NamedTuple{}

    caller_is_input = unit.uac == interface.target.uac ? true : false
    # ==true if interface is input of unit (caller puts energy in unit); 
    # ==false if interface is output of unit (caller gets energy from unit)

    return (
            balance = interface.balance,
            storage_potential = caller_is_input ? -(unit.capacity-unit.load) : unit.load,
            energy_potential = 0.0,
            temperature = interface.temperature
            )
end

function process(unit::Battery, parameters::Dict{String,Any})
    if unit.controller.state_machine.state != 2
        return
    end

    outface = unit.output_interfaces[unit.medium]
    exchange = balance_on(outface, outface.target)

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

function load(unit::Battery, parameters::Dict{String,Any})
    if unit.controller.state_machine.state != 1
        return
    end

    inface = unit.input_interfaces[unit.medium]
    exchange = balance_on(inface, inface.source)
    energy_available = exchange.balance

    if energy_available <= 0.0
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
    return ["IN", "OUT", "Load", "Load%", "Capacity"]
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
    end
    throw(KeyError(key.value_key))
end

export Battery