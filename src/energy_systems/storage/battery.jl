"""
Implementation of a battery energy system holding electric charge.

For the moment the implementation remains simple with only one state (its charge) and one
parameters (its capacity). However the default operation strategy is more complex and
toggles the production of the battery dependant on available PV power and its own charge.
"""
Base.@kwdef mutable struct Battery <: ControlledSystem
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

    caller_is_input = false   # ==true if interface is input of unit (caller puts energy in unit); 
                              # ==false if interface is output of unit (caller gets energy from unit)

    # check if caller is input or output of unit
    for (_, input_uac) in pairs(unit.input_interfaces)
        if input_uac == interface.source.uac
            caller_is_input = true
            break
        end
        if input_uac.source.uac == interface.source.uac
            caller_is_input = true
            break
        end
    end

    return (
            balance = interface.balance,
            storage_potential = caller_is_input ? -(unit.capacity-unit.load) : unit.load,
            energy_potential = 0.0,
            temperature = interface.temperature
            )
end

function produce(unit::Battery, parameters::Dict{String,Any}, watt_to_wh::Function)
    if unit.controller.state_machine.state != 2
        return
    end

    outface = unit.output_interfaces[unit.medium]
    InterfaceInfo = balance_on(outface, outface.target)

    if InterfaceInfo.balance >= 0.0
        return # produce is only concerned with moving energy to the target
    end

    if unit.load > abs(InterfaceInfo.balance)
        unit.load += InterfaceInfo.balance
        add!(outface, abs(InterfaceInfo.balance))
    else
        add!(outface, unit.load)
        unit.load = 0.0
    end
end

function load(unit::Battery, parameters::Dict{String,Any}, watt_to_wh::Function)
    if unit.controller.state_machine.state != 1
        return
    end

    inface = unit.input_interfaces[unit.medium]
    InterfaceInfo = balance_on(inface, inface.source)

    if InterfaceInfo.balance <= 0.0
        return # load is only concerned with receiving energy from the target
    end

    diff = unit.capacity - unit.load
    if diff > InterfaceInfo.balance
        unit.load += InterfaceInfo.balance
        sub!(inface, InterfaceInfo.balance)
    else
        unit.load = unit.capacity
        sub!(inface, diff)
    end
end

function output_values(unit::Battery)::Vector{String}
    return ["IN", "OUT", "Load", "Capacity"]
end

function output_value(unit::Battery, key::OutputKey)::Float64
    if key.value_key == "IN"
        return calculate_energy_flow(unit.input_interfaces[key.medium])
    elseif key.value_key == "OUT"
        return calculate_energy_flow(unit.output_interfaces[key.medium])
    elseif key.value_key == "Load"
        return unit.load
    elseif key.value_key == "Capacity"
        return unit.capacity
    end
    throw(KeyError(key.value_key))
end

export Battery