"""
Implementation of an energy system modeling the connection to a public grid of a certain medium.

Public grids are considered to have an unlimited amount of energy they can provide, but
might be limited in the power they can provide, although this behaviour is not yet
implemented. They exist to model real connections to a public grid that can provide any
remaining demand of energy or take in any excess of energy. To make it possible to model a
one-way connection they are split into two instances for providing or receiving energy and
must be handled as such in the input for constructing a project.
"""
mutable struct GridConnection <: ControlledSystem
    uac::String
    controller::Controller
    sys_function::SystemFunction
    medium::Symbol

    input_interfaces::InterfaceMap
    output_interfaces::InterfaceMap

    draw_sum::Float64
    load_sum::Float64

    function GridConnection(uac::String, config::Dict{String,Any})
        medium = Symbol(config["medium"])
        register_media([medium])

        return new(
            uac, # uac
            controller_for_strategy( # controller
                config["strategy"]["name"], config["strategy"]
            ),
            if Bool(config["is_source"])
                sf_dispatchable_source
            else
                sf_dispatchable_sink
            end, # sys_function
            medium, # medium
            InterfaceMap( # input_interfaces
                medium => nothing
            ),
            InterfaceMap( # output_interfaces
                medium => nothing
            ),
            0.0, # draw_sum,
            0.0, # load_sum
        )
    end
end

function control(
    unit::GridConnection,
    systems::Grouping,
    parameters::Dict{String,Any}
)
    move_state(unit, systems, parameters)
    if unit.sys_function === sf_dispatchable_source
        set_max_energy!(unit.output_interfaces[unit.medium], Inf)
    else
        set_max_energy!(unit.input_interfaces[unit.medium], -Inf)
    end
end

function produce(unit::GridConnection, parameters::Dict{String,Any}, watt_to_wh::Function)
    if unit.sys_function === sf_dispatchable_source
        outface = unit.output_interfaces[unit.medium]
        # @TODO: if grids should be allowed to load storage systems, then the potential
        # must be handled here instead of being ignored
        InterfaceInfo = balance_on(outface, outface.target)
        if InterfaceInfo.balance < 0.0
            unit.draw_sum += InterfaceInfo.balance
            add!(outface, abs(InterfaceInfo.balance))
        end
    else
        inface = unit.input_interfaces[unit.medium]
        InterfaceInfo = balance_on(inface, inface.source)
        if InterfaceInfo.balance > 0.0
            unit.load_sum += InterfaceInfo.balance
            sub!(inface, InterfaceInfo.balance)
        end
    end
end

function output_values(unit::GridConnection)::Vector{String}
    return ["IN", "OUT", "Draw sum", "Load sum"]
end

function output_value(unit::GridConnection, key::OutputKey)::Float64
    if key.value_key == "IN"
        return calculate_energy_flow(unit.input_interfaces[key.medium])
    elseif key.value_key == "OUT"
        return calculate_energy_flow(unit.output_interfaces[key.medium])
    elseif key.value_key == "Draw sum"
        return unit.draw_sum
    elseif key.value_key == "Load sum"
        return unit.load_sum
    end
    throw(KeyError(key.value_key))
end

export GridConnection