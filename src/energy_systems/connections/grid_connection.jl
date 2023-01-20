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
    uac :: String
    controller :: Controller
    sys_function :: SystemFunction
    medium :: MediumCategory

    input_interfaces :: InterfaceMap
    output_interfaces :: InterfaceMap

    draw_sum :: Float64
    load_sum :: Float64

    function GridConnection(uac :: String, config :: Dict{String, Any})
        medium = getproperty(EnergySystems, Symbol(config["medium"]))
        return new(
            uac, # uac
            controller_for_strategy( # controller
                config["strategy"]["name"], config["strategy"]
            ),
            if Bool(config["is_source"]) sf_dispatchable_source else sf_dispatchable_sink end, # sys_function
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

function produce(unit :: GridConnection, parameters :: Dict{String, Any}, watt_to_wh :: Function)
    if unit.sys_function === sf_dispatchable_source
        outface = unit.output_interfaces[unit.medium]
        # @TODO: if grids should be allowed to load storage systems, then the potential
        # must be handled here instead of being ignored
        balance, _, _ = balance_on(outface, outface.target)
        if balance < 0.0
            unit.draw_sum += balance
            add!(outface, abs(balance))
        end
    else
        inface = unit.input_interfaces[unit.medium]
        balance, _, _ = balance_on(inface, inface.source)
        if balance > 0.0
            unit.load_sum += balance
            sub!(inface, balance)
        end
    end
end

function output_values(unit :: GridConnection) :: Vector{String}
    return ["IN", "OUT", "Draw sum", "Load sum"]
end

function output_value(unit :: GridConnection, key :: OutputKey) :: Float64
    if key.value_key == "IN"
        return unit.input_interfaces[key.medium].sum_abs_change * 0.5
    elseif key.value_key == "OUT"
        return unit.output_interfaces[key.medium].sum_abs_change * 0.5
    elseif key.value_key == "Draw sum"
        return unit.draw_sum
    elseif key.value_key == "Load sum"
        return unit.load_sum
    end
    throw(KeyError(key.value_key))
end

export GridConnection