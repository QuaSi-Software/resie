Base.@kwdef mutable struct GridConnection <: ControlledSystem
    controller :: StateMachine
    sys_function :: SystemFunction
    medium :: MediumCategory

    input_interfaces :: InterfaceMap
    output_interfaces :: InterfaceMap

    draw_sum :: Float64
    load_sum :: Float64
end

function make_GridConnection(medium :: MediumCategory, is_source) :: GridConnection
    return GridConnection(
        StateMachine(), # controller
        if is_source infinite_source else infinite_sink end, # sys_function
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

function produce(unit :: GridConnection, parameters :: Dict{String, Any}, watt_to_wh :: Function)
    if unit.sys_function === infinite_source
        outface = unit.output_interfaces[unit.medium]
        # @TODO: if grids should be allowed to load storage systems, then the potential
        # must be handled here instead of being ignored
        balance, _ = balance_on(outface, outface.target)
        unit.draw_sum += balance
        set!(outface, 0.0)
    else
        inface = unit.input_interfaces[unit.medium]
        balance, _ = balance_on(inface, inface.source)
        unit.load_sum += balance
        set!(inface, 0.0)
    end
end

function specific_values(unit :: GridConnection, time :: Int) :: Vector{Tuple}
    return [
        ("Draw sum", "$(unit.draw_sum)"),
        ("Load sum", "$(unit.load_sum)")
    ]
end

export GridConnection, specific_values, make_GridConnection