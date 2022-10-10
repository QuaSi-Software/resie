Base.@kwdef mutable struct GridConnection <: ControlledSystem
    controller :: StateMachine
    is_storage :: Bool
    medium :: MediumCategory

    input_interfaces :: InterfaceMap
    output_interfaces :: InterfaceMap

    is_source :: Bool

    draw_sum :: Float64
    load_sum :: Float64
end

function make_GridConnection(medium :: MediumCategory, is_source) :: GridConnection
    return GridConnection(
        StateMachine(), # controller
        false, # is_storage
        medium, # medium
        InterfaceMap( # input_interfaces
            medium => nothing
        ),
        InterfaceMap( # output_interfaces
            medium => nothing
        ),
        is_source, # is_source
        0.0, # draw_sum,
        0.0, # load_sum
    )
end

function produce(unit :: GridConnection, parameters :: Dict{String, Any}, watt_to_wh :: Function)
    if unit.is_source
        outface = unit.output_interfaces[unit.medium]
        gather_from_all!(outface, outface.target)
        unit.draw_sum += outface.balance
        set!(outface, 0.0)
    else
        inface = unit.input_interfaces[unit.medium]
        gather_from_all!(inface, inface.source)
        unit.load_sum += inface.balance
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