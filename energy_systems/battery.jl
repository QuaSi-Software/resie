Base.@kwdef mutable struct Battery <: ControlledSystem
    controller :: StateMachine
    sys_function :: SystemFunction

    input_interfaces :: InterfaceMap
    output_interfaces :: InterfaceMap

    capacity :: Float64
    load :: Float64
end

function make_Battery(strategy :: String, capacity :: Float64, load :: Float64) :: Battery
    if strategy == "Economical discharge"
        controller = StateMachine(
            state=UInt(1),
            state_names=Dict{UInt, String}(
                1 => "Load",
                2 => "Discharge",
            ),
            time_in_state=UInt(0),
            transitions=Dict{UInt, TruthTable}(
                1 => TruthTable( # State: Load
                    conditions=[
                        Condition(
                            "Little PV power",
                            Dict{String, Any}(
                                "threshold" => capacity * 0.05
                            )
                        ),
                        Condition(
                            "Sufficient charge",
                            Dict{String, Any}(
                                "threshold" => 0.2
                            )
                        )
                    ],
                    table_data=Dict{Tuple, UInt}(
                        (false,false) => 1,
                        (false,true) => 1,
                        (true,false) => 1,
                        (true,true) => 2,
                    )
                ),

                2 => TruthTable( # State: Discharge
                    conditions=[
                        Condition(
                            "Much PV power",
                            Dict{String, Any}(
                                "threshold" => capacity * 0.05
                            )
                        ),
                        Condition(
                            "Insufficient charge",
                            Dict{String, Any}(
                                "threshold" => 0.05
                            )
                        )
                    ],
                    table_data=Dict{Tuple, UInt}(
                            (false,false) => 2,
                            (false,true) => 1,
                            (true,false) => 1,
                            (true,true) => 1,
                        )
                ),
            )
        )
    else
        controller = StateMachine()
    end

    return Battery(
        controller, # controller
        storage, # sys_function
        InterfaceMap( # input_interfaces
            m_e_ac_230v => nothing
        ),
        InterfaceMap( # output_interfaces
            m_e_ac_230v => nothing
        ),
        capacity, # capacity
        load # load
    )
end

function produce(unit :: Battery, parameters :: Dict{String, Any}, watt_to_wh :: Function)
    if unit.controller.state != 2
        return
    end

    outface = unit.output_interfaces[m_e_ac_230v]
    balance = check_balance(outface, outface.target)

    if balance >= 0.0
        return # produce is only concerned with moving energy to the target
    end

    # first we only check if there is a balance that produce needs to handle
    # without already gathering energy in the output interface. otherwise,
    # when the balance is positive, no energy is moved but the act of calling
    # gather_from_all! has incorrectly recorded a move of energy
    gather_from_all!(outface, outface.target)

    if unit.load > outface.balance
        unit.load += outface.balance
        set!(outface, 0.0)
    else
        add!(outface, unit.load)
        unit.load = 0.0
    end
end

function load(unit :: Battery, parameters :: Dict{String, Any}, watt_to_wh :: Function)
    if unit.controller.state != 1
        return
    end

    inface = unit.input_interfaces[m_e_ac_230v]
    gather_from_all!(inface, inface.source)

    if inface.balance <= 0.0
        return # load is only concerned with receiving energy from the target
    end

    unit.load += inface.balance # @TODO: check if loading exceeds capacity
    set!(inface, 0.0)
end

function specific_values(unit :: Battery, time :: Int) :: Vector{Tuple}
    return [
        ("Load", "$(unit.load)"),
        ("Capacity", "$(unit.capacity)")
    ]
end

export Battery, specific_values, make_Battery