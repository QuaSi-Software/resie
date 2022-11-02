"""
Implementation of a battery energy system holding electric charge.

For the moment the implementation remains simple with only one state (its charge) and one
parameters (its capacity). However the default operation strategy is more complex and
toggles the production of the battery dependant on available PV power and its own charge.
"""
Base.@kwdef mutable struct Battery <: ControlledSystem
    uac :: String
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
                                "threshold" => 0.15
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
                            "Little PV power",
                            Dict{String, Any}(
                                "threshold" => 0.15
                            )
                        ),
                        Condition(
                            "Sufficient charge",
                            Dict{String, Any}(
                                "threshold" => 0.05
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
    balance, _ = balance_on(outface, outface.target)

    if balance >= 0.0
        return # produce is only concerned with moving energy to the target
    end

    if unit.load > abs(balance)
        unit.load += balance
        add!(outface, abs(balance))
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
    balance, _ = balance_on(inface, inface.source)

    if balance <= 0.0
        return # load is only concerned with receiving energy from the target
    end

    diff = unit.capacity - unit.load
    if diff > balance
        unit.load += balance
        sub!(inface, balance)
    else
        unit.load = unit.capacity
        sub!(inface, diff)
    end
end

function output_values(unit :: Battery) :: Vector{String}
    return ["IN", "OUT", "Load", "Capacity"]
end

function output_value(unit :: Battery, key :: OutputKey) :: Float64
    if key.value_key == "IN"
        return unit.input_interfaces[key.medium].sum_abs_change * 0.5
    elseif key.value_key == "OUT"
        return unit.output_interfaces[key.medium].sum_abs_change * 0.5
    elseif key.value_key == "Load"
        return unit.load
    elseif key.value_key == "Capacity"
        return unit.capacity
    end
    raise(KeyError(key.value_key))
end

export Battery, make_Battery, output_values, output_value