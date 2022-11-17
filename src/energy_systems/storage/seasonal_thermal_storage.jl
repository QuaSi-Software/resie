"""
Implementation of a seasonal thermal storage system.

For the moment this remains a very simple implementation that only has a load (how much
energy is stored) and a capacity, with no temperatures being considered. Given how the
simulation engine works, there will likely always be the need to deal with energy being
transfered with water temperature being secondary input variables.
"""
mutable struct SeasonalThermalStorage <: ControlledSystem
    uac :: String
    controller :: Controller
    sys_function :: SystemFunction

    input_interfaces :: InterfaceMap
    output_interfaces :: InterfaceMap

    capacity :: Float64
    load :: Float64

    function SeasonalThermalStorage(uac :: String, config :: Dict{String, Any})
        if config["strategy"]["name"] == "No cycle loading"
            strategy = config["strategy"]["name"]

            machine = StateMachine(
                state=UInt(1),
                state_names=Dict{UInt, String}(
                    1 => "Load",
                    2 => "Produce",
                ),
                time_in_state=UInt(0),
                transitions=Dict{UInt, TruthTable}(
                    1 => TruthTable( # State: Load
                        conditions=[
                            Condition(
                                "Buffer >= X%",
                                Dict{String, Any}(
                                    "percentage" => config["strategy"]["percentage"]
                                )
                            ),
                            Condition(
                                "HP is running",
                                Dict{String, Any}()
                            )
                        ],
                        table_data=Dict{Tuple, UInt}(
                            (false,false) => 2,
                            (false,true) => 2,
                            (true,false) => 1,
                            (true,true) => 2,
                        )
                    ),

                    2 => TruthTable( # State: Produce
                        conditions=[
                            Condition(
                                "Buffer >= X%",
                                Dict{String, Any}(
                                    "percentage" => config["strategy"]["percentage"]
                                )
                            ),
                            Condition(
                                "HP is running",
                                Dict{String, Any}()
                            )
                        ],
                        table_data=Dict{Tuple, UInt}(
                                (false,false) => 2,
                                (false,true) => 2,
                                (true,false) => 1,
                                (true,true) => 2,
                            )
                    ),
                )
            )
        else
            strategy = "Default"
            machine = StateMachine()
        end

        return new(
            uac, # uac
            Controller(strategy, machine), # controller
            sf_storage, # sys_function
            InterfaceMap( # input_interfaces
                m_h_w_ht1 => nothing
            ),
            InterfaceMap( # output_interfaces
                m_h_w_lt1 => nothing
            ),
            config["capacity"], # capacity
            config["load"] # load
        )
    end
end

function produce(unit :: SeasonalThermalStorage, parameters :: Dict{String, Any}, watt_to_wh :: Function)
    outface = unit.output_interfaces[m_h_w_lt1]
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

function load(unit :: SeasonalThermalStorage, parameters :: Dict{String, Any}, watt_to_wh :: Function)
    inface = unit.input_interfaces[m_h_w_ht1]
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

function output_values(unit :: SeasonalThermalStorage) :: Vector{String}
    return ["IN", "OUT", "Load", "Capacity"]
end

function output_value(unit :: SeasonalThermalStorage, key :: OutputKey) :: Float64
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

export SeasonalThermalStorage