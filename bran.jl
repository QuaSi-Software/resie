const TIME_STEP = UInt(900)

abstract type EnergySystem end
abstract type ControlledSystem <: EnergySystem end

Base.@kwdef struct Condition
    name :: String
end

Base.@kwdef struct TruthTable
    conditions :: Vector{Condition}
    table_data :: Dict{Tuple, UInt}
end

Base.@kwdef mutable struct StateMachine
    state :: UInt
    state_names :: Dict{UInt, String}
    transitions :: Dict{UInt, TruthTable}
    time_in_state :: UInt
end

StateMachine() = StateMachine(
    1, # StateMachine.state
    Dict(1=>"Default"), # StateMachine.state_names
    Dict(1=>TruthTable( # StateMachine.transitions
        conditions=Vector(),
        table_data=Dict()
    )),
    0 # StateMachine.time_in_state
)

Base.@kwdef mutable struct BufferTank <: ControlledSystem
    controller :: StateMachine = StateMachine()

    capacity :: Float64
    load :: Float64
end

Base.@kwdef mutable struct PVPlant <: ControlledSystem
    controller :: StateMachine = StateMachine()

    amplitude :: Float64
end

Base.@kwdef mutable struct CHPP <: ControlledSystem
    controller :: StateMachine

    min_run_time :: UInt = 1800

    function CHPP(strategy::String)
        if strategy == "Heat-driven"
            return new(
                StateMachine( # CHPP.controller
                    state=UInt(1),
                    state_names=Dict{UInt, String}(
                        1 => "Off",
                        2 => "Load",
                        3 => "Produce"
                    ),
                    time_in_state=UInt(0),
                    transitions=Dict{UInt, TruthTable}(
                        1 => TruthTable( # State: Off
                            conditions=[
                                Condition("PS < 50%"),
                                Condition("Produce while space")
                            ],
                            table_data=Dict{Tuple, UInt}(
                                (false, false) => 1,
                                (false, true) => 3,
                                (true, false) => 2,
                                (true, true) => 1
                            )
                        ),

                        2 => TruthTable( # State: Load
                            conditions=[
                                Condition("PS < 100%")
                            ],
                            table_data=Dict{Tuple, UInt}(
                                (false,) => 1,
                                (true,) => 2
                            )
                        ),

                        3 => TruthTable( # State: Produce
                            conditions=[
                                Condition("Min time"),
                                Condition("PS < 100%"),
                                Condition("Profitability")
                            ],
                            table_data=Dict{Tuple, UInt}(
                                (false, false, false) => 1,
                                (false, false, true) => 1,
                                (false, true, false) => 3,
                                (true, false, false) => 1,
                                (false, true, true) => 3,
                                (true, false, true) => 1,
                                (true, true, false) => 1,
                                (true, true, true) => 3,
                            )
                        ),
                    )
                ),
                1800 # CHPP.min_run_time
            )
        end
    end
end

function check(
    condition :: Condition,
    unit :: ControlledSystem,
    system :: Vector{ControlledSystem},
    parameters :: Dict{String, Any}
) :: Bool
    buffer = [u for u in system if typeof(u) <: BufferTank][1]
    chpp = [u for u in system if typeof(u) <: CHPP][1]
    pv_plant = [u for u in system if typeof(u) <: PVPlant][1]

    if condition.name == "PS < 50%"
        return buffer.load < 0.5 * buffer.capacity
    elseif condition.name == "Produce while space"
        return buffer.load < 0.95 * buffer.capacity
    elseif condition.name == "PS < 100%"
        return buffer.load < 0.99 * buffer.capacity
    elseif condition.name == "Min time"
        return chpp.controller.time_in_state * TIME_STEP >= chpp.min_run_time
    elseif condition.name == "Profitability"
        return (20000 - production(pv_plant, parameters["time"])) / 20000 > parameters["price_factor"]
    end
end

function production(plant :: PVPlant, time :: Int) :: Float64
    seconds_in_day = 60 * 60 * 24
    return plant.amplitude * Base.Math.sin(
        Base.MathConstants.pi * (time % seconds_in_day) / seconds_in_day
    )
end

function move_state(
    unit :: ControlledSystem,
    system :: Vector{ControlledSystem},
    parameters :: Dict{String, Any}
)
    old_state = unit.controller.state
    table = unit.controller.transitions[unit.controller.state]

    if length(table.conditions) > 0
        evaluations = Tuple(
            check(condition, unit, system, parameters)
            for condition in table.conditions
        )
        new_state = table.table_data[evaluations]
        unit.controller.state = new_state
    else
        new_state = old_state
    end

    if old_state == new_state
        unit.controller.time_in_state += 1
    else
        unit.controller.time_in_state = 1
    end
end

function run_simulation()
    system = [
        CHPP("Heat-driven"),
        BufferTank(capacity=40000.0, load=21000.0),
        PVPlant(amplitude=30000.0)
    ]
    parameters = Dict{String, Any}(
        "time" => 0,
        "price_factor" => 0.5
    )
    plant = [u for u in system if typeof(u) <: CHPP][1]
    println("Starting state is ", plant.controller.state_names[plant.controller.state])

    for i = 1:96
        for unit in system
            move_state(unit, system, parameters)
            println("State is ", unit.controller.state_names[unit.controller.state])
        end
        parameters["time"] += Int(TIME_STEP)
    end
end

run_simulation()
