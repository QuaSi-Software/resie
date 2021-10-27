const TIME_STEP = UInt(900)

struct Condition
    name :: String
end

Base.@kwdef struct TruthTable
    conditions :: Vector{Condition}
    table_data :: Dict{Tuple, UInt}
end

Base.@kwdef mutable struct BufferTank
    capacity :: Float64
    load :: Float64
end

Base.@kwdef mutable struct PVPlant
    amplitude :: Float64
end

Base.@kwdef mutable struct CHPP
    state :: UInt = 1
    state_names :: Dict{UInt, String}
    transitions :: Dict{UInt, TruthTable}
    time_in_state :: UInt = 0

    min_run_time :: UInt = 1800

    function CHPP(strategy::String)
        if (strategy == "Heat-driven")
            return new(
                1, # starting state

                Dict{UInt, String}( # state names
                    1 => "Off",
                    2 => "Load",
                    3 => "Produce"
                ),

                Dict{UInt, TruthTable}( # transitions
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
            )
        end
    end
end

function check(
    condition :: Condition, chpp :: CHPP,
    buffer :: BufferTank,
    pv_plant :: PVPlant,
    time :: Int,
    price_factor :: Float64
) :: Bool
    if (condition.name == "PS < 50%")
        return buffer.load < 0.5 * buffer.capacity
    elseif (condition.name == "Produce while space")
        return buffer.load < 0.95 * buffer.capacity
    elseif (condition.name == "PS < 100%")
        return buffer.load < 0.99 * buffer.capacity
    elseif (condition.name == "Min time")
        return chpp.time_in_state * TIME_STEP >= chpp.min_run_time
    elseif (condition.name == "Profitability")
        return (20000 - production(pv_plant, time)) / 20000 > price_factor
    end
end

function production(plant :: PVPlant, time :: Int) :: Float64
    seconds_in_day = 60 * 60 * 24
    return plant.amplitude * Base.Math.sin(
        Base.MathConstants.pi * (time % seconds_in_day) / seconds_in_day
    )
end

function move_state(
    chpp :: CHPP,
    buffer :: BufferTank,
    pv_plant :: PVPlant,
    time :: Int
)
    old_state = chpp.state
    table = chpp.transitions[chpp.state]
    evaluations = Tuple(check(condition, chpp, buffer, pv_plant, time, 0.5) for condition in table.conditions)
    new_state = table.table_data[evaluations]
    chpp.state = new_state
    if (old_state == new_state)
        chpp.time_in_state += 1
    else
        chpp.time_in_state = 1
    end
end

function run_simulation()
    plant = CHPP("Heat-driven")
    buffer = BufferTank(capacity=40000, load=21000)
    pv_plant = PVPlant(amplitude=30000)
    time = Int(0)
    println("Starting state is ", plant.state_names[plant.state])

    for i = 1:96
        move_state(plant, buffer, pv_plant, time)
        println("State is ", plant.state_names[plant.state])
        time += Int(TIME_STEP)
    end
end

run_simulation()
