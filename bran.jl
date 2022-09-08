const TIME_STEP = UInt(900)

@enum MediumCategory m_e_ac_230v m_c_g_natgas m_h_w_60c

abstract type EnergySystem end
abstract type ControlledSystem <: EnergySystem end

function Wh(watts :: Float64) :: Float64
    return Float64(TIME_STEP) * watts / 3600.0
end

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

Base.@kwdef mutable struct Demand <: ControlledSystem
    controller :: StateMachine = StateMachine()
    medium :: MediumCategory

    load :: Float64
end

function specific_values(unit :: Demand) :: String
    return "$(unit.load)"
end

Base.@kwdef mutable struct GridConnection <: ControlledSystem
    controller :: StateMachine = StateMachine()
    medium :: MediumCategory

    draw_sum :: Float64 = 0.0
    load_sum :: Float64 = 0.0
end

function specific_values(unit :: GridConnection) :: String
    return "$(unit.draw_sum)/$(unit.load_sum)"
end

Base.@kwdef mutable struct Bus <: ControlledSystem
    controller :: StateMachine = StateMachine()
    medium :: MediumCategory

    balance :: Float64 = 0.0
end

function specific_values(unit :: Bus) :: String
    return "$(unit.balance)"
end

Base.@kwdef mutable struct BufferTank <: ControlledSystem
    controller :: StateMachine = StateMachine()

    capacity :: Float64
    load :: Float64
end

function specific_values(unit :: BufferTank) :: String
    return "$(unit.load)/$(unit.capacity)"
end

Base.@kwdef mutable struct PVPlant <: ControlledSystem
    controller :: StateMachine = StateMachine()

    last_produced_e :: Float64 = 0.0

    amplitude :: Float64
end

function specific_values(unit :: PVPlant) :: String
    return "$(unit.last_produced_e)"
end

Base.@kwdef mutable struct CHPP <: ControlledSystem
    controller :: StateMachine

    last_produced_e :: Float64 = 0.0
    last_produced_h :: Float64 = 0.0

    power :: Float64
    electricity_fraction :: Float64 = 0.4
    min_power_fraction :: Float64
    min_run_time :: UInt = 1800
end

function make_CHPP(strategy :: String, power :: Float64) :: CHPP
    if strategy == "Ensure storage"
        return CHPP(
            StateMachine( # CHPP.controller
                state=UInt(1),
                state_names=Dict{UInt, String}(
                    1 => "Off",
                    2 => "Load",
                ),
                time_in_state=UInt(0),
                transitions=Dict{UInt, TruthTable}(
                    1 => TruthTable( # State: Off
                        conditions=[
                            Condition("PS < 20%"),
                        ],
                        table_data=Dict{Tuple, UInt}(
                            (false,) => 1,
                            (true,) => 2,
                        )
                    ),

                    2 => TruthTable( # State: Load
                        conditions=[
                            Condition("PS >= 90%"),
                            Condition("Min time"),
                            Condition("Would overfill CHPP"),
                        ],
                        table_data=Dict{Tuple, UInt}(
                            (false, false, false) => 2,
                            (false, true, false) => 2,
                            (true, false, false) => 2,
                            (true, true, false) => 1,
                            (false, false, true) => 1,
                            (false, true, true) => 1,
                            (true, false, true) => 1,
                            (true, true, true) => 1,
                        )
                    ),
                )
            ),
            0.0, # CHPP.last_produced_e
            0.0, # CHPP.last_produced_h
            power, # CHPP.power
            0.4, # CHPP.electricity_fraction
            0.2, # CHPP.min_power_fraction
            1800 # CHPP.min_run_time
        )
    else
        return CHPP(controller=StateMachine(), power=power)
    end
end

function specific_values(unit :: CHPP) :: String
    return "$(unit.last_produced_e)/$(unit.last_produced_h)"
end

Base.@kwdef mutable struct HeatPump <: ControlledSystem
    controller :: StateMachine = StateMachine()

    last_consumed_e :: Float64 = 0.0
    last_produced_h :: Float64 = 0.0

    power :: Float64
    min_power_fraction :: Float64 = 0.2
    cop :: Float64
end

function make_HeatPump(strategy :: String, power :: Float64, cop :: Float64) :: HeatPump
    if strategy == "Ensure storage"
        return HeatPump(
            StateMachine( # HeatPump.controller
                state=UInt(1),
                state_names=Dict{UInt, String}(
                    1 => "Off",
                    2 => "Load"
                ),
                time_in_state=UInt(0),
                transitions=Dict{UInt, TruthTable}(
                    1 => TruthTable( # State: Off
                        conditions=[
                            Condition("PS < 10%")
                        ],
                        table_data=Dict{Tuple, UInt}(
                            (true,) => 2,
                            (false,) => 1
                        )
                    ),

                    1 => TruthTable( # State: Load
                        conditions=[
                            Condition("PS >= 50%"),
                            Condition("Would overfill HP")
                        ],
                        table_data=Dict{Tuple, UInt}(
                            (false, false) => 2,
                            (false, true) => 1,
                            (true, false) => 1,
                            (true, true) => 1,
                        )
                    ),
                )
            ),
            0.0, # HeatPump.last_consumed_e
            0.0, # HeatPump.last_produced_h
            power, # HeatPump.power
            0.2, # HeatPump.min_power_fraction
            cop, # HeatPump.electricity_fraction
        )
    else
        return HeatPump(controller=StateMachine(), power=power, cop=cop)
    end
end

function specific_values(unit :: HeatPump) :: String
    return "$(unit.last_consumed_e)/$(unit.last_produced_h)"
end

function represent(unit :: ControlledSystem) :: String
    return string(
        "$(typeof(unit)) ($(unit.controller.state_names[unit.controller.state])) ",
        "$(specific_values(unit))"
    )
end

pprint(unit :: ControlledSystem) = print(represent(unit))

function check(
    condition :: Condition,
    unit :: ControlledSystem,
    system :: Vector{ControlledSystem},
    parameters :: Dict{String, Any}
) :: Bool
    buffer = [u for u in system if typeof(u) <: BufferTank][1]
    chpp = [u for u in system if typeof(u) <: CHPP][1]
    hp = [u for u in system if typeof(u) <: HeatPump][1]
    pv_plant = [u for u in system if typeof(u) <: PVPlant][1]

    if condition.name == "PS < 20%"
        return buffer.load < 0.2 * buffer.capacity
    elseif condition.name == "PS >= 90%"
        return buffer.load >= 0.9 * buffer.capacity
    elseif condition.name == "PS < 10%"
        return buffer.load < 0.1 * buffer.capacity
    elseif condition.name == "PS >= 50%"
        return buffer.load >= 0.5 * buffer.capacity
    elseif condition.name == "Min time"
        return chpp.controller.time_in_state * TIME_STEP >= chpp.min_run_time
    elseif condition.name == "Would overfill CHPP"
        return buffer.capacity - buffer.load < chpp.power * chpp.min_power_fraction
    elseif condition.name == "Would overfill HP"
        return buffer.capacity - buffer.load < hp.power * hp.min_power_fraction
    end
end

function produce(
    system :: Vector{ControlledSystem},
    parameters :: Dict{String, Any}
)
    grid_e = [u for u in system if (typeof(u) <: GridConnection && u.medium == m_e_ac_230v)][1]
    demand_h = [u for u in system if (typeof(u) <: Demand && u.medium == m_h_w_60c)][1]
    demand_e = [u for u in system if (typeof(u) <: Demand && u.medium == m_e_ac_230v)][1]
    buffer = [u for u in system if typeof(u) <: BufferTank][1]
    e_bus = [u for u in system if (typeof(u) <: Bus && u.medium == m_e_ac_230v)][1]
    chpp = [u for u in system if typeof(u) <: CHPP][1]
    hp = [u for u in system if typeof(u) <: HeatPump][1]
    pv_plant = [u for u in system if typeof(u) <: PVPlant][1]

    # reset balances
    e_bus.balance = 0.0
    chpp.last_produced_e = 0.0
    chpp.last_produced_h = 0.0
    hp.last_consumed_e = 0.0
    hp.last_produced_h = 0.0
    pv_plant.last_produced_e = 0.0

    # unload buffer
    buffer.load -= Wh(demand_h.load)

    # run chpp
    if chpp.controller.state == 2
        space_in_buffer = buffer.capacity - buffer.load
        max_produce_h = Wh(chpp.power * (1.0 - chpp.electricity_fraction))
        max_produce_e = Wh(chpp.power * chpp.electricity_fraction)

        if space_in_buffer > max_produce_h
            usage_fraction = 1.0
        else
            usage_fraction = space_in_buffer / max_produce_h
        end

        buffer.load += max_produce_h * usage_fraction
        e_bus.balance += max_produce_e * usage_fraction
        chpp.last_produced_e = chpp.power * chpp.electricity_fraction * usage_fraction
        chpp.last_produced_h = chpp.power * (1.0 - chpp.electricity_fraction) * usage_fraction
    end

    # electricity demand and PV plant
    e_bus.balance -= Wh(demand_e.load)
    e_bus.balance += Wh(production(pv_plant, parameters["time"]))
    pv_plant.last_produced_e = Wh(production(pv_plant, parameters["time"]))

    # run heat pump
    if hp.controller.state == 2
        space_in_buffer = buffer.capacity - buffer.load
        max_produce_h = Wh(hp.power)

        if space_in_buffer > max_produce_h
            usage_fraction = 1.0
        else
            usage_fraction = space_in_buffer / max_produce_h
        end

        buffer.load += max_produce_h * usage_fraction
        e_bus.balance -= max_produce_h * usage_fraction / hp.cop
        hp.last_consumed_e = max_produce_h * usage_fraction / hp.cop
        hp.last_produced_h = max_produce_h * usage_fraction
    end

    # balance electricity bus
    if e_bus.balance >= 0
        grid_e.load_sum += e_bus.balance
    else
        grid_e.draw_sum += e_bus.balance
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

function print_system_state(system :: Vector{ControlledSystem}, time :: Int)
    println("Time is ", time)
    for unit in system
        pprint(unit)
        print(" | ")
    end
    print("\n")
end

function reset_file()
    open("./out.csv", "w") do file_handle
        write(file_handle, "")
    end
end

function write_to_file(system :: Vector{ControlledSystem}, time :: Int)
    open("./out.csv", "a") do file_handle
        write(file_handle, "$time")
        for unit in system
            write(file_handle, replace(
                replace(";$(specific_values(unit))", "/" => ";"),
                "." => ","
            ))
        end
        write(file_handle, "\n")
    end
end

function run_simulation()
    system = [
        GridConnection(medium=m_c_g_natgas),
        GridConnection(medium=m_e_ac_230v),
        make_CHPP("Ensure storage", 20000.0),
        make_HeatPump("Ensure storage", 20000.0, 3.0),
        BufferTank(capacity=40000.0, load=20000.0),
        PVPlant(amplitude=20000.0),
        Bus(medium=m_e_ac_230v),
        Demand(medium=m_h_w_60c, load=10000),
        Demand(medium=m_e_ac_230v, load=20000),
    ]
    parameters = Dict{String, Any}(
        "time" => 0,
        "price_factor" => 0.5
    )

    print_system_state(system, parameters["time"])
    reset_file()

    for i = 1:96
        for unit in system
            # control
            move_state(unit, system, parameters)
        end

        # production
        produce(system, parameters)

        # output and simulation update
        print_system_state(system, parameters["time"])
        write_to_file(system, parameters["time"])
        parameters["time"] += Int(TIME_STEP)
    end
end

run_simulation()
