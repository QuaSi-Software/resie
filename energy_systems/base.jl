module EnergySystems

export MediumCategory, EnergySystem, ControlledSystem, Condition, TruthTable, StateMachine,
    control, represent, pprint, check, produce, production, link_control_with, each, Grouping,
    link_production_with

const TIME_STEP = UInt(900)

@enum MediumCategory m_e_ac_230v m_c_g_natgas m_h_w_60c

abstract type EnergySystem end
abstract type ControlledSystem <: EnergySystem end

const Grouping = Dict{String, ControlledSystem}

function each(systems :: Grouping) :: Base.ValueIterator
    return values(systems)
end

Base.@kwdef mutable struct SystemInterface
    left :: Union{Nothing, ControlledSystem} = nothing
    right :: Union{Nothing, ControlledSystem} = nothing
    balance :: Float64 = 0.0
end

const InterfaceMap = Dict{MediumCategory, Union{Nothing, SystemInterface}}

function gather_from_all!(interface :: SystemInterface, unit :: ControlledSystem)
    return # the default implementation is to do nothing
end

function reset(unit :: ControlledSystem)
    for inface in values(unit.input_interfaces)
        if inface !== nothing
            inface.balance = 0.0
        end
    end
    for outface in values(unit.output_interfaces)
        if outface !== nothing
            outface.balance = 0.0
        end
    end
end

struct Condition
    name :: String
    parameters :: Dict{String, Any}
    required_systems :: Dict{String, Tuple{Type, MediumCategory}}
    linked_systems :: Grouping
end

function rel(condition :: Condition, name :: String) :: ControlledSystem
    return condition.linked_systems[name]
end

function Condition(
    name :: String,
    parameters :: Dict{String, Any}
) :: Condition
    required_systems = Dict{String, Tuple{Type, MediumCategory}}()
    default_params = Dict{String, Any}()

    if name == "Buffer < X%"
        required_systems["buffer"] = (BufferTank, m_h_w_60c)
        default_params["percentage"] = 0.5
    elseif name == "Buffer >= X%"
        required_systems["buffer"] = (BufferTank, m_h_w_60c)
        default_params["percentage"] = 0.5
    elseif name == "Min run time"
        # nothing to do
    elseif name == "Would overfill thermal buffer"
        required_systems["buffer"] = (BufferTank, m_h_w_60c)
    else
        throw(KeyError(name))
    end

    return Condition(
        name,
        merge(default_params, parameters),
        required_systems,
        Grouping()
    )
end

function link(condition :: Condition, systems :: Grouping)
    for (name, req_unit) in pairs(condition.required_systems)
        found_link = false
        for unit in each(systems)
            if isa(unit, req_unit[1]) && req_unit[2] in keys(unit.input_interfaces)
                condition.linked_systems[name] = unit
                found_link = true
            end
        end

        if !found_link
            throw(KeyError("Could not find match for required system $name "
                * "for condition $(condition.name)"))
        end
    end
end

function link_control_with(unit :: ControlledSystem, systems :: Grouping)
    for table in values(unit.controller.transitions)
        for condition in table.conditions
            link(condition, systems)
        end
    end
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

function represent(unit :: ControlledSystem, time :: Int) :: String
    repr = ""
    for val in specific_values(unit, time)
        repr = repr * " $(val[1]) $(val[2])"
    end
    return string(
        "$(typeof(unit)) ",
        "($(unit.controller.state_names[unit.controller.state]))",
        repr
    )
end

pprint(unit :: ControlledSystem, time :: Int) = print(represent(unit, time))

function move_state(
    unit :: ControlledSystem,
    systems :: Grouping,
    parameters :: Dict{String, Any}
)
    old_state = unit.controller.state
    table = unit.controller.transitions[unit.controller.state]

    if length(table.conditions) > 0
        evaluations = Tuple(
            check(condition, unit, parameters)
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

function control(
    unit :: ControlledSystem,
    systems :: Grouping,
    parameters :: Dict{String, Any}
)
    move_state(unit, systems, parameters)
end

function Wh(watts :: Float64) :: Float64
    return Float64(TIME_STEP) * watts / 3600.0
end

include("buffer_tank.jl")
include("bus.jl")
include("chpp.jl")
include("demand.jl")
include("grid_connection.jl")
include("heat_pump.jl")
include("pv_plant.jl")

function check(
    condition :: Condition,
    unit :: ControlledSystem,
    parameters :: Dict{String, Any}
) :: Bool
    if condition.name == "Buffer < X%"
        return (rel(condition, "buffer").load
            < condition.parameters["percentage"]
            * rel(condition, "buffer").capacity)
    elseif condition.name == "Buffer >= X%"
        return (rel(condition, "buffer").load
            >= condition.parameters["percentage"]
            * rel(condition, "buffer").capacity)
    elseif condition.name == "Min run time"
        return (unit.controller.time_in_state
            * parameters["time_step_seconds"]
            >= unit.min_run_time)
    elseif condition.name == "Would overfill thermal buffer"
        return (rel(condition, "buffer").capacity
            - rel(condition, "buffer").load
            < unit.power * unit.min_power_fraction)
    end
    throw(KeyError(condition.name))
end

function link_production_with(unit :: ControlledSystem, systems :: Grouping)
    if isa(unit, Bus)
        for system in each(systems)
            if isa(system, Bus)
                if out_medium == system.medium
                    connection = SystemInterface(left=unit, right=system)
                    push!(system.input_interfaces, connection)
                    push!(unit.output_interfaces, connection)
                end
            else
                for in_medium in keys(system.input_interfaces)
                    if in_medium == unit.medium
                        connection = SystemInterface(left=unit, right=system)
                        system.input_interfaces[in_medium] = connection
                        push!(unit.output_interfaces, connection)
                    end
                end
            end
        end
    else
        for out_medium in keys(unit.output_interfaces)
            for system in each(systems)
                if isa(system, Bus)
                    if out_medium == system.medium
                        connection = SystemInterface(left=unit, right=system)
                        push!(system.input_interfaces, connection)
                        unit.output_interfaces[out_medium] = connection
                    end
                else
                    for in_medium in keys(system.input_interfaces)
                        if out_medium == in_medium
                            connection = SystemInterface(left=unit, right=system)
                            unit.output_interfaces[out_medium] = connection
                            system.input_interfaces[in_medium] = connection
                        end
                    end
                end
            end
        end
    end
end

function produce(
    systems :: Grouping,
    parameters :: Dict{String, Any}
)
    grid_e = [u for u in each(systems) if (typeof(u) <: GridConnection && u.medium == m_e_ac_230v)][1]
    demand_h = [u for u in each(systems) if (typeof(u) <: Demand && u.medium == m_h_w_60c)][1]
    demand_e = [u for u in each(systems) if (typeof(u) <: Demand && u.medium == m_e_ac_230v)][1]
    buffer = [u for u in each(systems) if typeof(u) <: BufferTank][1]
    e_bus = [u for u in each(systems) if (typeof(u) <: Bus && u.medium == m_e_ac_230v)][1]
    chpp = [u for u in each(systems) if typeof(u) <: CHPP][1]
    hp = [u for u in each(systems) if typeof(u) <: HeatPump][1]
    pv_plant = [u for u in each(systems) if typeof(u) <: PVPlant][1]

    # reset balances
    e_bus.balance = 0.0
    chpp.last_produced_e = 0.0
    chpp.last_produced_h = 0.0
    hp.last_consumed_e = 0.0
    hp.last_produced_h = 0.0
    pv_plant.last_produced_e = 0.0

    # unload buffer
    buffer.load -= Wh(load_at_time(demand_h, parameters["time"]))

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
        chpp.last_produced_e = max_produce_e * usage_fraction
        chpp.last_produced_h = max_produce_h * usage_fraction
    end

    # electricity demand and PV plant
    e_bus.balance -= Wh(load_at_time(demand_e, parameters["time"]))
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
    base_sine = Base.Math.sin(Base.MathConstants.pi * (time % seconds_in_day) / seconds_in_day)
    return Base.Math.max(0.0, Base.Math.min(
        plant.amplitude,
        1.4 * plant.amplitude * base_sine * base_sine * base_sine - 0.2 * plant.amplitude
    ))
end

end