module EnergySystems

export MediumCategory, EnergySystem, ControlledSystem, Condition, TruthTable, StateMachine,
    control, represent, pprint, check, produce, production, link_control_with, each, Grouping,
    link_production_with, check_balances

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
    min :: Float64 = 0.0
    max :: Float64 = 0.0
end

function add!(interface :: SystemInterface, change :: Float64)
    interface.balance += change
    interface.min = min(interface.min, interface.balance)
    interface.max = max(interface.max, interface.balance)
end

function sub!(interface :: SystemInterface, change :: Float64)
    interface.balance -= change
    interface.min = min(interface.min, interface.balance)
    interface.max = max(interface.max, interface.balance)
end

function set!(interface :: SystemInterface, new_val :: Float64)
    interface.balance = new_val
    interface.min = min(interface.min, new_val)
    interface.max = max(interface.max, new_val)
end

function reset!(interface :: SystemInterface)
    interface.balance = 0.0
    interface.min = 0.0
    interface.max = 0.0
end

const InterfaceMap = Dict{MediumCategory, Union{Nothing, SystemInterface}}

function gather_from_all!(interface :: SystemInterface, unit :: ControlledSystem)
    return # the default implementation is to do nothing
end

function reset(unit :: ControlledSystem)
    for inface in values(unit.input_interfaces)
        if inface !== nothing reset!(inface) end
    end
    for outface in values(unit.output_interfaces)
        if outface !== nothing reset!(outface) end
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

# the order matters here as some energy systems require the definition of certain
# basic systems such as a bus or a grid connection
include("demand.jl")
include("grid_connection.jl")
include("bus.jl")
include("buffer_tank.jl")
include("chpp.jl")
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

function check_balances(
    systems :: Grouping,
    epsilon :: Float64
) :: Vector{Tuple{String, Float64}}
    warnings = []

    for (key, unit) in pairs(systems)
        balance = check_balance(unit)
        if balance > epsilon || balance < -epsilon
            push!(warnings, (key, balance))
        end
    end

    return warnings
end

function check_balance(unit :: ControlledSystem) :: Float64
    balance = 0.0

    for inface in values(unit.input_interfaces)
        if inface !== nothing balance += inface.balance end
    end

    for outface in values(unit.output_interfaces)
        if outface !== nothing balance += outface.balance end
    end

    return balance
end

function control(
    systems :: Grouping,
    order :: Vector{String},
    parameters :: Dict{String, Any}
)
    for key in order
        control(systems[key], systems, parameters)
    end
end

function produce(
    systems :: Grouping,
    order :: Vector{String},
    parameters :: Dict{String, Any}
)
    watt_to_wh = function (watts :: Float64)
        watts * Float64(parameters["time_step_seconds"]) / 3600.0
    end

    for unit in each(systems)
        reset(unit)
    end

    for key in order
        produce(systems[key], parameters, watt_to_wh)
    end

    for key in order
        unit = systems[key]
        if unit.is_storage
            load(unit, parameters, watt_to_wh)
        end
    end
end

end