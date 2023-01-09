"""Convenience type alias for requirements of energy systems."""
const EnSysRequirements = Dict{String, Tuple{Type, Union{Nothing, MediumCategory}}}

"""
A boolean decision variable for a transition in a state machine.

Because the implementation of conditions can be arbitrary but require the values of certain
energy systems, a condition must be parameterized with its requirements. This is part of
the input to the simulation and cannot be derived otherwise.
"""
struct Condition
    """An identifiable name."""
    name :: String

    """Parameters the condition requires and holds the values after loading."""
    parameters :: Dict{String, Any}

    """Defines which systems the condition requires, indexed by an internal name.

    For some systems a medium is required as they can take varying values.
    """
    required_systems :: EnSysRequirements

    """The systems linked to the condition indexed by an internal name."""
    linked_systems :: Grouping
end

"""
    rel(condition, name)

Get the linked system of the given name for a condition.
"""
function rel(condition :: Condition, name :: String) :: ControlledSystem
    return condition.linked_systems[name]
end

"""
Constructor for Condition.

# Arguments
- `name::String`: The name of the condition
- `parameters::Dict{String, Any]`: Parameters for the condition. Not to be confused with
    the project-wide parameters for the entire simulation.

# Returns
- `Condition`: A Condition instance with default parameter values and information on which
    energy systems are required, but systems have been linked yet
"""
function Condition(
    name :: String,
    parameters :: Dict{String, Any}
) :: Condition
    required_systems = EnSysRequirements()
    default_params = Dict{String, Any}()

    if name == "Buffer < X%"
        required_systems["buffer"] = (BufferTank, nothing)
        default_params["percentage"] = 0.5

    elseif name == "Buffer >= X%"
        required_systems["buffer"] = (BufferTank, nothing)
        default_params["percentage"] = 0.5

    elseif name == "Min run time"
        # nothing to do

    elseif name == "Would overfill thermal buffer"
        required_systems["buffer"] = (BufferTank, nothing)

    elseif name == "Little PV power"
        required_systems["pv_plant"] = (PVPlant, nothing)
        default_params["threshold"] = 1000

    elseif name == "Much PV power"
        required_systems["pv_plant"] = (PVPlant, nothing)
        default_params["threshold"] = 1000

    elseif name == "Sufficient charge"
        default_params["threshold"] = 0.2

    elseif name == "Insufficient charge"
        default_params["threshold"] = 0.05

    elseif name == "HP is running"
        required_systems["heat_pump"] = (HeatPump, nothing)

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

"""
    link(condition, systems)

Look for the condition's required systems in the given set and link the condition to them.

For example, if a condition required a system "grid_out" of type GridConnection and medium
m_e_ac_230v, it will look through the set of given systems and link to the first match.
"""
function link(condition :: Condition, systems :: Grouping)
    for (name, req_unit) in pairs(condition.required_systems)
        found_link = false
        for unit in each(systems)
            if isa(unit, req_unit[1])
                if (req_unit[2] !== nothing
                    && hasfield(typeof(unit), Symbol("medium"))
                    && unit.medium == req_unit[2]
                )
                    condition.linked_systems[name] = unit
                    found_link = true
                elseif req_unit[2] === nothing
                    condition.linked_systems[name] = unit
                    found_link = true
                end
            end
        end

        if !found_link
            throw(KeyError("Could not find match for required system $name "
                * "for condition $(condition.name)"))
        end
    end
end

"""
    link_control_with(unit, systems)

Link the given systems with all conditions of the given unit.

See also [`link`](@ref)
"""
function link_control_with(unit :: ControlledSystem, systems :: Grouping)
    for table in values(unit.controller.state_machine.transitions)
        for condition in table.conditions
            link(condition, systems)
        end
    end
end

"""
Maps a vector of boolean values to integers.

This is used to define the transitions of a StateMachine by defining which values of
conditions lead to which state.

# Examples
```
table = TruthTable(
    conditions=[Condition("foo is big"), Condition("bar is small")],
    table_data=Dict{Tuple, UInt}(
        (false, false) => 1,
        (false, true) => 1,
        (true, false) => 2,
        (true, true) => 1,
    )
)
```
This example defines a transitions for a state machine with two states. If the condition
"foo is big" is true and the condition "bar is small" is false, the new state should be the
second state, otherwise the first.
"""
Base.@kwdef struct TruthTable
    conditions :: Vector{Condition}
    table_data :: Dict{Tuple, UInt}
end

"""
Implementation of state machines with generalized conditions instead of an input alphabet.

Similar to the state machines used in regular languages, a state machine is always in one
of its states and certain conditions define how transitions between states occurs. Instead
of checking against characters in an input alphabet, these state machines are typically
checked once every time step of a simulation and the conditions can have arbitrary
implementations that require the simulation state as input.
"""
mutable struct StateMachine
    """The current state of the state machine."""
    state :: UInt

    """A map of state names indexes by their ID."""
    state_names :: Dict{UInt, String}

    """Maps states to a TruthTable that define the transitions in that state."""
    transitions :: Dict{UInt, TruthTable}

    """
    The number of steps the state machine has been in the current state.

    Starts counting at 1.
    """
    time_in_state :: UInt
end

"""
Constructor for non-default fields.
"""
StateMachine(
    state :: UInt,
    state_names :: Dict{UInt, String},
    transitions :: Dict{UInt, TruthTable}
) = StateMachine(
    state,
    state_names,
    transitions,
    UInt(0)
)

"""
Default constructor that creates a state machine with only one state called "Default".
"""
StateMachine() = StateMachine(
    UInt(1),
    Dict(UInt(1)=>"Default"),
    Dict(UInt(1)=>TruthTable(
        conditions=Vector(),
        table_data=Dict()
    ))
)

"""
Wraps around the mechanism of control for the operation strategy of an EnergySystem.

For now this merely the StateMachine handling the controller state and the name of the
operation strategy, but this can be easily extended.
"""
Base.@kwdef mutable struct Controller
    strategy :: String
    state_machine :: StateMachine
end

"""
    move_state(unit, systems, parameters)

Checks the controller of the given unit and moves the state machine to its new state.
"""
function move_state(
    unit :: ControlledSystem,
    systems :: Grouping,
    parameters :: Dict{String, Any}
)
    machine = unit.controller.state_machine
    old_state = machine.state
    table = machine.transitions[machine.state]

    if length(table.conditions) > 0
        evaluations = Tuple(
            check(condition, unit, parameters)
            for condition in table.conditions
        )
        new_state = table.table_data[evaluations]
        machine.state = new_state
    else
        new_state = old_state
    end

    if old_state == new_state
        machine.time_in_state += 1
    else
        machine.time_in_state = 1
    end
end

"""
    check(condition, unit, parameters)

Check the condition of the given unit and return the result.
"""
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
        return (unit.controller.state_machine.time_in_state
            * parameters["time_step_seconds"]
            >= unit.min_run_time)

    elseif condition.name == "Would overfill thermal buffer"
        return (rel(condition, "buffer").capacity
            - rel(condition, "buffer").load
            < unit.power * unit.min_power_fraction)

    elseif condition.name == "Little PV power"
        # by checking the current balance of the output interface we avoid the problem
        # that changes are counted doubled when the energy already has been consumed
        # but only once if the energy is still "within" the interface
        outface = rel(condition, "pv_plant").output_interfaces[m_e_ac_230v]
        return (if outface.balance != 0.0 outface.sum_abs_change else outface.sum_abs_change * 0.5 end
            < condition.parameters["threshold"] * rel(condition, "pv_plant").amplitude * 0.25)

    elseif condition.name == "Sufficient charge"
        return unit.load >= condition.parameters["threshold"] * unit.capacity

    elseif condition.name == "HP is running"
        return (rel(condition, "heat_pump").output_interfaces[m_h_w_ht1].sum_abs_change
            > parameters["epsilon"])

    end

    throw(KeyError(condition.name))
end

"""
A type of operational strategy that defines which parameters and systems a strategy requires.
"""
Base.@kwdef struct OperationalStrategyType
    """Machine-readable name of the strategy."""
    name :: String

    """Human-readable description that explains how to use the strategy."""
    description :: String

    """Constructor method for the state machine used by the strategy."""
    sm_constructor :: Function

    """A list of condition names that the strategy uses."""
    conditions :: Vector{String}

    """Required parameters for the strategy including those for the conditions."""
    strategy_parameters :: Dict{String, Any}

    """Energy systems that the strategy requires for the correct order of execution.

    This differs from the system the conditions of the strategy require.
    """
    required_systems :: EnSysRequirements
end

OP_STRATS = Dict{String, OperationalStrategyType}()

STRT_SM_PARAMS = Dict{String, Dict{String, Any}}()
STRT_SM_FUNCS = Dict{String, Function}()

include("strategies/economical_discharge.jl")
include("strategies/storage_driven.jl")
include("strategies/demand_driven.jl")
include("strategies/supply_driven.jl")
include("strategies/use_surplus_in_cycle.jl")

"""
    controller_for_strategy(strategy, parameters)

Construct the controller for the strategy of the given name using the given parameters.

# Arguments
- `strategy::String`: Must be an exact match to the name defined in the strategy's code file.
- `parameters::Dict{String, Any}`: Parameters for the configuration of the strategy. The
    names must match those in the default parameter values dictionary defined in the
    strategy's code file. Given values override default values.
# Returns
- `Controller`: The constructed controller for the given strategy.
"""
function controller_for_strategy(strategy :: String, parameters :: Dict{String, Any}) :: Controller
    if lowercase(strategy) == "default"
        return Controller("default", StateMachine())
    end

    if !(strategy in keys(STRT_SM_FUNCS) && strategy in keys(STRT_SM_PARAMS))
        throw(ArgumentError("Unknown strategy $strategy"))
    end

    params = merge(STRT_SM_PARAMS[strategy], parameters)
    machine = STRT_SM_FUNCS[strategy](params)
    return Controller(strategy, machine)
end

export Condition, TruthTable, StateMachine, link_control_with, controller_for_strategy