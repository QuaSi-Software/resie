"""
Maps a vector of boolean values to integers.

This is used to define the transitions of a StateMachine by defining which values of
conditions lead to which state.

# Examples
```
table = TruthTable(
    conditions=[x -> x > 5, x -> x <= 5],
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
    conditions::Vector{Function}
    table_data::Dict{Tuple,UInt}
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
    state::UInt

    """A map of state names indexes by their ID."""
    state_names::Dict{UInt,String}

    """Maps states to a TruthTable that define the transitions in that state."""
    transitions::Dict{UInt,TruthTable}

    """
    The number of steps the state machine has been in the current state.

    Starts counting at 1.
    """
    time_in_state::UInt
end

"""
Constructor for non-default fields.
"""
StateMachine(
    state::UInt,
    state_names::Dict{UInt,String},
    transitions::Dict{UInt,TruthTable}
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
    Dict(UInt(1) => "Default"),
    Dict(UInt(1) => TruthTable(
        conditions=Vector(),
        table_data=Dict()
    ))
)

"""
Base type for control modules.
"""
abstract type ControlModule end

"""
Wraps around the mechanism of control for the operation strategy of a Component.
"""
mutable struct Controller
    base_module::Union{Nothing,ControlModule}
    modules::Vector{ControlModule}
end

"""
Default constructor with empty fields.
"""
Controller() = Controller(nothing, [])

function load_storages(controller::Controller, medium::Symbol)::Bool
    return default(
        controller.base_module.parameters, "load_storages " * String(medium), true
    )
end

function unload_storages(controller::Controller, medium::Symbol)::Bool
    return default(
        controller.base_module.parameters, "unload_storages " * String(medium), true
    )
end

function update(controller::Controller)
    for control_module in controller.modules
        update(control_module)
    end
end

function upper_plr_limit(
    controller::Controller,
    sim_params::Dict{String,Any}
)::Float64
    if length(controller.modules) > 0
        return Base.maximum(upper_plr_limit(mod, sim_params) for mod in controller.modules)
    else
        return 1.0
    end
end

function charge_is_allowed(controller::Controller, sim_params::Dict{String,Any})
    return all(charge_is_allowed(mod, sim_params) for mod in controller.modules)
end

function discharge_is_allowed(controller::Controller, sim_params::Dict{String,Any})
    return all(discharge_is_allowed(mod, sim_params) for mod in controller.modules)
end

"""
    move_state(unit, components, sim_params)
"""
function move_state(
    machine::StateMachine,
)
    old_state = machine.state
    table = machine.transitions[machine.state]

    if length(table.conditions) > 0
        evaluations = Tuple(condition(machine) for condition in table.conditions)
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
