"""
Maps a vector of boolean values to integers.

This is used to define the transitions of a StateMachine by defining which values of
conditions lead to which state.

# Examples
```
table = TruthTable(
    conditions=[sm -> sm.x > 5, sm -> sm.y <= 3],
    table_data=Dict{Tuple, UInt}(
        (false, false) => 1,
        (false, true) => 1,
        (true, false) => 2,
        (true, true) => 1,
    )
)
```
This example defines transitions from the first state of a state machine with two states to
the second based on the results of two conditions. Only if the first condition is true and
the second condition is false does the state machine transition to the second state
otherwise it remains in the first.

# Conditions
The conditions given in the constructor must be functions that take only one argument,
namely the state machine itself. All other information required for calculation should be
given in a closure around calling the constructor. An example for a condition checking if a
linked component has a value for "foo" greater than 5:
```
parameters["component"] = SomeComponent(...)
table = TruthTable(
    conditions=[
        function(sm)
            return parameters["component"].foo > 5
        end
    ],
    ...
)
```
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
Constructor of StateMachine for non-default fields.
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
Default constructor of StateMachine that creates a state machine with only one state called
"Default" and no transitions (as there no other states).
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
Advances the given state machine by checking the conditions in its current state.

This also sets or increments the number of timesteps the state machine has been in the
current state. If the machine switched to a different state, this number is then 1 as the
current timestep is considered to be part of the new state.

# Arguments
- `machine::StateMachine`: The state machine to advance
"""
function move_state(machine::StateMachine)
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

"""
Base type for control modules.

Implementing types are expected to have the fields:
  - name::String
  - parameters::Dict{String,Any}
"""
abstract type ControlModule end

"""
Encodes for which functions a control module type implements methods for, which is required
to filter modules to a selection of modules that have methods for a specific function.
"""
@enum ControlModuleFunction begin
    cmf_upper_plr_limit
    cmf_charge_is_allowed
    cmf_discharge_is_allowed
end

"""
Wraps around the mechanism of control for the operational strategy of a Component.

Holds general parameters of control and acts as container for control modules.
"""
mutable struct Controller
    parameters::Dict{String,Any}
    modules::Vector{ControlModule}

    function Controller(config::Union{Nothing,Dict{String,Any}})::Controller
        return new(
            Base.merge( # parameters
                Dict{String,Any}(
                    "aggregation_plr_limit" => "max",
                    "aggregation_charge" => "all",
                    "aggregation_discharge" => "all",
                    "consider_m_el_in" => true,
                    "consider_m_el_out" => true,
                    "consider_m_gas_in" => true,
                    "consider_m_fuel_in" => true,
                    "consider_m_h2_out" => true,
                    "consider_m_o2_out" => true,
                    "consider_m_heat_out" => true,
                    "consider_m_heat_ht_out" => true,
                    "consider_m_heat_lt_out" => true,
                    "consider_m_heat_in" => true
                ),
                config === nothing ? Dict{String,Any}() : config
            ),
            [] # modules
        )
    end
end

"""
Default constructor with empty modules fields and only default parameters.
"""
Controller() = Controller(nothing)

"""
Returns if the given medium is configured to be allowed to load storages.

This checks if the corresponding parameter is set and is set to true. The default behaviour
is to allow storage loading.

# Arguments
- `controller::Controller`: The controller of the component
- `medium::Symbol`: Which medium to check
# Returns
- `Bool`: True if the parameter is set and is false, true otherwise
"""
function load_storages(controller::Controller, medium::Symbol)::Bool
    return default(
        controller.parameters, "load_storages " * String(medium), true
    )
end

"""
Returns if the given medium is configured to be allowed to unload storages.

This checks if the corresponding parameter is set and is set to true. The default behaviour
is to allow storage unloading.

# Arguments
- `controller::Controller`: The controller of the component
- `medium::Symbol`: Which medium to check
# Returns
- `Bool`: True if the parameter is set and is false, true otherwise
"""
function unload_storages(controller::Controller, medium::Symbol)::Bool
    return default(
        controller.parameters, "unload_storages " * String(medium), true
    )
end

"""
Update the controller for a timestep, more specifically the control modules.

# Arguments
- `controller::Controller`: The controller to update
"""
function update(controller::Controller)
    for control_module in controller.modules
        update(control_module)
    end
end

"""
Returns if the given control module type has a method for the given control module function.

The default is false, as the abstract type ControlModule has no methods for any of the
control module function (although it would be possible to define one).

# Arguments
- `mod::ControlModule`: The module to check
- `func::ControlModuleFunction`: The control module function to check
# Returns
- `Bool`: True if the type has a method for the function, false otherwise
"""
function has_method_for(mod::ControlModule, func::ControlModuleFunction)::Bool
    return false # no default implementation
end

"""
Aggregate the result for the upper PLR limit of matching control modules.

The aggregation can either be the maximum or minimum over the limits, which leads to two
behaviours. For example if two modules impose a limit and the aggregation `max` is used,
then the larger of the two numbers is used, which means the component operates as long any
of the modules has a limit larger than zero. Conversly if `min` is used all modules must
have a value larger zero for the component to run.

# Arguments
- `controller::Controller`: The controller containing modules
- `sim_params::Dict{String,Any}`: Simulation parameters
# Returns
- `Float64`: The aggregated upper PLR limit
"""
function upper_plr_limit(controller::Controller, sim_params::Dict{String,Any})::Float64
    limits = collect(
        upper_plr_limit(mod, sim_params)
        for mod in controller.modules
        if has_method_for(mod, cmf_upper_plr_limit)
    )
    if length(limits) == 0
        return 1.0
    end

    if controller.parameters["aggregation_plr_limit"] == "max"
        return Base.maximum(limits)
    elseif controller.parameters["aggregation_plr_limit"] == "min"
        return Base.minimum(limits)
    end

    return 1.0
end

"""
Aggregate the result for charge allowance of matching control modules.

The aggregation can be either the union or intersection of the result flags of all modules.
For example if two modules control the charging of a battery, both modules must return true
in order for charging to be allowed, assuming aggregation "all" is used. For the aggregation
method "any" it is sufficient if any one module returns true.

# Arguments
- `controller::Controller`: The controller containing modules
- `sim_params::Dict{String,Any}`: Simulation parameters
# Returns
- `Bool`: The aggregated charging flag with true meaning charging is allowed
"""
function charge_is_allowed(controller::Controller, sim_params::Dict{String,Any})::Bool
    flags = collect(
        charge_is_allowed(mod, sim_params)
        for mod in controller.modules
        if has_method_for(mod, cmf_charge_is_allowed)
    )
    if length(flags) == 0
        return true
    end

    if controller.parameters["aggregation_charge"] == "all"
        return all(flags)
    elseif controller.parameters["aggregation_charge"] == "any"
        return any(flags)
    end

    return true
end

"""
Aggregate the result for discharge allowance of matching control modules.

The aggregation can be either the union or intersection of the result flags of all modules.
For example if two modules control the discharging of a battery, both modules must return
true in order for discharging to be allowed, assuming aggregation "all" is used. For the
aggregation method "any" it is sufficient if any one module returns true.

# Arguments
- `controller::Controller`: The controller containing modules
- `sim_params::Dict{String,Any}`: Simulation parameters
# Returns
- `Bool`: The aggregated charging flag with true meaning discharging is allowed
"""
function discharge_is_allowed(controller::Controller, sim_params::Dict{String,Any})::Bool
    flags = collect(
        discharge_is_allowed(mod, sim_params)
        for mod in controller.modules
        if has_method_for(mod, cmf_discharge_is_allowed)
    )
    if length(flags) == 0
        return true
    end

    if controller.parameters["aggregation_discharge"] == "all"
        return all(flags)
    elseif controller.parameters["aggregation_discharge"] == "any"
        return any(flags)
    end

    return true
end
