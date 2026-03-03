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
function StateMachine(state::UInt,
                      state_names::Dict{UInt,String},
                      transitions::Dict{UInt,TruthTable})
    StateMachine(state,
                 state_names,
                 transitions,
                 UInt(0))
end

"""
Default constructor of StateMachine that creates a state machine with only one state called
"Default" and no transitions (as there no other states).
"""
function StateMachine()
    StateMachine(UInt(1),
                 Dict(UInt(1) => "Default"),
                 Dict(UInt(1) => TruthTable(; conditions=Vector(),
                                            table_data=Dict())))
end

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
    cmf_reorder_inputs
    cmf_reorder_outputs
    cmf_negotiate_temperature
    cmf_limit_cooling_input_temperature
    cmf_check_src_to_snk
    cmf_change_bus_priorities
    cmf_reorder_operations
end

"""
Default method for function control_module_name. Control module files should provide a
method that returns the name, as used in the project config.
"""
control_module_name(::Type{ControlModule})::String = "control_module"

"""
Wraps around the mechanism of control for the operational strategy of a Component.

Holds general parameters of control and acts as container for control modules.
"""
mutable struct Controller
    parameters::Dict{String,Any}
    modules::Vector{ControlModule}

    function Controller(config::Union{Nothing,Dict{String,Any}})::Controller
        return new(Base.merge(Dict{String,Any}(
                                  "aggregation_plr_limit" => "max",
                                  "aggregation_charge" => "all",
                                  "aggregation_discharge" => "all",
                                  "aggregation_check_src_to_snk" => "all",
                                  "consider_m_el_in" => true,
                                  "consider_m_el_out" => true,
                                  "consider_m_gas_in" => true,
                                  "consider_m_fuel_in" => true,
                                  "consider_m_h2_out" => true,
                                  "consider_m_o2_out" => true,
                                  "consider_m_heat_out" => true,
                                  "consider_m_heat_ht_out" => true,
                                  "consider_m_heat_lt_out" => true,
                                  "consider_m_heat_in" => true,
                              ),
                              config === nothing ? Dict{String,Any}() : config),
                   [])
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
    return default(controller.parameters, "load_storages " * String(medium), true)
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
    return default(controller.parameters, "unload_storages " * String(medium), true)
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
    limits = collect(upper_plr_limit(mod, sim_params)
                     for mod in controller.modules
                     if has_method_for(mod, cmf_upper_plr_limit))
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
    flags = collect(charge_is_allowed(mod, sim_params)
                    for mod in controller.modules
                    if has_method_for(mod, cmf_charge_is_allowed))
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
    flags = collect(discharge_is_allowed(mod, sim_params)
                    for mod in controller.modules
                    if has_method_for(mod, cmf_discharge_is_allowed))
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

"""
Callback for reordering inputs of a component according to the temperatures of the inputs.

As there is no clear way to aggregate multiple reorderings the last module to apply a
reordering "wins" and its index permutation is returned. If no module performs any
reordering the indices of the `temps_max` argument is returned, resulting in no changes.

It is assumed, but not checked, that all vectors, to which the permutation is applied, have
the same length.

# Arguments
- `controller::Controller`: The controller containing modules
- `temps_min::Vector{<:Temperature}`: The minimum temperature vector
- `temps_max::Vector{<:Temperature}`: The maximum temperature vector
# Returns
- `Vector{Integer}`: The index permutation
"""
function reorder_inputs(controller::Controller,
                        temps_min::Vector{<:Temperature},
                        temps_max::Vector{<:Temperature})::Vector{Integer}
    reordering = [i for i in 1:length(temps_max)]

    for mod in controller.modules
        if !has_method_for(mod, cmf_reorder_inputs)
            continue
        end

        reordering = reorder_inputs(mod, temps_min, temps_max)
    end

    return reordering
end

"""
Callback for reordering outputs of a component according to the temperatures of the outputs.

As there is no clear way to aggregate multiple reorderings the last module to apply a
reordering "wins" and its index permutation is returned. If no module performs any
reordering the indices of the `temps_max` argument is returned, resulting in no changes.

It is assumed, but not checked, that all vectors, to which the permutation is applied, have
the same length.

# Arguments
- `controller::Controller`: The controller containing modules
- `temps_min::Vector{<:Temperature}`: The minimum temperature vector
- `temps_max::Vector{<:Temperature}`: The maximum temperature vector
# Returns
- `Vector{Integer}`: The index permutation
"""
function reorder_outputs(controller::Controller,
                         temps_min::Vector{<:Temperature},
                         temps_max::Vector{<:Temperature})::Vector{Integer}
    reordering = [i for i in 1:length(temps_min)]

    for mod in controller.modules
        if !has_method_for(mod, cmf_reorder_outputs)
            continue
        end

        reordering = reorder_outputs(mod, temps_min, temps_max)
    end

    return reordering
end

"""
Callback for determine_temperature_and_energy of a component according to the temperatures of the outputs.

This function only checks if a control module exists for the chosen interconnection from source
to target. No aggregation with other control modules. 

# Arguments
- `controller::Controller`: The controller containing modules
- `components::Grouping`: All components of the energy system.
- `source_uac::String`: The source unit uac.
- `target_uac::String`: The target unit uac.
- `sim_params::Dict{String,Any}`: Simulation parameters.

# Returns
- `success::Bool`: A bool indicating if a control module exists between source and target (true) or not (false).
- `temperature::Temperature`: The negotiated temperature.
- `max_energy::Float64`: The energy corresponding to the negotiated temperature.
"""
function determine_temperature_and_energy(controller::Controller,
                                          components::Grouping,
                                          source_uac::String,
                                          target_uac::Stringing,
                                          sim_params::Dict{String,Any})::Tuple{Bool,Temperature,Float64}
    for mod in controller.modules
        if !has_method_for(mod, cmf_negotiate_temperature)
            continue
        end
        if mod.parameters["target_uac"] !== target_uac
            continue
        end

        results = determine_temperature_and_energy(mod, components, source_uac, target_uac, sim_params)
        return true, results[1], results[2]
    end

    if sim_params["time"] == 0
        if components[source_uac] isa TemperatureNegotiateSource &&
           target_uac !== nothing && components[target_uac] isa TemperatureNegotiateTarget
            @warn "From $(source_uac) to $(target_uac), no control module is activated. This can lead to unexpected " *
                  "results. Add a `negotiate_temperature` control module at $(source_uac)!"
        end
    end

    return false, nothing, 0.0
end

"""
Callback for cooling_input_temperature_exceeded


# Arguments
- `controller::Controller`: The controller containing modules
- `components::Grouping`: All components of the energy system.
- `target_uac::String`: The target unit uac.

# Returns
- `Bool`: A bool indicating if the input temperature is exceeded (true, meaning no energy flow is allowed) 
          or not (false, energy flow is allowed)

"""
function cooling_input_temperature_exceeded(controller::Controller,
                                            target_uac::Stringing,
                                            sim_params::Dict{String,Any})::Bool
    for mod in controller.modules
        if !has_method_for(mod, cmf_limit_cooling_input_temperature)
            continue
        end
        if mod.parameters["target_uac"] !== target_uac
            continue
        end
        components = Dict{String,Component}(get_run(sim_params["run_ID"]).components)
        return cooling_input_temperature_exceeded(mod, components, target_uac)
    end

    return false
end

"""
Callback for checking if a specific source is allowed to be used to supply a specific sink.

If multiple modules exist that implement this callback, the results of all modules are
aggregated as defined by the control parameter `aggregation_check_src_to_snk` with possible
options being `all` (all modules have to return true) or `any` (at least one module has to
return true).

# Arguments
- `controller::Controller`: The controller containing modules
- `in_uac::Stringing`: The UAC of the source.
- `out_uac::Stringing`: The UAC of the sink.
# Returns
- `Bool`: True, if the source is allowed to be used to supply the sink. False otherwise.
"""
function check_src_to_snk(controller::Controller,
                          in_uac::Stringing,
                          out_uac::Stringing)::Bool
    flags = collect(check_src_to_snk(mod, in_uac, out_uac)
                    for mod in controller.modules
                    if has_method_for(mod, cmf_check_src_to_snk))
    if length(flags) == 0
        return true
    end

    if controller.parameters["aggregation_check_src_to_snk"] == "all"
        return all(flags)
    elseif controller.parameters["aggregation_check_src_to_snk"] == "any"
        return any(flags)
    end

    return true
end

"""
Callback for reordering the operations of the simulation for one time step.

If multiple control modules exist that implement this callback, they will be called one
after the other with the result of the previous one. The order is not deterministic. As
such it is assumed that they can be called in any order for the same result.

This callback is also not specific to a controller (and thus a component) as the callback
is used once at the beginning of a new time step.

# Arguments
- `components::Grouping`: All components of the energy system.
- `order_of_operations::OrderOfOperations`: The unmodified order of operations determined
    at the start of the simulation
- `sim_params::Dict{String,Any})`: Simulation parameters
# Returns
- `OrderOfOperations`: The modified order of operations
"""
function reorder_operations(components::Grouping,
                            order_of_operations::OrderOfOperations,
                            sim_params::Dict{String,Any})::OrderOfOperations
    for unit in each(components)
        for mod in unit.controller.modules
            if !has_method_for(mod, cmf_reorder_operations)
                continue
            end

            order_of_operations = reorder_operations(mod, order_of_operations, sim_params)
        end
    end

    return order_of_operations
end

"""
Callback for changing a bus' priorities for one time step.

This is particularly relevant for control modules that also change the order of operations
as both is required to ensure correct calculations when the priorities on a bus is changed.

This callback is also not specific to a controller (and thus a component) as the callback
is used once at the beginning of a new time step.

# Arguments
- `components::Grouping`: All components of the energy system.
- `sim_params::Dict{String,Any})`: Simulation parameters
"""
function change_bus_priorities!(components::Grouping,
                                sim_params::Dict{String,Any})
    for unit in each(components)
        for mod in unit.controller.modules
            if !has_method_for(mod, cmf_change_bus_priorities)
                continue
            end

            change_bus_priorities!(mod, components, sim_params)
        end
    end
end
