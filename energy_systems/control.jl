
struct Condition
    name :: String
    parameters :: Dict{String, Any}
    required_systems :: Dict{String, Tuple{Type, Union{Nothing, MediumCategory}}}
    linked_systems :: Grouping
end

function rel(condition :: Condition, name :: String) :: ControlledSystem
    return condition.linked_systems[name]
end

function Condition(
    name :: String,
    parameters :: Dict{String, Any}
) :: Condition
    required_systems = Dict{String, Tuple{Type, Union{Nothing, MediumCategory}}}()
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

    elseif condition.name == "Little PV power"
        # by checking the current balance of the output interface we avoid the problem
        # that changes are counted doubled when the energy already has been consumed
        # but only once if the energy is still "within" the interface
        outface = rel(condition, "pv_plant").output_interfaces[m_e_ac_230v]
        return (if outface.balance != 0.0 outface.sum_abs_change else outface.sum_abs_change * 0.5 end
            < condition.parameters["threshold"])

    elseif condition.name == "Much PV power"
        # see condition "Little PV power" for how this works
        outface = rel(condition, "pv_plant").output_interfaces[m_e_ac_230v]
        return (if outface.balance != 0.0 outface.sum_abs_change else outface.sum_abs_change * 0.5 end
            >= condition.parameters["threshold"])

    elseif condition.name == "Sufficient charge"
        return unit.load >= condition.parameters["threshold"]

    elseif condition.name == "Insufficient charge"
        return unit.load < condition.parameters["threshold"]
    end

    throw(KeyError(condition.name))
end
