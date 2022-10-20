module EnergySystems

export MediumCategory, EnergySystem, ControlledSystem, Condition, TruthTable, StateMachine,
    represent, pprint, link_control_with, each, Grouping, link_production_with, check_balances,
    perform_steps

@enum MediumCategory m_e_ac_230v m_c_g_natgas m_h_w_60c

@enum SystemFunction infinite_sink infinite_source limited_sink limited_source transformer storage bus

@enum Step s_reset s_control s_produce s_load s_distribute

abstract type EnergySystem end
abstract type ControlledSystem <: EnergySystem end

const Grouping = Dict{String, ControlledSystem}

function each(systems :: Grouping) :: Base.ValueIterator
    return values(systems)
end

Base.@kwdef mutable struct SystemInterface
    source :: Union{Nothing, ControlledSystem} = nothing
    target :: Union{Nothing, ControlledSystem} = nothing
    balance :: Float64 = 0.0
    sum_abs_change :: Float64 = 0.0
end

function add!(interface :: SystemInterface, change :: Float64)
    interface.balance += change
    interface.sum_abs_change += abs(change)
end

function sub!(interface :: SystemInterface, change :: Float64)
    interface.balance -= change
    interface.sum_abs_change += abs(change)
end

function set!(interface :: SystemInterface, new_val :: Float64)
    interface.sum_abs_change += abs(interface.balance - new_val)
    interface.balance = new_val
end

function reset!(interface :: SystemInterface)
    interface.balance = 0.0
    interface.sum_abs_change = 0.0
end

const InterfaceMap = Dict{MediumCategory, Union{Nothing, SystemInterface}}

function balance_on(
    interface :: SystemInterface,
    unit :: ControlledSystem
) :: Tuple{Float64, Float64}
    return interface.balance, 0.0
end

function balance(unit :: ControlledSystem) :: Float64
    balance = 0.0

    for inface in values(unit.input_interfaces)
        if inface !== nothing balance += inface.balance end
    end

    for outface in values(unit.output_interfaces)
        if outface !== nothing balance += outface.balance end
    end

    return balance
end

function reset(unit :: ControlledSystem)
    for inface in values(unit.input_interfaces)
        if inface !== nothing reset!(inface) end
    end
    for outface in values(unit.output_interfaces)
        if outface !== nothing reset!(outface) end
    end
end

function control(
    unit :: ControlledSystem,
    systems :: Grouping,
    parameters :: Dict{String, Any}
)
    move_state(unit, systems, parameters)
end

function produce(
    unit :: ControlledSystem,
    parameters :: Dict{String, Any},
    watt_to_wh :: Function
)
    # default implementation is to do nothing
end

function load(
    unit :: ControlledSystem,
    parameters :: Dict{String, Any},
    watt_to_wh :: Function
)
    # default implementation is to do nothing
end

function distribute!(unit :: ControlledSystem)
    # default implementation is to do nothing
end

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

# for the moment control must be an include as it contains circular dependencies
# of definitions
include("control.jl")

# the order of includes of the indicivual systems matters here as some energy systems require the
# definition of certain basic systems such as a bus or a grid connection
include("demand.jl")
include("grid_connection.jl")
include("bus.jl")
include("battery.jl")
include("buffer_tank.jl")
include("chpp.jl")
include("heat_pump.jl")
include("pv_plant.jl")

function link_production_with(unit :: ControlledSystem, systems :: Grouping)
    if isa(unit, Bus)
        for system in each(systems)
            if isa(system, Bus)
                if out_medium == system.medium
                    connection = SystemInterface(source=unit, target=system)
                    push!(system.input_interfaces, connection)
                    push!(unit.output_interfaces, connection)
                end
            else
                for in_medium in keys(system.input_interfaces)
                    if in_medium == unit.medium
                        connection = SystemInterface(source=unit, target=system)
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
                        connection = SystemInterface(source=unit, target=system)
                        push!(system.input_interfaces, connection)
                        unit.output_interfaces[out_medium] = connection
                    end
                else
                    for in_medium in keys(system.input_interfaces)
                        if out_medium == in_medium
                            connection = SystemInterface(source=unit, target=system)
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
        unit_balance = balance(unit)
        if unit_balance > epsilon || unit_balance < -epsilon
            push!(warnings, (key, unit_balance))
        end
    end

    return warnings
end

function perform_steps(
    systems :: Grouping,
    order_of_steps :: Vector{Vector{Any}},
    parameters :: Dict{String, Any}
)
    watt_to_wh = function (watts :: Float64)
        watts * Float64(parameters["time_step_seconds"]) / 3600.0
    end

    for entry in order_of_steps
        if length(entry) < 2
            continue
        end
        unit = systems[entry[1]]

        for step in entry[2:lastindex(entry)]
            if step == s_reset
                reset(unit)
            elseif step == s_control
                control(unit, systems, parameters)
            elseif step == s_produce
                produce(unit, parameters, watt_to_wh)
            elseif step == s_load
                load(unit, parameters, watt_to_wh)
            elseif step == s_distribute
                distribute!(unit)
            end
        end
    end
end

end