"""
Implementations of energy systems and required functions to construct a network of systems.

Code in this module attempts to match the model description as close as possible as well as
adds utility features and makes abstract descriptions work with actual code.

# Notes on this file in particular

Functionality that is common to all energy systems should go here, while code handling
specific features should be placed in the corresponding file for the energy system. These
individual files are included in the middle of the module. This has been proven necessary
as the implementations do require the fundamental types, definitions and functions in this
module while also extending the module. To avoid circular dependencies, they cannot be put
into a module block of the same module in a different file.

Utility and interface functions should be placed after the part that includes the
individual energy systems. These interact with the previous parts to provide an interface
to the simulation as a whole as well as provide functionality on groups of energy systems.
"""
module EnergySystems

export MediumCategory, EnergySystem, ControlledSystem, Condition, TruthTable, StateMachine,
    represent, pprint, link_control_with, each, Grouping, link_production_with, check_balances,
    perform_steps

"""
Categories that each represent a physical medium in conjunction with additional attributes,
such as temperature or voltage. These attributes are not necessarily unchanging, but are
considered the nominal range. For example, a heating system might circulate water anywhere
from 30°C to 60°C, but the nominal temperature is considered to be 60°C. This is intended
so it becomes possible to prevent linking systems that do not work together because they
work on different nominal temperatures, while both work with the same physical medium,
for example water.

The names are structured in a composite of segments. For example, these are:
    m_e_ac_230v

    1. m: This segment is used to distinguish its symbols from the symbols of other types
    2. e: The energy type, in this case electricity
    3. ac: The physical medium, in this case AC current
    4. 230v: Additional attributes of nominal value, in this case 230V
"""
@enum MediumCategory m_e_ac_230v m_c_g_natgas m_h_w_60c

"""
Enumerations of the archetype of an energy system describing its general function.

These are described in more detail in the accompanying documentation of the simulation
model.
"""
@enum SystemFunction infinite_sink infinite_source limited_sink limited_source transformer storage bus

"""
Enumerations of a simulation step that can be performed on an energy system.

The names are prefixed with `s` to avoid shadowing functions of the same name.
"""
@enum Step s_reset s_control s_produce s_load s_distribute

"""
The basic type of all energy systems.
"""
abstract type EnergySystem end

"""
A type describing an energy system that has control functionality.

Because Julia does not have interface->implementation like OOP languages such as Java,
types implementing this abstract type are further to be assumed to have the fields required
by all energy systems, in particular the field `controller` of type `StateMachine`.
"""
abstract type ControlledSystem <: EnergySystem end

"""
Convenience alias to a dict mapping UAC keys to an energy system.
"""
const Grouping = Dict{String, ControlledSystem}

"""
    each(systems :: Grouping)

Generator over each of the energy systems in the given grouping.
"""
function each(systems :: Grouping) :: Base.ValueIterator
    return values(systems)
end

"""
Handles the tracking of energy being transfered from one energy system to another.

This abstraction is useful to avoid energy system having to "know" other types. Instead of
calling functions to transfer energy, a system can draw from or load into a SystemInterface
instance instead.

Energy is considered to always flow from the source to the target. A negative balance is a
lack of energy that needs to be covered in order for the energy to be balanced.

A system interface keeps track of how much energy (in absolute terms) was transfered
via the interface. Assuming the energy balance holds, at the end of a time step a system
interface's field "sum_abs_change" will have a value of twice the total energy transfered.
"""
Base.@kwdef mutable struct SystemInterface
    """The source system providing energy"""
    source :: Union{Nothing, ControlledSystem} = nothing

    """The target system receiving energy"""
    target :: Union{Nothing, ControlledSystem} = nothing

    """The current balance of the interface"""
    balance :: Float64 = 0.0

    """The sum of absolute changes to the interface's balance"""
    sum_abs_change :: Float64 = 0.0
end

"""
    add!(interface, change)

Add the given amount of energy (in Wh) to the balance of the interface.
"""
function add!(interface :: SystemInterface, change :: Float64)
    interface.balance += change
    interface.sum_abs_change += abs(change)
end

"""
    sub!(interface, change)

Subtract the given amount of energy (in Wh) from the balance of the interface.
"""
function sub!(interface :: SystemInterface, change :: Float64)
    interface.balance -= change
    interface.sum_abs_change += abs(change)
end

"""
    set!(interface, new_val)

Set the balance of the interface to the given new value.
"""
function set!(interface :: SystemInterface, new_val :: Float64)
    interface.sum_abs_change += abs(interface.balance - new_val)
    interface.balance = new_val
end

"""
    reset!(interface)

Reset the interface back to zero.
"""
function reset!(interface :: SystemInterface)
    interface.balance = 0.0
    interface.sum_abs_change = 0.0
end

"""
Convenience type used define the required system interfaces of an energy system.

To simultaneously define what is required and then hold references to instances after the
whole system has been loaded, it maps a medium category to either nothing (before systems
are linked) or a SystemInterface instance.
"""
const InterfaceMap = Dict{MediumCategory, Union{Nothing, SystemInterface}}

"""
    balance_on(interface, unit)

Return the balance of a unit in respect to the given interface.

For most energy systems this is simply the balance of the interface itself, but for Bus
instances this is handled differently. This function helps to implement an energy system
without having to check if its connected to a Bus or directly to a system.

# Arguments
- `interface::SystemInterface`: The interface "on which" the balance is calculated. This
    also defines which system is the source.
- `unit::ControlledSystem`: The receiving system

# Returns
- `Float64`: The balance of the target system that can be considered a demand on the source
    system
- `Float64`: An additional demand that covers the free storage space connected to the
    target system
"""
function balance_on(
    interface :: SystemInterface,
    unit :: ControlledSystem
) :: Tuple{Float64, Float64}
    return interface.balance, 0.0
end

"""
    balance(unit)

Calculate the energy balance of the given unit as a whole.

This is expected to start at zero at the beginning of a time step and return to zero at
the end of it. If it is not zero, either the simulation failed to correctly calculate the
energy balance of the entire system or the simulated network was not able to ensure the
balance on the current time step. In either case, something went wrong.
"""
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

"""
    reset(unit)

Reset the given energy system back to zero.

For most energy systems this only resets the balances on the system interfaces but some
systems might require more complex reset handling.
"""
function reset(unit :: ControlledSystem)
    for inface in values(unit.input_interfaces)
        if inface !== nothing reset!(inface) end
    end
    for outface in values(unit.output_interfaces)
        if outface !== nothing reset!(outface) end
    end
end

"""
    control(unit, systems, parameters)

Perform the control calculations for the given energy system.

# Arguments
- `unit::ControlledSystem`: The system for which control is handled
- `systems::Grouping`: A reference dict to all energy systems in the project
- `parameters::Dict{String, Any}`: Project-wide parameters
"""
function control(
    unit :: ControlledSystem,
    systems :: Grouping,
    parameters :: Dict{String, Any}
)
    move_state(unit, systems, parameters)
end

"""
    produce(unit, parameters, watt_to_wh)

Perform the production calculations for the given energy system.

# Arguments
- `unit::ControlledSystem`: The system for which production is calculated
- `parameters::Dict{String, Any}`: Project-wide parameters
- `watt_to_wh::Function`: Utility function to calculate work from a given power
"""
function produce(
    unit :: ControlledSystem,
    parameters :: Dict{String, Any},
    watt_to_wh :: Function
)
    # default implementation is to do nothing
end

"""
    load(unit, parameters, watt_to_wh)

Load excess energy into storage energy systems.

For non-storage systems this function does nothing.

# Arguments
- `unit::ControlledSystem`: The system for which production is calculated
- `parameters::Dict{String, Any}`: Project-wide parameters
- `watt_to_wh::Function`: Utility function to calculate work from a given power
"""
function load(
    unit :: ControlledSystem,
    parameters :: Dict{String, Any},
    watt_to_wh :: Function
)
    # default implementation is to do nothing
end

"""
    distribute(unit)

Distribute the energy inputs and outputs of a bus.

For non-bus systems this function does nothing.
"""
function distribute!(unit :: ControlledSystem)
    # default implementation is to do nothing
end

"""
    output_values(unit)

Specify which outputs an energy system can provide.

For the special values "IN" and "OUT" a medium category is required for fetching the actual
value, while this method only specifies that there is an input or output.
"""
function output_values(unit :: EnergySystem) :: Vector{String}
    return ["IN", "OUT"]
end

"""
    output_value(unit, key)

Return the value for the output with the given output key.

Note that for the "IN" and "OUT" output values, the value corresponds to the sum of
absolute changes of the system interfaces and divided by 2. This behaviour is part of the
expected use of the method.

Args:
- `unit::EnergySystem`: The energy system for which to fetch the output
- `key::OutputKey`: An OutputKey specifying which output to return. This should be one of
    the options provided by `output_values()` as well as "IN" or "OUT"
Returns:
- `Float64`: The value of the desired output
Raises:
- `KeyError`: The key value requested must be one the energy system can provide
"""
function output_value(unit :: EnergySystem, key :: OutputKey) :: Float64
    if key.key_value == "IN"
        return unit.input_interfaces[key.medium].sum_abs_change * 0.5
    elseif key.key_value == "OUT"
        return unit.output_interfaces[key.medium].sum_abs_change * 0.5
    end
    raise(KeyError(key.key_value))
end

"""
    represent(unit, time)

Represent the state of a unit and any values apart from the energy balance.

# Arguments
- `unit::ControlledSystem`: The energy system to represent
- `time::Int`: The time of the simulation in second. Used mostly for time-dependant values

# Returns
- `String`: A short representation of the given unit
"""
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

"""
    pprint(unit, time)

Print a representation of the given unit to the console output.

# Arguments
- `unit::ControlledSystem`: The energy system to print
- `time::Int`: The time of the simulation in second. Used mostly for time-dependant values
"""
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

"""
    link_production_with(unit, systems)

Set the production targets of the given unit to the given energy systems.

This function is used to construct the network of energy system from a graph input that
determines which systems provide energy to which other systems.

# Arguments
- `unit::ControlledSystem`: The unit providing energy
- `systems::Grouping`: A set of systems receiving energy. As systems might have multiple
    outputs, this is used to set them all at once.
"""
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

"""
    check_balances(systems, epsilon)

Check the energy balance of the given systems and return warnings of any violations.

# Arguments
- `system::Grouping`: The systems to check
- `epsilon::Float64`: A balance is only considered violated if the absolute value of the
    sum is larger than this value. This helps with spurious floating point issues

# Returns
- `Vector{Tuple{String, Float64}}`: A list of tuples, where each tuple is the key of the
    system that has a non-zero energy balance and the value of that balance.
"""
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

"""
    perform_steps(systems, order_of_steps, parameters)

Perform the simulation steps of one time step for the given systems in the given order.

# Arguments
- `systems::Grouping`: The entirety of the energy systems
- `order_of_steps::Vector{Vector{Any}}`: Defines which steps are performed in which order.
    Each system must go through the simulation steps defined in EnergySystems.Step, but the
    order is not the same for all simulations. Determining the order must be handled
    elsewhere, as this function only goes through and calls the appropriate functions. The
    first item of each entry must be the key of the system for which the following steps
    are performed.
- `parameters::Dict{String, Any}`: Project-wide parameters

# Examples
```
    systems = Grouping(
        "system_a" => EnergySystem(),
        "system_b" => EnergySystem(),
    )
    order = [
        ["system_a", EnergySystems.s_control]
        ["system_b", EnergySystems.s_control, EnergySystems.s_produce]
        ["system_a", EnergySystems.s_produce]
    ]
    parameters = Dict{String, Any}("time" => 0)
    perform_steps(systems, order, parameters)
```
In this example the control of system A is performed first, then control and production of
system B and finally production of system A.
"""
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