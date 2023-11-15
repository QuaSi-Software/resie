"""
Implementations of energy system components and required functions to construct a network
of components.

Code in this module attempts to match the model description as close as possible as well as
adds utility features and makes abstract descriptions work with actual code.

# Notes on this file in particular

Functionality that is common to all components should go here, while code handling
specific features should be placed in the corresponding file for the component. These
individual files are included in the middle of the module. This has been proven necessary
as the implementations do require the fundamental types, definitions and functions in this
module while also extending the module. To avoid circular dependencies, they cannot be put
into a module block of the same module in a different file.

Utility and interface functions should be placed after the part that includes the
individual components. These interact with the previous parts to provide an interface
to the simulation as a whole as well as provide functionality on groups of components.
"""
module EnergySystems

export check_balances, Component, each, Grouping, link_output_with, perform_steps,
    output_values, output_value, StepInstruction, StepInstructions, calculate_energy_flow,
    highest_temperature, collect_all_energy_potentials_of_interface_balance, 
    collect_all_storage_potentials_of_interface_balance,
    collect_all_temperatures_of_interface_balance

"""
Convenience function to get the value of a key from a config dict using a default value.
"""
default(config::Dict{String,Any}, name::String, default_val::Any)::Any =
    return name in keys(config) ? config[name] : default_val


"""
The number of hours per time step, used by various utility functions.
"""
HOURS_PER_TIME_STEP::Float64 = 0.25

"""
Update the time step, in seconds.
"""
function set_timestep(step_seconds::Integer)
    global HOURS_PER_TIME_STEP = Float64(step_seconds) / 3600.0
end

"""
Calculate energy from power by using the simulation time step.

This function must be assigned a method after simulation parameters (such as the timestep)
have been loaded.
"""
function watt_to_wh(watts::Float64)
    return watts * HOURS_PER_TIME_STEP
end

"""
Categories that each represent a physical medium in conjunction with additional attributes,
such as temperature or voltage. These attributes are not necessarily unchanging, but are
considered the nominal range. For example, a heating component might circulate water anywhere
from 30°C to 60°C, but the nominal temperature is considered to be 60°C. This is intended
so it becomes possible to prevent linking components that do not work together because they
work on different nominal temperatures, while both work with the same physical medium,
for example water.

The names are structured in a composite of segments. For example, these are:
    m_e_ac_230v

    1. m: This segment is used to distinguish its symbols from the symbols of other types
    2. e: The energy type, in this case electricity
    3. ac: The physical medium, in this case AC current
    4. 230v: Additional attributes of nominal value or ranges numbered through
"""
medium_categories = Set{Symbol}([
    # electricity
    :m_e_ac_230v,

    # chemicals - gasses
    :m_c_g_natgas,
    :m_c_g_h2,
    :m_c_g_o2,

    # heat - low temperature water
    :m_h_w_lt1,
    :m_h_w_lt2,
    :m_h_w_lt3,
    :m_h_w_lt4,
    :m_h_w_lt5,

    # heat - high temperature water
    :m_h_w_ht1,
    :m_h_w_ht2,
    :m_h_w_ht3,
    :m_h_w_ht4,
    :m_h_w_ht5,
])

"""
register_media(categories)

Add the given medium categories to the set of all medium categories.
"""
function register_media(categories::Vector{Symbol})
    for cat in categories
        push!(medium_categories, cat)
    end
end

register_media(category::Symbol) = register_media([category,])

"""
Enumerations of the archetype of a component describing its general function.

These are described in more detail in the accompanying documentation of the simulation
model.
"""
@enum(SystemFunction, sf_bounded_sink, sf_bounded_source, sf_fixed_sink,
    sf_fixed_source, sf_transformer, sf_storage, sf_bus)

"""
Enumerations of a simulation step that can be performed on a component.

The names are prefixed with `s` to avoid shadowing functions of the same name.
"""
@enum Step s_reset s_control s_process s_load s_distribute s_potential

"""
Convenience type for holding the instruction for one component and one step.
"""
const StepInstruction = Tuple{String,Step}

"""
Holds the order of steps as instructions for how to perform the simulation.
"""
const StepInstructions = Vector{StepInstruction}

"""
The basic type of all energy system components.

Because Julia does not have field inheritance as OOP languages such as Java do, types
implementing this abstract type are further to be assumed to have the fields required by all
components:
* uac::String: The user address code (UAC) of the component
* controller::Controller: Handles control functionality
* sys_function::SystemFunction: The system function the component has within the energy
    system (storage, transformer, etc.)

Some fields are expected, but behave differently for Bus components. For components of other
system function these are:
* input_interfaces::InterfaceMap: The input interfaces to the component indexed by medium
* output_interfaces::InterfaceMap: The output interfaces to the component indexed by medium
"""
abstract type Component end

"""
Convenience alias to a dict mapping UAC keys to a component.
"""
const Grouping = Dict{String,Component}

"""
Convenience alias for temperatures which can be a number or "nothing".
"""
const Temperature = Union{Nothing,Float64}

"""
Holds the options which output values should be recorded.

This is a specific data structure intended to speed up recording output by avoiding the
need to parse the user-submitted config options for every time step.
"""
Base.@kwdef struct OutputKey
    unit::Component
    medium::Union{Nothing,Symbol}
    value_key::String
end

"""
    each(components :: Grouping)

Generator over each of the components in the given grouping.
"""
function each(components::Grouping)::Base.ValueIterator
    return values(components)
end

"""
Handles the tracking of energy being transfered from one component to another.

This abstraction is useful to avoid components having to "know" other types. Instead of
calling functions to transfer energy, a component can draw from or load into a SystemInterface
instance instead.

Energy is considered to always flow from the source to the target. A negative balance is a
lack of energy that needs to be covered in order for the energy to be balanced.

A system interface keeps track of how much energy (in absolute terms) was transfered
via the interface. Assuming the energy balance holds, at the end of a time step a system
interface's field "sum_abs_change" will have a value of twice the total energy transfered.
"""
Base.@kwdef mutable struct SystemInterface
    """The source component providing energy"""
    source::Union{Nothing,Component} = nothing

    """The target component receiving energy"""
    target::Union{Nothing,Component} = nothing

    """The current balance of the interface"""
    balance::Float64 = 0.0

    """The sum of absolute changes to the interface's balance"""
    sum_abs_change::Float64 = 0.0

    """Current temperature of the medium on this interface"""
    temperature::Temperature = nothing

    """Maximum energy the source can provide in the current timestep"""
    max_energy::Union{Nothing,Float64} = nothing
end

"""
    add!(interface, change, temperature)

Add the given amount of energy (in Wh) to the balance of the interface.
"""
function add!(
    interface::SystemInterface,
    change::Float64,
    temperature::Temperature=nothing
)
    interface.balance += change
    interface.sum_abs_change += abs(change)

    if temperature !== nothing
        if interface.temperature !== nothing && interface.temperature < temperature
            println("Warning: Temperatures are chilled on interface $(interface.source.uac)"
            * " -> $(interface.target.uac) from $temperature to $(interface.temperature) °C")
        elseif interface.temperature !== nothing && interface.temperature > temperature
            println("Warning: Temperature in interface is higher than delivered temperature on interface $(interface.source.uac)"
            * " -> $(interface.target.uac). Requested: $(interface.temperature) °C; Delivered: $temperature °C")
        end
        interface.temperature = temperature
    end
end

"""
    sub!(interface, change, temperature)

Subtract the given amount of energy (in Wh) from the balance of the interface.
"""
function sub!(
    interface::SystemInterface,
    change::Float64,
    temperature::Temperature=nothing
)
    interface.balance -= change
    interface.sum_abs_change += abs(change)

    if temperature !== nothing
        if interface.temperature !== nothing && interface.temperature < temperature
            println("Warning: Temperature in interface is lower than requested temperature on interface $(interface.source.uac)"
                    * " -> $(interface.target.uac). Requested: $temperature °C; Delivered: $(interface.temperature) °C")
        elseif interface.temperature !== nothing && interface.temperature > temperature
            println("Warning: Temperatures are chilled on interface $(interface.source.uac)"
            * " -> $(interface.target.uac) from $(interface.temperature) to $temperature °C")
        end
        interface.temperature = temperature
    end
end

"""
    set!(interface, new_val, temperature)

Set the balance of the interface to the given new value.
"""
function set!(
    interface::SystemInterface,
    new_val::Float64,
    temperature::Temperature=nothing
)
    interface.sum_abs_change += abs(interface.balance - new_val)
    interface.balance = new_val
    interface.temperature = temperature
end

"""
    set_max_energy!(interface, value)

Set the maximum power that can be delivered to the given value.
"""
function set_max_energy!(
    interface::SystemInterface,
    value::Union{Nothing,Float64}
)
    interface.max_energy = value
end

"""
    reset!(interface)

Reset the interface back to zero.
"""
function reset!(interface::SystemInterface)
    interface.balance = 0.0
    interface.sum_abs_change = 0.0
    interface.temperature = nothing
    interface.max_energy = nothing
end


"""
    highest_temperature(temperature_1, temperature_2)::Temperature

Returns the highest temperature::Temperature and handles nothing-values:
- If both of the inputs are floats, the maximum will be returned.
- If one of the inputs is nothing and one a float, the float will be returned.
- If both of the inputs are nothing, nothing will be returned.
"""
function highest_temperature(temperature_1, temperature_2)::Temperature
    if temperature_1 !== nothing && temperature_2 !== nothing
        return max(temperature_1, temperature_2)
    elseif temperature_1 === nothing && temperature_2 !== nothing
        return temperature_2
    elseif temperature_1 !== nothing && temperature_2 === nothing
        return temperature_1
    end
end


"""
    collect_all_energy_potentials_of_interface_balance(input::Any)

Input is exchange.energy, an array of touples or a single touple as
return value of the balance_on().

Returns a vector containing all energy_potentials listed in the touples.
"""
function collect_all_energy_potentials_of_interface_balance(input::Any)
    energy_tuples = isa(input, NamedTuple) ? [input] : input
    return [t.energy_potential for t in energy_tuples]
end


"""
    collect_all_storage_potentials_of_interface_balance(input::Any)

Input is exchange.energy, an array of touples or a single touple as
return value of the balance_on().

Returns a vector containing all storage_potentials listed in the touples.
"""
function collect_all_storage_potentials_of_interface_balance(input::Any)
    energy_tuples = isa(input, NamedTuple) ? [input] : input
    return [t.storage_potential for t in energy_tuples]
end


"""
collect_all_temperatures_of_interface_balance(input::Any)

Input is exchange.energy, an array of touples or a single touple as
return value of the balance_on().

Returns a vector containing all temperatures listed in the touples.
"""
function collect_all_temperatures_of_interface_balance(input::Any)
    energy_tuples = isa(input, NamedTuple) ? [input] : input
    return [t.temperature for t in energy_tuples]
end


"""
Convenience type used define the required system interfaces of a component.

To simultaneously define what is required and then hold references to instances after the
whole energy system has been loaded, it maps a medium category to either nothing (before
components are linked) or a SystemInterface instance.
"""
const InterfaceMap = Dict{Symbol,Union{Nothing,SystemInterface}}

"""
    balance_on(interface, unit)

Return the balance of a unit in respect to the given interface.

For most components this is simply the balance of the interface itself, but for Bus
instances this is handled differently. This function helps to implement a component
without having to check if its connected to a Bus or directly to a component.

# Arguments
- `interface::SystemInterface`: The interface "on which" the balance is calculated. This
    also defines which component is the source.
- `unit::Component`: The receiving component

# Returns NamedTuple with
- "balance"::Float64:           The balance of the target component that can be considered a 
                                demand on the source component
- "storage_potential"::Float64: An additional demand that covers the free storage space 
                                connected to the target component
- "energy_potential"::Float64:  The maximum enery an interface can provide or consume
- "temperature"::Temperature:   The temperature of the interface
"""
function balance_on(
    interface::SystemInterface,
    unit::Component
)::NamedTuple{}
    balance_written = interface.max_energy === nothing || interface.sum_abs_change > 0.0
    input_sign = unit.uac == interface.target.uac ? -1 : +1

    return (
            balance = interface.balance,
            energy = (  uac=unit.uac,
                        energy_potential=(balance_written ? 0.0 : input_sign * interface.max_energy),
                        storage_potential=0.0,
                        temperature=interface.temperature,
                        pressure=nothing,
                        voltage=nothing),
            )
end

"""
    balance(unit)

Calculate the energy balance of the given unit as a whole.

This is expected to start at zero at the beginning of a time step and return to zero at
the end of it. If it is not zero, either the simulation failed to correctly calculate the
energy balance of the energy system or the simulated network was not able to ensure the
balance on the current time step. In either case, something went wrong.
"""
function balance(unit::Component)::Float64
    balance = 0.0

    for inface in values(unit.input_interfaces)
        if inface !== nothing
            balance += inface.balance
        end
    end

    for outface in values(unit.output_interfaces)
        if outface !== nothing
            balance += outface.balance
        end
    end

    return balance    
end

"""
    reset(unit)

Reset the given component back to zero.

For most components this only resets the balances on the system interfaces but some
components might require more complex reset handling.
"""
function reset(unit::Component)
    for inface in values(unit.input_interfaces)
        if inface !== nothing
            reset!(inface)
        end
    end
    for outface in values(unit.output_interfaces)
        if outface !== nothing
            reset!(outface)
        end
    end
end

"""
    control(unit, components, parameters)

Perform the control calculations for the given component.

# Arguments
- `unit::Component`: The component for which control is handled
- `components::Grouping`: A reference dict to all components in the project
- `parameters::Dict{String, Any}`: Project-wide parameters
"""
function control(
    unit::Component,
    components::Grouping,
    parameters::Dict{String,Any}
)
    move_state(unit, components, parameters)
end

"""
    potential(unit, parameters)

Calculate potential energy processing for the given component.

# Arguments
- `unit::Component`: The component for which potentials are calculated
- `parameters::Dict{String, Any}`: Project-wide parameters
"""
function potential(
    unit::Component,
    parameters::Dict{String,Any}
)
    # default implementation is to do nothing
end

"""
    process(unit, parameters)

Perform the processing calculations for the given component.

# Arguments
- `unit::Component`: The component that is processed
- `parameters::Dict{String, Any}`: Project-wide parameters
"""
function process(
    unit::Component,
    parameters::Dict{String,Any}
)
    # default implementation is to do nothing
end

"""
    load(unit, parameters)

Load excess energy into storage components.

For non-storage components this function does nothing.

# Arguments
- `unit::Component`: The storage loading excess energy
- `parameters::Dict{String, Any}`: Project-wide parameters
"""
function load(
    unit::Component,
    parameters::Dict{String,Any}
)
    # default implementation is to do nothing
end

"""
    distribute(unit)

Distribute the energy inputs and outputs of a bus.

For non-bus components this function does nothing.
"""
function distribute!(unit::Component)
    # default implementation is to do nothing
end

"""
    output_values(unit)

Specify which outputs a component can provide.

For the special values "IN" and "OUT" a medium category is required for fetching the actual
value, while this method only specifies that there is an input or output.
"""
function output_values(unit::Component)::Vector{String}
    return ["IN", "OUT"]
end

"""
    calculate_energy_flow(interface)

Calculates the energy flow in an interface and returns the energy.
If the balance in an interface was not zero, the requested amount of energy and not
the acutal energy that has been provided is returned! 
This is conformily of the output interfaces of busses to non-busses which are set
to the requested energy, not regarding the energy that has been delivered!

Replacing the - below by a +, this can be changed, then the delivered energy and not
the requested energy is returned. Comparing this to the bahaviour of busses, this can 
be non-intuitive.
"""
function calculate_energy_flow(interface::SystemInterface)::Float64
    return (interface.sum_abs_change - interface.balance) / 2
end

"""
    output_value(unit, key)

Return the value for the output with the given output key.

Note that for the "IN" and "OUT" output values, the value corresponds to the sum of
absolute changes of the system interfaces and divided by 2. This behaviour is part of the
expected use of the method.

Args:
- `unit::Component`: The component for which to fetch the output
- `key::OutputKey`: An OutputKey specifying which output to return. This should be one of
    the options provided by `output_values()` as well as "IN" or "OUT"
Returns:
- `Float64`: The value of the desired output
Throws:
- `KeyError`: The key value requested must be one the component can provide
"""
function output_value(unit::Component, key::OutputKey)::Float64
    if key.value_key == "IN"
        return calculate_energy_flow(unit.input_interfaces[key.medium])
    elseif key.value_key == "OUT"
        return calculate_energy_flow(unit.output_interfaces[key.medium])
    end
    throw(KeyError(key.value_key))
end

# for the moment control must be an include as it contains circular dependencies
# of definitions
include("control.jl")

using ..Profiles

# the order of includes of the individual components matters here as some components
# require the definition of certain basic components such as a bus or a grid connection
include("general/fixed_sink.jl")
include("general/fixed_supply.jl")
include("general/bounded_supply.jl")
include("general/bounded_sink.jl")
include("connections/grid_connection.jl")
include("connections/bus.jl")
include("storage/battery.jl")
include("storage/buffer_tank.jl")
include("storage/seasonal_thermal_storage.jl")
include("heat_sources/geothermal_probes.jl")
include("heat_sources/geothermal_heat_collectors.jl")
include("electric_producers/chpp.jl")
include("others/electrolyser.jl")
include("heat_producers/fuel_boiler.jl")
include("heat_producers/heat_pump.jl")
include("electric_producers/pv_plant.jl")

load_condition_prototypes()

"""
    link_output_with(unit, components)

Set the output targets of the given unit to the given components.

This function is used to construct the network of components from a graph input that
determines which components provide energy to which other components.

# Arguments
- `unit::Component`: The unit providing energy
- `components::Grouping`: A set of components receiving energy. As components might have multiple
    outputs, this is used to set them all at once.
"""
function link_output_with(unit::Component, components::Grouping)
    if isa(unit, Bus)
        for component in each(components)
            if isa(component, Bus)
                if unit.medium == component.medium
                    connection = SystemInterface(source=unit, target=component)
                    push!(component.input_interfaces, connection)
                    push!(unit.output_interfaces, connection)
                end
            else
                for in_medium in keys(component.input_interfaces)
                    if in_medium == unit.medium
                        connection = SystemInterface(source=unit, target=component)
                        component.input_interfaces[in_medium] = connection
                        push!(unit.output_interfaces, connection)
                    end
                end
            end
        end
    else
        for out_medium in keys(unit.output_interfaces)
            for component in each(components)
                if isa(component, Bus)
                    if out_medium == component.medium
                        connection = SystemInterface(source=unit, target=component)
                        push!(component.input_interfaces, connection)
                        unit.output_interfaces[out_medium] = connection
                    end
                else
                    for in_medium in keys(component.input_interfaces)
                        if out_medium == in_medium
                            connection = SystemInterface(source=unit, target=component)
                            unit.output_interfaces[out_medium] = connection
                            component.input_interfaces[in_medium] = connection
                        end
                    end
                end
            end
        end
    end
end

"""
    check_balances(components, epsilon)

Check the energy balance of the given components and return warnings of any violations.

# Arguments
- `component::Grouping`: The components to check
- `epsilon::Float64`: A balance is only considered violated if the absolute value of the
    sum is larger than this value. This helps with spurious floating point issues

# Returns
- `Vector{Tuple{String, Float64}}`: A list of tuples, where each tuple is the key of the
    component that has a non-zero energy balance and the value of that balance.
"""
function check_balances(
    components::Grouping,
    epsilon::Float64
)::Vector{Tuple{String,Float64}}
    warnings = []

    for (key, unit) in pairs(components)
        unit_balance = balance(unit)
        if unit_balance > epsilon || unit_balance < -epsilon
            push!(warnings, (key, unit_balance))
        end
    end

    return warnings
end

"""
    perform_steps(components, order_of_operations, parameters)

Perform the simulation steps of one time step for the given components in the given order.

# Arguments
- `components::Grouping`: The entirety of the components
- `order_of_operations::Vector{Vector{Any}}`: Defines which steps are performed in which order.
    Each component must go through the simulation steps defined in EnergySystems.Step, but the
    order is not the same for all simulations. Determining the order must be handled
    elsewhere, as this function only goes through and calls the appropriate functions. The
    first item of each entry must be the key of the component for which the following steps
    are performed.
- `parameters::Dict{String, Any}`: Project-wide parameters

# Examples
```
    components = Grouping(
        "component_a" => Component(),
        "component_b" => Component(),
    )
    order = [
        ["component_a", EnergySystems.s_control]
        ["component_b", EnergySystems.s_control, EnergySystems.s_process]
        ["component_a", EnergySystems.s_process]
    ]
    parameters = Dict{String, Any}("time" => 0)
    perform_steps(components, order, parameters)
```
In this example the control of component A is performed first, then control and processing of
component B and finally processing of component A.
"""
function perform_steps(
    components::Grouping,
    order_of_operations::StepInstructions,
    parameters::Dict{String,Any}
)
    for entry in order_of_operations
        unit = components[entry[1]]
        step = entry[2]

        if step == s_reset
            reset(unit)
        elseif step == s_control
            control(unit, components, parameters)
        elseif step == s_potential
            potential(unit, parameters)
        elseif step == s_process
            process(unit, parameters)
        elseif step == s_load
            load(unit, parameters)
        elseif step == s_distribute
            distribute!(unit)
        end
    end
end

end