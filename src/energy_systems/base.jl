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
       highest, default, plot_optional_figures_begin, plot_optional_figures_end

using ..Profiles

"""
Convenience function to get the value of a key from a config dict using a default value.
"""
default(config::Dict{String,Any}, name::String, default_val::Any)::Any = return name in keys(config) ?
                                                                                config[name] :
                                                                                default_val

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
                                 :m_h_w_ht5])

"""
register_media(categories)

Add the given medium categories to the set of all medium categories.
"""
function register_media(categories::Vector{Symbol})
    for cat in categories
        push!(medium_categories, cat)
    end
end

register_media(category::Symbol) = register_media([category])

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
Convenience alias for a float that can also have a value of "nothing".
"""
const Floathing = Union{Nothing,Float64}

"""
Convenience alias for a string that can also have a value of "nothing".
"""
const Stringing = Union{Nothing,String}

"""
Struct for the max_energy in the SystemInterface of in the bus balance_table. 
Can handle multiple maximum energies at different temperatures or for different components, 
also in the case when for each max_energy the maximum energy is calculated by the components. 
Then, has_calculated_all_maxima is set to true.
"""
mutable struct MaxEnergy
    max_energy::Vector{Floathing}
    has_calculated_all_maxima::Bool
    purpose_uac::Vector{Stringing}
end

MaxEnergy() = MaxEnergy([nothing], false, [])

"""
Copy function for custom type MaxEnergy
"""
function Base.copy(original::MaxEnergy)
    return MaxEnergy(copy(original.max_energy),
                     copy(original.has_calculated_all_maxima),
                     copy(original.purpose_uac))
end

"""
Copy function for type Nothing
"""
function Base.copy(value::Nothing)
    return value
end

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

    """Minimum temperature of the medium on this interface"""
    temperature_min::Temperature = nothing

    """Maximum temperature of the medium on this interface"""
    temperature_max::Temperature = nothing

    """Maximum energy the source can provide in the current timestep"""
    max_energy::MaxEnergy = MaxEnergy()

    """Flag to decide if storage potentials are transferred over the interface."""
    do_storage_transfer::Bool = true
end

"""
Custom error handler for exception "InputError".
Call with throw(InputError)
"""
struct InputError <: Exception end

"""
    set_storage_transfer!(interface, value)

Sets the flag to decide over storage potential transfer to the given boolean value. Note
that if the flag is already set to false, no further changes will set it back to true. This
is to prevent components overwriting the values of others.
"""
function set_storage_transfer!(interface::SystemInterface, flag::Bool)
    interface.do_storage_transfer = interface.do_storage_transfer && flag
end

"""
The following functions are effectively arithmetic operations on floats where one or both
operands may be nothing. It is possible to overwrite the typical operators + and -,
however this was deemed too dangerous.
"""
function _sub(first::Float64, second::Float64)
    return first - second
end
function _sub(first::Nothing, second::Float64)
    return -second
end
function _sub(first::Float64, second::Nothing)
    return first
end
function _sub(first::Nothing, second::Nothing)
    return nothing
end

function _add(first::Float64, second::Float64)
    return first + second
end
function _add(first::Nothing, second::Float64)
    return second
end
function _add(first::Float64, second::Nothing)
    return first
end
function _add(first::Nothing, second::Nothing)
    return nothing
end

function _abs(val::Union{Floathing,Vector{<:Floathing}})
    if !isa(val, AbstractVector)
        val = [val]
    end
    abs_val = Vector{Floathing}(collect(x === nothing ? nothing : abs(x) for x in val))
    return abs_val
end

function _sum(vector::Union{Floathing,Vector{<:Floathing}})
    sum = nothing
    for entry in vector
        sum = _add(sum, entry)
    end
    return sum
end

function _mean(vector::Union{Floathing,Vector{<:Floathing}})
    number_of_values = count(x -> x !== nothing, vector)
    sum = _sum(vector)
    if number_of_values > 0 && sum !== nothing
        return sum / number_of_values
    else
        return nothing
    end
end

function _weighted_mean(values::Union{Floathing,Vector{<:Floathing}},
                        weights::Union{Floathing,Vector{<:Floathing}})
    if isempty(values) || isempty(weights)
        return 0.0
    end

    if !isa(values, Vector)
        values = [values]
    end
    if !isa(weights, Vector)
        weights = [weights]
    end

    valid_weights = filter(!isnothing, weights)
    normalized_weights = valid_weights ./ _sum(valid_weights)
    valid_values = filter(!isnothing, values)

    if length(valid_values) != length(normalized_weights)
        return 0.0
    end

    return _sum(valid_values .* normalized_weights)
end

function _isless(first::Nothing, second::Nothing)
    return false
end
function _isless(first::Float64, second::Nothing)
    return false
end
function _isless(first::Nothing, second::Float64)
    return true
end
function _isless(first::Float64, second::Float64)
    return first < second
end

"""
    check_epsilon(value, sim_params)

Function to check if a value is withing +/- epsilon.
Returns 0.0 if value is within the epsilon boundaries and
the original value if value is out of the boundaries.
"""
function check_epsilon(value::Float64, sim_params::Dict{String,Any})
    if abs(value) < sim_params["epsilon"]
        return 0.0
    else
        return value
    end
end

"""
    add!(interface, change, temperature)

Add the given amount of energy (in Wh) to the balance of the interface.

If the source or target of the interface is a bus, communicates the change to the bus.
"""
function add!(interface::SystemInterface,
              change::Union{Floathing,Vector{<:Floathing}},
              temperature::Temperature=nothing,
              purpose_uac::Union{Stringing,Vector{Stringing}}=nothing)
    interface.balance += sum(change)
    interface.sum_abs_change += sum(abs.(change))

    if temperature !== nothing
        if interface.temperature_min !== nothing && temperature < interface.temperature_min
            @warn ("Given temperature $temperature on interface $(interface.source.uac) " *
                   "-> $(interface.target.uac) lower than minimum $(interface.temperature_min)")
        elseif interface.temperature_max !== nothing && temperature > interface.temperature_max
            @warn ("Given temperature $temperature on interface $(interface.source.uac) " *
                   "-> $(interface.target.uac) higher than maximum $(interface.temperature_max)")
        end
    end

    if interface.source.sys_function == sf_bus
        add_balance!(interface.source, interface.target, false, change, purpose_uac)
    elseif interface.target.sys_function == sf_bus
        add_balance!(interface.target, interface.source, true, change, purpose_uac)
    end
end

"""
    sub!(interface, change, temperature)

Subtract the given amount of energy (in Wh) from the balance of the interface.

If the source or target of the interface is a bus, communicates the change to the bus.
"""
function sub!(interface::SystemInterface,
              change::Union{Floathing,Vector{<:Floathing}},
              temperature::Temperature=nothing,
              purpose_uac::Union{Stringing,Vector{Stringing}}=nothing)
    interface.balance -= sum(change)
    interface.sum_abs_change += sum(abs.(change))

    if temperature !== nothing
        if interface.temperature_min !== nothing && temperature < interface.temperature_min
            @warn ("Given temperature $temperature on interface $(interface.source.uac) " *
                   "-> $(interface.target.uac) lower than minimum $(interface.temperature_min)")
        elseif interface.temperature_max !== nothing && temperature > interface.temperature_max
            @warn ("Given temperature $temperature on interface $(interface.source.uac) " *
                   "-> $(interface.target.uac) higher than maximum $(interface.temperature_max)")
        end
    end

    if interface.source.sys_function == sf_bus
        sub_balance!(interface.source, interface.target, false, change, purpose_uac)
    elseif interface.target.sys_function == sf_bus
        sub_balance!(interface.target, interface.source, true, change, purpose_uac)
    end
end

"""
    set!(interface, new_val, temperature)

Set the balance of the interface to the given new value.
"""
function set!(interface::SystemInterface,
              new_val::Float64)
    interface.sum_abs_change += abs(interface.balance - new_val)
    interface.balance = new_val
end

"""
    set_temperature!(interface, temp_min, temp_max)

Set the minimum and maximum temperature for the interface.

If the source or target of the interface is a bus, also calls the set_temperature!
function on that bus.
"""
function set_temperature!(interface::SystemInterface,
                          temperature_min::Temperature=nothing,
                          temperature_max::Temperature=nothing)
    interface.temperature_min = highest(interface.temperature_min, temperature_min)
    interface.temperature_max = lowest(interface.temperature_max, temperature_max)

    if interface.source.sys_function == sf_bus
        set_temperatures!(interface.source, interface.target, false,
                          temperature_min, temperature_max)
    elseif interface.target.sys_function == sf_bus
        set_temperatures!(interface.target, interface.source, true,
                          temperature_min, temperature_max)
    end
end

"""
    set_max_energy!(interface, value, purpose_uac, has_calculated_all_maxima)

Set the maximum energy in the `interface` to the given `value` representing the maximum 
energy that the calling component can deliver.
For 1-to-1 connection between two components, the minimum value of the current max_energy 
and the given `value` is written to the interface.

If the source or target of the interface is a bus, the set_max_energy! function
on that bus is called as well to write the max_energy into the balance table of the bus.

The optional `purpose_uac` can be set by the calling component to define the component 
that is supposed to receive or deliver the given energy. This makes not much sense on 1-to-1 
interfaces, but is set for uniformity as well.

Both `value` and `purpose_uac` can be given either as scalars or as vectors, defining multiple
maximum energies for multiple other components (purposes). The given `purpose_uac` do not have
to be necessarily unique. 

The flag `has_calculated_all_maxima` can be set to true if the calling component has not known
the combination of temperature and energy of a source or a sink during its pre-calculation. 
Then, the calling component can calculate the maximum energy that could potentially be taken 
or delivered for each source or for each sink seperately, ignoring all other inputs or outputs. 
Setting the flag, the max_energy is then later handeled accordingly by reduce_max_energy!(),
considering that the single values can not be summed up but they all have to be linearly
decreased if one of them is reduced. An analogy would be to divide the current time step 
into smaller time steps, in each of which a different component is supplied or drawn with
a subset of the maximum possible energy of the current time step.
"""
function set_max_energy!(interface::SystemInterface,
                         value::Union{Floathing,Vector{<:Floathing}},
                         purpose_uac::Union{Stringing,Vector{Stringing}}=nothing,
                         has_calculated_all_maxima::Bool=false)
    if interface.source.sys_function == sf_bus
        set_max_energy!(interface.max_energy, value, purpose_uac, has_calculated_all_maxima)
        set_max_energy!(interface.source,
                        interface.target,
                        false,
                        value,
                        purpose_uac,
                        has_calculated_all_maxima)
    elseif interface.target.sys_function == sf_bus
        set_max_energy!(interface.max_energy, value, purpose_uac, has_calculated_all_maxima)
        set_max_energy!(interface.target,
                        interface.source,
                        true,
                        value,
                        purpose_uac,
                        has_calculated_all_maxima)
    else
        # 1-to-1 interface between two components 
        # --> check if value is smaller than current max_energy
        if is_max_energy_nothing(interface.max_energy) || _isless(_sum(value), get_max_energy(interface.max_energy))
            set_max_energy!(interface.max_energy, value, purpose_uac, has_calculated_all_maxima)
        end
    end
end

"""
    get_max_energy(max_energy, purpose_uac)

This function extracts the `max_energy` values from a given `MaxEnergy` struct.
If `purpose_uac` is provided and found within `max_energy.purpose_uac`, it returns the sum
of all `max_energy` values corresponding to the specific `purpose_uac`.
If `purpose_uac` is provided but not found, the function returns 0.0.
If `purpose_uac` is not provided or if `max_energy.purpose_uac` is empty, the sum of 
all `max_energy` values is returned.
"""
function get_max_energy(max_energy::EnergySystems.MaxEnergy, purpose_uac::Stringing=nothing)
    if purpose_uac === nothing || is_purpose_uac_nothing(max_energy)
        return max_energy_sum(max_energy)
    else
        idx = findall(==(purpose_uac), max_energy.purpose_uac)
        if isempty(idx)
            return 0.0
        else
            return sum(max_energy.max_energy[i] for i in idx)
        end
    end
end

"""
    reduce_max_energy!(max_energy, amount, uac_to_reduce)

This function decreases the maximum energy values stored in `max_energy` by the given `amount`.
If `uac_to_reduce` is specified and found in `max_energy.purpose_uac`, only the corresponding
max_energy value will be reduced. 
f the maximum possible energy for all max_energy have been calculated (`has_calculated_all_maxima` 
is true), the other max_energy values will be reduced proportionally by a linear percentage. 
If `uac_to_reduce` could not be found, an error arises.

If no `uac_to_reduce` is specified of if `max_energy.purpose_uac` is empty, the amount 
decreases the entries in max_energy.max_energy, starting with the first one, until the amount 
is fully distributed.
"""
function reduce_max_energy!(max_energy::EnergySystems.MaxEnergy,
                            amount::Union{Float64,Vector{Float64}},
                            uac_to_reduce::Union{Stringing,Vector{Stringing}}=nothing)
    for amount_idx in eachindex(amount)
        current_amount = length(amount) == 1 ? amount : amount[amount_idx]
        if uac_to_reduce === nothing || is_purpose_uac_nothing(max_energy)
            for i in eachindex(max_energy.max_energy)
                to_reduce = min(current_amount, max_energy.max_energy[i])
                max_energy.max_energy[i] -= to_reduce
                current_amount -= to_reduce
                if current_amount <= 0.0
                    break
                end
            end
        else
            currrent_uac_to_reduce = isa(uac_to_reduce, AbstractVector) ? uac_to_reduce[amount_idx] : uac_to_reduce
            idx = findfirst(==(currrent_uac_to_reduce), max_energy.purpose_uac)
            if idx !== nothing
                if max_energy.has_calculated_all_maxima
                    max_energy.max_energy .*= 1 - (current_amount / max_energy.max_energy[idx])
                else
                    max_energy.max_energy[idx] -= current_amount
                end
            else
                @error "The uac could not be found in the max_energy."
            end
        end
    end
end

"""
    increase_max_energy!(max_energy, amount, uac_to_reduce)

This function increases the maximum energy values stored in `max_energy` by the given `amount`.
If `uac_to_increase` is specified and found in `max_energy.purpose_uac`, the corresponding
max_energy value will be increased. 
If `uac_to_increase` could not be found, a new value in max_energy is added.
If no `uac_to_increase` is given of if `max_energy.purpose_uac` is empty, the given amount is added
to the first max_energy.max_energy as in this case, only one value should be present (e.g 1-to-1 
connection or components that don't care).
"""
function increase_max_energy!(max_energy::EnergySystems.MaxEnergy,
                              amount::Union{Float64,Vector{Float64}},
                              uac_to_increase::Union{Stringing,Vector{Stringing}}=nothing)
    for amount_idx in eachindex(amount)
        current_amount = length(amount) == 1 ? amount : amount[amount_idx]
        currrent_uac_to_increase = isa(uac_to_increase, AbstractVector) ? uac_to_increase[amount_idx] : uac_to_increase
        if uac_to_increase === nothing || is_purpose_uac_nothing(max_energy)
            if is_max_energy_nothing(max_energy) || length(max_energy.max_energy) == 0
                set_max_energy!(max_energy, current_amount, currrent_uac_to_increase, false)
            else
                max_energy.max_energy[1] += current_amount
            end
        else
            idx = findfirst(==(currrent_uac_to_increase), max_energy.purpose_uac)
            if idx === nothing
                set_max_energy!(max_energy, current_amount, currrent_uac_to_increase, false, true)
            else
                max_energy.max_energy[idx] += current_amount
            end
        end
    end
end

"""
    max_energy_sum(max_energy)

Returns the sum of all `max_energy` in a MaxEnergy struct while ignoring nothing values.
"""
function max_energy_sum(max_energy::EnergySystems.MaxEnergy)
    if is_max_energy_nothing(max_energy)
        return 0.0
    else
        return sum(e for e in max_energy.max_energy; init=0.0)
    end
end

"""
    set_max_energy!(max_energy,  values, purpose_uac, has_calculated_all_maxima)

Fills a MaxEnergy struct with given values.
If the given values are scalars, they are converted into vectors.
If `values` is empty, a zero is added to ensure accessibility.
"""
function set_max_energy!(max_energy::EnergySystems.MaxEnergy,
                         values::Union{Floathing,Vector{<:Floathing}},
                         purpose_uac::Union{Stringing,Vector{Stringing}},
                         has_calculated_all_maxima::Bool,
                         append::Bool=false)
    if !isa(values, AbstractVector)
        values = [values]
    end
    if !isa(purpose_uac, AbstractVector)
        purpose_uac = [purpose_uac]
    end

    if length(values) == 0
        push!(values, 0.0)
    end

    """
    Overwriting is possible here as either the interface has no max energy yet, or a component
    or the set_max_energy!(interface, ...) has already read out the max energy and has considered 
    it (it can only be equal or smaller) or a component is overwriting its own max_energy.
    """
    if append
        append!(max_energy.max_energy, values)
        append!(max_energy.purpose_uac, purpose_uac)
        max_energy.has_calculated_all_maxima = has_calculated_all_maxima
    else
        max_energy.max_energy = values
        max_energy.purpose_uac = purpose_uac
        max_energy.has_calculated_all_maxima = has_calculated_all_maxima
    end
end

"""
    is_max_energy_nothing(max_energy)

Checks a MaxEnergy struct if the `max_energy` is nothing (true), meaning that no potential has
beed performed (yet).
"""
function is_max_energy_nothing(max_energy::EnergySystems.MaxEnergy)
    return max_energy.max_energy[1] === nothing && length(max_energy.max_energy) == 1
end

"""
    is_purpose_uac_nothing(max_energy)

Checks a MaxEnergy struct if the `purpose_uac` is nothing/empty (true), meaning that the 
max_energy is not intended for a specific component.
"""
function is_purpose_uac_nothing(max_energy::EnergySystems.MaxEnergy)
    return (isempty(max_energy.purpose_uac) ||
            max_energy.purpose_uac[1] === nothing && length(max_energy.purpose_uac) == 1)
end

"""
    reset!(interface)

Reset the interface back to zero.
"""
function reset!(interface::SystemInterface)
    interface.balance = 0.0
    interface.sum_abs_change = 0.0
    interface.temperature_min = nothing
    interface.temperature_max = nothing
    interface.max_energy = MaxEnergy()
end

"""
    highest(temperature_1, temperature_2)::Temperature

Returns the highest temperature of the two inputs and handles nothing-values:
- If both of the inputs are floats, the maximum will be returned.
- If one of the inputs is nothing and one a float, the float will be returned.
- If both of the inputs are nothing, nothing will be returned.
"""
function highest(temperature_1::Temperature,
                 temperature_2::Temperature)::Temperature
    if temperature_1 !== nothing && temperature_2 !== nothing
        return max(temperature_1, temperature_2)
    elseif temperature_1 === nothing && temperature_2 !== nothing
        return temperature_2
    elseif temperature_1 !== nothing && temperature_2 === nothing
        return temperature_1
    end
end

"""
    highest(temperature_vector::Vector{Temperature})::Temperature

Returns the highest temperature of a vector of temperatures and handles nothing-values:
- If all of the inputs are floats, the maximum will be returned.
- If some of the inputs are nothing and the others are float, the highest float will be returned.
- If all of the inputs are nothing, nothing will be returned.
"""
function highest(temperature_vector::Vector{Temperature})::Temperature
    n_temperatures = length(temperature_vector)
    if n_temperatures == 0
        return nothing
    elseif n_temperatures == 1
        return temperature_vector[1]
    else
        current_highest = highest(temperature_vector[1], temperature_vector[2])
        if n_temperatures > 2
            for idx in 3:n_temperatures
                current_highest = highest(current_highest, temperature_vector[idx])
            end
        end
        return current_highest
    end
end

"""
    lowest(temperature_1, temperature_2)::Temperature

Returns the lowest temperature of the two inputs and handles nothing-values:
- If both of the inputs are floats, the minimum will be returned.
- If one of the inputs is nothing and one a float, the float will be returned.
- If both of the inputs are nothing, nothing will be returned.
"""
function lowest(temperature_1::Temperature,
                temperature_2::Temperature)::Temperature
    if temperature_1 !== nothing && temperature_2 !== nothing
        return min(temperature_1, temperature_2)
    elseif temperature_1 === nothing && temperature_2 !== nothing
        return temperature_2
    elseif temperature_1 !== nothing && temperature_2 === nothing
        return temperature_1
    end
end

"""
    lowest(temperature_vector::Vector{Temperature})::Temperature

Returns the lowest temperature of a vector of temperatures and handles nothing-values:
- If all of the inputs are floats, the minimum will be returned.
- If some of the inputs are nothing and the others are float, the lowest float will be returned.
- If all of the inputs are nothing, nothing will be returned.
"""
function lowest(temperature_vector::Vector{Temperature})::Temperature
    n_temperatures = length(temperature_vector)
    if n_temperatures == 0
        return nothing
    elseif n_temperatures == 1
        return temperature_vector[1]
    else
        current_lowest = lowest(temperature_vector[1], temperature_vector[2])
        if n_temperatures > 2
            for idx in 3:n_temperatures
                current_lowest = lowest(current_lowest, temperature_vector[idx])
            end
        end
        return current_lowest
    end
end

"""
Convenience type used define the required system interfaces of a component.

To simultaneously define what is required and then hold references to instances after the
whole energy system has been loaded, it maps a medium category to either nothing (before
components are linked) or a SystemInterface instance.
"""
const InterfaceMap = Dict{Symbol,Union{Nothing,SystemInterface}}

"""
Contains the data on the energy exchance (and related information) on an interface.
"""
Base.@kwdef mutable struct EnergyExchange
    balance::Float64
    energy_potential::Float64
    purpose_uac::Stringing
    temperature_min::Temperature
    temperature_max::Temperature
    pressure::Floathing
    voltage::Floathing
end

"""
Convenience alias to EnergyExchange.
"""
const EnEx = EnergyExchange

"""
    balance(exchanges)

Sum of balances over the given list of energy exchanges.

Args:
    `entries::Vector{EnergyExchange}`: The exchanges to sum over
Returns:
    `Float64`: Sum of balances
"""
function balance(entries::Vector{EnergyExchange})::Float64
    return sum(e.balance for e in entries; init=0.0)
end

"""
    energy_potential(exchanges)

Sum of energy potentials over the given list of energy exchanges.

Args:
    `entries::Vector{EnergyExchange}`: The exchanges to sum over
Returns:
    `Float64`: Sum of energy potentials
"""
function energy_potential(entries::Vector{EnergyExchange})::Float64
    return sum(e.energy_potential for e in entries; init=0.0)
end

"""
    temp_min_first(exchanges)

First not-nothing temperature of the given list of energy exchanges. If no not-nothing
temperature can be found, returns nothing.

Args:
    `entries::Vector{EnergyExchange}`: The exchanges over which to search for a temperature
Returns:
    `Temperature`: First not-nothing temperature found or nothing if no such exists
"""
function temp_min_first(entries::Vector{EnergyExchange})::Temperature
    temps = [e.temperature_min for e in entries if e.temperature_min !== nothing]
    return length(temps) > 0 ? first(temps) : nothing
end

"""
    temp_min_highest(exchanges)

Highest not-nothing minimum temperature of the given list of energy exchanges. If no
not-nothing temperature can be found, returns nothing.

Args:
    `entries::Vector{EnergyExchange}`: The exchanges over which to search for a temperature
Returns:
    `Temperature`: First not-nothing minimum temperature found or nothing if no such exists
"""
function temp_min_highest(entries::Vector{EnergyExchange})::Temperature
    temps = [e.temperature_min for e in entries if e.temperature_min !== nothing]
    return length(temps) > 0 ? maximum(temps) : nothing
end

"""
    temp_max_highest(exchanges)

Highest not-nothing maximum temperature of the given list of energy exchanges. If no
not-nothing temperature can be found, returns nothing.

Args:
    `entries::Vector{EnergyExchange}`: The exchanges over which to search for a temperature
Returns:
    `Temperature`: First not-nothing maximum temperature found or nothing if no such exists
"""
function temp_max_highest(entries::Vector{EnergyExchange})::Temperature
    temps = [e.temperature_max for e in entries if e.temperature_max !== nothing]
    return length(temps) > 0 ? maximum(temps) : nothing
end

"""
    temp_min_all(exchanges)

A list of all minimum temperatures of the given list of energy exchanges.

Args:
    `entries::Vector{EnergyExchange}`: The exchanges over which to list temperatures
Returns:
    `Vector{Temperature}`: A list of minimum temperatures
"""
function temp_min_all(entries::Vector{EnergyExchange})::Vector{Temperature}
    return [e.temperature_min for e in entries]
end

"""
    temp_max_all(exchanges)

A list of all maximum temperatures of the given list of energy exchanges.

Args:
    `entries::Vector{EnergyExchange}`: The exchanges over which to list temperatures
Returns:
    `Vector{Temperature}`: A list of maximum temperatures
"""
function temp_max_all(entries::Vector{EnergyExchange})::Vector{Temperature}
    return [e.temperature_max for e in entries]
end

"""
    temp_min_all_non_empty(exchanges)

A list of all not-nothing minimum temperatures of the given list of energy exchanges.

Args:
    `entries::Vector{EnergyExchange}`: The exchanges over which to list temperatures
Returns:
    `Vector{Temperature}`: A (possibly empty) list of not-nothing minimum temperatures
"""
function temp_min_all_non_empty(entries::Vector{EnergyExchange})::Vector{Temperature}
    return [e.temperature_min for e in entries if e.temperature_min !== nothing]
end

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

# Returns
- `Vector{EnergyExchange}`: A list of energy exchange items, which contain the information
    required to perform the energy flow calculations.
"""
function balance_on(interface::SystemInterface, unit::Component)::Vector{EnergyExchange}
    balance_written = interface.max_energy.max_energy[1] === nothing || interface.sum_abs_change > 0.0
    input_sign = unit.uac == interface.target.uac ? -1 : +1
    purpose_uac = unit.uac == interface.target.uac ? interface.target.uac : interface.source.uac

    return [EnEx(; balance=interface.balance,
                 energy_potential=(balance_written ? 0.0 : input_sign * get_max_energy(interface.max_energy)),
                 purpose_uac=purpose_uac,
                 temperature_min=interface.temperature_min,
                 temperature_max=interface.temperature_max,
                 pressure=nothing,
                 voltage=nothing)]
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
    initialise!(unit)

Perform steps to initialise the component in addition to the constructor.

# Arguments
- `unit::Component`: The component to initialise
- `sim_params::Dict{String,Any}`: Simulation parameters
"""
function initialise!(unit::Component, sim_params::Dict{String,Any})
    # default implementation is to do nothing
end

"""
    reset(unit)

Reset the given component back to zero.

For most components this only resets the losses and the balances on the system interfaces,
but some components might require more complex reset handling like for electrolysers due to
several different losses present.
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
    if hasfield(typeof(unit), Symbol("losses"))
        unit.losses = 0.0
    end
end

"""
    control(unit, components, sim_params)

Perform the control calculations for the given component.

# Arguments
- `unit::Component`: The component for which control is handled
- `components::Grouping`: A reference dict to all components in the project
- `sim_params::Dict{String, Any}`: Project-wide simulation parameters
"""
function control(unit::Component,
                 components::Grouping,
                 sim_params::Dict{String,Any})
    update(unit.controller)
end

"""
    potential(unit, sim_params)

Calculate potential energy processing for the given component.

# Arguments
- `unit::Component`: The component for which potentials are calculated
- `sim_params::Dict{String, Any}`: Project-wide simulation parameters
"""
function potential(unit::Component,
                   sim_params::Dict{String,Any})
    # default implementation is to do nothing
end

"""
    process(unit, sim_params)

Perform the processing calculations for the given component.

# Arguments
- `unit::Component`: The component that is processed
- `sim_params::Dict{String, Any}`: Project-wide simulation parameters
"""
function process(unit::Component,
                 sim_params::Dict{String,Any})
    # default implementation is to do nothing
end

"""
    load(unit, sim_params)

Load excess energy into storage components.

For non-storage components this function does nothing.

# Arguments
- `unit::Component`: The storage loading excess energy
- `sim_params::Dict{String, Any}`: Project-wide simulation parameters
"""
function load(unit::Component,
              sim_params::Dict{String,Any})
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
    plot_optional_figures_begin(unit, output_path, output_formats, sim_params)

Plot optional figures that are potentially created after initialisation of each 
component. Saves all figures to "output_path" for each specified "output_formats".
Possible output formats are:
    - html
    - pdf
    - png
    - ps
    - svg

# Arguments
- `unit::Component`: The unit that plots additional figures
- `output_path::String`: The output folder as string (absolute/relative) for the additional plots
- `output_formats::Vector{Any}`: A Vector of output file formats, each as string without dot
- `sim_params::Dict{String,Any}`: Simulation parameters of ReSiE
 
# Returns:
- Bool: True if a figure was created, false if no figure was created
"""
function plot_optional_figures_begin(unit::Component,
                                     output_path::String,
                                     output_formats::Vector{String},
                                     sim_params::Dict{String,Any})
    # default implementation is to do nothing
    return false
end

"""
    plot_optional_figures_end(unit, sim_params)

Plot optional figures that are potentially created at the end of the simulation for each 
component.

# Arguments
- `unit::Component`: The unit that plots additional figures
- `sim_params::Dict{String,Any}`: Simulation parameters of ReSiE
 
# Returns:
- Bool: True if a figure was created, false if no figure was created
"""
function plot_optional_figures_end(unit::Component,
                                   sim_params::Dict{String,Any})
    # default implementation is to do nothing
    return false
end

"""
    output_values(unit)

Specify which data outputs a component can provide, including the medium of each output.
Output at this point means data output in every timestep. This can include the energy on 
the input- and output interfaces of a unit or its current state (like "LOAD").

This methods provides the actual output type (like "IN" or "OUT") and the corresponding
media of the data outputs. A medium is only needed for inputs and outputs, not for states.
"""
function output_values(unit::Component)::Vector{String}
    return []  # base implementation returns an empty output vector as the output values 
    # have to be specified in every component.
end

"""
    calculate_energy_flow(interface)

Calculates the energy flow in an interface and returns the energy.
If the balance in an interface was not zero, the actual transferred energy is returned.

"""
function calculate_energy_flow(interface::SystemInterface)::Float64
    return (interface.sum_abs_change - abs(interface.balance)) / 2
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

# the control functionality cannot be its own module as it depends on definitions of
# components and other members of the EnergySystems module. this file contains the general
# definitions, the modules are loaded later
include("control.jl")

# the order of includes of the individual components matters here as some components
# require the definition of certain basic components such as a bus or a grid connection
include("general/fixed_sink.jl")
include("general/fixed_supply.jl")
include("general/bounded_supply.jl")
include("general/bounded_sink.jl")
include("general/storage.jl")
include("connections/grid_connection.jl")
include("connections/bus.jl")
include("storage/battery.jl")
include("storage/buffer_tank.jl")
include("storage/seasonal_thermal_storage.jl")
include("heat_sources/geothermal_probes.jl")
include("heat_sources/geothermal_heat_collectors.jl")
include("heat_sources/generic_heat_source.jl")
include("electric_producers/chpp.jl")
include("others/electrolyser.jl")
include("heat_producers/fuel_boiler.jl")
include("heat_producers/heat_pump.jl")
include("electric_producers/pv_plant.jl")

# additional functionality applicable to multiple component types, that belongs in the
# base module and has been moved into seperate files for less clutter
include("efficiency.jl")

# now that the components are defined we can load the control modules, which depend on their
# definitions
const StorageComponent = Union{Battery,BufferTank,SeasonalThermalStorage,Storage}

include("control_modules/economical_discharge.jl")
include("control_modules/profile_limited.jl")
include("control_modules/storage_driven.jl")
include("control_modules/temperature_sorting.jl")

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
                    connection = SystemInterface(; source=unit, target=component)
                    push!(component.input_interfaces, connection)
                    push!(unit.output_interfaces, connection)
                end
            else
                for in_medium in keys(component.input_interfaces)
                    if in_medium == unit.medium
                        connection = SystemInterface(; source=unit, target=component)
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
                        connection = SystemInterface(; source=unit, target=component)
                        push!(component.input_interfaces, connection)
                        unit.output_interfaces[out_medium] = connection
                    end
                else
                    for in_medium in keys(component.input_interfaces)
                        if out_medium == in_medium
                            connection = SystemInterface(; source=unit, target=component)
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
    initialise_components(component, sim_params)

Initialises the given components.

# Arguments
`components::Grouping`: The components
`sim_params::Dict{String,Any}`: Simulation parameters
"""
function initialise_components(components::Grouping, sim_params::Dict{String,Any})
    for component in [c for c in values(components) if c.sys_function !== sf_bus]
        initialise!(component, sim_params)
    end
    # initialise non-bus components first, so they write do_storage_transfer flags on the
    # interfaces before the busses are initialised, using these flags for their own input
    for component in [c for c in values(components) if c.sys_function === sf_bus]
        initialise!(component, sim_params)
    end
end

"""
    merge_bus_chains(chains, components, sim_params)

For each of the given bus chains, creates a proxy bus as a merger of all busses in the
chain and sets the proxy.

Each chain is a set of interconnected busses of the same medium. For an example
@See EnergySystems.find_chains which can identify the chains to merge.

# Arguments
`chains::Vector{Set{Component}}`: A list of chains of busses
`components::Grouping`: All components of the energy systems
`sim_params::Dict{String,Any}`: Simulation parameters
"""
function merge_bus_chains(chains::Vector{Set{Component}},
                          components::Grouping,
                          sim_params::Dict{String,Any})
    for chain in chains
        if length(chain) < 2
            continue
        end

        comp_as_grouping = Grouping(comp.uac => comp for comp in chain)
        merged = merge_busses(comp_as_grouping, components)
        for bus in chain
            bus.proxy = merged
        end
        components[merged.uac] = merged
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
function check_balances(components::Grouping,
                        epsilon::Float64)::Vector{Tuple{String,Float64}}
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
    perform_steps(components, order_of_operations, sim_params)

Perform the simulation steps of one time step for the given components in the given order.

# Arguments
- `components::Grouping`: The entirety of the components
- `order_of_operations::Vector{Vector{Any}}`: Defines which steps are performed in which order.
    Each component must go through the simulation steps defined in EnergySystems.Step, but the
    order is not the same for all simulations. Determining the order must be handled
    elsewhere, as this function only goes through and calls the appropriate functions. The
    first item of each entry must be the key of the component for which the following steps
    are performed.
- `sim_params::Dict{String, Any}`: Project-wide simulation parameters

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
    sim_params = Dict{String, Any}("time" => 0)
    perform_steps(components, order, sim_params)
```
In this example the control of component A is performed first, then control and processing of
component B and finally processing of component A.
"""
function perform_steps(components::Grouping,
                       order_of_operations::StepInstructions,
                       sim_params::Dict{String,Any})
    for entry in order_of_operations
        unit = components[entry[1]]
        step = entry[2]

        if step == s_reset
            reset(unit)
        elseif step == s_control
            control(unit, components, sim_params)
        elseif step == s_potential
            potential(unit, sim_params)
        elseif step == s_process
            process(unit, sim_params)
        elseif step == s_load
            load(unit, sim_params)
        elseif step == s_distribute
            distribute!(unit)
        end
    end
end

"""
    get_parameter_profile_from_config(config,
                                      sim_params,
                                      param_symbol,
                                      profile_file_key,
                                      from_global_file_key,
                                      constant_key,
                                      uac;
                                      required=false)

Function to get specified data, either from a given constant value, a user-defined profile or from the global weather file:
* If no information is given: If required=false, nothing will be returned and an info message is logged. 
                              If required=true, an error message will be thrown                              
* If a constant value is set, this will be used as first argument
* If a profile_file_key is given, the data will be read from the user-defined profile.
* If 'from_global_file_key' is set to a valid entry of the global weather file, this will be used.

Inputs:
- config::Dict{String, Any}: input file config
- sim_params::Dict{String, Any}: Project-wide simulation parameters
- param_symbol::Union{Symbol, String}: will be used for info/warning messages (e.g., `:temperature`, `"temperature"`, etc.).
- profile_file_key::String: The configuration key in `config` that points to a `.prf` file path (e.g. `"temperature_profile_file_path"`).
- from_global_file_key::String: The configuration key in `config` that names a field in `weather_data` (e.g. `"temperature_from_global_file"`).
- constant_key::String: The configuration key in `config` for a constant parameter value (e.g. `"constant_temperature"`).
- uac::String: The calling unit uac as string for log/info/error messages
- required:Bool=false: flag if an error (true) or info (false) should be thrown in case no data is given

Returns:
- constant_value::Floathing: If a constant value is set, this will be returned as first return value
- profile_values::Union{Profile,Nothing}: If a profile is read out, it will be returned as second argument.

"""
function get_parameter_profile_from_config(config::Dict{String,Any},
                                           sim_params::Dict{String,Any},
                                           param_symbol::Union{Symbol,String},
                                           profile_file_key::String,
                                           from_global_file_key::String,
                                           constant_key::String,
                                           uac::String;
                                           required::Bool=false)
    # Count how many sources are specified 
    sources_specified = (haskey(config, profile_file_key) ? 1 : 0) +
                        (haskey(config, from_global_file_key) ? 1 : 0) +
                        (haskey(config, constant_key) ? 1 : 0)
    if sources_specified > 1
        @warn "For '$uac', the '$param_symbol' is specified with two or more sources in the input file!"
    end

    # 1. If a constant value is specified
    if haskey(config, constant_key) && config[constant_key] isa Number
        val = config[constant_key]
        @info "For '$uac', the '$param_symbol' is set to the constant value of $val."
        return val, nothing
    end

    # 2. If a `.prf` file path is specified
    if haskey(config, profile_file_key)
        path = config[profile_file_key]
        @info "For '$uac', the '$param_symbol' is taken from the user-defined .prf file located at: $path"
        return nothing, Profile(path, sim_params)
    end

    # 3. If the config says to read from the global weather data
    if haskey(config, from_global_file_key) && haskey(sim_params, "weather_data")
        field_name_str = config[from_global_file_key]
        # Check if it matches a field in weather_data
        wd = sim_params["weather_data"]
        field_symbols = fieldnames(typeof(wd))
        if any(occursin(field_name_str, string(sym)) for sym in field_symbols)
            @info "For '$uac', the '$param_symbol' is taken from the project-wide weather file: $field_name_str"
            return nothing, getfield(wd, Symbol(field_name_str))
        else
            @error "In '$uac', the '$param_symbol' given as '$field_name_str' must be one of: " *
                   "$(join(string.(field_symbols), ", "))."
            throw(InputError)
        end
    end

    # 4. Otherwise, nothing is set
    if required
        @error "For '$uac', no '$param_symbol' is set, but it is required."
        throw(InputError)
    else
        @info "For '$uac', no '$param_symbol' is set."
        return nothing, nothing
    end
end

end
