using ..Resie: get_run

"""
Implementation of a heat pump component, elevating heat to a higher temperature using
electricity.

This models includes several effects that complicate the calculation of efficiency (as COP),
some of which involve the temperatures at which the heat input is provided and heat output
is requested. While this is similar to the calculation of efficiencies in other transformers
such as a CHPP, the implementation does not numerically inverse the given input functions
and determine energies from those. Instead the efficiencies are calculated "forward", which
has problems in reaching the equilibrium point if kappa_opt != 1.0. This might be improved
in the future, but was deemed to difficult to implement at time of writing.

The heat pump can also be configured to skip most complicated calculations and use a
constant or Carnot-based COP and no PLRDE, etc. This too can be improved for heat pumps with
exactly one input and exactly one output, since the slicing algorithm is not required then.
"""
mutable struct HeatPump <: Component
    uac::String
    controller::Controller
    sys_function::SystemFunction

    input_interfaces::InterfaceMap
    output_interfaces::InterfaceMap

    m_el_in::Symbol
    m_heat_out::Symbol
    m_heat_in::Symbol

    design_power_th::Float64
    max_power_function::Function
    min_power_function::Function
    constant_cop::Floathing
    dynamic_cop::Function
    bypass_cop::Float64

    consider_icing::Bool
    icing_coefficients::Vector{Float64}
    optimise_slice_dispatch::Bool
    nr_optimisation_passes::UInt
    optimal_plr::Float64
    output_temperature::Temperature
    input_temperature::Temperature

    cop::Float64
    mix_temp_input::Float64
    mix_temp_output::Float64

    function HeatPump(uac::String, config::Dict{String,Any}, sim_params::Dict{String,Any})
        m_el_in = Symbol(default(config, "m_el_in", "m_e_ac_230v"))
        m_heat_out = Symbol(default(config, "m_heat_out", "m_h_w_ht1"))
        m_heat_in = Symbol(default(config, "m_heat_in", "m_h_w_lt1"))
        register_media([m_el_in, m_heat_out, m_heat_in])

        func_def = default(config, "cop_function", "carnot:0.4")
        constant_cop, cop_function = parse_cop_function(func_def)

        func_def = default(config, "max_power_function", "const:1.0")
        max_power_function = parse_2dim_function(func_def)

        func_def = default(config, "min_power_function", "const:0.2")
        min_power_function = parse_2dim_function(func_def)

        coeff_def = default(config, "icing_coefficients", "3,-0.42,15,2,30")
        icing_coefficients = parse.(Float64, split(coeff_def, ","))

        optimise_slice_dispatch = default(config, "optimise_slice_dispatch", false)
        optimal_plr = default(config, "optimal_plr", 1.0)
        if optimise_slice_dispatch && 1.0 - optimal_plr < sim_params["epsilon"]
            @warn "Heat pump $(uac) is configured to optimise slice dispatch but has an " *
                  "optimal PLR of 1.0. Consider toggling optimisation off as it has no effect."
        end
        if optimise_slice_dispatch && constant_cop !== nothing
            @error "Heat pump $(uac) is configured to optimise slice dispatch but has " *
                   "a constant COP. Toggle optimisation off as the algorithm is unstable " *
                   "in this case."
        end

        return new(uac,
                   Controller(default(config, "control_parameters", nothing)),
                   sf_transformer,
                   InterfaceMap(m_heat_in => nothing,
                                m_el_in => nothing),
                   InterfaceMap(m_heat_out => nothing),
                   m_el_in,
                   m_heat_out,
                   m_heat_in,
                   config["power_th"],
                   max_power_function,
                   min_power_function,
                   constant_cop,
                   cop_function,
                   default(config, "bypass_cop", 15.0),
                   default(config, "consider_icing", false),
                   icing_coefficients,
                   optimise_slice_dispatch,
                   UInt(default(config, "nr_optimisation_passes", 10)),
                   optimal_plr,
                   default(config, "output_temperature", nothing),
                   default(config, "input_temperature", nothing),
                   0.0, # cop
                   0.0, # mixing temperature in the input interface
                   0.0) # mixing temperature in the output interface
    end
end

"""
A struct containing lots of temporary values for internal calculations of the heat pump.
This is used to avoid excessively long argument and return lists and to improve readability.
Unless you're changing the code for heat pumps, this can be ignored entirely.
"""
mutable struct HPEnergies
    potential_energy_el::Floathing
    potentials_energies_heat_in::Vector{<:Floathing}
    potentials_energies_heat_out::Vector{<:Floathing}
    available_el_in::Float64
    available_heat_in::Vector{<:Floathing}
    available_heat_out::Vector{<:Floathing}
    max_usage_fraction::Float64
    in_indices::Vector{Integer}
    in_temps_min::Vector{<:Temperature}
    in_temps_max::Vector{<:Temperature}
    in_uacs::Vector{<:Stringing}
    out_indices::Vector{Integer}
    out_temps_min::Vector{<:Temperature}
    out_temps_max::Vector{<:Temperature}
    out_uacs::Vector{<:Stringing}
    heat_in_has_inf_energy::Bool
    heat_out_has_inf_energy::Bool
    slices_idx_to_plr::Dict{Integer,Tuple{Integer,Integer}}
    slices_el_in::Vector{Floathing}
    slices_heat_in::Vector{Floathing}
    slices_heat_in_temperature::Vector{Temperature}
    slices_heat_in_uac::Vector{Stringing}
    slices_heat_out::Vector{Floathing}
    slices_heat_out_temperature::Vector{Temperature}
    slices_heat_out_uac::Vector{Stringing}
    slices_el_in_temp::Vector{Floathing}
    slices_heat_in_temp::Vector{Floathing}
    slices_heat_in_temperature_temp::Vector{Temperature}
    slices_heat_in_uac_temp::Vector{Stringing}
    slices_heat_out_temp::Vector{Floathing}
    slices_heat_out_temperature_temp::Vector{Temperature}
    slices_heat_out_uac_temp::Vector{Stringing}

    function HPEnergies()
        return new(0.0,
                   Vector{Floathing}(),
                   Vector{Floathing}(),
                   0.0,
                   Vector{Floathing}(),
                   Vector{Floathing}(),
                   0.0,
                   Vector{Integer}(),
                   Vector{Temperature}(),
                   Vector{Temperature}(),
                   Vector{Stringing}(),
                   Vector{Integer}(),
                   Vector{Temperature}(),
                   Vector{Temperature}(),
                   Vector{Stringing}(),
                   false,
                   false,
                   Dict{Integer,Tuple{Integer,Integer}}(),
                   Vector{Floathing}(),
                   Vector{Floathing}(),
                   Vector{Temperature}(),
                   Vector{Stringing}(),
                   Vector{Floathing}(),
                   Vector{Temperature}(),
                   Vector{Stringing}(),
                   Vector{Floathing}(),
                   Vector{Floathing}(),
                   Vector{Temperature}(),
                   Vector{Stringing}(),
                   Vector{Floathing}(),
                   Vector{Temperature}(),
                   Vector{Stringing}())
    end
end

"""
Reset the temp slices of an energies container.
"""
function reset_temp_slices!(energies::HPEnergies)::HPEnergies
    energies.slices_el_in_temp = Vector{Floathing}()
    energies.slices_heat_in_temp = Vector{Floathing}()
    energies.slices_heat_in_temperature_temp = Vector{Temperature}()
    energies.slices_heat_in_uac_temp = Vector{Stringing}()
    energies.slices_heat_out_temp = Vector{Floathing}()
    energies.slices_heat_out_temperature_temp = Vector{Temperature}()
    energies.slices_heat_out_uac_temp = Vector{Stringing}()
    return energies
end

"""
Copies the temp slices of an energies container to the slices output fields and overwrites
existing values.
"""
function copy_temp_to_slices!(energies::HPEnergies)::HPEnergies
    energies.slices_el_in = energies.slices_el_in_temp
    energies.slices_heat_in = energies.slices_heat_in_temp
    energies.slices_heat_in_temperature = energies.slices_heat_in_temperature_temp
    energies.slices_heat_in_uac = energies.slices_heat_in_uac_temp
    energies.slices_heat_out = energies.slices_heat_out_temp
    energies.slices_heat_out_temperature = energies.slices_heat_out_temperature_temp
    energies.slices_heat_out_uac = energies.slices_heat_out_uac_temp
    return energies
end

"""
Resets the available energies of an energies container based on the potentials.
If an specific input or output layer index is given, that layer will be set as infinite and
the calculation will be done with only that layer on the input or output side.

# Arguments
- `energies::HPEnergies`: The energies container
- `as_inf::Integer`: If 1, sets the given input layer as infinite. If 2, the outputs.
    Defaults to 0, which does nothing.
- `idx::Integer`: The index of the layer to set as infinte. Defaults to 0, which does
    nothing.
# Returns
- `HPEnergies`: The energies container
"""
function reset_available!(energies::HPEnergies; as_inf::Integer=0, idx::Integer=0)::HPEnergies
    energies.available_el_in = copy(energies.potential_energy_el)
    energies.available_heat_in = copy(energies.potentials_energies_heat_in)
    energies.available_heat_out = copy(energies.potentials_energies_heat_out)
    if as_inf == 1
        energies.available_heat_in[idx] = Inf
    elseif as_inf == 2
        energies.available_heat_out[idx] = Inf
    end
    energies.in_indices = as_inf == 1 ? [idx] : [Int(i) for i in eachindex(energies.available_heat_in)]
    energies.out_indices = as_inf == 2 ? [idx] : [Int(i) for i in eachindex(energies.available_heat_out)]
    return energies
end

"""
Adds the heat_in temp slices to the existing slices output fields.
"""
function add_heat_in_temp_to_slices(energies::HPEnergies)::HPEnergies
    append!(energies.slices_heat_in, energies.slices_heat_in_temp)
    append!(energies.slices_heat_in_temperature, energies.slices_heat_in_temperature_temp)
    append!(energies.slices_heat_in_uac, energies.slices_heat_in_uac_temp)
    return energies
end

"""
Adds the heat_out temp slices to the existing slices output fields.
"""
function add_heat_out_temp_to_slices(energies::HPEnergies)::HPEnergies
    append!(energies.slices_heat_out, energies.slices_heat_out_temp)
    append!(energies.slices_heat_out_temperature, energies.slices_heat_out_temperature_temp)
    append!(energies.slices_heat_out_uac, energies.slices_heat_out_uac_temp)
    return energies
end

"""
Copies the electricity and heat_in temp slices to the slices output fields and overwrites
existing values based on whether the sums over the temp slices are smaller than the existing
values.
"""
function copy_heat_in_temp_to_slices(energies::HPEnergies)::HPEnergies
    if length(energies.slices_heat_in) == 0
        energies.slices_el_in = copy(energies.slices_el_in_temp)
        energies.slices_heat_in = copy(energies.slices_heat_in_temp)
        energies.slices_heat_in_temperature = copy(energies.slices_heat_in_temperature_temp)
        energies.slices_heat_in_uac = copy(energies.slices_heat_in_uac_temp)
    else
        if _isless(_sum(energies.slices_el_in), _sum(energies.slices_el_in_temp))
            energies.slices_el_in = energies.slices_el_in_temp
        end
        if _isless(_sum(energies.slices_heat_in), _sum(energies.slices_heat_in_temp))
            energies.slices_heat_in = energies.slices_heat_in_temp
            energies.slices_heat_in_temperature = energies.slices_heat_in_temperature_temp
            energies.slices_heat_in_uac = energies.slices_heat_in_uac_temp
        end
    end
    return energies
end

"""
Copies the electricity and heat_out temp slices to the slices output fields and overwrites
existing values based on whether the sums over the temp slices are smaller than the existing
values.
"""
function copy_heat_out_temp_to_slices(energies::HPEnergies)::HPEnergies
    if length(energies.slices_heat_out) == 0
        energies.slices_el_in = copy(energies.slices_el_in_temp)
        energies.slices_heat_out = copy(energies.slices_heat_out_temp)
        energies.slices_heat_out_temperature = copy(energies.slices_heat_out_temperature_temp)
        energies.slices_heat_out_uac = copy(energies.slices_heat_out_uac_temp)
    else
        if _isless(_sum(energies.slices_el_in), _sum(energies.slices_el_in_temp))
            energies.slices_el_in = energies.slices_el_in_temp
        end
        if _isless(_sum(energies.slices_heat_out), _sum(energies.slices_heat_out_temp))
            energies.slices_heat_out = energies.slices_heat_out_temp
            energies.slices_heat_out_temperature = energies.slices_heat_out_temperature_temp
            energies.slices_heat_out_uac = energies.slices_heat_out_uac_temp
        end
    end
    return energies
end

function initialise!(unit::HeatPump, sim_params::Dict{String,Any})
    set_storage_transfer!(unit.input_interfaces[unit.m_heat_in],
                          unload_storages(unit.controller, unit.m_heat_in))
    set_storage_transfer!(unit.input_interfaces[unit.m_el_in],
                          unload_storages(unit.controller, unit.m_el_in))
    set_storage_transfer!(unit.output_interfaces[unit.m_heat_out],
                          load_storages(unit.controller, unit.m_heat_out))
end

function control(unit::HeatPump, components::Grouping, sim_params::Dict{String,Any})
    update(unit.controller)
    # for fixed input/output temperatures, overwrite the interface with those. otherwise
    # highest will choose the interface's temperature (including nothing)
    if unit.output_temperature !== nothing
        set_temperature!(unit.output_interfaces[unit.m_heat_out],
                         unit.output_temperature,
                         unit.output_temperature)
    end
    if unit.input_temperature !== nothing
        set_temperature!(unit.input_interfaces[unit.m_heat_in],
                         unit.input_temperature,
                         unit.input_temperature)
    end
end

function set_max_energies!(unit::HeatPump,
                           el_in::Union{Floathing,Vector{<:Floathing}},
                           heat_in::Union{Floathing,Vector{<:Floathing}},
                           heat_out::Union{Floathing,Vector{<:Floathing}},
                           purpose_uac_heat_in::Union{Stringing,Vector{Stringing}}=nothing,
                           purpose_uac_heat_out::Union{Stringing,Vector{Stringing}}=nothing,
                           has_calculated_all_maxima_heat_in::Bool=false,
                           has_calculated_all_maxima_heat_out::Bool=false)
    set_max_energy!(unit.input_interfaces[unit.m_el_in], el_in)
    set_max_energy!(unit.input_interfaces[unit.m_heat_in], heat_in, purpose_uac_heat_in,
                    has_calculated_all_maxima_heat_in)
    set_max_energy!(unit.output_interfaces[unit.m_heat_out], heat_out, purpose_uac_heat_out,
                    has_calculated_all_maxima_heat_out)
end

"""
Determines the temperature of the input/output layer and if it should be skipped.

Layers are skipped if they have a temperature, but it does fall into the temperature band
defined by the heat pump input/output temperatures (if any).

# Arguments
- `unit::HeatPump`: The heat pump.
- `current_idx::Integer`: The index of the layer.
- `temps_min::Vector{<:Temperature}`: The minimum temperatures of all layers.
- `temps_max::Vector{<:Temperature}`: The maximum temperatures of all layers.
- `input::Bool`: If the layer is an input or not. Defaults to true.
# Returns
- `Bool`: If the layer should be skipped.
- `Temperature`: The layer temperature.
"""
function get_layer_temperature(unit::HeatPump,
                               current_idx::Integer,
                               temps_min::Vector{<:Temperature},
                               temps_max::Vector{<:Temperature};
                               input::Bool=true)::Tuple{Bool,Temperature}
    static_temperature = input ? unit.input_temperature : unit.output_temperature
    if static_temperature === nothing
        layer_temp = highest(temps_min[current_idx], temps_max[current_idx])
        if layer_temp === nothing && unit.constant_cop === nothing
            term = input ? "input" : "output"
            @error "Error: The $(term) temperature for $(unit.uac) could not be detected. " *
                   "Please specify one with the parameter '$(term)_temperature' or check " *
                   "the connected components."
            throw(InputError)
        end
        return false, layer_temp
    end

    # skip slice if given static temperature is not within the temperature band of the layer
    if (temps_min[current_idx] !== nothing && temps_min[current_idx] > static_temperature
        ||
        temps_max[current_idx] !== nothing && temps_max[current_idx] < static_temperature)
        # end of condition
        return true, nothing
    else
        return false, static_temperature
    end
end

"""
Calculates the new COP due to icing losses.

# Arguments
`unit::HeatPump`: The heat pump.
`cop::Floathing`: The COP without icing losses.
`in_temp::Temperature`: The input temperature of the current layer.
# Returns
- `Floathing`: The COP corrected due to icing losses. Returns nothing if the input COP is
    nothing.
"""
function icing_correction(unit::HeatPump, cop::Floathing, in_temp::Temperature)::Floathing
    if cop === nothing
        return nothing
    end

    lin_factor = unit.icing_coefficients[1] + unit.icing_coefficients[2] * in_temp
    exp_factor = unit.icing_coefficients[3] * exp(-1.0
                                                  * (in_temp - unit.icing_coefficients[4])^2
                                                  /
                                                  unit.icing_coefficients[5])
    return cop * (1 - 0.01 * (max(0.0, lin_factor) + exp_factor))
end

"""
Determines if the optimisation pass is accepted.

The optimisation might calculate a pass that would lead to implausible results, which is why
this check must be done. The pass is accepted if the "signature" of the given pass is the
same as that of the best pass and the timestep hasn't been more than used up. The signature
is a tuple with flags of which remaining energies of in- or outputs are non-zero. This is
checked because the optimisation should not change how the heat pump operates, it should
only improve the efficiency of the selected slices.

# Arguments
- `energies::HPEnergies`: The energies container.
- `times::Vector{Float64}`: The times of the slices.
- `sim_params::Dict{String,Any}`: Simulation parameters.
# Returns
- `Bool`: If the pass is accepted or not.
"""
function accept_pass(energies::HPEnergies, times::Vector{Float64}, sim_params::Dict{String,Any})::Bool
    if sum(times; init=0.0) > sim_params["time_step_seconds"]
        return false
    end

    EPS = sim_params["epsilon"]
    signature_best = (energies.potential_energy_el - sum(energies.slices_el_in; init=0.0) > EPS,
                      sum(energies.potentials_energies_heat_in; init=0.0)
                      -
                      sum(energies.slices_heat_in; init=0.0)
                      >
                      EPS,
                      sum(energies.potentials_energies_heat_out; init=0.0)
                      -
                      sum(energies.slices_heat_out; init=0.0)
                      >
                      EPS)
    signature_given = (energies.available_el_in > EPS,
                       sum(energies.available_heat_in; init=0.0) > EPS,
                       sum(energies.available_heat_out; init=0.0) > EPS)
    return signature_given == signature_best
end

"""
Determines if the optimisation pass is better than the currently best pass.

# Arguments
- `energies::HPEnergies`: The energies container.
- `sim_params::Dict{String,Any}`: Simulation parameters.
# Returns
- `Bool`: If the pass is better or not.
"""
function pass_is_better(energies::HPEnergies, sim_params::Dict{String,Any})::Bool
    return sum(energies.slices_el_in_temp; init=0.0) < sum(energies.slices_el_in; init=0.0)
end

"""
Updates the PLRs chosen for each slices as input for the next optimisation pass.

This is essentially how the optimisation is performed, by iteratively adjusting the PLRs for
slices according to a heuristic, which approaches a local optimum if iterated long enough.

The heuristic determines the most promising slice by its distance to the optimal PLR of the
heat pump multiplies with the amount of electricity that the slice requires, then halves its
distance to the optimal PLR.

Can be given a slice to ignore, which is used to avoid trying to update the same slice twice
in a row, which avoids some but not all problems with the algorithm getting "stuck".

# Arguments
`unit::HeatPump`: The heat pump.
`energies::HPEnergies`: The energies container.
`plrs::Array{<:Floathing,2}`: Current selection of PLRs of the slices.
`ignore_idx::Tuple{Integer,Integer}`: Index of the slice to ignore.
# Returns
- `Tuple{Integer,Integer}`: The index of the updated slice.
"""
function update_plrs!(unit::HeatPump,
                      energies::HPEnergies,
                      plrs::Array{<:Floathing,2},
                      ignore_idx::Tuple{Integer,Integer})::Tuple{Integer,Integer}
    best_idx = (0, 0)
    best_val = 0.0

    for (slice_idx, idxs) in pairs(energies.slices_idx_to_plr)
        idx_src, idx_snk = idxs
        if idx_src == ignore_idx[1] && idx_snk == ignore_idx[2]
            continue
        end

        val = (plrs[idx_src, idx_snk] - unit.optimal_plr) * energies.slices_el_in[slice_idx]
        if val > best_val
            best_idx = (idx_src, idx_snk)
            best_val = val
        end
    end

    if best_idx != (0, 0)
        plrs[best_idx[1], best_idx[2]] -= 0.5 * (plrs[best_idx[1], best_idx[2]] - unit.optimal_plr)
    end

    return best_idx
end

"""
    filter_by_transformer(energies, sim_param; heat_in=true)

Filters the heat input or heat output potentials by system function so only transformers are
listed.

# Args
- `energie::HPEnergies`: The energies container
- `sim_params::Dict{String,Any}`: Simulation parameters
- `heat_in::Bool`: (Optional) If true, filters the inputs. Otherwise the outputs. Defaults
    to true.
# Returns
- `Vector{Floathing}`: The filtered input or output potentials.
"""
function filter_by_transformer(energies::HPEnergies,
                               sim_params::Dict{String,Any};
                               heat_in::Bool=true)::Vector{Floathing}
    components = get_run(sim_params["run_ID"]).components
    filtered = []
    for (idx, uac) in enumerate(heat_in ? energies.in_uacs : energies.out_uacs)
        if uac !== nothing && components[uac].sys_function === sf_transformer
            if heat_in
                push!(filtered, energies.potentials_energies_heat_in[idx])
            else
                push!(filtered, energies.potentials_energies_heat_out[idx])
            end
        end
    end
    return filtered
end

"""
Calculate the energies for the given slice.

This is essentially how the heat pump would be calculated in the 1:1 case. The COP is
determined first, then the available energies on all three inputs/outputs are checked as
limiting factors. Bypass operation and icing losses are also considered. The function
specifically does not determine the PLR. For the 1:1 case this could be done here too, but
for the overarching slicing algorithm the PLR must be a required input.

# Arguments
- `unit::HeatPump`: The heat pump.
- `available_el_in::Float64`: Available electricity.
- `available_heat_in::Floathing`: Available heat input.
- `available_heat_out::Floathing`: Available heat output.
- `in_temp::Temperature`: Input temperature.
- `out_temp::Temperature`: Output temperature.
- `plr::Float64`: The PLR of the slice.
# Returns
- `Floathing`: Used input heat.
- `Floathing`: Used electricity.
- `Floathing`: Used output heat.
- `Temperature`: Temperature of input heat. Same as input argument.
- `Temperature`: Temperature of output heat. Same as output argument unless in bypass
    operation, when it's equal to the input temperature.
"""
function handle_slice(unit::HeatPump,
                      available_el_in::Float64,
                      available_heat_in::Floathing,
                      available_heat_out::Floathing,
                      in_temp::Temperature,
                      out_temp::Temperature,
                      plr::Float64)::Tuple{Floathing,Floathing,Floathing,Temperature,Temperature}
    # determine COP depending on three cases. a constant COP precludes the use of a bypass
    do_bypass = false
    if unit.constant_cop !== nothing
        cop = unit.constant_cop
    elseif in_temp >= out_temp
        cop = unit.bypass_cop
        do_bypass = true
    else
        cop = unit.dynamic_cop(in_temp, out_temp)(plr)
        if unit.consider_icing
            cop = icing_correction(unit, cop, in_temp)
        end
    end
    if cop === nothing
        @error ("Input and/or output temperature for heatpump $(unit.uac) is not " *
                "given. Provide temperatures or fixed cop.")
        throw(InputError)
    end
    if cop < 1.0
        cop = 1.0
        @warn ("Calculated COP of heat pump $(unit.uac) was below 1.0. Please check " *
               "input for mistakes as this should not happen. COP was set to 1.0")
    end

    # calculate energies with the current cop
    # energies for current slice with potential heat in as basis
    used_heat_in = available_heat_in
    used_el_in = used_heat_in / (cop - 1.0)
    used_heat_out = used_heat_in + used_el_in

    # check heat out as limiter
    if used_heat_out > available_heat_out
        used_heat_out = available_heat_out
        used_el_in = used_heat_out / cop
        used_heat_in = used_el_in * (cop - 1.0)
    end

    # check electricity in as limiter
    if used_el_in > available_el_in
        used_el_in = available_el_in
        used_heat_in = used_el_in * (cop - 1.0)
        used_heat_out = used_heat_in + used_el_in
    end

    return (used_heat_in,
            used_el_in,
            used_heat_out,
            in_temp,
            do_bypass ? in_temp : out_temp)
end

"""
Performs the slicing algorithm given a set list of PLRs for each possible slice (each
combination input and output layer). If the PLR for a slice is nothing, sets it to 1.0 and
calculates the slice based on that.

Note the difference in indexing between the PLRs, where all possible slices are listed via
their input layer and output layer indices, and the return values, where the selected slices
are listed linearly and doesn't necessarily list all slices.

# Arguments
- `unit::HeatPump`: The heat pump.
- `sim_params::Dict{String,Any}`: Simulation parameters.
- `energies::HPEnergies`: The energies container.
- `plrs::Array{<:Floathing,2}`: The PLRs for each slice (may be nothing). Will be modified.
# Returns
- `HPEnergies`: The energies container.
- `Vector{Float64}`: The time each selected slice takes up assuming the heat pump runs with
    minimal power for that slice.
- `Vector{Float64}`: The time each selected slice takes up assuming the heat pump runs with
    the selected PLR for that slice.
- `Array{<:Floathing,2}`: The PLRs for each slice.
"""
function calculate_slices(unit::HeatPump,
                          sim_params::Dict{String,Any},
                          energies::HPEnergies,
                          plrs::Array{<:Floathing,2})::Tuple{HPEnergies,Vector{Float64},Vector{Float64},
                                                             Array{<:Floathing,2}}
    # times are the timespans each slice takes. times_min means the time slice assuming
    # minimum power. times_min are actually larger than times
    times_min = Vector{Float64}()
    times = Vector{Float64}()

    sum_usage = 0.0
    slice_idx::Int = 1
    current_in_idx::Int = 0
    current_out_idx::Int = 0
    EPS = sim_params["epsilon"]

    while (true)
        try
            current_in_idx = first(energies.in_indices)
            current_out_idx = first(energies.out_indices)
        catch BoundsError
            break
        end

        # we continue for as long as there is energy to be distributed and max power has
        # not been reached. during slice dispatch optimisation it might also be the case
        # that extra slices appear, with the same input and output as the previous slice,
        # because the target PLR is too low for full operation. we need to catch that too
        if energies.available_el_in < EPS ||
           sum(energies.available_heat_in; init=0.0) < EPS ||
           sum(energies.available_heat_out; init=0.0) < EPS ||
           energies.max_usage_fraction - sum_usage < EPS ||
           slice_idx > length(energies.in_indices) * length(energies.out_indices)
            # end of condition
            break
        end

        # skip layers with zero energy remaining
        if energies.available_heat_in[current_in_idx] < EPS
            current_in_idx = popfirst!(energies.in_indices)
            continue
        end
        if energies.available_heat_out[current_out_idx] < EPS
            current_out_idx = popfirst!(energies.out_indices)
            continue
        end

        # check and determine input temperature of layer
        skip_slice, current_in_temp = get_layer_temperature(unit, current_in_idx,
                                                            energies.in_temps_min,
                                                            energies.in_temps_max;
                                                            input=true)
        if skip_slice
            energies.available_heat_in[current_in_idx] = 0.0
            current_in_idx = popfirst!(energies.in_indices)
            continue
        end

        # check and determine output temperature of layer
        skip_slice, current_out_temp = get_layer_temperature(unit, current_out_idx,
                                                             energies.out_temps_min,
                                                             energies.out_temps_max;
                                                             input=false)
        if skip_slice
            energies.available_heat_out[current_out_idx] = 0.0
            current_out_idx = popfirst!(energies.out_indices)
            continue
        end

        if plrs[current_in_idx, current_out_idx] === nothing
            plrs[current_in_idx, current_out_idx] = 1.0
        end

        min_power_frac = min(1.0, unit.min_power_function(current_in_temp, current_out_temp))
        min_power = max(EPS, unit.design_power_th * min_power_frac)
        used_power_frac = min(plrs[current_in_idx, current_out_idx],
                              unit.max_power_function(current_in_temp, current_out_temp))
        used_power = max(min_power, unit.design_power_th * used_power_frac)

        remaining_heat_out = min(energies.available_heat_out[current_out_idx],
                                 (energies.max_usage_fraction - sum_usage) * sim_params["watt_to_wh"](used_power))

        used_heat_in,
        used_el_in,
        used_heat_out,
        current_in_temp,
        current_out_temp = handle_slice(unit,
                                        energies.available_el_in,
                                        energies.available_heat_in[current_in_idx],
                                        remaining_heat_out,
                                        current_in_temp,
                                        current_out_temp,
                                        plrs[current_in_idx, current_out_idx])

        # finally all checks done, we add the slice and update remaining energies
        push!(energies.slices_el_in_temp, used_el_in)
        push!(energies.slices_heat_in_temp, used_heat_in)
        push!(energies.slices_heat_in_temperature_temp, current_in_temp)
        push!(energies.slices_heat_in_uac_temp, energies.in_uacs[current_in_idx])
        push!(energies.slices_heat_out_temp, used_heat_out)
        push!(energies.slices_heat_out_temperature_temp, current_out_temp)
        push!(energies.slices_heat_out_uac_temp, energies.out_uacs[current_out_idx])

        energies.available_el_in -= used_el_in
        energies.available_heat_in[current_in_idx] -= used_heat_in
        energies.available_heat_out[current_out_idx] -= used_heat_out

        sum_usage += used_heat_out / sim_params["watt_to_wh"](unit.design_power_th)
        push!(times_min, used_heat_out * 3600 / min_power)
        push!(times, used_heat_out * 3600 / used_power)

        energies.slices_idx_to_plr[slice_idx] = (current_in_idx, current_out_idx)
        slice_idx += 1
    end

    return energies, times_min, times, plrs
end

"""
The inner part of the calculate_energies function.

The steps are:
1. Initialises the PLRs for all possible slices as nothing (they will be set as 1.0 upon
    selection)
2. Performs the slicing algorithm with the assumption of full power
3. Checks if minimum power restrictions were observed.
4. Optionally performs optimisation of the PLRs

# Arguments
- `unit::HeatPump`: The heat pump.
- `sim_params::Dict{String,Any}`: Simulation parameters.
- `energies::HPEnergies`: The energies container.
- `optimise_slice_dispatch::Bool`: Defaults to false.
# Returns
- `HPEnergies`: The energies container.
"""
function calculate_energies_heatpump(unit::HeatPump,
                                     sim_params::Dict{String,Any},
                                     energies::HPEnergies,
                                     optimise_slice_dispatch::Bool=false)::HPEnergies
    energies = reset_temp_slices!(energies)

    # initial PLRs are nothing, as we do not yet know which slices are used and which aren't
    plrs = Array{Floathing,2}(nothing,
                              length(energies.available_heat_in),
                              length(energies.available_heat_out))

    # first pass with full power (or only pass if no optimisation is used)
    energies, times_min, times, plrs = calculate_slices(unit, sim_params, energies, plrs)

    # as long as the sum of times with minimum power per slice is larger than the time step
    # and the sum of times with chosen power per slice is smaller or equal to the time step
    # (which is implicit in the calculation), we can find a dispatch of each slice by
    # "slowing them down" from the minimum power up to the chosen power, so the sum works
    # out. we don't have to calculate this dispatch, it is theoretical. if the minimum power
    # cannot be observed, the heat pump should not run at all
    if sum(times_min; init=0.0) < sim_params["time_step_seconds"]
        energies = reset_temp_slices!(energies)
        return energies
    end

    # warning: this optimisation is janky and doesn't work as well as it could
    if optimise_slice_dispatch
        energies = copy_temp_to_slices!(energies)
        last_updated = (0, 0)

        for _ in 1:(unit.nr_optimisation_passes)
            energies = reset_temp_slices!(energies)
            energies = reset_available!(energies)
            last_updated = update_plrs!(unit, energies, plrs, last_updated)
            energies, times_min, times, plrs = calculate_slices(unit, sim_params, energies, plrs)

            if accept_pass(energies, times, sim_params) && pass_is_better(energies, sim_params)
                energies = copy_temp_to_slices!(energies)
            end
        end
    end

    return energies
end

"""
Calculates the energies the heat pump can process in the timestep.

This is split in an inner function (see calculate_energies_heatpump) and this outer
function, because the calculation happens differently depending on whether there are other
transformers in the input exchanges, outputs, both or neither.

# Arguments
- `unit::HeatPump`: The heat pump.
- `sim_params::Dict{String,Any}`: Simulation parameters.
# Returns
- `Bool`: If false the heat pump cannot process any energy for one of several reasons. The
    values might still be zero even if this flag is true, but a value of false definitely
    means no operation.
- `HPEnergies`: The energies container.
"""
function calculate_energies(unit::HeatPump, sim_params::Dict{String,Any})
    energies = HPEnergies()

    # get usage fraction from control modules
    energies.max_usage_fraction = upper_plr_limit(unit.controller, sim_params)
    if energies.max_usage_fraction <= 0.0
        return false, (nothing)
    end

    # get potentials from inputs/outputs. The heat input and output are calculated as 
    # vectors to allow for temperature layers, while the electricity input is a scalar
    energies.potential_energy_el = check_el_in(unit, sim_params)
    if energies.potential_energy_el === nothing
        # shortcut if we're limited by zero input electricity
        return false, (nothing)
    end

    energies.potentials_energies_heat_in,
    energies.in_temps_min,
    energies.in_temps_max,
    energies.in_uacs = check_heat_in_layered(unit, sim_params)

    energies.potentials_energies_heat_out,
    energies.out_temps_min,
    energies.out_temps_max,
    energies.out_uacs = check_heat_out_layered(unit, sim_params)

    # in the following we want to work with positive values as it is easier
    energies.potentials_energies_heat_in = abs.(energies.potentials_energies_heat_in)
    energies.potentials_energies_heat_out = abs.(energies.potentials_energies_heat_out)

    # reorder inputs and outputs according to control modules
    index = reorder_inputs(unit.controller, energies.in_temps_min, energies.in_temps_max)
    energies.in_temps_min = energies.in_temps_min[index]
    energies.in_temps_max = energies.in_temps_max[index]
    energies.in_uacs = energies.in_uacs[index]
    energies.potentials_energies_heat_in = energies.potentials_energies_heat_in[index]

    index = reorder_outputs(unit.controller, energies.out_temps_min, energies.out_temps_max)
    energies.out_temps_min = energies.out_temps_min[index]
    energies.out_temps_max = energies.out_temps_max[index]
    energies.out_uacs = energies.out_uacs[index]
    energies.potentials_energies_heat_out = energies.potentials_energies_heat_out[index]

    # there are three different cases of how to handle the layered approach of operating the
    # heat pump, depending on wether or not any input heat or output heat transformer has a
    # value of infinite as the potential. if this is the case for both the input and output,
    # the calculation cannot be resolved.
    energies.heat_in_has_inf_energy = any(isinf, filter_by_transformer(energies, sim_params; heat_in=true))
    energies.heat_out_has_inf_energy = any(isinf, filter_by_transformer(energies, sim_params; heat_in=false))

    if energies.heat_in_has_inf_energy && energies.heat_out_has_inf_energy
        @warn "The heat pump $(unit.uac) has unknown energies in both its inputs and " *
              "outputs. This cannot be resolved. Please check the order of operation and " *
              "make sure that either the inputs or the outputs have been fully calculated " *
              "before the heat pump $(unit.uac) has its potential step."
        return false, energies
    end

    if energies.heat_in_has_inf_energy
        for heat_in_idx in eachindex(energies.potentials_energies_heat_in)
            energies = reset_available!(energies; as_inf=1, idx=heat_in_idx)
            energies = calculate_energies_heatpump(unit, sim_params, energies)
            energies = add_heat_in_temp_to_slices(energies)
            energies = copy_heat_out_temp_to_slices(energies)
        end
    elseif energies.heat_out_has_inf_energy
        for heat_out_idx in eachindex(energies.potentials_energies_heat_out)
            energies = reset_available!(energies; as_inf=2, idx=heat_out_idx)
            energies = calculate_energies_heatpump(unit, sim_params, energies)
            energies = add_heat_out_temp_to_slices(energies)
            energies = copy_heat_in_temp_to_slices(energies)
        end
    else
        energies = reset_available!(energies)
        energies = calculate_energies_heatpump(unit, sim_params, energies, unit.optimise_slice_dispatch)
        energies = copy_temp_to_slices!(energies)
    end

    return true, energies
end

function potential(unit::HeatPump, sim_params::Dict{String,Any})
    success, energies = calculate_energies(unit, sim_params)

    if !success
        set_max_energies!(unit, 0.0, 0.0, 0.0)
    else
        set_max_energies!(unit,
                          energies.slices_el_in,
                          energies.slices_heat_in,
                          energies.slices_heat_out,
                          energies.slices_heat_in_uac,
                          energies.slices_heat_out_uac,
                          energies.heat_in_has_inf_energy,
                          energies.heat_out_has_inf_energy)
        set_temperature!(unit.input_interfaces[unit.m_heat_in],
                         lowest(energies.slices_heat_in_temperature),
                         nothing)
        set_temperature!(unit.output_interfaces[unit.m_heat_out],
                         nothing,
                         highest(energies.slices_heat_out_temperature))
    end
end

function process(unit::HeatPump, sim_params::Dict{String,Any})
    success, energies = calculate_energies(unit, sim_params)

    if !success
        set_max_energies!(unit, 0.0, 0.0, 0.0)
        return
    end

    el_in = sum(energies.slices_el_in; init=0.0)
    heat_out = sum(energies.slices_heat_out; init=0.0)

    if heat_out < sim_params["epsilon"]
        set_max_energies!(unit, 0.0, 0.0, 0.0)
        return
    end

    # calculate average cop of current time step
    if el_in > sim_params["epsilon"]
        unit.cop = heat_out / el_in
    end
    unit.mix_temp_input = _weighted_mean(energies.slices_heat_in_temperature, energies.slices_heat_in)
    unit.mix_temp_output = _weighted_mean(energies.slices_heat_out_temperature, energies.slices_heat_out)

    sub!(unit.input_interfaces[unit.m_el_in], el_in)
    sub!(unit.input_interfaces[unit.m_heat_in],
         energies.slices_heat_in,
         lowest(energies.slices_heat_in_temperature),
         energies.slices_heat_in_uac)
    add!(unit.output_interfaces[unit.m_heat_out],
         energies.slices_heat_out,
         highest(energies.slices_heat_out_temperature),
         energies.slices_heat_out_uac)
end

# has its own reset function as here more parameters are present that need to be reset in
# every timestep
function reset(unit::HeatPump)
    invoke(reset, Tuple{Component}, unit)

    # reset other parameter
    unit.cop = 0.0
    unit.mix_temp_input = 0.0
    unit.mix_temp_output = 0.0
end

function output_values(unit::HeatPump)::Vector{String}
    return [string(unit.m_el_in) * " IN",
            string(unit.m_heat_in) * " IN",
            string(unit.m_heat_out) * " OUT",
            "COP",
            "MixingTemperature_Input",
            "MixingTemperature_Output"]
end

function output_value(unit::HeatPump, key::OutputKey)::Float64
    if key.value_key == "IN"
        return calculate_energy_flow(unit.input_interfaces[key.medium])
    elseif key.value_key == "OUT"
        return calculate_energy_flow(unit.output_interfaces[key.medium])
    elseif key.value_key == "COP"
        return unit.cop
    elseif key.value_key == "MixingTemperature_Input"
        return unit.mix_temp_input
    elseif key.value_key == "MixingTemperature_Output"
        return unit.mix_temp_output
    end
    throw(KeyError(key.value_key))
end

export HeatPump
