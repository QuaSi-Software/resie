using ..Resie: get_run

"""
Thermal Booster

Takes thermal energy of low temperature to preheats a demand and adds additional energy to
reach the desired output temperature of the demand.

"""
mutable struct ThermalBooster <: Component
    uac::String
    controller::Controller
    sys_function::SystemFunction

    input_interfaces::InterfaceMap
    output_interfaces::InterfaceMap

    m_el_in::Symbol
    m_heat_out::Symbol
    m_heat_in::Symbol

    output_temperature::Temperature
    input_temperature::Temperature

    power_losses_factor::Float64
    heat_losses_factor::Float64

    power_el::Float64
    cp_medium_out::Float64
    terminal_dT::Float64
    demand_input_temperature_profile::Union{Profile,Nothing}
    demand_input_temperature::Temperature
    allow_boost_solely::Bool
    allow_boost_additional::Bool
    mean_intermediate_temperature::Temperature
    intermediate_temperatures::Vector{Temperature}
    intermediate_temperature_energies::Vector{Float64}

    losses::Float64
    losses_power::Float64
    losses_heat::Float64

    function ThermalBooster(uac::String, config::Dict{String,Any}, sim_params::Dict{String,Any})
        m_el_in = Symbol(default(config, "m_el_in", "m_e_ac_230v"))
        m_heat_out = Symbol(default(config, "m_heat_out", "m_h_w_ht1"))
        m_heat_in = Symbol(default(config, "m_heat_in", "m_h_w_lt1"))
        register_media([m_el_in, m_heat_out, m_heat_in])

        constant_demand_input_temperature,
        demand_input_temperature_profile = get_parameter_profile_from_config(config,
                                                                             sim_params,
                                                                             "demand_input_temperature",
                                                                             "demand_input_temperature_profile_file_path",
                                                                             "",
                                                                             "constant_demand_input_temperature",
                                                                             uac;
                                                                             required=true)

        return new(uac,
                   Controller(default(config, "control_parameters", nothing)),
                   sf_transformer,
                   InterfaceMap(m_heat_in => nothing,
                                m_el_in => nothing),
                   InterfaceMap(m_heat_out => nothing),
                   m_el_in,
                   m_heat_out,
                   m_heat_in,
                   default(config, "output_temperature", nothing),    # [°C]
                   default(config, "input_temperature", nothing),     # [°C]
                   default(config, "power_losses_factor", 0.97),      # [-]
                   default(config, "heat_losses_factor", 0.95),       # [-]
                   config["power_el"],                                # [W]
                   default(config, "cp_medium_out", 4180.0),          # [J/(kgK)]
                   default(config, "terminal_dT", 2.0),               # [K]
                   demand_input_temperature_profile,                  # [°C] demand_input_temperature_profile
                   constant_demand_input_temperature,                 # [°C] demand_input_temperature
                   default(config, "allow_boost_solely", true),       # allow_boost_solely: allow use of energy if no heat_in is available at all
                   default(config, "allow_boost_additional", true),   # allow_boost_additional: allow use of energy to increase mass flow if heat_in is not sufficient --> allow intermediate temperature to be lower than heat_in temperature
                   0.0,                                               # mean_intermediate_temperature
                   Temperature[],                                     # intermediate_temperatures
                   Float64[],                                         # intermediate_temperature_energies
                   0.0,                                               # losses
                   0.0,                                               # losses_power
                   0.0)                                               # losses_heat
    end
end

mutable struct TBEnergies
    potential_el_in::Floathing
    potential_el_in_layered::Vector{<:Floathing}
    potentials_heat_in::Vector{<:Floathing}
    potentials_heat_out::Vector{<:Floathing}
    available_el_in::Float64
    in_uacs_el::Vector{<:Stringing}
    available_heat_in::Vector{<:Floathing}
    available_heat_out::Vector{<:Floathing}
    in_temps_min::Vector{<:Temperature}
    in_temps_max::Vector{<:Temperature}
    in_uacs_heat::Vector{<:Stringing}
    out_temps_min::Vector{<:Temperature}
    out_temps_max::Vector{<:Temperature}
    out_uacs::Vector{<:Stringing}
    heat_in_has_inf_energy::Bool
    heat_out_has_inf_energy::Bool
    slices_temp_el_in::Vector{Floathing}
    slices_temp_heat_in::Vector{Floathing}
    slices_temp_heat_in_temperature::Vector{Temperature}
    slices_temp_heat_in_uac::Vector{Stringing}
    slices_temp_heat_out::Vector{Floathing}
    slices_temp_heat_out_temperature::Vector{Temperature}
    slices_temp_heat_out_uac::Vector{Stringing}

    slices_el_in::Vector{Floathing}
    slices_heat_in::Vector{Floathing}
    slices_heat_in_temperature::Vector{Temperature}
    slices_heat_in_uac::Vector{Stringing}
    slices_heat_out::Vector{Floathing}
    slices_heat_out_temperature::Vector{Temperature}
    slices_heat_out_uac::Vector{Stringing}
    in_indices::Vector{Integer}
    out_indices::Vector{Integer}

    function TBEnergies()
        return new(0.0,
                   Vector{Floathing}(),
                   Vector{Floathing}(),
                   Vector{Floathing}(),
                   0.0,
                   Vector{Stringing}(),
                   Vector{Floathing}(),
                   Vector{Floathing}(),
                   Vector{Temperature}(),
                   Vector{Temperature}(),
                   Vector{Stringing}(),
                   Vector{Temperature}(),
                   Vector{Temperature}(),
                   Vector{Stringing}(),
                   false,
                   false,
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
                   Vector{Stringing}(),
                   Vector{Integer}(),
                   Vector{Integer}())
    end
end

function initialise!(unit::ThermalBooster, sim_params::Dict{String,Any})
    set_storage_transfer!(unit.input_interfaces[unit.m_heat_in],
                          unload_storages(unit.controller, unit.m_heat_in))
    set_storage_transfer!(unit.input_interfaces[unit.m_el_in],
                          unload_storages(unit.controller, unit.m_el_in))
    set_storage_transfer!(unit.output_interfaces[unit.m_heat_out],
                          load_storages(unit.controller, unit.m_heat_out))
end

function control(unit::ThermalBooster, components::Grouping, sim_params::Dict{String,Any})
    update(unit.controller)

    # update current demand input temperature
    if unit.demand_input_temperature_profile !== nothing
        unit.demand_input_temperature = Profiles.value_at_time(unit.demand_input_temperature_profile, sim_params)
    end

    if unit.output_temperature !== nothing
        set_max_energy!(unit.output_interfaces[unit.m_heat_out],
                        nothing,
                        nothing,
                        unit.output_temperature)
    end
    if unit.input_temperature !== nothing
        set_max_energy!(unit.input_interfaces[unit.m_heat_in],
                        nothing,
                        unit.input_temperature,
                        nothing)
    end
end

function set_max_energies!(unit::ThermalBooster,
                           is_transformer_potential::Bool,
                           el_in::Union{Floathing,Vector{<:Floathing}},
                           heat_in::Union{Floathing,Vector{<:Floathing}},
                           heat_out::Union{Floathing,Vector{<:Floathing}},
                           slices_heat_in_temperature::Union{Temperature,Vector{<:Temperature}}=nothing,
                           slices_heat_out_temperature::Union{Temperature,Vector{<:Temperature}}=nothing,
                           purpose_uac_heat_in::Union{Stringing,Vector{Stringing}}=nothing,
                           purpose_uac_heat_out::Union{Stringing,Vector{Stringing}}=nothing,
                           has_calculated_all_maxima_heat_in::Bool=false,
                           has_calculated_all_maxima_heat_out::Bool=false)
    set_max_energy!(unit.input_interfaces[unit.m_el_in], el_in; is_transformer_potential=is_transformer_potential)
    set_max_energy!(unit.input_interfaces[unit.m_heat_in], heat_in, slices_heat_in_temperature, nothing,
                    purpose_uac_heat_in, has_calculated_all_maxima_heat_in;
                    is_transformer_potential=is_transformer_potential)
    set_max_energy!(unit.output_interfaces[unit.m_heat_out], heat_out, nothing, slices_heat_out_temperature,
                    purpose_uac_heat_out, has_calculated_all_maxima_heat_out;
                    is_transformer_potential=is_transformer_potential)
end

"""
    filter_by_transformer(energies, sim_param; heat_in=true)

Filters the heat input or heat output potentials by system function so only transformers are
listed.

# Args
- `energie::TBEnergies`: The energies container
- `sim_params::Dict{String,Any}`: Simulation parameters
- `heat_in::Bool`: (Optional) If true, filters the inputs. Otherwise the outputs. Defaults
    to true.
# Returns
- `Vector{Floathing}`: The filtered input or output potentials.
"""
function filter_by_transformer(energies::TBEnergies,
                               sim_params::Dict{String,Any};
                               heat_in::Bool=true)::Vector{Floathing}
    components = get_run(sim_params["run_ID"]).components
    filtered = []
    for (idx, uac) in enumerate(heat_in ? energies.in_uacs_heat : energies.out_uacs)
        if uac !== nothing && components[uac].sys_function === sf_transformer
            if heat_in
                push!(filtered, energies.potentials_heat_in[idx])
            else
                push!(filtered, energies.potentials_heat_out[idx])
            end
        end
    end
    return filtered
end

"""
Reset the temporary slices of an energies container.
"""
function reset_temp_slices!(energies::TBEnergies)::TBEnergies
    energies.slices_temp_el_in = Vector{Floathing}()
    energies.slices_temp_heat_in = Vector{Floathing}()
    energies.slices_temp_heat_in_temperature = Vector{Temperature}()
    energies.slices_temp_heat_in_uac = Vector{Stringing}()
    energies.slices_temp_heat_out = Vector{Floathing}()
    energies.slices_temp_heat_out_temperature = Vector{Temperature}()
    energies.slices_temp_heat_out_uac = Vector{Stringing}()
    return energies
end

"""
Adds the heat_in temp slices to the existing slices output fields.
"""
function add_heat_in_temp_to_slices(energies::TBEnergies)::TBEnergies
    append!(energies.slices_heat_in, energies.slices_temp_heat_in)
    append!(energies.slices_heat_in_temperature, energies.slices_temp_heat_in_temperature)
    append!(energies.slices_heat_in_uac, energies.slices_temp_heat_in_uac)
    return energies
end

"""
Adds the heat_out temp slices to the existing slices output fields.
"""
function add_heat_out_temp_to_slices(energies::TBEnergies)::TBEnergies
    append!(energies.slices_heat_out, energies.slices_temp_heat_out)
    append!(energies.slices_heat_out_temperature, energies.slices_temp_heat_out_temperature)
    append!(energies.slices_heat_out_uac, energies.slices_temp_heat_out_uac)
    return energies
end

"""
Copies the electricity and heat_in temp slices to the slices output fields and overwrites
existing values based on whether the sums over the temp slices are smaller than the existing
values.
"""
function copy_heat_in_temp_to_slices(energies::TBEnergies)::TBEnergies
    if length(energies.slices_heat_in) == 0
        energies.slices_el_in = copy(energies.slices_temp_el_in)
        energies.slices_heat_in = copy(energies.slices_temp_heat_in)
        energies.slices_heat_in_temperature = copy(energies.slices_temp_heat_in_temperature)
        energies.slices_heat_in_uac = copy(energies.slices_temp_heat_in_uac)
    else
        if _isless(_sum(energies.slices_el_in), _sum(energies.slices_temp_el_in))
            energies.slices_el_in = energies.slices_temp_el_in
        end
        if _isless(_sum(energies.slices_heat_in), _sum(energies.slices_temp_heat_in))
            energies.slices_heat_in = energies.slices_temp_heat_in
            energies.slices_heat_in_temperature = energies.slices_temp_heat_in_temperature
            energies.slices_heat_in_uac = energies.slices_temp_heat_in_uac
        end
    end
    return energies
end

"""
Copies the electricity and heat_out temp slices to the slices output fields and overwrites
existing values based on whether the sums over the temp slices are smaller than the existing
values.
"""
function copy_heat_out_temp_to_slices(energies::TBEnergies)::TBEnergies
    if length(energies.slices_heat_out) == 0
        energies.slices_el_in = copy(energies.slices_temp_el_in)
        energies.slices_heat_out = copy(energies.slices_temp_heat_out)
        energies.slices_heat_out_temperature = copy(energies.slices_temp_heat_out_temperature)
        energies.slices_heat_out_uac = copy(energies.slices_temp_heat_out_uac)
    else
        if _isless(_sum(energies.slices_el_in), _sum(energies.slices_temp_el_in))
            energies.slices_el_in = energies.slices_temp_el_in
        end
        if _isless(_sum(energies.slices_heat_out), _sum(energies.slices_temp_heat_out))
            energies.slices_heat_out = energies.slices_temp_heat_out
            energies.slices_heat_out_temperature = energies.slices_temp_heat_out_temperature
            energies.slices_heat_out_uac = energies.slices_temp_heat_out_uac
        end
    end
    return energies
end

"""
Copies the temporary slices of the given energy container to the actual slices.

This overwrites any existing actual slices.
"""
function copy_temp_to_slices!(energies::TBEnergies)::TBEnergies
    energies.slices_el_in = copy(energies.slices_temp_el_in)
    energies.slices_heat_in = copy(energies.slices_temp_heat_in)
    energies.slices_heat_in_temperature = copy(energies.slices_temp_heat_in_temperature)
    energies.slices_heat_in_uac = copy(energies.slices_temp_heat_in_uac)
    energies.slices_heat_out = copy(energies.slices_temp_heat_out)
    energies.slices_heat_out_temperature = copy(energies.slices_temp_heat_out_temperature)
    energies.slices_heat_out_uac = copy(energies.slices_temp_heat_out_uac)
    return energies
end

"""
Determines the temperature of the input/output layer and if it should be skipped.

Layers are skipped if they have a temperature that are out of the temperature band
defined by the thermal booster input/output temperatures (if any).

# Arguments
- `unit::ThermalBooster`: The thermal booster.
- `current_idx::Integer`: The index of the layer.
- `temps_min::Vector{<:Temperature}`: The minimum temperatures of all layers.
- `temps_max::Vector{<:Temperature}`: The maximum temperatures of all layers.
- `input::Bool`: If the layer is an input or not. Defaults to true.
# Returns
- `Bool`: If the layer should be skipped.
- `Temperature`: The layer temperature.
"""
function get_layer_temperature(unit::ThermalBooster,
                               current_idx::Integer,
                               temps_min::Vector{<:Temperature},
                               temps_max::Vector{<:Temperature};
                               input::Bool=true)::Tuple{Bool,Temperature}
    static_temperature = input ? unit.input_temperature : unit.output_temperature
    if static_temperature === nothing
        layer_temp = highest(temps_min[current_idx], temps_max[current_idx])
        if layer_temp === nothing
            term = input ? "input" : "output"
            @error "Error: The $(term) temperature for $(unit.uac) could not be detected. " *
                   "Please specify one with the parameter '$(term)_temperature' or check " *
                   "the connected components."
            throw(InputError())
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
Resets the available energies of an energies container based on the potentials.

If a specific input or output layer index is given, that layer will be set as infinite and
the calculation will be done with only that layer on the input or output side.

# Arguments
- `energies::TBEnergies`: The energies container
- `fixed_heat_in::Union{Nothing,Integer}`: If a layer index is given, sets that input layer
    to infinite and selects it exclusively for sources.
- `fixed_heat_out::Union{Nothing,Integer}`: If a layer index is given, sets that output
    layer to infinite and selects it exclusively for sinks.
# Returns
- `TBEnergies`: The energies container
"""
function reset_available!(energies::TBEnergies;
                          fixed_heat_in::Union{Nothing,Integer}=nothing,
                          fixed_heat_out::Union{Nothing,Integer}=nothing)::TBEnergies
    energies.available_el_in = copy(energies.potential_el_in)
    energies.available_heat_in = copy(energies.potentials_heat_in)
    energies.available_heat_out = copy(energies.potentials_heat_out)
    if fixed_heat_in !== nothing
        energies.available_heat_in[fixed_heat_in] = Inf
    elseif fixed_heat_out !== nothing
        energies.available_heat_out[fixed_heat_out] = Inf
    end
    energies.in_indices = fixed_heat_in !== nothing ? [fixed_heat_in] :
                          [Int(i) for i in eachindex(energies.available_heat_in)]
    energies.out_indices = fixed_heat_out !== nothing ? [fixed_heat_out] :
                           [Int(i) for i in eachindex(energies.available_heat_out)]
    return energies
end

# Calculate inner physics of Thermal Booster
function calculate_booster(unit::ThermalBooster,
                           sim_params::Dict{String,Any},
                           energies::TBEnergies;
                           fixed_heat_in::Union{Nothing,Integer}=nothing,
                           fixed_heat_out::Union{Nothing,Integer}=nothing)::TBEnergies
    energies = reset_available!(energies; fixed_heat_in, fixed_heat_out)
    energies = reset_temp_slices!(energies)

    src_idx::Int = 1
    snk_idx::Int = 1
    EPS = sim_params["epsilon"]

    ### Calculate Thermal Booster ###
    # loop over input and output layer
    while (src_idx <= length(energies.available_heat_in) && snk_idx <= length(energies.available_heat_out))
        if sum(energies.available_heat_out; init=0.0) < EPS ||
           (sum(energies.available_heat_in; init=0.0) < EPS && !unit.allow_boost_solely)
            # end of condition
            break
        end

        # apply restrictions of control modules for a slice
        if !check_src_to_snk(unit.controller, energies.in_uacs_heat[src_idx], energies.out_uacs[snk_idx])
            snk_idx += 1
            continue
        end

        # check temperatures of output
        skip_slice, demand_temperature = get_layer_temperature(unit, snk_idx,
                                                               energies.out_temps_min,
                                                               energies.out_temps_max;
                                                               input=false)
        if skip_slice || !(snk_idx in energies.out_indices)
            snk_idx += 1
            continue
        end

        demand_delta_temperature = demand_temperature - unit.demand_input_temperature
        if demand_delta_temperature < EPS
            snk_idx += 1
            continue
        end

        # defaults for this slice
        out_energy = 0.0
        low_temp_heat = 0.0
        required_boost_heat = 0.0
        intermediate_temperature = unit.demand_input_temperature
        source_temperature = nothing

        # no heat_in available at all -> solely boost (if allowed)
        boost_only = unit.allow_boost_solely && sum(energies.available_heat_in; init=0.0) < EPS
        no_heat_in_in_next_slices = src_idx == length(energies.available_heat_in) ||
                                    sum(energies.available_heat_in[(src_idx + 1):end]; init=0.0) < EPS
        if boost_only
            out_energy = min(energies.available_heat_out[snk_idx], energies.available_el_in)
            required_boost_heat = out_energy
            low_temp_heat = 0.0
            intermediate_temperature = unit.demand_input_temperature
        else
            # check input slice temperature
            skip_slice, source_temperature = get_layer_temperature(unit, src_idx,
                                                                   energies.in_temps_min,
                                                                   energies.in_temps_max;
                                                                   input=true)
            if skip_slice || !(src_idx in energies.in_indices)
                src_idx += 1
                continue
            end

            # cap intermediate temperature by source and terminal_dT
            intermediate_temperature_cap = min(demand_temperature, source_temperature - unit.terminal_dT)

            # if the source is too cold to contribute (after terminal_dT)
            if intermediate_temperature_cap <= unit.demand_input_temperature + EPS
                if unit.allow_boost_solely && no_heat_in_in_next_slices
                    out_energy = min(energies.available_heat_out[snk_idx], energies.available_el_in)
                    required_boost_heat = out_energy
                    low_temp_heat = 0.0
                    intermediate_temperature = unit.demand_input_temperature
                end
            else
                # max thermal fraction (at intermediate_temperature_cap)
                low_temp_fraction_max = (intermediate_temperature_cap - unit.demand_input_temperature) /
                                        demand_delta_temperature
                low_temp_fraction_max = clamp(low_temp_fraction_max, 0.0, 1.0)
                boost_fraction_min = 1.0 - low_temp_fraction_max

                if unit.allow_boost_additional && no_heat_in_in_next_slices
                    # variable split: Qlow <= low_temp_fraction_max * Qout, Qel = Qout - Qlow
                    electric_limited_out_energy = energies.available_heat_out[snk_idx]
                    if boost_fraction_min > EPS
                        electric_limited_out_energy = energies.available_el_in / boost_fraction_min
                    end

                    out_energy = min(energies.available_heat_out[snk_idx],
                                     energies.available_heat_in[src_idx] + energies.available_el_in,
                                     electric_limited_out_energy)

                    low_temp_heat = min(energies.available_heat_in[src_idx], low_temp_fraction_max * out_energy)
                    required_boost_heat = out_energy - low_temp_heat
                else
                    # fixed split at intermediate_temperature_cap (if any low heat is used)
                    if energies.available_heat_in[src_idx] < EPS
                        if unit.allow_boost_solely && no_heat_in_in_next_slices
                            out_energy = min(energies.available_heat_out[snk_idx], energies.available_el_in)
                            required_boost_heat = out_energy
                            low_temp_heat = 0.0
                            intermediate_temperature = unit.demand_input_temperature
                        end
                    else
                        thermal_limited_out_energy = energies.available_heat_in[src_idx] / low_temp_fraction_max

                        electric_limited_out_energy = energies.available_heat_out[snk_idx]
                        if boost_fraction_min > EPS
                            electric_limited_out_energy = energies.available_el_in / boost_fraction_min
                        end

                        out_energy = min(energies.available_heat_out[snk_idx],
                                         thermal_limited_out_energy,
                                         electric_limited_out_energy)

                        low_temp_heat = low_temp_fraction_max * out_energy
                        required_boost_heat = out_energy - low_temp_heat
                    end
                end
            end
        end

        # if nothing can be delivered for this pairing, move on (avoid stalling)
        if out_energy < EPS
            if source_temperature !== nothing &&
               (source_temperature - unit.terminal_dT) <= unit.demand_input_temperature + EPS &&
               !unit.allow_boost_solely
                src_idx += 1
            else
                snk_idx += 1
            end
            continue
        end

        # compute mass and intermediate temperature from chosen split
        mass_out_current_layer = out_energy / (convert_J_in_Wh(unit.cp_medium_out) * demand_delta_temperature)

        if low_temp_heat > EPS
            intermediate_temperature = unit.demand_input_temperature +
                                       low_temp_heat / (mass_out_current_layer * convert_J_in_Wh(unit.cp_medium_out))
        else
            intermediate_temperature = unit.demand_input_temperature
        end

        # update available energies
        energies.available_el_in -= required_boost_heat
        energies.available_heat_in[src_idx] -= low_temp_heat
        energies.available_heat_out[snk_idx] -= out_energy

        # map to temp slices
        push!(energies.slices_temp_el_in, required_boost_heat)
        push!(energies.slices_temp_heat_in, low_temp_heat)
        push!(energies.slices_temp_heat_in_temperature, source_temperature)
        push!(energies.slices_temp_heat_in_uac, energies.in_uacs_heat[src_idx])
        push!(energies.slices_temp_heat_out, out_energy)
        push!(energies.slices_temp_heat_out_temperature, demand_temperature)
        push!(energies.slices_temp_heat_out_uac, energies.out_uacs[snk_idx])

        push!(unit.intermediate_temperatures, intermediate_temperature)
        push!(unit.intermediate_temperature_energies, out_energy)

        # go to next layer(s)
        if !boost_only && energies.available_heat_in[src_idx] < EPS
            src_idx += 1
        end
        if energies.available_heat_out[snk_idx] < EPS
            snk_idx += 1
        end
    end

    return energies
end

# Calculate Energies
function calculate_energies(unit::ThermalBooster, sim_params::Dict{String,Any})::TBEnergies
    energies = TBEnergies()
    unit.intermediate_temperatures = Temperature[]
    unit.intermediate_temperature_energies = Float64[]

    # get electricity potential and reduce it by constant power draw (or however much
    # is available)
    energies.potential_el_in_layered, energies.in_uacs_el = check_el_in_layered(unit, sim_params)
    energies.potential_el_in = sum(energies.potential_el_in_layered)

    # limit el_in to max usage fraction from control modules and to maximum power_el
    max_usage_fraction = upper_plr_limit(unit.controller, sim_params)
    energies.potential_el_in = min(energies.potential_el_in, unit.power_el * max_usage_fraction)

    if energies.potential_el_in > 0.0     # shortcut if we're limited by zero input electricity
        # get vectored values for the input and output heat potentials
        energies.potentials_heat_in,
        energies.in_temps_min,
        energies.in_temps_max,
        energies.in_uacs_heat = check_heat_in_layered(unit, sim_params)

        energies.potentials_heat_out,
        energies.out_temps_min,
        energies.out_temps_max,
        energies.out_uacs = check_heat_out_layered(unit, sim_params)

        # in the following we want to work with positive values as it is easier
        energies.potentials_heat_in = abs.(energies.potentials_heat_in)
        energies.potentials_heat_out = abs.(energies.potentials_heat_out)

        # reduce available input energies by the power/heat losses that would occur if the
        # sources would be fully utilised. since the actual usage is equal or less than that,
        # it works out even if the PLR for that source is not 1.0
        energies.potentials_heat_in .*= unit.heat_losses_factor
        energies.potential_el_in *= unit.power_losses_factor

        # reorder inputs and outputs according to control modules
        index = reorder_inputs(unit.controller, energies.in_temps_min, energies.in_temps_max)
        energies.in_temps_min = energies.in_temps_min[index]
        energies.in_temps_max = energies.in_temps_max[index]
        energies.in_uacs_heat = energies.in_uacs_heat[index]
        energies.potentials_heat_in = energies.potentials_heat_in[index]

        index = reorder_outputs(unit.controller, energies.out_temps_min, energies.out_temps_max)
        energies.out_temps_min = energies.out_temps_min[index]
        energies.out_temps_max = energies.out_temps_max[index]
        energies.out_uacs = energies.out_uacs[index]
        energies.potentials_heat_out = energies.potentials_heat_out[index]

        # there are three different cases of how to handle the layered approach of operating the
        # thermal booster, depending on wether or not any input heat or output heat transformer has a
        # value of infinite as the potential. if this is the case for both the input and output,
        # the calculation cannot be resolved.
        energies.heat_in_has_inf_energy = any(isinf, filter_by_transformer(energies, sim_params; heat_in=true))
        energies.heat_out_has_inf_energy = any(isinf, filter_by_transformer(energies, sim_params; heat_in=false))

        if energies.heat_in_has_inf_energy && energies.heat_out_has_inf_energy
            @warn "The thermal booster $(unit.uac) has unknown energies in both its inputs and " *
                  "outputs. This cannot be resolved. Please check the order of operation and " *
                  "make sure that either the inputs or the outputs have been fully calculated " *
                  "before the thermal booster $(unit.uac) has its potential step."
            return energies
        end

        if energies.heat_in_has_inf_energy
            for heat_in_idx in eachindex(energies.potentials_heat_in)
                energies = calculate_booster(unit, sim_params, energies; fixed_heat_in=heat_in_idx)
                energies = add_heat_in_temp_to_slices(energies)
                energies = copy_heat_out_temp_to_slices(energies)
            end
        elseif energies.heat_out_has_inf_energy
            for heat_out_idx in eachindex(energies.potentials_heat_out)
                energies = calculate_booster(unit, sim_params, energies; fixed_heat_out=heat_out_idx)
                energies = add_heat_out_temp_to_slices(energies)
                energies = copy_heat_in_temp_to_slices(energies)
            end
        else
            energies = calculate_booster(unit, sim_params, energies)
            energies = copy_temp_to_slices!(energies)
        end
    end

    # now set losses of the thermal booster and add the losses to the actually consumed
    # power / heat for the slices
    el_in = sum(energies.slices_el_in; init=0.0)
    heat_in = sum(energies.slices_heat_in; init=0.0)
    unit.losses_power = -1.0 * el_in
    unit.losses_heat = -1.0 * heat_in
    energies.slices_el_in ./= unit.power_losses_factor
    energies.slices_heat_in ./= unit.heat_losses_factor

    el_in = sum(energies.slices_el_in; init=0.0)
    heat_in = sum(energies.slices_heat_in; init=0.0)
    unit.losses_power += el_in
    unit.losses_heat += heat_in

    return energies
end

function potential(unit::ThermalBooster, sim_params::Dict{String,Any})
    energies = calculate_energies(unit, sim_params)

    if sum(energies.slices_heat_out; init=0.0) < sim_params["epsilon"]
        set_max_energies!(unit, true, sum(energies.slices_el_in; init=0.0), 0.0, 0.0)
        return
    end

    set_max_energies!(unit,
                      true,
                      energies.slices_el_in,
                      energies.slices_heat_in,
                      energies.slices_heat_out,
                      energies.slices_heat_in_temperature,
                      energies.slices_heat_out_temperature,
                      energies.slices_heat_in_uac,
                      energies.slices_heat_out_uac,
                      energies.heat_in_has_inf_energy,
                      energies.heat_out_has_inf_energy)
end

function process(unit::ThermalBooster, sim_params::Dict{String,Any})
    energies = calculate_energies(unit, sim_params)

    el_in = sum(energies.slices_el_in; init=0.0)
    heat_out = sum(energies.slices_heat_out; init=0.0)

    if heat_out < sim_params["epsilon"]
        set_max_energies!(unit, false, el_in, 0.0, 0.0, 0.0)
        sub!(unit.input_interfaces[unit.m_el_in], el_in)
    else
        sub!(unit.input_interfaces[unit.m_el_in], el_in)
        sub!(unit.input_interfaces[unit.m_heat_in],
             energies.slices_heat_in,
             energies.slices_heat_in_temperature,
             [nothing for _ in energies.slices_heat_in_temperature],
             energies.slices_heat_in_uac)
        add!(unit.output_interfaces[unit.m_heat_out],
             energies.slices_heat_out,
             [nothing for _ in energies.slices_heat_out_temperature],
             energies.slices_heat_out_temperature,
             energies.slices_heat_out_uac)
    end
    # calculate total losses
    unit.losses = unit.losses_power + unit.losses_heat

    # calculate mean intermediate temperature
    unit.mean_intermediate_temperature = sum(unit.intermediate_temperatures .* unit.intermediate_temperature_energies) /
                                         sum(unit.intermediate_temperature_energies)
end

# ThermalBooster has its own reset function as here more parameters are present 
# that need to be reset in every timestep
function reset(unit::ThermalBooster)
    invoke(reset, Tuple{Component}, unit)

    # reset other parameter
    unit.losses = 0.0
    unit.losses_heat = 0.0
    unit.losses_power = 0.0
end

function output_values(unit::ThermalBooster)::Vector{String}
    output_vals = [string(unit.m_el_in) * ":IN",
                   string(unit.m_heat_in) * ":IN",
                   string(unit.m_heat_out) * ":OUT"]
    append!(output_vals, ["LossesGains",
                          "Losses_power",
                          "Losses_heat",
                          "mean_intermediate_temperature"])
    return output_vals
end

function output_value(unit::ThermalBooster, key::OutputKey)::Float64
    if key.value_key == "IN"
        return calculate_energy_flow(unit.input_interfaces[key.medium])
    elseif key.value_key == "OUT"
        return calculate_energy_flow(unit.output_interfaces[key.medium])
    elseif key.value_key == "Losses_power"
        return -unit.losses_power
    elseif key.value_key == "Losses_heat"
        return -unit.losses_heat
    elseif key.value_key == "LossesGains"
        return -unit.losses
    elseif key.value_key == "mean_intermediate_temperature"
        return unit.mean_intermediate_temperature
    end
    throw(KeyError(key.value_key))
end

export ThermalBooster
