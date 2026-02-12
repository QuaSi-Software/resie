using ..Resie: get_run

"""
Thermal Booster

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
    demand_input_temperature::Float64
    allow_boost_only::Bool

    losses::Float64
    losses_power::Float64
    losses_heat::Float64

    function ThermalBooster(uac::String, config::Dict{String,Any}, sim_params::Dict{String,Any})
        m_el_in = Symbol(default(config, "m_el_in", "m_e_ac_230v"))
        m_heat_out = Symbol(default(config, "m_heat_out", "m_h_w_ht1"))
        m_heat_in = Symbol(default(config, "m_heat_in", "m_h_w_lt1"))
        register_media([m_el_in, m_heat_out, m_heat_in])

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
                   default(config, "demand_input_temperature", 12.0), # [°C] TODO Could also be a profile
                   default(config, "allow_boost_only", false))
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
    slices_temp_times::Vector{Floathing}
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
                           el_in::Union{Floathing,Vector{<:Floathing}},
                           heat_in::Union{Floathing,Vector{<:Floathing}},
                           heat_out::Union{Floathing,Vector{<:Floathing}},
                           slices_heat_in_temperature::Union{Temperature,Vector{<:Temperature}}=nothing,
                           slices_heat_out_temperature::Union{Temperature,Vector{<:Temperature}}=nothing,
                           purpose_uac_heat_in::Union{Stringing,Vector{Stringing}}=nothing,
                           purpose_uac_heat_out::Union{Stringing,Vector{Stringing}}=nothing,
                           has_calculated_all_maxima_heat_in::Bool=false,
                           has_calculated_all_maxima_heat_out::Bool=false)
    set_max_energy!(unit.input_interfaces[unit.m_el_in], el_in)
    set_max_energy!(unit.input_interfaces[unit.m_heat_in], heat_in, slices_heat_in_temperature, nothing,
                    purpose_uac_heat_in, has_calculated_all_maxima_heat_in)
    set_max_energy!(unit.output_interfaces[unit.m_heat_out], heat_out, nothing, slices_heat_out_temperature,
                    purpose_uac_heat_out, has_calculated_all_maxima_heat_out)
end

function calculate_energies(unit::ThermalBooster, sim_params::Dict{String,Any})::TBEnergies
    energies = TBEnergies()

    # get electricity potential and reduce it by constant power draw (or however much
    # is available)
    energies.potential_el_in_layered, energies.in_uacs_el = check_el_in_layered(unit, sim_params)
    energies.potential_el_in = sum(energies.potential_el_in_layered)
    energies.potential_el_in = min(unit.power_el, energies.potential_el_in)

    # shortcut if we're limited by zero input electricity
    if energies.potential_el_in <= 0.0
        do_calculation = false
    else
        do_calculation = true
    end

    if do_calculation
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
            end
        elseif energies.heat_out_has_inf_energy
            for heat_out_idx in eachindex(energies.potentials_heat_out)
                energies = calculate_booster(unit, sim_params, energies; fixed_heat_out=heat_out_idx)
            end
        else
            energies = calculate_booster(unit, sim_params, energies)
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

function calculate_booster(unit::ThermalBooster,
                           sim_params::Dict{String,Any},
                           energies::TBEnergies;
                           fixed_heat_in::Union{Nothing,Integer}=nothing,
                           fixed_heat_out::Union{Nothing,Integer}=nothing)::TBEnergies
    energies = reset_available!(energies; fixed_heat_in, fixed_heat_out)
    src_idx::Int = 1
    snk_idx::Int = 1
    EPS = sim_params["epsilon"]

    # loop over input and output layer
    while (src_idx <= length(energies.available_heat_in) && snk_idx <= length(energies.available_heat_out))
        if energies.available_el_in < EPS || sum(energies.available_heat_in; init=0.0) < EPS ||
           sum(energies.available_heat_out; init=0.0) < EPS
            # end of condition
            break
        end

        # apply restrictions of control modules for a slice
        if !check_src_to_snk(unit.controller, energies.in_uacs_heat[src_idx], energies.out_uacs[snk_idx])
            snk_idx += 1
            continue
        end

        # check temperatures
        skip_slice, src_temperature = get_layer_temperature(unit, src_idx,
                                                            energies.in_temps_min,
                                                            energies.in_temps_max;
                                                            input=true)
        if skip_slice || !(src_idx in energies.in_indices)
            src_idx += 1
            continue
        end

        # same for output
        skip_slice, snk_temperature = get_layer_temperature(unit, snk_idx,
                                                            energies.out_temps_min,
                                                            energies.out_temps_max;
                                                            input=false)
        if skip_slice || !(snk_idx in energies.out_indices)
            snk_idx += 1
            continue
        end

        # reduce source temperature by terminal_dT
        src_temperature_reduced = src_temperature - unit.terminal_dT
        # TODO check for negative src_temperature_reduced

        # calculate required mass in output for current slice
        mass_out_current_layer = energies.available_heat_out[snk_idx] / (convert_J_in_Wh(unit.cp_medium_out) *
                                                                         (snk_temperature - unit.demand_input_temperature))

        # calculate required energy from thermal input for current slice
        required_low_temp_heat = mass_out_current_layer * convert_J_in_Wh(unit.cp_medium_out) *
                                 (src_temperature_reduced - unit.demand_input_temperature)
        # limit to maximum available energy
        required_low_temp_heat = min(required_low_temp_heat, energies.available_heat_in[src_idx])
        if unit.allow_boost_only
            # calculate intermediate temperature to start from boosting while keeping the mass_out_current_layer
            src_temperature_reduced = required_low_temp_heat /
                                      (mass_out_current_layer * convert_J_in_Wh(unit.cp_medium_out)) +
                                      unit.demand_input_temperature
        else
            # recalculate output layer mass that can be heated
            mass_out_current_layer = required_low_temp_heat / (convert_J_in_Wh(unit.cp_medium_out) *
                                                               (src_temperature_reduced - unit.demand_input_temperature))
        end

        # calculate the remaining energy to boost the input temperature to the output temperature in current slice
        # handle cases where output is fully provided by power with toggle!
        required_boost_heat = mass_out_current_layer * convert_J_in_Wh(unit.cp_medium_out) *
                              (snk_temperature - src_temperature_reduced)
        # limit to maximum available boost energy
        if energies.available_el_in < required_boost_heat
            # limit
            required_boost_heat = min(required_boost_heat, energies.available_el_in)
            # recalculate mass than can be heated
            mass_out_current_layer = required_boost_heat /
                                     (convert_J_in_Wh(unit.cp_medium_out) * (snk_temperature - src_temperature_reduced))
            # update also low temp heat
        end

        # calculate output energy
        out_energy = mass_out_current_layer * convert_J_in_Wh(unit.cp_medium_out) *
                     (snk_temperature - unit.demand_input_temperature)
        if (required_low_temp_heat + required_boost_heat) - out_energy > EPS
            @error "something went wrong in the thermal booster... \nrequired_low_temp_heat: $(required_low_temp_heat) \nrequired_boost_heat: $(required_boost_heat) \nout_energy: $(out_energy)"
            # TODO remove later...
        end

        # update available energies
        energies.available_el_in -= required_boost_heat
        energies.available_heat_in[src_idx] -= required_low_temp_heat
        energies.available_heat_out[snk_idx] -= out_energy

        # map to slices
        push!(energies.slices_el_in, required_boost_heat)
        push!(energies.slices_heat_in, required_low_temp_heat)
        push!(energies.slices_heat_in_temperature, src_temperature)
        push!(energies.slices_heat_in_uac, energies.in_uacs_heat[src_idx])
        push!(energies.slices_heat_out, out_energy)
        push!(energies.slices_heat_out_temperature, snk_temperature)
        push!(energies.slices_heat_out_uac, energies.out_uacs[snk_idx])

        # go to next layer(s)
        if energies.available_heat_in[src_idx] < EPS
            src_idx += 1
        end
        if energies.available_heat_out[snk_idx] < EPS
            snk_idx += 1
        end
    end

    return energies
end

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

function potential(unit::ThermalBooster, sim_params::Dict{String,Any})
    energies = calculate_energies(unit, sim_params)

    if sum(energies.slices_heat_out; init=0.0) < sim_params["epsilon"]
        set_max_energies!(unit, sum(energies.slices_el_in; init=0.0), 0.0, 0.0)
        return
    end

    set_max_energies!(unit,
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
        set_max_energies!(unit, el_in, 0.0, 0.0, 0.0)
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
end

# has its own reset function as here more parameters are present that need to be reset in
# every timestep
function reset(unit::ThermalBooster)
    invoke(reset, Tuple{Component}, unit)

    # reset other parameter
    unit.losses = 0.0
end

function output_values(unit::ThermalBooster)::Vector{String}
    output_vals = [string(unit.m_el_in) * ":IN",
                   string(unit.m_heat_in) * ":IN",
                   string(unit.m_heat_out) * ":OUT"]
    append!(output_vals, ["LossesGains",
                          "Losses_power",
                          "Losses_heat"])
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
    end
    throw(KeyError(key.value_key))
end

export ThermalBooster
