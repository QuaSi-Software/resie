"""
Implementation of a heat pump component.

Can be given fixed input and output temperatures instead of checking the temperature
potentials of other components. Can also be given a fixed COP instead of dynamically
calculating it as a Carnot-COP by input/output temperatures.
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
    min_power_fraction::Float64
    constant_cop::Floathing
    dynamic_cop::Function
    bypass_cop::Float64

    consider_icing::Bool
    icing_coefficients::Vector{Float64}
    output_temperature::Temperature
    input_temperature::Temperature

    cop::Float64
    losses::Float64
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
                   default(config, "min_power_fraction", 0.2),
                   constant_cop,
                   cop_function,
                   default(config, "bypass_cop", 15.0),
                   default(config, "considering_icing", false),
                   icing_coefficients,
                   default(config, "output_temperature", nothing),
                   default(config, "input_temperature", nothing),
                   0.0, # cop
                   0.0, # losses
                   0.0, # mixing temperature in the input interface
                   0.0) # mixing temperature in the output interface
    end
end

mutable struct HPEnergies
    potential_energy_el::Float64
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

function calculate_slices(unit::HeatPump,
                          sim_params::Dict{String,Any},
                          energies::HPEnergies,
                          plrs::Array{<:Floathing,2})::Tuple{HPEnergies,Vector{Float64},Vector{Float64},
                                                             Array{<:Floathing,2}}

    # times_min means the time slice assuming minimum power. times_min are actually larger
    # than times_max
    times_min = Vector{Float64}()
    times_max = Vector{Float64}()
    sum_usage = 0.0
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
        # not been reached
        if energies.available_el_in < EPS ||
           sum(energies.available_heat_in; init=0.0) < EPS ||
           sum(energies.available_heat_out; init=0.0) < EPS ||
           energies.max_usage_fraction - sum_usage < EPS
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

        min_power_frac = min(1.0, unit.min_power_function(current_in_temp, current_out_temp))
        min_power = max(0.0, unit.design_power_th * min_power_frac)
        max_power_frac = min(1.0, unit.max_power_function(current_in_temp, current_out_temp))
        max_power = max(min_power, unit.design_power_th * max_power_frac)

        remaining_heat_out = min(energies.available_heat_out[current_out_idx],
                                 (energies.max_usage_fraction - sum_usage) * watt_to_wh(max_power))

        if plrs[current_in_idx, current_out_idx] === nothing
            plrs[current_in_idx, current_out_idx] = 1.0
        end

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

        sum_usage += used_heat_out / watt_to_wh(unit.design_power_th)
        push!(times_min, used_heat_out * 3600 / min_power)
        push!(times_max, used_heat_out * 3600 / max_power)
    end

    return energies, times_min, times_max, plrs
end

function calculate_energies_heatpump(unit::HeatPump,
                                     sim_params::Dict{String,Any},
                                     energies::HPEnergies)::HPEnergies
    energies.slices_el_in_temp = Vector{Floathing}()
    energies.slices_heat_in_temp = Vector{Floathing}()
    energies.slices_heat_in_temperature_temp = Vector{Temperature}()
    energies.slices_heat_in_uac_temp = Vector{Stringing}()
    energies.slices_heat_out_temp = Vector{Floathing}()
    energies.slices_heat_out_temperature_temp = Vector{Temperature}()
    energies.slices_heat_out_uac_temp = Vector{Stringing}()

    # initial PLRs are nothing, as we do not yet know which slices are used and which aren't
    plrs = Array{Floathing,2}(nothing,
                              length(energies.available_heat_in),
                              length(energies.available_heat_out))

    energies, times_min, times_max, plrs = calculate_slices(unit, sim_params, energies, plrs)

    # as long as sum of times with minimum power per slice is larger than the time step
    # multiplied with the min power fraction, we can find a dispatch of each slice by
    # "slowing them down", so the sum exceeds the minimum power fraction. if not, the heat
    # pump should not run at all
    if sum(times_min; init=0.0) < unit.min_power_fraction * sim_params["time_step_seconds"]
        energies.slices_el_in_temp = []
        energies.slices_heat_in_temp = []
        energies.slices_heat_in_temperature_temp = []
        energies.slices_heat_in_uac_temp = []
        energies.slices_heat_out_temp = []
        energies.slices_heat_out_temperature_temp = []
        energies.slices_heat_out_uac_temp = []
    end

    return energies
end

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
    # heat pump, depending on wether or not the input or output heat has infinite values. if
    # both have infinite values the calculation cannot be resolved.
    energies.heat_in_has_inf_energy = any(isinf, energies.potentials_energies_heat_in)
    energies.heat_out_has_inf_energy = any(isinf, energies.potentials_energies_heat_out)

    if energies.heat_in_has_inf_energy && energies.heat_out_has_inf_energy
        @warn "The heat pump $(unit.uac) has unknown energies in both its inputs and " *
              "outputs. This cannot be resolved. Please check the order of operation and " *
              "make sure that either the inputs or the outputs have been fully calculated " *
              "before the heat pump $(unit.uac) has its potential step."
        return false, energies
    end

    if energies.heat_in_has_inf_energy
        for heat_in_idx in eachindex(energies.potentials_energies_heat_in)
            energies.available_el_in = copy(energies.potential_energy_el)
            energies.available_heat_in = copy(energies.potentials_energies_heat_in)
            energies.available_heat_out = copy(energies.potentials_energies_heat_out)
            energies.available_heat_in[heat_in_idx] = Inf
            energies.in_indices = [heat_in_idx]
            energies.out_indices = [Int(i) for i in eachindex(energies.available_heat_out)]
            energies = calculate_energies_heatpump(unit, sim_params, energies)

            append!(energies.slices_heat_in, energies.slices_heat_in_temp)
            append!(energies.slices_heat_in_temperature, energies.slices_heat_in_temperature_temp)
            append!(energies.slices_heat_in_uac, energies.slices_heat_in_uac_temp)
            # Is this correct? Using the highest energy as worst case should be good for now...
            if heat_in_idx == 1
                energies.slices_el_in = energies.slices_el_in_temp
                energies.slices_heat_out = energies.slices_heat_out_temp
                energies.slices_heat_out_temperature = energies.slices_heat_out_temperature_temp
                energies.slices_heat_out_uac = energies.slices_heat_out_uac_temp
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
        end
    elseif energies.heat_out_has_inf_energy
        for heat_out_idx in eachindex(energies.potentials_energies_heat_out)
            energies.available_el_in = copy(energies.potential_energy_el)
            energies.available_heat_in = copy(energies.potentials_energies_heat_in)
            energies.available_heat_out = copy(energies.potentials_energies_heat_out)
            energies.available_heat_out[heat_out_idx] = Inf
            energies.in_indices = [Int(i) for i in eachindex(energies.available_heat_in)]
            energies.out_indices = [heat_out_idx]
            energies = calculate_energies_heatpump(unit, sim_params, energies)

            append!(energies.slices_heat_out, energies.slices_heat_out_temp)
            append!(energies.slices_heat_out_temperature, energies.slices_heat_out_temperature_temp)
            append!(energies.slices_heat_out_uac, energies.slices_heat_out_uac_temp)
            # Is this correct? Using the highest energy as worst case should be good for now...
            if heat_out_idx == 1
                energies.slices_el_in = energies.slices_el_in_temp
                energies.slices_heat_in = energies.slices_heat_in_temp
                energies.slices_heat_in_temperature = energies.slices_heat_in_temperature_temp
                energies.slices_heat_in_uac = energies.slices_heat_in_uac_temp
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
        end
    else # fully known energies and temperatures in inputs and outputs
        energies.available_el_in = copy(energies.potential_energy_el)
        energies.available_heat_in = copy(energies.potentials_energies_heat_in)
        energies.available_heat_out = copy(energies.potentials_energies_heat_out)
        energies.in_indices = [Int(i) for i in eachindex(energies.available_heat_in)]
        energies.out_indices = [Int(i) for i in eachindex(energies.available_heat_out)]
        energies = calculate_energies_heatpump(unit, sim_params, energies)

        energies.slices_el_in = energies.slices_el_in_temp
        energies.slices_heat_in = energies.slices_heat_in_temp
        energies.slices_heat_in_temperature = energies.slices_heat_in_temperature_temp
        energies.slices_heat_in_uac = energies.slices_heat_in_uac_temp
        energies.slices_heat_out = energies.slices_heat_out_temp
        energies.slices_heat_out_temperature = energies.slices_heat_out_temperature_temp
        energies.slices_heat_out_uac = energies.slices_heat_out_uac_temp
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

    sub!(unit.input_interfaces[unit.m_el_in], energies.slices_el_in)
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
            "Losses",
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
    elseif key.value_key == "Losses"
        return unit.losses
    elseif key.value_key == "MixingTemperature_Input"
        return unit.mix_temp_input
    elseif key.value_key == "MixingTemperature_Output"
        return unit.mix_temp_output
    end
    throw(KeyError(key.value_key))
end

export HeatPump
