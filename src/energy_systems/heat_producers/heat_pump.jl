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

    power_th::Float64
    min_power_fraction::Float64
    min_run_time::UInt
    constant_cop::Any
    output_temperature::Temperature
    input_temperature::Temperature
    cop::Float64

    losses::Float64

    function HeatPump(uac::String, config::Dict{String,Any}, sim_params::Dict{String,Any})
        m_el_in = Symbol(default(config, "m_el_in", "m_e_ac_230v"))
        m_heat_out = Symbol(default(config, "m_heat_out", "m_h_w_ht1"))
        m_heat_in = Symbol(default(config, "m_heat_in", "m_h_w_lt1"))
        register_media([m_el_in, m_heat_out, m_heat_in])

        return new(
            uac, # uac
            controller_for_strategy( # controller
                config["strategy"]["name"], config["strategy"], sim_params
            ),
            sf_transformer, # sys_function
            InterfaceMap( # input_interfaces
                m_heat_in => nothing,
                m_el_in => nothing
            ),
            InterfaceMap( # output_interfaces
                m_heat_out => nothing
            ),
            m_el_in,
            m_heat_out,
            m_heat_in,
            config["power_th"], # power_th
            default(config, "min_power_fraction", 0.2),
            default(config, "min_run_time", 0),
            default(config, "constant_cop", nothing),
            default(config, "output_temperature", nothing),
            default(config, "input_temperature", nothing),
            0.0, # cop
            0.0, # losses
        )
    end
end

function initialise!(unit::HeatPump, sim_params::Dict{String,Any})
    set_storage_transfer!(
        unit.input_interfaces[unit.m_heat_in],
        default(
            unit.controller.parameter, "unload_storages " * String(unit.m_heat_in), true
        )
    )
    set_storage_transfer!(
        unit.input_interfaces[unit.m_el_in],
        default(
            unit.controller.parameter, "unload_storages " * String(unit.m_el_in), true
        )
    )
    set_storage_transfer!(
        unit.output_interfaces[unit.m_heat_out],
        default(
            unit.controller.parameter, "load_storages " * String(unit.m_heat_out), true
        )
    )
end

function control(
    unit::HeatPump,
    components::Grouping,
    sim_params::Dict{String,Any}
)
    move_state(unit, components, sim_params)

    # for fixed input/output temperatures, overwrite the interface with those. otherwise
    # highest will choose the interface's temperature (including nothing)
    if unit.output_temperature !== nothing
        set_temperature!(
            unit.output_interfaces[unit.m_heat_out],
            unit.output_temperature,
            unit.output_temperature
        )
    end
    if unit.input_temperature !== nothing
        set_temperature!(
            unit.input_interfaces[unit.m_heat_in],
            unit.input_temperature,
            unit.input_temperature
        )
    end
end

function set_max_energies!(
    unit::HeatPump,
    el_in::Union{Floathing, Vector{Floathing}},
    heat_in::Union{Floathing, Vector{Floathing}},
    heat_out::Union{Floathing, Vector{Floathing}},
    purpose_uac_heat_in::Union{Stringing, Vector{Stringing}}=nothing,
    purpose_uac_heat_out::Union{Stringing, Vector{Stringing}}=nothing,
    temperature_heat_in::Union{Temperature, Vector{Temperature}}=nothing,
    temperature_heat_out::Union{Temperature, Vector{Temperature}}=nothing,
    has_calculated_all_maxima_heat_in::Bool=false,
    has_calculated_all_maxima_heat_out::Bool=false
)
    set_max_energy!(unit.input_interfaces[unit.m_el_in], el_in)
    set_max_energy!(unit.input_interfaces[unit.m_heat_in], heat_in, purpose_uac_heat_in, temperature_heat_in, has_calculated_all_maxima_heat_in)
    set_max_energy!(unit.output_interfaces[unit.m_heat_out], heat_out, purpose_uac_heat_out, temperature_heat_out, has_calculated_all_maxima_heat_out)
end

function dynamic_cop(in_temp::Temperature, out_temp::Temperature)::Union{Nothing,Float64}
    if (in_temp === nothing || out_temp === nothing)
        return nothing
    end

    # Carnot-COP with 40 % efficiency
    return 0.4 * (273.15 + out_temp) / (out_temp - in_temp)
end

function check_el_in(
    unit::HeatPump,
    sim_params::Dict{String,Any}
)
    if unit.controller.parameter["consider_m_el_in"] == true
        if (unit.input_interfaces[unit.m_el_in].source.sys_function == sf_transformer  # HP has direct connection to a transfomer...
            && is_max_energy_nothing(unit.input_interfaces[unit.m_el_in].max_energy)   # ...and none of them have had their potential step
        )
            return (Inf)
        else
            exchanges = balance_on(
                unit.input_interfaces[unit.m_el_in],
                unit.input_interfaces[unit.m_el_in].source
            )
            potential_energy_el = balance(exchanges) + energy_potential(exchanges)
            if potential_energy_el <= sim_params["epsilon"]
                return (0.0)
            end
            return (potential_energy_el)
        end
    else
        return (Inf)
    end
end

function check_heat_in(
    unit::HeatPump,
    sim_params::Dict{String,Any}
)
    if unit.controller.parameter["consider_m_heat_in"] == true
        if (unit.input_interfaces[unit.m_heat_in].source.sys_function == sf_transformer  # HP has direct connection to a transfomer...
             && is_max_energy_nothing(unit.input_interfaces[unit.m_heat_in].max_energy)  # ...and none of them have had their potential step
        )
            return ([Inf],
                    [unit.input_interfaces[unit.m_heat_in].temperature_min], 
                    [unit.input_interfaces[unit.m_heat_in].temperature_max],
                    [unit.input_interfaces[unit.m_heat_in].source.uac]
                )
        else
            exchanges = balance_on(
                unit.input_interfaces[unit.m_heat_in],
                unit.input_interfaces[unit.m_heat_in].source
            )
            return (
                [e.balance + e.energy_potential for e in exchanges],
                temp_min_all(exchanges),
                temp_max_all(exchanges),
                [e.purpose_uac for e in exchanges]
            )
        end
    else
        return ([Inf],
                [unit.input_interfaces[unit.m_heat_in].temperature_min], 
                [unit.input_interfaces[unit.m_heat_in].temperature_max],
                [unit.input_interfaces[unit.m_heat_in].source.uac]
            )
    end
end

function check_heat_out(
    unit::HeatPump,
    sim_params::Dict{String,Any}
)
    if unit.controller.parameter["consider_m_heat_out"] == true
        if (unit.output_interfaces[unit.m_heat_out].target.sys_function == sf_transformer   # HP has direct connection to a transfomer...
            && is_max_energy_nothing(unit.output_interfaces[unit.m_heat_out].max_energy)    # ...and none of them have had their potential step
        )
            return ([-Inf],
                    [unit.output_interfaces[unit.m_heat_out].temperature_min], 
                    [unit.output_interfaces[unit.m_heat_out].temperature_max],
                    [unit.output_interfaces[unit.m_heat_out].target.uac]
                )
        else
            exchanges = balance_on(
                unit.output_interfaces[unit.m_heat_out],
                unit.output_interfaces[unit.m_heat_out].target
            )
            return (
                [e.balance + e.energy_potential for e in exchanges],
                temp_min_all(exchanges),
                temp_max_all(exchanges),
                [e.purpose_uac for e in exchanges]
            )
        end
    else
        return ([-Inf],
                [unit.output_interfaces[unit.m_heat_out].temperature_min], 
                [unit.output_interfaces[unit.m_heat_out].temperature_max],
                [unit.output_interfaces[unit.m_heat_out].target.uac]
                )
    end
end


function calculate_energies_heatpump(unit::HeatPump,
                                     sim_params::Dict{String,Any},
                                     available_el_in,
                                     available_heat_in,
                                     available_heat_out,
                                     max_usage_fraction,
                                     in_temps_min,
                                     in_temps_max,
                                     in_uacs,
                                     out_temps_min,
                                     out_temps_max,
                                     out_uacs)
    layers_el_in = Vector{Floathing}()
    layers_heat_in = Vector{Floathing}()
    layers_heat_in_temperature = Vector{Temperature}()
    layers_heat_in_uac = Vector{Stringing}()
    layers_heat_out = Vector{Floathing}()
    layers_heat_out_temperature = Vector{Temperature}()
    layers_heat_out_uac = Vector{Stringing}()

    current_in_idx = 1
    current_out_idx = 1
    while (
           sum(available_el_in; init=0.0) > sim_params["epsilon"]
           && sum(available_heat_in; init=0.0) > sim_params["epsilon"]
           && sum(available_heat_out; init=0.0) > sim_params["epsilon"]
           && max_usage_fraction * watt_to_wh(unit.power_th) - sum(layers_heat_out; init=0.0) > sim_params["epsilon"]
       )
        # find first non-zero index
        while available_heat_in[current_in_idx] <= sim_params["epsilon"]
            current_in_idx += 1
        end
        while available_heat_out[current_out_idx] <= sim_params["epsilon"]
            current_out_idx += 1
        end

        # detect and check temperatures. If a input or an output temperature is given, this will be set!
        if unit.input_temperature === nothing
            current_in_temp = highest(in_temps_min[current_in_idx], in_temps_max[current_in_idx])
            if current_in_temp === nothing && unit.constant_cop === nothing
                @error "Error: The input temperature for $(unit.uac) could not be detected. Please specify one with the parameter 'input_temperature' or check the connected components."
                throw(InputError)
            end
        else
            # skip layer if given fixed input temperature is not within the temperature band of the layer
            if (
                in_temps_min[current_in_idx] !== nothing
                   && in_temps_min[current_in_idx] > unit.input_temperature
                || in_temps_max[current_in_idx] !== nothing
                   && in_temps_max[current_in_idx] < unit.input_temperature
            )
                available_heat_in[current_in_idx] = 0.0
                current_in_idx += 1
                continue
            else
                current_in_temp = unit.input_temperature
            end
        end

        if unit.output_temperature === nothing
            current_out_temp = lowest(out_temps_min[current_out_idx], out_temps_max[current_out_idx])
            if current_out_temp === nothing && unit.constant_cop === nothing
                @error "Error: The output temperature for $(unit.uac) could not be detected. Please specify one with the parameter 'output_temperature' or check the connected components."
                throw(InputError)
            end
        else
            # skip layer if given fixed output temperature is not within the temperature band of the layer
            if (
                out_temps_min[current_out_idx] !== nothing
                    && out_temps_min[current_out_idx] > unit.output_temperature
                || out_temps_max[current_out_idx] !== nothing
                    && out_temps_max[current_out_idx] < unit.output_temperature
            )
                available_heat_out[current_out_idx] = 0.0
                current_out_idx += 1
                continue
            else
                current_out_temp = unit.output_temperature
            end
        end 

        if current_in_temp >= current_out_temp 
            # bypass
            heat_transfer = minimum([
                available_heat_in[current_in_idx],
                available_heat_out[current_out_idx],
                max_usage_fraction * watt_to_wh(unit.power_th) - sum(layers_heat_out; init=0.0)
            ])
            used_heat_in = heat_transfer
            used_heat_out = heat_transfer
            # The el. energy can not be zero as this would lead to problems in process step
            # as then it seems that no electrical energy is available!
            used_el_in = 0.01 # TODO
            current_out_temp = current_in_temp
            # TODO: Whats about the usage_fraction, should the bypass be kept included in the usage_fraction calculation?
            # Or should we deny the bypass and force users to implement heat pumps in parallel if they want a bypass?
        else
            # calculate cop
            cop = unit.constant_cop === nothing ?
                dynamic_cop(current_in_temp, current_out_temp) :
                unit.constant_cop
            if cop === nothing
                @error ("Input and/or output temperature for heatpump $(unit.uac) is not given. Provide temperatures or fixed cop.")
                throw(InputError)
            end

            # calculate energies with the current cop
            # energies for current layer with potential heat in as basis
            used_heat_in = copy(available_heat_in[current_in_idx])
            used_el_in = used_heat_in / (cop - 1.0)
            used_heat_out = used_heat_in + used_el_in

            # check heat out as limiter, also checking for limit of heat pump
            max_heat_out = min(available_heat_out[current_out_idx], 
                               max_usage_fraction * watt_to_wh(unit.power_th) - sum(layers_heat_out; init=0.0))
            if used_heat_out > max_heat_out
                used_heat_out = max_heat_out
                used_el_in = used_heat_out / cop
                used_heat_in = used_el_in * (cop - 1.0)
            end

            # check electricity in as limiter
            if used_el_in > available_el_in
                used_el_in = available_el_in
                used_heat_in = used_el_in * (cop - 1.0)
                used_heat_out = used_heat_in + used_el_in
            end
        end

        # finally all checks done, we add the layer and update remaining energies
        push!(layers_el_in, used_el_in)
        push!(layers_heat_in, used_heat_in)
        push!(layers_heat_in_temperature, current_in_temp)
        push!(layers_heat_in_uac, in_uacs[current_in_idx])
        push!(layers_heat_out, used_heat_out)
        push!(layers_heat_out_temperature, current_out_temp)
        push!(layers_heat_out_uac, out_uacs[current_out_idx])

        available_el_in -= used_el_in
        available_heat_in[current_in_idx] -= used_heat_in
        available_heat_out[current_out_idx] -= used_heat_out
    end

    return layers_el_in,
           layers_heat_in,
           layers_heat_in_temperature,
           layers_heat_in_uac,
           layers_heat_out,
           layers_heat_out_temperature,
           layers_heat_out_uac
end


function calculate_energies(
    unit::HeatPump,
    sim_params::Dict{String,Any}
)
    # check operational strategy, specifically storage_driven
    if (
        unit.controller.strategy == "storage_driven"
        && unit.controller.state_machine.state != 2
    )
        return false, (nothing, nothing, nothing)
    end

    # get usage fraction of external profile (normalized from 0 to 1, related to unit.power_th)
    max_usage_fraction = (
        unit.controller.parameter["operation_profile_path"] === nothing
        ? 1.0
        : value_at_time(unit.controller.parameter["operation_profile"], sim_params["time"])
    )
    if max_usage_fraction <= 0.0
        return false, (nothing, nothing, nothing)
    end

    # get potentials from inputs/outputs. only the heat input is calculated as vector,
    # the electricity input and heat output are calculated as scalars
    potential_energy_el = check_el_in(unit, sim_params)
    potentials_energies_heat_in,
        in_temps_min,
        in_temps_max,
        in_uacs = check_heat_in(unit, sim_params)
    potential_energies_heat_out, 
        out_temps_min,
        out_temps_max,
        out_uacs = check_heat_out(unit, sim_params)

    # in the following we want to work with positive values as it is easier
    available_el_in = copy(potential_energy_el)
    available_heat_in = copy(abs.(potentials_energies_heat_in))
    available_heat_out = copy(abs.(potential_energies_heat_out))

    heat_in_has_inf_energy = any(isinf, available_heat_in)
    heat_out_has_inf_energy = any(isinf, available_heat_out)

    if heat_in_has_inf_energy && heat_out_has_inf_energy
        # can not perform calculation if both inputs and outputs have inf energies
        @warn "The heat pump $(unit.uac) has unknown energies in both its inputs and outputs. This cannot be resolved. 
Please check the order of operation and make sure that either the inputs or the outputs have been fully calculated 
before the heat pump $(unit.uac) has its potential step."
    elseif heat_in_has_inf_energy 
        layers_el_in = Vector{Floathing}()
        layers_heat_in = Vector{Floathing}()
        layers_heat_in_temperature = Vector{Temperature}()
        layers_heat_in_uac = Vector{Stringing}()
        layers_heat_out = Vector{Floathing}()
        layers_heat_out_temperature = Vector{Temperature}()
        layers_heat_out_uac = Vector{Stringing}()

        for heat_in_idx in eachindex(available_heat_in)
            layers_el_in_temp,
            layers_heat_in_temp,
            layers_heat_in_temperature_temp,
            layers_heat_in_uac_temp,
            layers_heat_out_temp,
            layers_heat_out_temperature_temp,
            layers_heat_out_uac_temp  = calculate_energies_heatpump(unit,
                                                            sim_params,
                                                            copy(available_el_in),
                                                            [Inf],
                                                            copy(available_heat_out),
                                                            max_usage_fraction,
                                                            [in_temps_min[heat_in_idx]],
                                                            [in_temps_max[heat_in_idx]],
                                                            [in_uacs[heat_in_idx]],
                                                            out_temps_min,
                                                            out_temps_max,
                                                            out_uacs)
                                                        

            append!(layers_heat_in, layers_heat_in_temp)
            append!(layers_heat_in_temperature, layers_heat_in_temperature_temp)
            append!(layers_heat_in_uac, layers_heat_in_uac_temp)
            if heat_in_idx == 1 #|| layers_heat_out[1] < layers_heat_out_temp[1]   # TODO: Think about if this is right...
                layers_el_in = layers_el_in_temp
                layers_heat_out = layers_heat_out_temp
                layers_heat_out_temperature = layers_heat_out_temperature_temp
                layers_heat_out_uac = layers_heat_out_uac_temp
            end
        end
    elseif heat_out_has_inf_energy
        layers_el_in = Vector{Floathing}()
        layers_heat_in = Vector{Floathing}()
        layers_heat_in_temperature = Vector{Temperature}()
        layers_heat_in_uac = Vector{Stringing}()
        layers_heat_out = Vector{Floathing}()
        layers_heat_out_temperature = Vector{Temperature}()
        layers_heat_out_uac = Vector{Stringing}()

        for heat_out_idx in eachindex(available_heat_out)
            layers_el_in_temp,
            layers_heat_in_temp,
            layers_heat_in_temperature_temp,
            layers_heat_in_uac_temp,
            layers_heat_out_temp,
            layers_heat_out_temperature_temp,
            layers_heat_out_uac_temp  = calculate_energies_heatpump(unit,
                                                            sim_params,
                                                            copy(available_el_in),
                                                            copy(available_heat_in),
                                                            [Inf],
                                                            max_usage_fraction,
                                                            in_temps_min,
                                                            in_temps_max,
                                                            in_uacs,
                                                            [out_temps_min[heat_out_idx]],
                                                            [out_temps_max[heat_out_idx]],
                                                            [out_uacs[heat_out_idx]])
                                                        

            append!(layers_heat_out, layers_heat_out_temp)
            append!(layers_heat_out_temperature, layers_heat_out_temperature_temp)
            append!(layers_heat_out_uac, layers_heat_out_uac_temp)
            if heat_out_idx == 1 #|| layers_heat_in[1] < layers_heat_in_temp[1]   # TODO: Think about if this is right...
                layers_el_in = layers_el_in_temp
                layers_heat_in = layers_heat_in_temp
                layers_heat_in_temperature = layers_heat_in_temperature_temp
                layers_heat_in_uac = layers_heat_in_uac_temp
            end
        end
    else
        layers_el_in,
        layers_heat_in,
        layers_heat_in_temperature,
        layers_heat_in_uac,
        layers_heat_out,
        layers_heat_out_temperature,
        layers_heat_out_uac  = calculate_energies_heatpump(unit,
                                                           sim_params,
                                                           available_el_in,
                                                           available_heat_in,
                                                           available_heat_out,
                                                           max_usage_fraction,
                                                           in_temps_min,
                                                           in_temps_max,
                                                           in_uacs,
                                                           out_temps_min,
                                                           out_temps_max,
                                                           out_uacs)
    end

    # if all chosen heat layers combined are not enough to meet minimum power fraction,
    # the heat pump doesn't run at all
    usage_fraction = (sum(layers_heat_out; init=0.0)) / watt_to_wh(unit.power_th)
    if usage_fraction < unit.min_power_fraction
        return false, (nothing, nothing, nothing)
    end

    return true,                          
           (layers_el_in,                 # 1
           layers_heat_in,                # 2 
           layers_heat_in_temperature,    # 3
           layers_heat_in_uac,            # 4
           layers_heat_out,               # 5
           layers_heat_out_temperature,   # 6
           layers_heat_out_uac,           # 7
           heat_in_has_inf_energy,        # 8
           heat_out_has_inf_energy        # 9 
       )
end

function potential(unit::HeatPump, sim_params::Dict{String,Any})
    success, energies = calculate_energies(unit, sim_params)

    if !success
        set_max_energies!(unit, 0.0, 0.0, 0.0)
    else
        set_max_energies!(
            unit,
            energies[1],
            energies[2],
            energies[5],
            energies[4],
            energies[7],
            energies[3],
            energies[6],
            energies[8],
            energies[9]
        )
        set_temperature!(
            unit.input_interfaces[unit.m_heat_in],
            lowest(energies[3]),
            nothing
        )
        set_temperature!(
            unit.output_interfaces[unit.m_heat_out],
            nothing,
            highest(energies[6])
        )
    end
end

function process(unit::HeatPump, sim_params::Dict{String,Any})
    success, energies = calculate_energies(unit, sim_params)

    if !success
        set_max_energies!(unit, 0.0, 0.0, 0.0)
        return
    end

    el_in = sum(energies[1]; init=0.0)
    heat_in = sum(energies[2]; init=0.0)
    heat_out = sum(energies[5]; init=0.0)

    if heat_out < sim_params["epsilon"]
        set_max_energies!(unit, 0.0, 0.0, 0.0)
        return
    end

    # calculate average cop of current time step
    if el_in > sim_params["epsilon"]
        unit.cop = heat_out / el_in
    end

    sub!(unit.input_interfaces[unit.m_el_in], el_in)
    sub!(unit.input_interfaces[unit.m_heat_in], heat_in, lowest(energies[3]))
    add!(unit.output_interfaces[unit.m_heat_out], heat_out, highest(energies[6]))
end

function output_values(unit::HeatPump)::Vector{String}
    return [string(unit.m_el_in)*" IN", 
            string(unit.m_heat_in)*" IN",
            string(unit.m_heat_out)*" OUT",
            "COP", 
            "Losses"]
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
    end
    throw(KeyError(key.value_key))
end

export HeatPump