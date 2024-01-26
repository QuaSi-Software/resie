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
    unit::HeatPump, el_in::Float64,
    heat_in::Float64, heat_out::Float64
)
    set_max_energy!(unit.input_interfaces[unit.m_el_in], el_in)
    set_max_energy!(unit.input_interfaces[unit.m_heat_in], heat_in)
    set_max_energy!(unit.output_interfaces[unit.m_heat_out], heat_out)
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
        if (
            unit.input_interfaces[unit.m_el_in].source.sys_function == sf_transformer
            && unit.input_interfaces[unit.m_el_in].max_energy === nothing
        )
            return (Inf, Inf)
        else
            exchanges = balance_on(
                unit.input_interfaces[unit.m_el_in],
                unit.input_interfaces[unit.m_el_in].source
            )
            potential_energy_el = balance(exchanges) + energy_potential(exchanges)
            potential_storage_el = storage_potential(exchanges)
            if (
                unit.controller.parameter["unload_storages"]
                ? potential_energy_el + potential_storage_el
                : potential_energy_el
            ) <= sim_params["epsilon"]
                return (0.0, 0.0)
            end
            return (potential_energy_el, potential_storage_el)
        end
    else
        return (Inf, Inf)
    end
end

function check_heat_in(
    unit::HeatPump,
    sim_params::Dict{String,Any}
)
    if unit.controller.parameter["consider_m_heat_in"] == true
        if (
            unit.input_interfaces[unit.m_heat_in].source.sys_function == sf_transformer
            && unit.input_interfaces[unit.m_heat_in].max_energy === nothing
        )
            return ([Inf], 
                    [Inf], 
                    [unit.input_interfaces[unit.m_heat_in].temperature_min], 
                    [unit.input_interfaces[unit.m_heat_in].temperature_max]
                )
        else
            exchanges = balance_on(
                unit.input_interfaces[unit.m_heat_in],
                unit.input_interfaces[unit.m_heat_in].source
            )
            return (
                [e.balance + e.energy_potential for e in exchanges],
                [e.storage_potential for e in exchanges],
                temp_min_all(exchanges),
                temp_max_all(exchanges)
            )
        end
    else
        return ([Inf], 
                [Inf], 
                [unit.input_interfaces[unit.m_heat_in].temperature_min], 
                [unit.input_interfaces[unit.m_heat_in].temperature_max]
            )
    end
end

function check_heat_out(
    unit::HeatPump,
    sim_params::Dict{String,Any}
)
    if unit.controller.parameter["consider_m_heat_out"] == true
        exchanges = balance_on(
            unit.output_interfaces[unit.m_heat_out],
            unit.output_interfaces[unit.m_heat_out].target
        )
        potential_energy_heat_out = balance(exchanges) + energy_potential(exchanges)
        potential_storage_heat_out = storage_potential(exchanges)
        temperature = temp_min_highest(exchanges)
        if (
            unit.controller.parameter["load_storages"]
            ? potential_energy_heat_out + potential_storage_heat_out
            : potential_energy_heat_out
        ) >= -sim_params["epsilon"]
            return (0.0, 0.0, temperature)
        end
        return (potential_energy_heat_out, potential_storage_heat_out, temperature)
    else
        return (-Inf, -Inf, nothing)
    end
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
        return (false, nothing, nothing, nothing)
    end

    # get usage fraction of external profile (normalized from 0 to 1)
    max_usage_fraction = (
        unit.controller.parameter["operation_profile_path"] === nothing
        ? 1.0
        : value_at_time(unit.controller.parameter["operation_profile"], sim_params["time"])
    )
    if max_usage_fraction <= 0.0
        return (false, nothing, nothing, nothing)
    end

    # get potentials from inputs/outputs. only the heat input is calculated as vector,
    # the electricity input and heat output are calculated as scalars
    potential_energy_el, potential_storage_el = check_el_in(unit, sim_params)
    potentials_energy_heat_in,
        potentials_storage_heat_in,
        in_temps_min,
        in_temps_max = check_heat_in(unit, sim_params)
    potential_energy_heat_out,
        potential_storage_heat_out,
        out_temp = check_heat_out(unit, sim_params)

    available_el_in = unit.controller.parameter["unload_storages"] ?
                      potential_energy_el + potential_storage_el :
                      potential_energy_el

    available_heat_out = unit.controller.parameter["load_storages"] ?
                         potential_energy_heat_out + potential_storage_heat_out :
                         potential_energy_heat_out

    # in the following we want to work with positive values as it is easier
    available_heat_out = abs(available_heat_out)

    # limit heat output to design power
    available_heat_out = min(available_heat_out, watt_to_wh(unit.power_th))

    # shortcut if we're limited by electricity input or heat output or the requested
    # temperature is higher than a set limit of the heat pump
    if (
        available_el_in <= sim_params["epsilon"]
        || available_heat_out <= sim_params["epsilon"]
        || unit.output_temperature !== nothing && out_temp > unit.output_temperature
    )
        return (false, nothing, nothing, nothing)
    end

    layers_el_in = []
    layers_heat_in = []
    layers_heat_in_temperature = []
    layers_heat_out = []

    for (idx_layer, pot_heat_in) in pairs(potentials_energy_heat_in)
        # if all electricity input or heat output was used up, skip through the rest
        # of the heat input layers
        if (
            available_el_in <= sim_params["epsilon"]
            || available_heat_out <= sim_params["epsilon"]
        )
            continue
        end

        # check if it is an "empty" layer, usually from other outputs on a bus, which are
        # included for balance calculations but cannot offer energy
        if (
            unit.controller.parameter["unload_storages"]
            ? pot_heat_in + potentials_storage_heat_in[idx_layer]
            : pot_heat_in
        ) <= sim_params["epsilon"]
            continue
        end

        # skip layer if optionally given fixed input temperature is not within
        # temperature band of the layer
        if (
            unit.input_temperature !== nothing
                && in_temps_min[idx_layer] !== nothing
                && in_temps_min[idx_layer] > unit.input_temperature
            || unit.input_temperature !== nothing
                && in_temps_max[idx_layer] !== nothing
                && in_temps_max[idx_layer] < unit.input_temperature
        )
            continue
        end
        in_temp = unit.input_temperature !== nothing ?
            unit.input_temperature :
            highest(in_temps_min[idx_layer], in_temps_max[idx_layer])

        # a constant COP has priority. if it's not given the dynamic cop requires temperatures
        cop = unit.constant_cop === nothing ?
              dynamic_cop(in_temp, out_temp) :
              unit.constant_cop
        if cop === nothing
            @error ("Input and/or output temperature for heatpump $(unit.uac) is not given. Provide temperatures or fixed cop.")
            exit()
        end

        # energies for current layer with potential (+storage) heat in as basis
        used_heat_in = unit.controller.parameter["unload_storages"] ?
                       pot_heat_in + potentials_storage_heat_in[idx_layer] :
                       pot_heat_in
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

        # check if usage fraction went over the maximum, in which case the last layer added
        # can't be fully utilised and is added with the remaining fraction to the max
        old_usage_fraction = (sum(layers_heat_out; init=0.0)) / watt_to_wh(unit.power_th)
        new_usage_fraction = (sum(layers_heat_out; init=0.0) + used_heat_out) /
                             watt_to_wh(unit.power_th)
        if new_usage_fraction > max_usage_fraction
            used_heat_out *= (max_usage_fraction - old_usage_fraction)
            used_el_in = used_heat_out * cop
            used_heat_in = used_el_in * (cop - 1.0)
        end

        # finally all checks done, we add the layer and update remaining energies
        push!(layers_el_in, used_el_in)
        push!(layers_heat_in, used_heat_in)
        push!(layers_heat_in_temperature, in_temp)
        push!(layers_heat_out, used_heat_out)
        available_el_in -= used_el_in
        available_heat_out -= used_heat_out
    end

    # if all chosen heat layers combined are not enough to meet minimum power fraction,
    # the heat pump doesn't run at all
    usage_fraction = (sum(layers_heat_out; init=0.0)) / watt_to_wh(unit.power_th)
    if usage_fraction < unit.min_power_fraction
        return (false, nothing, nothing, nothing)
    end

    return (
        true,
        layers_el_in,
        layers_heat_in,
        layers_heat_in_temperature,
        layers_heat_out,
        out_temp
    )
end

function potential(unit::HeatPump, sim_params::Dict{String,Any})
    energies = calculate_energies(unit, sim_params)

    if !energies[1]
        set_max_energies!(unit, 0.0, 0.0, 0.0)
    else
        set_max_energies!(
            unit,
            sum(energies[2]; init=0.0),
            sum(energies[3]; init=0.0),
            sum(energies[5]; init=0.0)
        )
    end
end

function process(unit::HeatPump, sim_params::Dict{String,Any})
    energies = calculate_energies(unit, sim_params)

    if !energies[1]
        set_max_energies!(unit, 0.0, 0.0, 0.0)
        return
    end

    el_in = sum(energies[2]; init=0.0)
    heat_in = sum(energies[3]; init=0.0)
    heat_out = sum(energies[5]; init=0.0)

    if heat_out < sim_params["epsilon"]
        set_max_energies!(unit, 0.0, 0.0, 0.0)
        return
    end

    if el_in > sim_params["epsilon"]
        unit.cop = heat_out / el_in
    end

    sub!(unit.input_interfaces[unit.m_el_in], el_in)
    sub!(unit.input_interfaces[unit.m_heat_in], heat_in, nothing)
    add!(unit.output_interfaces[unit.m_heat_out], heat_out, energies[6])
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