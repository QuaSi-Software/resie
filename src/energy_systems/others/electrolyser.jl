"""
Implementation of an electrolyser, turning electricity and water into H2, O2 and heat.

For the moment this remains a simple implementation that converts electricity into
the gases and heat (as medium m_h_w_ht1) at a defined ratio (default 1:0.6:0.4). Has a
minimum run time taken into consideration in its control behaviour and a minimum power
fraction in its processing. The power is considered the maximum amount of electricity that
the electrolyser can consume.
"""
mutable struct Electrolyser <: Component
    uac::String
    controller::Controller
    sys_function::SystemFunction

    input_interfaces::InterfaceMap
    output_interfaces::InterfaceMap

    m_el_in::Symbol
    m_heat_out::Symbol
    m_h2_out::Symbol
    m_o2_out::Symbol

    power::Float64
    heat_fraction::Float64
    min_power_fraction::Float64
    min_run_time::UInt
    output_temperature::Temperature

    function Electrolyser(uac::String, config::Dict{String,Any})
        m_el_in = Symbol(default(config, "m_el_in", "m_e_ac_230v"))
        m_heat_out = Symbol(default(config, "m_heat_out", "m_h_w_lt1"))
        m_h2_out = Symbol(default(config, "m_h2_out", "m_c_g_h2"))
        m_o2_out = Symbol(default(config, "m_o2_out", "m_c_g_o2"))
        register_media([m_el_in, m_heat_out, m_h2_out, m_o2_out])

        return new(
            uac, # uac
            controller_for_strategy( # controller
                config["strategy"]["name"], config["strategy"]
            ),
            sf_transformer, # sys_function
            InterfaceMap( # input_interfaces
                m_el_in => nothing
            ),
            InterfaceMap( # output_interfaces
                m_heat_out => nothing,
                m_h2_out => nothing,
                m_o2_out => nothing,
            ),
            m_el_in,
            m_heat_out,
            m_h2_out,
            m_o2_out,
            config["power"], # power
            default(config, "heat_fraction", 0.4),
            default(config, "min_power_fraction", 0.2),
            default(config, "min_run_time", 3600),
            default(config, "output_temperature", 55.0),
        )
    end
end

function control(
    unit::Electrolyser,
    components::Grouping,
    parameters::Dict{String,Any}
)
    move_state(unit, components, parameters)
    if unit.output_interfaces[unit.m_heat_out].temperature === nothing
        unit.output_interfaces[unit.m_heat_out].temperature = unit.output_temperature
    end
end

function set_max_energies!(
    unit::Electrolyser, el_in::Float64, heat_out::Float64,
    h2_out::Float64, o2_out::Float64
)
    set_max_energy!(unit.input_interfaces[unit.m_el_in], el_in)
    set_max_energy!(unit.output_interfaces[unit.m_heat_out], heat_out)
    set_max_energy!(unit.output_interfaces[unit.m_h2_out], h2_out)
    set_max_energy!(unit.output_interfaces[unit.m_o2_out], o2_out)
end

function check_el_in(
    unit::Electrolyser,
    parameters::Dict{String,Any}
)
    if unit.controller.parameter["m_el_in"] == true
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
            ) <= parameters["epsilon"]
                return (0.0, 0.0)
            end
            return (potential_energy_el, potential_storage_el)
        end
    else
        return (Inf, Inf)
    end
end

function check_heat_out(
    unit::Electrolyser,
    parameters::Dict{String,Any}
)
    if unit.controller.parameter["m_heat_out"] == true
        exchanges = balance_on(
            unit.output_interfaces[unit.m_heat_out],
            unit.output_interfaces[unit.m_heat_out].target
        )
        return (
            [e.balance + e.energy_potential for e in exchanges],
            [e.storage_potential for e in exchanges],
            temperature_all(exchanges)
        )
    else
        return (-Inf, -Inf, unit.output_interfaces[unit.m_heat_out].temperature)
    end
end

function check_h2_out(
    unit::Electrolyser,
    parameters::Dict{String,Any}
)
    if unit.controller.parameter["m_h2_out"] == true
        exchanges = balance_on(
            unit.output_interfaces[unit.m_h2_out],
            unit.output_interfaces[unit.m_h2_out].target
        )
        potential_energy_h2 = balance(exchanges) + energy_potential(exchanges)
        potential_storage_h2 = storage_potential(exchanges)
        if (
            unit.controller.parameter["load_storages"]
            ? potential_energy_h2 + potential_storage_h2
            : potential_energy_h2
        ) >= -parameters["epsilon"]
            return (0.0, 0.0)
        end
        return (potential_energy_h2, potential_storage_h2)
    else
        return (-Inf, -Inf)
    end
end

function check_o2_out(
    unit::Electrolyser,
    parameters::Dict{String,Any}
)
    if unit.controller.parameter["m_o2_out"] == true
        exchanges = balance_on(
            unit.output_interfaces[unit.m_o2_out],
            unit.output_interfaces[unit.m_o2_out].target
        )
        potential_energy_o2 = balance(exchanges) + energy_potential(exchanges)
        potential_storage_o2 = storage_potential(exchanges)
        if (
            unit.controller.parameter["load_storages"]
            ? potential_energy_o2 + potential_storage_o2
            : potential_energy_o2
        ) >= -parameters["epsilon"]
            return (0.0, 0.0)
        end
        return (potential_energy_o2, potential_storage_o2)
    else
        return (-Inf, -Inf)
    end
end

function calculate_energies(
    unit::Electrolyser,
    parameters::Dict{String,Any}
)
    # get usage fraction of external profile (normalized from 0 to 1)
    max_usage_fraction = (
        unit.controller.parameter["operation_profile_path"] === nothing
        ? 1.0
        : value_at_time(unit.controller.parameter["operation_profile"], parameters["time"])
    )
    if max_usage_fraction <= 0.0
        return (false, nothing, nothing, nothing, nothing, nothing)
    end

    # get potentials from inputs/outputs. only the heat output is calculated as vector,
    # the electricity input and h2/o2 outputs are calculated as scalars
    potential_energy_el, potential_storage_el = check_el_in(unit, parameters)
    potentials_energy_heat_out,
        potentials_storage_heat_out,
        out_temps = check_heat_out(unit, parameters)
    potential_energy_h2_out, potential_storage_h2_out = check_h2_out(unit, parameters)
    potential_energy_o2_out, potential_storage_o2_out = check_o2_out(unit, parameters)

    available_el_in = unit.controller.parameter["unload_storages"] ?
                      potential_energy_el + potential_storage_el :
                      potential_energy_el

    available_h2_out = unit.controller.parameter["load_storages"] ?
                         potential_energy_h2_out + potential_storage_h2_out :
                         potential_energy_h2_out

    available_o2_out = unit.controller.parameter["load_storages"] ?
                        potential_energy_o2_out + potential_storage_o2_out :
                        potential_energy_o2_out

    # in the following we want to work with positive values as it is easier
    available_h2_out = abs(available_h2_out)
    available_o2_out = abs(available_o2_out)

    # limit electricity input to design power
    available_el_in = min(available_el_in, watt_to_wh(unit.power))

    # shortcut if we're limited by electricity input or h2/o2 output
    if (
        available_el_in <= parameters["epsilon"]
        || available_h2_out <= parameters["epsilon"]
        || available_o2_out <= parameters["epsilon"]
    )
        return (false, nothing, nothing, nothing, nothing, nothing)
    end

    layers_el_in = []
    layers_heat_out = []
    layers_heat_out_temperature = []
    layers_h2_out = []
    layers_o2_out = []

    for (idx_layer, pot_heat_out) in pairs(potentials_energy_heat_out)
        # if the entire amount of one of the limiting inputs/outputs was used up, skip
        # through the rest of the heat output layers
        if (
            available_el_in <= parameters["epsilon"]
            || available_h2_out <= parameters["epsilon"]
            || available_o2_out <= parameters["epsilon"]
        )
            continue
        end

        # check if it is an "empty" layer, usually from other inputs on a bus, which are
        # included for balance calculations but cannot take in energy from the electorlyser
        if (
            unit.controller.parameter["load_storages"]
            ? pot_heat_out + potentials_storage_heat_out[idx_layer]
            : pot_heat_out
        ) <= parameters["epsilon"]
            continue
        end

        # skip layer if requested temperature is higher than the output temperature
        if out_temps[idx_layer] > unit.output_temperature
            continue
        end

        # energies for current layer with potential (+storage) heat out as basis
        used_heat_out = unit.controller.parameter["load_storages"] ?
            pot_heat_out + potentials_storage_heat_out[idx_layer] :
            pot_heat_out
        used_el_in = used_heat_out / unit.heat_fraction
        used_h2_out = used_el_in * (1.0 - unit.heat_fraction)
        used_o2_out = used_h2_out * 0.5

        # check electricity in as limiter
        if used_el_in > available_el_in
            used_el_in = available_el_in
            used_heat_out = used_el_in * unit.heat_fraction
            used_h2_out = used_el_in * (1.0 - unit.heat_fraction)
            used_o2_out = used_h2_out * 0.5
        end

        # check h2 out as limiter
        if used_h2_out > available_h2_out
            used_h2_out = available_h2_out
            used_o2_out = used_h2_out * 0.5
            used_el_in = used_h2_out / (1.0 - unit.heat_fraction)
            used_heat_out = used_el_in * unit.heat_fraction
        end

        # check o2 out as limiter
        if used_o2_out > available_o2_out
            used_o2_out = available_o2_out
            used_h2_out = used_ho_out * 2.0
            used_el_in = used_h2_out / (1.0 - unit.heat_fraction)
            used_heat_out = used_el_in * unit.heat_fraction
        end

        # check if usage fraction went over the maximum, in which case the last layer added
        # can't be fully utilised and is added with the remaining fraction to the max
        old_usage_fraction = (sum(layers_el_in; init=0.0)) / watt_to_wh(unit.power)
        new_usage_fraction = (sum(layers_el_in; init=0.0) + used_el_in) /
                             watt_to_wh(unit.power)
        if new_usage_fraction > max_usage_fraction
            used_el_in *= (max_usage_fraction - old_usage_fraction)
            used_heat_out = used_el_in * unit.heat_fraction
            used_h2_out = used_el_in * (1.0 - unit.heat_fraction)
            used_o2_out = used_h2_out * 0.5
        end

        # finally all checks done, we add the layer and update remaining energies
        push!(layers_el_in, used_el_in)
        push!(layers_heat_out, used_heat_out)
        push!(layers_heat_out_temperature, out_temps[idx_layer])
        push!(layers_h2_out, used_h2_out)
        push!(layers_o2_out, used_o2_out)
        available_el_in -= used_el_in
        available_h2_out -= used_h2_out
        available_o2_out -= used_o2_out
    end

    # if all chosen heat layers combined are not enough to induce enough electricity demand
    # to meet minimum power fraction, the electorlyser doesn't run at all
    usage_fraction = (sum(layers_el_in; init=0.0)) / watt_to_wh(unit.power)
    if usage_fraction < unit.min_power_fraction
        return (false, nothing, nothing, nothing, nothing, nothing)
    end

    return (
        true,
        layers_el_in,
        layers_heat_out,
        layers_heat_out_temperature,
        layers_h2_out,
        layers_o2_out
    )
end

function potential(
    unit::Electrolyser,
    parameters::Dict{String,Any}
)
    energies = calculate_energies(unit, parameters)

    if !energies[1]
        set_max_energies!(unit, 0.0, 0.0, 0.0, 0.0)
    else
        set_max_energies!(
            unit,
            -sum(energies[2]; init=0.0),
            sum(energies[3]; init=0.0),
            sum(energies[5]; init=0.0),
            sum(energies[6]; init=0.0)
        )
    end
end

function process(unit::Electrolyser, parameters::Dict{String,Any})
    energies = calculate_energies(unit, parameters)

    if !energies[1]
        return
    end

    el_in = sum(energies[2]; init=0.0)
    heat_out = sum(energies[3]; init=0.0)
    h2_out = sum(energies[5]; init=0.0)
    o2_out = sum(energies[6]; init=0.0)

    if el_in < parameters["epsilon"]
        return
    end

    # calculate mixed temperature for heat output, as interfaces do not support vectorized
    # energy balances (yet)
    mixed_temperature = 0.0
    for (layer_idx, layer_heat_out) in pairs(energies[3])
        mixed_temperature += energies[4][layer_idx] * (layer_heat_out / heat_out)
    end

    sub!(unit.input_interfaces[unit.m_el_in], el_in)
    add!(unit.output_interfaces[unit.m_heat_out], heat_out, mixed_temperature)
    add!(unit.output_interfaces[unit.m_h2_out], h2_out)
    add!(unit.output_interfaces[unit.m_o2_out], o2_out)
end

export Electrolyser