"""
Implementation of a heat pump component.

For the moment this remains a simple implementation that requires a low temperature heat
and electricity input and produces high temperature heat. Has a fixed coefficient of
performance (COP) of 3 and a minimum power fraction of 20%. The power parameter is
considered the maximum power of heat output the heat pump can produce.

The only currently implemented operation strategy involves checking the load of a linked
buffer tank and en-/disabling the heat pump when a threshold is reached, in addition to an
overfill shutoff condition.
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

    power::Float64
    min_power_fraction::Float64
    min_run_time::UInt
    fixed_cop::Any
    output_temperature::Temperature
    input_temperature::Temperature
    cop::Float64

    function HeatPump(uac::String, config::Dict{String,Any})
        m_el_in = Symbol(default(config, "m_el_in", "m_e_ac_230v"))
        m_heat_out = Symbol(default(config, "m_heat_out", "m_h_w_ht1"))
        m_heat_in = Symbol(default(config, "m_heat_in", "m_h_w_lt1"))
        register_media([m_el_in, m_heat_out, m_heat_in])

        return new(
            uac, # uac
            controller_for_strategy( # controller
                config["strategy"]["name"], config["strategy"]
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
            config["power"], # power
            default(config, "min_power_fraction", 0.2),
            default(config, "min_run_time", 0),
            default(config, "fixed_cop", nothing),
            default(config, "output_temperature", nothing),
            default(config, "input_temperature", nothing),
            0.0, # cop
        )
    end
end

function control(
    unit::HeatPump,
    components::Grouping,
    parameters::Dict{String,Any}
)
    move_state(unit, components, parameters)
    unit.output_interfaces[unit.m_heat_out].temperature = highest(unit.output_temperature, unit.output_interfaces[unit.m_heat_out].temperature)
    unit.input_interfaces[unit.m_heat_in].temperature = highest(unit.input_temperature, unit.input_interfaces[unit.m_heat_in].temperature)

end

function output_values(unit::HeatPump)::Vector{String}
    return ["IN", "OUT", "COP"]
end

function output_value(unit::HeatPump, key::OutputKey)::Float64
    if key.value_key == "IN"
        return calculate_energy_flow(unit.input_interfaces[key.medium])
    elseif key.value_key == "OUT"
        return calculate_energy_flow(unit.output_interfaces[key.medium])
    elseif key.value_key == "COP"
        return unit.cop
    end
    throw(KeyError(key.value_key))
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
    
    return 0.4 * (273.15 + out_temp) / (out_temp - in_temp) # Carnot-COP with 40 % efficiency
end


function check_el_in(
    unit::HeatPump,
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

function check_heat_in(
    unit::HeatPump,
    parameters::Dict{String,Any}
)
    if unit.controller.parameter["m_heat_in"] == true
        if (
            unit.input_interfaces[unit.m_heat_in].source.sys_function == sf_transformer
            && unit.input_interfaces[unit.m_heat_in].max_energy === nothing
        )
            return ([Inf], [Inf], [unit.input_interfaces[unit.m_heat_in].temperature])
        else
            exchanges = balance_on(
                unit.input_interfaces[unit.m_heat_in],
                unit.input_interfaces[unit.m_heat_in].source
            )
            return (
                [e.balance + e.energy_potential for e in exchanges],
                [e.storage_potential for e in exchanges],
                temperature_all(exchanges)
            )
        end
    else
        return ([Inf], [Inf], [unit.input_interfaces[unit.m_heat_in].temperature])
    end
end

function check_heat_out(
    unit::HeatPump,
    parameters::Dict{String,Any}
)
    if unit.controller.parameter["m_heat_out"] == true
        exchanges = balance_on(
            unit.output_interfaces[unit.m_heat_out],
            unit.output_interfaces[unit.m_heat_out].target
        )
        potential_energy_heat_out = balance(exchanges)
        potential_storage_heat_out = storage_potential(exchanges)
        temperature = temperature_highest(exchanges)
        if (
            unit.controller.parameter["load_storages"]
            ? potential_energy_heat_out + potential_storage_heat_out
            : potential_energy_heat_out
        ) >= -parameters["epsilon"]
            return (0.0, 0.0, temperature)
        end
        return (potential_energy_heat_out, potential_storage_heat_out, temperature)
    else
        return (-Inf, -Inf, unit.output_interfaces[unit.m_heat_out].temperature)
    end
end

function calculate_energies(
    unit::HeatPump,
    parameters::Dict{String,Any}
)
    # get usage fraction of external profile (normalized from 0 to 1)
    max_usage_fraction = (
        unit.controller.parameter["operation_profile_path"] === nothing
        ? 1.0
        : value_at_time(unit.controller.parameter["operation_profile"], parameters["time"])
    )
    if max_usage_fraction <= 0.0
        return (false, nothing, nothing, nothing)
    end

    # get potentials from inputs/outputs. only the heat input is calculated as vector,
    # the electricity input and heat output are calculated as scalars
    potential_energy_el, potential_storage_el = check_el_in(unit, parameters)
    potentials_energy_heat_in,
        potentials_storage_heat_in,
        in_temps = check_heat_in(unit, parameters)
    potential_energy_heat_out,
        potential_storage_heat_out,
        out_temp = check_heat_out(unit, parameters)

    available_el_in = unit.controller.parameter["unload_storages"] ?
                      potential_energy_el + potential_storage_el :
                      potential_energy_el

    available_heat_out = unit.controller.parameter["load_storages"] ?
                         potential_energy_heat_out + potential_storage_heat_out :
                         potential_energy_heat_out

    # in the following we want to work with positive values as it is easier
    available_heat_out = abs(available_heat_out)

    # limit heat output to design power
    available_heat_out = min(available_heat_out, watt_to_wh(unit.power))

    # shortcut if we're limited by electricity input or heat output or the requested
    # temperature is higher than a set limit of the heat pump
    if (
        available_el_in <= parameters["epsilon"]
        || available_heat_out <= parameters["epsilon"]
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
            available_el_in <= parameters["epsilon"]
            || available_heat_out <= parameters["epsilon"]
        )
            continue
        end

        # check if it is an "empty" layer, usually from other outputs on a bus, which are
        # included for balance calculations but cannot offer energy
        if (
            unit.controller.parameter["unload_storages"]
            ? pot_heat_in + potentials_storage_heat_in[idx_layer]
            : pot_heat_in
        ) <= parameters["epsilon"]
            continue
        end

        # skip layer if available temperature is lower than optionally given
        # minimum input_temperature
        if (
            unit.input_temperature !== nothing
            && in_temps[idx_layer] < unit.input_temperature
        )
            continue
        end

        # a fixed COP has priority. if it's not given the dynamic cop requires temperatures
        cop = unit.fixed_cop === nothing ?
              dynamic_cop(in_temps[idx_layer], out_temp) :
              unit.fixed_cop
        if cop === nothing
            throw(ArgumentError("Input and/or output temperature for heatpump $(unit.uac) is not given. Provide temperatures or fixed cop."))
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
        old_usage_fraction = (sum(layers_heat_in; init=0.0)) / watt_to_wh(unit.power)
        new_usage_fraction = (sum(layers_heat_in; init=0.0) + used_heat_out) /
                             watt_to_wh(unit.power)
        if new_usage_fraction > max_usage_fraction
            used_heat_out *= (max_usage_fraction - old_usage_fraction)
            used_el_in = used_heat_out * cop
            used_heat_in = used_el_in * (cop - 1.0)
        end

        # finally all checks done, we add the layer and update remaining energies
        push!(layers_el_in, used_el_in)
        push!(layers_heat_in, used_heat_in)
        push!(layers_heat_in_temperature, in_temps[idx_layer])
        push!(layers_heat_out, used_heat_out)
        available_el_in -= used_el_in
        available_heat_out -= used_heat_out
    end

    # if all chosen heat layers combined are not enough to meet minimum power fraction,
    # the heat pump doesn't run at all
    usage_fraction = (sum(layers_heat_out; init=0.0)) / watt_to_wh(unit.power)
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

function potential(unit::HeatPump, parameters::Dict{String,Any})
    energies = calculate_energies(unit, parameters)

    if !energies[1]
        set_max_energies!(unit, 0.0, 0.0, 0.0)
    else
        set_max_energies!(
            unit,
            -sum(energies[2]; init=0.0),
            -sum(energies[3]; init=0.0),
            sum(energies[5]; init=0.0)
        )
    end
end

function process(unit::HeatPump, parameters::Dict{String,Any})
    energies = calculate_energies(unit, parameters)

    el_in = sum(energies[2]; init=0.0)
    heat_in = sum(energies[3]; init=0.0)
    heat_out = sum(energies[5]; init=0.0)

    if heat_out < parameters["epsilon"]
        return
    end

    if el_in > parameters["epsilon"]
        unit.cop = heat_out / el_in
    end

    # calculate mixed temperature for heat input, as interfaces do not support vectorized
    # energy balances (yet)
    mixed_temperature = 0.0
    for (layer_idx, layer_heat_in) in pairs(energies[3])
        mixed_temperature += energies[4][layer_idx] * (layer_heat_in / heat_in)
    end

    if energies[1]
        sub!(unit.input_interfaces[unit.m_el_in], el_in)
        sub!(unit.input_interfaces[unit.m_heat_in], heat_in, mixed_temperature)
        add!(unit.output_interfaces[unit.m_heat_out], heat_out, energies[6])
    end
end

export HeatPump