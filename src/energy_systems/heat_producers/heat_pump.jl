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
    unit.output_interfaces[unit.m_heat_out].temperature = highest_temperature(unit.output_temperature, unit.output_interfaces[unit.m_heat_out].temperature)
    unit.input_interfaces[unit.m_heat_in].temperature = highest_temperature(unit.input_temperature, unit.input_interfaces[unit.m_heat_in].temperature)

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
            &&
            unit.input_interfaces[unit.m_el_in].max_energy === nothing
        )
            return (Inf, Inf)
        else
            exchange = balance_on(
                unit.input_interfaces[unit.m_el_in],
                unit.input_interfaces[unit.m_el_in].source
            )
            potential_energy_el = exchange.balance == 0.0 ? collect_all_energy_potentials_of_interface_balance(exchange.energy) : [exchange.balance]
            potential_storage_el = collect_all_storage_potentials_of_interface_balance(exchange.energy)
            if (unit.controller.parameter["unload_storages"] ? sum(potential_energy_el) + sum(potential_storage_el) : sum(potential_energy_el)) <= parameters["epsilon"]
                return (nothing, nothing)
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
            &&
            unit.input_interfaces[unit.m_heat_in].max_energy === nothing
        )
            return (Inf, Inf, unit.input_interfaces[unit.m_heat_in].temperature)
        else
            exchange = balance_on(
                unit.input_interfaces[unit.m_heat_in],
                unit.input_interfaces[unit.m_heat_in].source
            )
            potential_energy_heat_in = exchange.balance == 0.0 ? collect_all_energy_potentials_of_interface_balance(exchange.energy) : [exchange.balance]
            potential_storage_heat_in = collect_all_storage_potentials_of_interface_balance(exchange.energy)
            temperatures = collect_all_temperatures_of_interface_balance(exchange.energy)
            if (unit.controller.parameter["unload_storages"] ? sum(potential_energy_heat_in) + sum(potential_storage_heat_in) : sum(potential_energy_heat_in)) <= parameters["epsilon"]
                return (nothing, nothing, temperatures)
            end
            return (potential_energy_heat_in, potential_storage_heat_in, temperatures)
        end
    else
        return (Inf, Inf, unit.input_interfaces[unit.m_heat_in].temperature)  #ToDo not sure if that is correct here. Where is this written to the interface?
    end

end

function check_heat_out(
    unit::HeatPump,
    parameters::Dict{String,Any}
)
    if unit.controller.parameter["m_heat_out"] == true
        exchange = balance_on(
            unit.output_interfaces[unit.m_heat_out],
            unit.output_interfaces[unit.m_heat_out].target
        )
        potential_energy_heat_out = exchange.balance == 0.0 ? collect_all_energy_potentials_of_interface_balance(exchange.energy) : [exchange.balance]
        potential_storage_heat_out = collect_all_storage_potentials_of_interface_balance(exchange.energy)
        temperatures = collect_all_temperatures_of_interface_balance(exchange.energy)
        if (unit.controller.parameter["load_storages"] ? sum(potential_energy_heat_out) + sum(potential_storage_heat_out) : sum(potential_energy_heat_out)) >= -parameters["epsilon"]
            return (nothing, nothing, temperatures)
        end
        return (potential_energy_heat_out, potential_storage_heat_out, temperatures)
    else
        return (-Inf, -Inf, unit.output_interfaces[unit.m_heat_out].temperature) #ToDo not sure if that is correct here. Where is this written to the interface?
    end
end

function calculate_energies(
    unit::HeatPump,
    parameters::Dict{String,Any},
    temperatures::Tuple{Vector{Float64}, Vector{Float64}},
    potentials::Vector{Vector{Float64}}
)
    # rename variable for better readability
    potential_energy_el = sum(potentials[1])        # no differenciation of power sources for now
    potential_storage_el = sum(potentials[2])       # no differenciation of power sources for now
    potential_energy_heat_in = potentials[3]
    potential_storage_heat_in = potentials[4]
    potential_energy_heat_out = sum(potentials[5])  # no differenciation of heat output layers for now
    potential_storage_heat_out = sum(potentials[6]) # no differenciation of heat output layers for now

    temperatures_in = temperatures[1]
    temperatures_out = maximum(temperatures[2])     # no differenciation of heat output layers for now, maximum temperature is assumed

    n_heat_in = length(potential_energy_heat_in)    # number of heat layers in input

    # get usage fraction of external profile (normalized from 0 to 1)
    usage_fraction_operation_profile = unit.controller.parameter["operation_profile_path"] === nothing ? 1.0 : value_at_time(unit.controller.parameter["operation_profile"], parameters["time"])
    if usage_fraction_operation_profile <= 0.0
        return (false, nothing, nothing, nothing)
    end

    # a fixed COP has priority. if it's not given the dynamic cop requires temperatures
    cop = zeros(n_heat_in)
    for n=1:n_heat_in
        cop[n] = unit.fixed_cop === nothing ? dynamic_cop(temperatures_in[n], temperatures_out) : unit.fixed_cop
        if cop[n] === nothing
            throw(ArgumentError("Input and/or output temperature for heatpump $(unit.uac) is not given. Provide temperatures or fixed cop."))
        end
    end

    # extermal limits in possible in- and outputs of heat pump, not regarding the design power
    # Note: el_in and heat_out are scalars while heat_in is a vector representing the different sources
    el_in_potential = potential_energy_el + (unit.controller.parameter["unload_storages"] ? potential_storage_el : 0)
    heat_in_potential = potential_energy_heat_in + (unit.controller.parameter["unload_storages"] ? potential_storage_heat_in : 0)
    heat_out_potential = -(potential_energy_heat_out + (unit.controller.parameter["load_storages"] ? potential_storage_heat_out : 0))
    
    # limit heat_out by design power of heat pump and extermal usage_factor from profile
    heat_out_potential = min(heat_out_potential, usage_fraction_operation_profile * watt_to_wh(unit.power) )

    # calculate heat_out and el_in for given heat_in_potential vector while respecting the limits by heat_in_potential 
    heat_out_max = cop ./ (cop .- 1) .* heat_in_potential

    # check for limits by heat_out_potential 
    for i in 1:n_heat_in
        if heat_out_potential >= heat_out_max[i]
            heat_out_potential -= heat_out_max[i]
        else
            heat_out_max[i] = heat_out_potential
            heat_out_potential = 0.0
        end
    end
    # # faster but harder to read
    # for i in 1:n_heat_in
    #     available_heat_out = min(heat_out_potential, heat_out_max[i])
    #     heat_out_potential -= available_heat_out
    #     heat_out_max[i] = available_heat_out
    # end

    # calculate el_in_max
    heat_in_max = heat_out_max .* (cop .- 1) ./ cop
    el_in_max = heat_out_max .- heat_in_max

    # check for limits by el_in_potential
    for i in 1:n_heat_in
        if el_in_potential >= el_in_max[i]
            el_in_potential -= el_in_max[i]
        else
            el_in_max[i] = el_in_potential
            el_in_potential = 0.0
        end
    end

    # recalculate energies to match all limits
    heat_out_max = el_in_max .* cop
    heat_in_max = heat_out_max .- el_in_max

    # calculate sum_abs_change
    heat_out = sum(heat_out_max)
    heat_in = sum(heat_in_max)
    el_in = sum(el_in_max)

    # calculate overall cop
    unit.cop = heat_out / (heat_out-heat_in)

    # calculate overall input temperature
    input_temperature = sum(heat_in_max .* temperatures_in) / heat_in

    # calculate usage fraction related to he heat output and check for limits
    usage_fraction = heat_out / watt_to_wh(unit.power)
    if usage_fraction < unit.min_power_fraction
        return (false, nothing, nothing, nothing)
    end

    return (
        true,
        el_in,
        heat_in,
        heat_out
    )
end

function potential(
    unit::HeatPump,
    parameters::Dict{String,Any}
)
    potential_energy_el, potential_storage_el = check_el_in(unit, parameters)
    if potential_energy_el === nothing && potential_storage_el === nothing
        set_max_energies!(unit, 0.0, 0.0, 0.0)
        return
    end

    potential_energy_heat_in, potential_storage_heat_in, in_temp = check_heat_in(unit, parameters)
    if potential_energy_heat_in === nothing && potential_storage_heat_in === nothing
        set_max_energies!(unit, 0.0, 0.0, 0.0)
        return
    end

    potential_energy_heat_out, potential_storage_heat_out, out_temp = check_heat_out(unit, parameters)
    if potential_energy_heat_out === nothing && potential_storage_heat_out === nothing
        set_max_energies!(unit, 0.0, 0.0, 0.0)
        return
    end

    # check if temperature has already been read from input and output interface, if not
    # check by balance_on
    if in_temp === nothing
        exchange = balance_on(
            unit.input_interfaces[unit.m_heat_in],
            unit
        )
        in_temp = exchange.temperature
    end

    if out_temp === nothing
        exchange = balance_on(
            unit.output_interfaces[unit.m_heat_out],
            unit.output_interfaces[unit.m_heat_out].target
        )
        out_temp = exchange.temperature
    end

    # quit if requested temperature is higher than optionally given output_temperature
    if unit.output_temperature !== nothing && out_temp > unit.output_temperature
        set_max_energies!(unit, 0.0, 0.0, 0.0)
        return
    end
        
    # quit if available temperature is lower than optionally given input_temperature
    if unit.input_temperature !== nothing && in_temp < unit.input_temperature
        set_max_energies!(unit, 0.0, 0.0, 0.0)
        return
    end

    if out_temp === nothing && unit.fixed_cop === nothing
        throw(ArgumentError("Heat pump $(unit.uac) has no requested output temperature. Please provide a temperature."))
    elseif in_temp === nothing && unit.fixed_cop === nothing
        throw(ArgumentError("Heat pump $(unit.uac) has no requested input temperature. Please provide a temperature."))
    end

    energies = calculate_energies(
        unit, parameters,
        (in_temp, out_temp),
        [
            potential_energy_el, potential_storage_el,
            potential_energy_heat_in, potential_storage_heat_in,
            potential_energy_heat_out, potential_storage_heat_out
        ]
    )

    if !energies[1]
        set_max_energies!(unit, 0.0, 0.0, 0.0)
    else
        set_max_energies!(unit, energies[2], energies[3], energies[4])
    end
end

function process(unit::HeatPump, parameters::Dict{String,Any})
    potential_energy_el, potential_storage_el = check_el_in(unit, parameters)
    if potential_energy_el === nothing && potential_storage_el === nothing
        set_max_energies!(unit, 0.0, 0.0, 0.0)
        return
    end

    potential_energy_heat_in, potential_storage_heat_in, in_temp = check_heat_in(unit, parameters)
    if potential_energy_heat_in === nothing && potential_storage_heat_in === nothing
        set_max_energies!(unit, 0.0, 0.0, 0.0)
        return
    end

    potential_energy_heat_out, potential_storage_heat_out, out_temp = check_heat_out(unit, parameters)
    if potential_energy_heat_out === nothing && potential_storage_heat_out === nothing
        set_max_energies!(unit, 0.0, 0.0, 0.0)
        return
    end

    # check if temperature has already been read from input and output interface, if not
    # check by balance_on
    if in_temp === nothing
        exchange = balance_on(
            unit.input_interfaces[unit.m_heat_in],
            unit
        )
        in_temp = exchange.temperature
    end

    if out_temp === nothing
        exchange = balance_on(
            unit.output_interfaces[unit.m_heat_out],
            unit.output_interfaces[unit.m_heat_out].target
        )
        out_temp = exchange.temperature
    end

    # quit if requested temperature is higher than optionally given output_temperature
    if unit.output_temperature !== nothing && out_temp > unit.output_temperature
        set_max_energies!(unit, 0.0, 0.0, 0.0)
        return
    end
        
    # quit if available temperature is higher than optionally given input_temperature
    if unit.input_temperature !== nothing && in_temp < unit.input_temperature
        set_max_energies!(unit, 0.0, 0.0, 0.0)
        return
    end

    if out_temp === nothing && unit.fixed_cop === nothing
        throw(ArgumentError("Heat pump $(unit.uac) has no requested output temperature. Please provide a temperature."))
    elseif in_temp === nothing && unit.fixed_cop === nothing
        throw(ArgumentError("Heat pump $(unit.uac) has no requested input temperature. Please provide a temperature."))
    end

    energies = calculate_energies(
        unit, parameters,
        (in_temp, out_temp),
        [
            potential_energy_el, potential_storage_el,
            potential_energy_heat_in, potential_storage_heat_in,
            potential_energy_heat_out, potential_storage_heat_out
        ]
    )

    if energies[1]
        sub!(unit.input_interfaces[unit.m_el_in], energies[2])
        sub!(unit.input_interfaces[unit.m_heat_in], energies[3]) # revise! Should here the temperature (vector) be returned? ToDo 
        add!(unit.output_interfaces[unit.m_heat_out], energies[4], maximum(out_temp))  # revise! This is not valid for cooling? ToDo
    end
end

export HeatPump