"""
Implementation of an electrolyser, turning electricity and water into H2, O2 and heat.

For the moment this remains a simple implementation that converts electricity into
the gases and heat (as medium m_h_w_ht1) at a defined ratio of 1:0.6:0.4. Has a minimum
run time of 3600s taken into consideration in its control behaviour and a minimum power
fraction of 20%. The power is considered the maximum amount of electricity that the
electrolyser can consume.

At the moment there is no operation strategy is implemented and the production of the
electrolyser is controlled by the demand it is linked to requires.
"""
mutable struct Electrolyser <: ControlledSystem
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
            &&
            unit.input_interfaces[unit.m_el_in].max_energy === nothing
        )
            return (Inf, Inf)
        else
            exchange = balance_on(
                unit.input_interfaces[unit.m_el_in],
                unit.input_interfaces[unit.m_el_in].source
            )
            potential_energy_el_in = exchange.balance + exchange.energy_potential
            potential_storage_el_in = exchange.storage_potential
            if (unit.controller.parameter["unload_storages"] ? potential_energy_el_in + potential_storage_el_in : potential_energy_el_in) <= parameters["epsilon"]
                return (nothing, nothing)
            end
            return (potential_energy_el_in, potential_storage_el_in)
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
        exchange = balance_on(
            unit.output_interfaces[unit.m_heat_out],
            unit.output_interfaces[unit.m_heat_out].target
        )
        potential_energy_heat_out = exchange.balance + exchange.energy_potential
        potential_storage_heat_out = exchange.storage_potential
        if (unit.controller.parameter["load_storages"] ? potential_energy_heat_out + potential_storage_heat_out : potential_energy_heat_out) >= -parameters["epsilon"]
            return (nothing, nothing)
        end
        return (potential_energy_heat_out, potential_storage_heat_out)
    else
        return (-Inf, -Inf)
    end
end

function check_h2_out(
    unit::Electrolyser,
    parameters::Dict{String,Any}
)
    if unit.controller.parameter["m_h2_out"] == true
        exchange = balance_on(
            unit.output_interfaces[unit.m_h2_out],
            unit.output_interfaces[unit.m_h2_out].target
        )
        potential_energy_h2_out = exchange.balance + exchange.energy_potential
        potential_storage_h2_out = exchange.storage_potential
        if (unit.controller.parameter["load_storages"] ? potential_energy_h2_out + potential_storage_h2_out : potential_energy_h2_out) >= -parameters["epsilon"]
            return (nothing, nothing)
        end
        return (potential_energy_h2_out, potential_storage_h2_out)
    else
        return (-Inf, -Inf)
    end
end

function check_o2_out(
    unit::Electrolyser,
    parameters::Dict{String,Any}
)
    if unit.controller.parameter["m_o2_out"] == true
        exchange = balance_on(
            unit.output_interfaces[unit.m_o2_out],
            unit.output_interfaces[unit.m_o2_out].target
        )
        potential_energy_o2_out = exchange.balance + exchange.energy_potential
        potential_storage_o2_out = exchange.storage_potential
        if (unit.controller.parameter["load_storages"] ? potential_energy_o2_out + potential_storage_o2_out : potential_energy_o2_out) >= -parameters["epsilon"]
            return (nothing, nothing)
        end
        return (potential_energy_o2_out, potential_storage_o2_out)
    else
        return (-Inf, -Inf)
    end
end

function calculate_energies(
    unit::Electrolyser,
    parameters::Dict{String,Any},
    watt_to_wh::Function,
    potentials::Vector{Float64}
)
    potential_energy_el_in = potentials[1]
    potential_storage_el_in = potentials[2]
    potential_energy_heat_out = potentials[3]
    potential_storage_heat_out = potentials[4]
    potential_energy_h2_out = potentials[5]
    potential_storage_h2_out = potentials[6]
    potential_energy_o2_out = potentials[7]
    potential_storage_o2_out = potentials[8]

    max_produce_heat = watt_to_wh(unit.power * unit.heat_fraction)
    max_produce_h2 = watt_to_wh(unit.power * (1.0 - unit.heat_fraction))
    # @TODO: handle O2 calculation if it ever becomes relevant. for now use molar ratio
    max_produce_o2 = 0.5 * max_produce_h2
    max_consume_el = watt_to_wh(unit.power)

    # get usage fraction of external profile (normalized from 0 to 1)
    usage_fraction_operation_profile = unit.controller.parameter["operation_profile_path"] === nothing ? 1.0 : value_at_time(unit.controller.parameter["operation_profile"], parameters["time"])
    if usage_fraction_operation_profile <= 0.0
        return (false, nothing, nothing, nothing, nothing)
    end

    # all three standard operating strategies behave the same, but it is better to be
    # explicit about the behaviour rather than grouping all together
    if unit.controller.strategy == "storage_driven" && unit.controller.state_machine.state == 2
        usage_fraction_el = +(unit.controller.parameter["unload_storages"] ? potential_energy_el_in + potential_storage_el_in : potential_energy_el_in) / max_consume_el
        usage_fraction_heat = -(unit.controller.parameter["load_storages"] ? potential_energy_heat_out + potential_storage_heat_out : potential_energy_heat_out) / max_produce_heat
        usage_fraction_h2 = -(unit.controller.parameter["load_storages"] ? potential_energy_h2_out + potential_storage_h2_out : potential_energy_h2_out) / max_produce_h2
        usage_fraction_o2 = -(unit.controller.parameter["load_storages"] ? potential_energy_o2_out + potential_storage_o2_out : potential_energy_o2_out) / max_produce_o2

    elseif unit.controller.strategy == "storage_driven"
        return (false, nothing, nothing, nothing, nothing)

    elseif unit.controller.strategy == "supply_driven"
        usage_fraction_el = +(unit.controller.parameter["unload_storages"] ? potential_energy_el_in + potential_storage_el_in : potential_energy_el_in) / max_consume_el
        usage_fraction_heat = -(unit.controller.parameter["load_storages"] ? potential_energy_heat_out + potential_storage_heat_out : potential_energy_heat_out) / max_produce_heat
        usage_fraction_h2 = -(unit.controller.parameter["load_storages"] ? potential_energy_h2_out + potential_storage_h2_out : potential_energy_h2_out) / max_produce_h2
        usage_fraction_o2 = -(unit.controller.parameter["load_storages"] ? potential_energy_o2_out + potential_storage_o2_out : potential_energy_o2_out) / max_produce_o2

    elseif unit.controller.strategy == "demand_driven"
        usage_fraction_el = +(unit.controller.parameter["unload_storages"] ? potential_energy_el_in + potential_storage_el_in : potential_energy_el_in) / max_consume_el
        usage_fraction_heat = -(unit.controller.parameter["load_storages"] ? potential_energy_heat_out + potential_storage_heat_out : potential_energy_heat_out) / max_produce_heat
        usage_fraction_h2 = -(unit.controller.parameter["load_storages"] ? potential_energy_h2_out + potential_storage_h2_out : potential_energy_h2_out) / max_produce_h2
        usage_fraction_o2 = -(unit.controller.parameter["load_storages"] ? potential_energy_o2_out + potential_storage_o2_out : potential_energy_o2_out) / max_produce_o2
    end

    # limit actual usage by limits of inputs, outputs and profile
    usage_fraction = min(
        1.0,
        usage_fraction_el,
        usage_fraction_heat,
        usage_fraction_h2,
        usage_fraction_o2,
        usage_fraction_operation_profile
    )

    if usage_fraction < unit.min_power_fraction
        return (false, nothing, nothing, nothing, nothing)
    end

    return (
        true,
        max_consume_el * usage_fraction,
        max_produce_heat * usage_fraction,
        max_produce_h2 * usage_fraction,
        max_produce_o2 * usage_fraction,
    )
end

function potential(
    unit::Electrolyser,
    parameters::Dict{String,Any},
    watt_to_wh::Function
)
    potential_energy_el_in, potential_storage_el_in = check_el_in(unit, parameters)
    if potential_energy_el_in === nothing && potential_storage_el_in === nothing
        set_max_energies!(unit, 0.0, 0.0, 0.0, 0.0)
        return
    end

    potential_energy_heat_out, potential_storage_heat_out = check_heat_out(unit, parameters)
    if potential_energy_heat_out === nothing && potential_storage_heat_out === nothing
        set_max_energies!(unit, 0.0, 0.0, 0.0, 0.0)
        return
    end

    potential_energy_h2_out, potential_storage_h2_out = check_h2_out(unit, parameters)
    if potential_energy_h2_out === nothing && potential_storage_h2_out === nothing
        set_max_energies!(unit, 0.0, 0.0, 0.0, 0.0)
        return
    end

    potential_energy_o2_out, potential_storage_o2_out = check_o2_out(unit, parameters)
    if potential_energy_o2_out === nothing && potential_storage_o2_out === nothing
        set_max_energies!(unit, 0.0, 0.0, 0.0, 0.0)
        return
    end

    energies = calculate_energies(
        unit, parameters, watt_to_wh,
        [
            potential_energy_el_in, potential_storage_el_in,
            potential_energy_heat_out, potential_storage_heat_out,
            potential_energy_h2_out, potential_storage_h2_out,
            potential_energy_o2_out, potential_storage_o2_out
        ]
    )

    if !energies[1]
        set_max_energies!(unit, 0.0, 0.0, 0.0, 0.0)
    else
        set_max_energies!(unit, energies[2],energies[3], energies[4], energies[5])
    end
end

function produce(unit::Electrolyser, parameters::Dict{String,Any}, watt_to_wh::Function)
    potential_energy_el_in, potential_storage_el_in = check_el_in(unit, parameters)
    if potential_energy_el_in === nothing && potential_storage_el_in === nothing
        set_max_energies!(unit, 0.0, 0.0, 0.0, 0.0)
        return
    end

    potential_energy_heat_out, potential_storage_heat_out = check_heat_out(unit, parameters)
    if potential_energy_heat_out === nothing && potential_storage_heat_out === nothing
        set_max_energies!(unit, 0.0, 0.0, 0.0, 0.0)
        return
    end

    potential_energy_h2_out, potential_storage_h2_out = check_h2_out(unit, parameters)
    if potential_energy_h2_out === nothing && potential_storage_h2_out === nothing
        set_max_energies!(unit, 0.0, 0.0, 0.0, 0.0)
        return
    end

    potential_energy_o2_out, potential_storage_o2_out = check_o2_out(unit, parameters)
    if potential_energy_o2_out === nothing && potential_storage_o2_out === nothing
        set_max_energies!(unit, 0.0, 0.0, 0.0, 0.0)
        return
    end

    energies = calculate_energies(
        unit, parameters, watt_to_wh,
        [
            potential_energy_el_in, potential_storage_el_in,
            potential_energy_heat_out, potential_storage_heat_out,
            potential_energy_h2_out, potential_storage_h2_out,
            potential_energy_o2_out, potential_storage_o2_out
        ]
    )

    if energies[1]
        sub!(unit.input_interfaces[unit.m_el_in], energies[2])
        add!(unit.output_interfaces[unit.m_heat_out], energies[3], unit.output_temperature)
        add!(unit.output_interfaces[unit.m_h2_out], energies[4])
        add!(unit.output_interfaces[unit.m_o2_out], energies[5])
    end
end

export Electrolyser