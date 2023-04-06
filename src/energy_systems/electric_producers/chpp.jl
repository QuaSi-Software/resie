"""
Implementation of a combined-heat-power-plant (CHPP) energy system.

For the moment this remains a simple implementation that converts natural gas into
electricity and heat (as medium m_h_w_ht1) at a defined ratio of 1:0.4:0.6. Has a minimum
run time of 1800s taken into consideration in its control behaviour and a minimum power
fraction of 20%. The power is considered the maximum amount of both heat and electricity
that the CHPP can produce.

The only currently implemented operation strategy involves checking the load of a linked
buffer tank and en-/disabling the CHPP when a threshold is reached, in addition to an
overfill shutoff condition.
"""
mutable struct CHPP <: ControlledSystem
    uac::String
    controller::Controller
    sys_function::SystemFunction

    input_interfaces::InterfaceMap
    output_interfaces::InterfaceMap

    m_gas_in::Symbol
    m_heat_out::Symbol
    m_el_out::Symbol

    power::Float64
    electricity_fraction::Float64
    min_power_fraction::Float64
    min_run_time::UInt
    output_temperature::Temperature

    function CHPP(uac::String, config::Dict{String,Any})
        m_gas_in = Symbol(default(config, "m_gas_in", "m_c_g_natgas"))
        m_heat_out = Symbol(default(config, "m_heat_out", "m_h_w_ht1"))
        m_el_out = Symbol(default(config, "m_el_out", "m_e_ac_230v"))
        register_media([m_gas_in, m_heat_out, m_el_out])

        return new(
            uac, # uac
            controller_for_strategy( # controller
                config["strategy"]["name"], config["strategy"]
            ),
            sf_transformer, # sys_function
            InterfaceMap( # input_interfaces
                m_gas_in => nothing
            ),
            InterfaceMap( # output_interfaces
                m_heat_out => nothing,
                m_el_out => nothing
            ),
            m_gas_in,
            m_heat_out,
            m_el_out,
            config["power"], # power
            default(config, "electricity_fraction", 0.4),
            default(config, "min_power_fraction", 0.2),
            default(config, "min_run_time", 1800),
            default(config, "output_temperature", nothing),

        )
    end
end

function control(
    unit::CHPP,
    systems::Grouping,
    parameters::Dict{String,Any}
)
    move_state(unit, systems, parameters)
    unit.output_interfaces[unit.m_heat_out].temperature = highest_temperature(unit.output_temperature, unit.output_interfaces[unit.m_heat_out].temperature)
end

function set_max_energies!(unit::CHPP, gas_in::Float64, el_out::Float64, heat_out::Float64)
    set_max_energy!(unit.input_interfaces[unit.m_gas_in], gas_in)
    set_max_energy!(unit.output_interfaces[unit.m_el_out], el_out)
    set_max_energy!(unit.output_interfaces[unit.m_heat_out], heat_out)
end

function check_gas_in(
    unit::CHPP,
    parameters::Dict{String,Any}
)
    if unit.controller.parameter["m_gas_in"] == true
        if (
            unit.input_interfaces[unit.m_gas_in].source.sys_function == sf_transformer
            &&
            unit.input_interfaces[unit.m_gas_in].max_energy === nothing
        )
            return (Inf, Inf)
        else
            exchange = balance_on(
                unit.input_interfaces[unit.m_gas_in],
                unit.input_interfaces[unit.m_gas_in].source
            )
            potential_energy_gas = exchange.balance + exchange.energy_potential
            potential_storage_gas = exchange.storage_potential
            if (unit.controller.parameter["unload_storages"] ? potential_energy_gas + potential_storage_gas : potential_energy_gas) <= parameters["epsilon"]
                return (nothing, nothing)
            end
            return (potential_energy_gas, potential_storage_gas, exchange.temperature)
        end
    else
        return (Inf, Inf)
    end
end

function check_el_out(
    unit::CHPP,
    parameters::Dict{String,Any}
)
    if unit.controller.parameter["m_el_out"] == true
        exchange = balance_on(
            unit.output_interfaces[unit.m_el_out],
            unit.output_interfaces[unit.m_el_out].target
        )
        potential_energy_el = exchange.balance + exchange.energy_potential
        potential_storage_el = exchange.storage_potential
        if (unit.controller.parameter["unload_storages"] ? potential_energy_el + potential_storage_el : potential_energy_el) >= -parameters["epsilon"]
            return (nothing, nothing)
        end
        return (potential_energy_el, potential_storage_el)
    else
        return (-Inf, -Inf)
    end
end

function check_heat_out(
    unit::CHPP,
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

function calculate_energies(
    unit::CHPP,
    parameters::Dict{String,Any},
    watt_to_wh::Function,
    potentials::Vector{Float64}
)
    potential_energy_gas_in = potentials[1]
    potential_storage_gas_in = potentials[2]
    potential_energy_el_out = potentials[3]
    potential_storage_el_out = potentials[4]
    potential_energy_heat_out = potentials[5]
    potential_storage_heat_out = potentials[6]

    max_produce_heat = watt_to_wh(unit.power * (1.0 - unit.electricity_fraction))
    max_produce_el = watt_to_wh(unit.power * unit.electricity_fraction)
    max_consume_gas = max_produce_heat + max_produce_el

    # get usage fraction of external profile (normalized from 0 to 1)
    usage_fraction_operation_profile = unit.controller.parameter["operation_profile_path"] === nothing ? 1.0 : value_at_time(unit.controller.parameter["operation_profile"], parameters["time"])
    if usage_fraction_operation_profile <= 0.0
        return # no operation allowed from external profile
    end

    # all three standard operating strategies behave the same, but it is better to be
    # explicit about the behaviour rather than grouping all together
    if unit.controller.strategy == "storage_driven" && unit.controller.state_machine.state == 2
        usage_fraction_heat_out = -((unit.controller.parameter["load_storages"] ? potential_energy_heat_out + potential_storage_heat_out : potential_energy_heat_out) / max_produce_heat)
        usage_fraction_el_out = -((unit.controller.parameter["load_storages"] ? potential_energy_el_out + potential_storage_el_out : potential_energy_el_out) / max_produce_el)
        usage_fraction_gas_in = +((unit.controller.parameter["unload_storages"] ? potential_energy_gas_in + potential_storage_gas_in : potential_energy_gas_in) / max_consume_gas)

    elseif unit.controller.strategy == "storage_driven"
        return (false, nothing, nothing, nothing)

    elseif unit.controller.strategy == "supply_driven"
        usage_fraction_heat_out = -((unit.controller.parameter["load_storages"] ? potential_energy_heat_out + potential_storage_heat_out : potential_energy_heat_out) / max_produce_heat)
        usage_fraction_el_out = -((unit.controller.parameter["load_storages"] ? potential_energy_el_out + potential_storage_el_out : potential_energy_el_out) / max_produce_el)
        usage_fraction_gas_in = +((unit.controller.parameter["unload_storages"] ? potential_energy_gas_in + potential_storage_gas_in : potential_energy_gas_in) / max_consume_gas)

    elseif unit.controller.strategy == "demand_driven"
        usage_fraction_heat_out = -((unit.controller.parameter["load_storages"] ? potential_energy_heat_out + potential_storage_heat_out : potential_energy_heat_out) / max_produce_heat)
        usage_fraction_el_out = -((unit.controller.parameter["load_storages"] ? potential_energy_el_out + potential_storage_el_out : potential_energy_el_out) / max_produce_el)
        usage_fraction_gas_in = +((unit.controller.parameter["unload_storages"] ? potential_energy_gas_in + potential_storage_gas_in : potential_energy_gas_in) / max_consume_gas)

    end

    # limit actual usage by limits of inputs, outputs and profile
    usage_fraction = min(
        1.0,
        usage_fraction_heat_out,
        usage_fraction_el_out,
        usage_fraction_gas_in,
        usage_fraction_operation_profile
    )

    if usage_fraction < unit.min_power_fraction
        return (false, nothing, nothing, nothing)
    end

    return (
        true,
        max_consume_gas * usage_fraction,
        max_produce_el * usage_fraction,
        max_produce_heat * usage_fraction
    )
end

function potential(
    unit::CHPP,
    parameters::Dict{String,Any},
    watt_to_wh::Function
)
    potential_energy_gas_in, potential_storage_gas_in = check_gas_in(unit, parameters)
    if potential_energy_gas_in === nothing && potential_storage_gas_in === nothing
        set_max_energies!(unit, 0.0, 0.0, 0.0)
        return
    end

    potential_energy_el_out, potential_storage_el_out = check_el_out(unit, parameters)
    if potential_energy_el_out === nothing && potential_storage_el_out === nothing
        set_max_energies!(unit, 0.0, 0.0, 0.0)
        return
    end

    potential_energy_heat_out, potential_storage_heat_out = check_heat_out(unit, parameters)
    if potential_energy_heat_out === nothing && potential_storage_heat_out === nothing
        set_max_energies!(unit, 0.0, 0.0, 0.0)
        return
    end

    energies = calculate_energies(
        unit, parameters, watt_to_wh,
        [
            potential_energy_gas_in, potential_storage_gas_in,
            potential_energy_el_out, potential_storage_el_out,
            potential_energy_heat_out, potential_storage_heat_out
        ]
    )

    if !energies[1]
        set_max_energies!(unit, 0.0, 0.0, 0.0)
    else
        set_max_energies!(unit, energies[2], energies[3], energies[4])
    end
end

function produce(unit::CHPP, parameters::Dict{String,Any}, watt_to_wh::Function)
    potential_energy_gas_in, potential_storage_gas_in = check_gas_in(unit, parameters)
    if potential_energy_gas_in === nothing && potential_storage_gas_in === nothing
        set_max_energies!(unit, 0.0, 0.0, 0.0)
        return
    end

    potential_energy_el_out, potential_storage_el_out = check_el_out(unit, parameters)
    if potential_energy_el_out === nothing && potential_storage_el_out === nothing
        set_max_energies!(unit, 0.0, 0.0, 0.0)
        return
    end

    potential_energy_heat_out, potential_storage_heat_out = check_heat_out(unit, parameters)
    if potential_energy_heat_out === nothing && potential_storage_heat_out === nothing
        set_max_energies!(unit, 0.0, 0.0, 0.0)
        return
    end

    energies = calculate_energies(
        unit, parameters, watt_to_wh,
        [
            potential_energy_gas_in, potential_storage_gas_in,
            potential_energy_el_out, potential_storage_el_out,
            potential_energy_heat_out, potential_storage_heat_out
        ]
    )

    if energies[1]
        sub!(unit.input_interfaces[unit.m_gas_in], energies[2])
        add!(unit.output_interfaces[unit.m_el_out], energies[3])
        add!(unit.output_interfaces[unit.m_heat_out], energies[4])
    end
end

export CHPP