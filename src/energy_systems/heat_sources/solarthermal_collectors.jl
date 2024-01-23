"""
Implementation of a solarthermal collector.
Works for flat plate collectors, vacuum tube collectors and PVT modules.

Stagnation is ignored under the assumption that the system has either measures to prevent 
stagnation harming the collectors or the designed size is small enough for stagnation not 
to become a problem.

## ATTENTION: Geothermal heat collector is currently work in progress and not completed!!
"""
mutable struct SolarthermalCollector <: Component
    uac::String
    controller::Controller
    sys_function::SystemFunction

    input_interfaces::InterfaceMap
    output_interfaces::InterfaceMap

    m_heat_out::Symbol

    ambient_temperature_profile::Union{Profile,Nothing}

    current_output_temperature::Temperature
    current_input_temperature::Temperature
    ambient_temperature::Temperature

    function SolarthermalCollector(uac::String, config::Dict{String,Any})

        m_heat_out = Symbol(default(config, "m_heat_out", "m_h_w_ht1"))
        register_media([m_heat_out])

        return new(
            uac, # uac
            controller_for_strategy( # controller
                config["strategy"]["name"], config["strategy"]
            ),
            sf_bounded_source, # sys_function
            InterfaceMap( # input_interfaces

            ),
            InterfaceMap( # output_interfaces
                m_heat_out => nothing
            ),
            m_heat_out, # medium name of output interface
            ambient_temperature_profile, # ambient temperature profile

            nothing, # output temperature in current time step, calculated in control()
            nothing, # input temperature in current time step, calculated in control()
            nothing, # ambient temperature in current time step, calculated in control()
        )
    end
end

function output_values(unit::SolarthermalCollector)::Vector{String}

    return [string(unit.m_heat_out)*" OUT",
            "Temperature"]
end

function output_value(unit::SolarthermalCollector, key::OutputKey)::Float64
    if key.value_key == "OUT"
        return calculate_energy_flow(unit.output_interfaces[key.m_heat_out])
    elseif key.value_key == "Temperature"
        return unit.temperature
    end
    throw(KeyError(key.value_key))
end

function control(
    unit::SolarthermalCollector,
    components::Grouping,
    parameters::Dict{String,Any}
)
    move_state(unit, components, parameters)

    # set_max_energy!(unit.output_interfaces[unit.m_heat_out], unit.max_energy)

    unit.output_interfaces[unit.m_heat_out].temperature = highest_temperature(
        unit.temperature,
        unit.output_interfaces[unit.m_heat_out].temperature
    )
end

function process(unit::SolarthermalCollector, parameters::Dict{String,Any})
    outface = unit.output_interfaces[unit.m_heat_out]
    exchange = balance_on(outface, outface.target)
    if exchange.balance < 0.0
        add!(
            outface,
            min(abs(exchange.balance), unit.max_energy),
            unit.temperature
        )
    end
end

export SolarthermalCollector