"""
Implementation of a component modeling the connection to a public grid of a certain medium.

Public grids are considered to have an unlimited amount of energy they can provide, but
might be limited in the power they can provide, although this behaviour is not yet
implemented. They exist to model real connections to a public grid that can provide any
remaining demand of energy or take in any excess of energy. To make it possible to model a
one-way connection they are split into two instances for providing or receiving energy and
must be handled as such in the input for constructing a project.
"""
mutable struct GridConnection <: Component
    uac::String
    controller::Controller
    sys_function::SystemFunction
    medium::Symbol

    input_interfaces::InterfaceMap
    output_interfaces::InterfaceMap

    output_sum::Float64
    input_sum::Float64

    function GridConnection(uac::String, config::Dict{String,Any}, sim_params::Dict{String,Any})
        medium = Symbol(config["medium"])
        register_media([medium])

        return new(
            uac, # uac
            Controller(default(config, "control_parameters", nothing)),
            if Bool(config["is_source"])
                sf_bounded_source
            else
                sf_bounded_sink
            end, # sys_function
            medium, # medium
            InterfaceMap( # input_interfaces
                medium => nothing
            ),
            InterfaceMap( # output_interfaces
                medium => nothing
            ),
            0.0, # output_sum,
            0.0, # input_sum
        )
    end
end

function initialise!(unit::GridConnection, sim_params::Dict{String,Any})
    if unit.sys_function === sf_bounded_source
        set_storage_transfer!(
            unit.output_interfaces[unit.medium],
            load_storages(unit.controller, unit.medium)
        )
    else
        set_storage_transfer!(
            unit.input_interfaces[unit.medium],
            unload_storages(unit.controller, unit.medium)
        )
    end
end

function control(
    unit::GridConnection,
    components::Grouping,
    sim_params::Dict{String,Any}
)
    update(unit.controller)
    if unit.sys_function === sf_bounded_source
        set_max_energy!(unit.output_interfaces[unit.medium], Inf)
    else
        set_max_energy!(unit.input_interfaces[unit.medium], Inf)
    end
end

function process(unit::GridConnection, sim_params::Dict{String,Any})
    if unit.sys_function === sf_bounded_source
        outface = unit.output_interfaces[unit.medium]
        exchanges = balance_on(outface, outface.target)
        blnc = balance(exchanges) + energy_potential(exchanges)
        if blnc < 0.0
            unit.output_sum += blnc
            add!(outface, abs(blnc))
        end
    else
        inface = unit.input_interfaces[unit.medium]
        exchanges = balance_on(inface, inface.source)
        blnc = balance(exchanges) + energy_potential(exchanges)
        if blnc > 0.0
            unit.input_sum += blnc
            sub!(inface, abs(blnc))
        end
    end
end

function output_values(unit::GridConnection)::Vector{String}
    if unit.sys_function == sf_bounded_source
        return [string(unit.medium)*" OUT",
                "Output_sum"]
    elseif unit.sys_function == sf_bounded_sink
        return [string(unit.medium)*" IN",
                "Input_sum"]
    end
end

function output_value(unit::GridConnection, key::OutputKey)::Float64
    if key.value_key == "IN"
        return calculate_energy_flow(unit.input_interfaces[key.medium])
    elseif key.value_key == "OUT"
        return calculate_energy_flow(unit.output_interfaces[key.medium])
    elseif key.value_key == "Output_sum"
        return unit.output_sum
    elseif key.value_key == "Input_sum"
        return unit.input_sum
    end
    throw(KeyError(key.value_key))
end

export GridConnection