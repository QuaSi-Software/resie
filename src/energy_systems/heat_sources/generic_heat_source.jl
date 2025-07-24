"""
A generic heat source that models a number of different technologies used to supply heat.

Typically these are sources of low temperature heat for heat pumps, that only supply heat
and do not store it. This makes them different from other sources e.g. geothermal sources.
Examples of technologies: river water, lake water, waste water, atmosphere, industrial
processes or one-way connections to a heat network.
"""
mutable struct GenericHeatSource <: Component
    uac::String
    controller::Controller
    sys_function::SystemFunction
    medium::Symbol

    input_interfaces::InterfaceMap
    output_interfaces::InterfaceMap

    max_power_profile::Union{Profile,Nothing}
    temperature_profile::Union{Profile,Nothing}
    scaling_factor::Float64
    constant_power::Union{Nothing,Float64}
    constant_temperature::Temperature

    temperature_reduction_model::String
    min_source_in_temperature::Temperature
    max_source_in_temperature::Temperature
    avg_source_in_temperature::Temperature
    lmtd_min::Float64

    max_energy::Float64
    temperature_src_in::Temperature
    temperature_snk_out::Temperature

    function GenericHeatSource(uac::String,
                               config::Dict{String,Any},
                               sim_params::Dict{String,Any})
        max_power_profile = "max_power_profile_file_path" in keys(config) ?
                            Profile(config["max_power_profile_file_path"], sim_params) :
                            nothing

        constant_temperature,
        temperature_profile = get_parameter_profile_from_config(config,
                                                                sim_params,
                                                                "temperature",
                                                                "temperature_profile_file_path",
                                                                "temperature_from_global_file",
                                                                "constant_temperature",
                                                                uac; required=true)

        medium = Symbol(config["medium"])
        register_media([medium])

        return new(uac, # uac
                   Controller(default(config, "control_parameters", nothing)),
                   sf_bounded_source,               # sys_function
                   medium,                          # medium
                   InterfaceMap(medium => nothing), # input_interfaces
                   InterfaceMap(medium => nothing), # output_interfaces
                   max_power_profile,
                   temperature_profile,
                   default(config, "scale", 1.0),   # scaling_factor
                   default(config, "constant_power", nothing),
                   constant_temperature,
                   default(config, "temperature_reduction_model", "none"),
                   default(config, "min_source_in_temperature", nothing),
                   default(config, "max_source_in_temperature", nothing),
                   nothing,                                   # avg_source_in_temperature
                   default(config, "minimal_reduction", 2.0), # lmtd_min
                   0.0,     # max_energy
                   nothing, # temperature_src_in
                   nothing) # temperature_snk_out
    end
end

function initialise!(unit::GenericHeatSource, sim_params::Dict{String,Any})
    set_storage_transfer!(unit.output_interfaces[unit.medium],
                          load_storages(unit.controller, unit.medium))

    if unit.temperature_reduction_model == "lmtd"
        if unit.min_source_in_temperature === nothing && unit.temperature_profile !== nothing
            unit.min_source_in_temperature = Profiles.minimum(unit.temperature_profile)
        end
        if unit.max_source_in_temperature === nothing && unit.temperature_profile !== nothing
            unit.max_source_in_temperature = Profiles.maximum(unit.temperature_profile)
        end
        unit.avg_source_in_temperature = 0.5 * (unit.min_source_in_temperature + unit.max_source_in_temperature)
    end
end

function control(unit::GenericHeatSource,
                 components::Grouping,
                 sim_params::Dict{String,Any})
    update(unit.controller)

    if unit.constant_power !== nothing
        unit.max_energy = sim_params["watt_to_wh"](unit.constant_power)
    elseif unit.max_power_profile !== nothing
        unit.max_energy = unit.scaling_factor * Profiles.work_at_time(unit.max_power_profile, sim_params)
    else
        unit.max_energy = 0.0
    end

    if unit.constant_temperature !== nothing
        unit.temperature_src_in = unit.constant_temperature
    elseif unit.temperature_profile !== nothing
        unit.temperature_src_in = Profiles.value_at_time(unit.temperature_profile, sim_params)
    end

    if unit.temperature_reduction_model == "constant" && unit.temperature_src_in !== nothing
        unit.temperature_snk_out = unit.temperature_src_in - unit.lmtd_min

    elseif unit.temperature_reduction_model == "lmtd" && unit.temperature_src_in !== nothing
        alpha = min(0.95, # alpha_max
                    1 -
                    abs(unit.avg_source_in_temperature - unit.temperature_src_in)
                    /
                    (unit.max_source_in_temperature - unit.min_source_in_temperature))
        unit.temperature_snk_out = unit.temperature_src_in - unit.lmtd_min * log(1 / alpha) / (1 - alpha)

    else
        unit.temperature_snk_out = unit.temperature_src_in
    end
    set_max_energy!(unit.output_interfaces[unit.medium], unit.max_energy, nothing, unit.temperature_snk_out)
end

function process(unit::GenericHeatSource, sim_params::Dict{String,Any})
    outface = unit.output_interfaces[unit.medium]
    exchanges = balance_on(outface, outface.target)

    # if we get the exchanges from a bus, the temperature check has already been performed
    if outface.target.sys_function == EnergySystems.sf_bus
        energy_demand = [e.balance + e.energy_potential for e in exchanges]
        temperature_min = [nothing for _ in exchanges]
        temperature_max = [unit.temperature_snk_out for _ in exchanges]
    else # check temperatures
        energy_demand,
        temperature_min,
        temperature_max = check_temperatures_source(exchanges, unit.temperature_snk_out, unit.max_energy)
    end

    if sum(energy_demand; init=0.0) < 0.0
        add!(outface, abs.(energy_demand), temperature_min, temperature_max)
    end
end

function output_values(unit::GenericHeatSource)::Vector{String}
    if unit.temperature_profile === nothing && unit.constant_temperature === nothing
        return [string(unit.medium) * " OUT",
                "Max_Energy"]
    else
        return [string(unit.medium) * " OUT",
                "Max_Energy",
                "Temperature_src_in",
                "Temperature_snk_out"]
    end
end

function output_value(unit::GenericHeatSource, key::OutputKey)::Float64
    if key.value_key == "OUT"
        return calculate_energy_flow(unit.output_interfaces[key.medium])
    elseif key.value_key == "Max_Energy"
        return unit.max_energy
    elseif key.value_key == "Temperature_src_in"
        return unit.temperature_src_in
    elseif key.value_key == "Temperature_snk_out"
        return unit.temperature_snk_out
    end
    throw(KeyError(key.value_key))
end

export GenericHeatSource
