#! format: off
const GENERIC_HEAT_SOURCE_PARAMETERS = Dict(
    "medium" => (
        description="Medium of the heat source",
        display_name="Medium",
        required=true,
        type=String,
        json_type="string",
        unit="-"
    ),
    "temperature_profile_file_path" => (
        default=nothing,
        description="Path to a temperature profile file",
        display_name="Temperature profile file",
        required=false,
        conditionals=[
            ("temperature_from_global_file", "mutex"),
            ("constant_temperature", "mutex")
        ],
        type=String,
        json_type="string",
        unit="-"
    ),
    "temperature_from_global_file" => (
        default=false,
        description="If true, take the temperature profile from the global weather data file",
        display_name="Temperature from global file",
        required=false,
        conditionals=[
            ("temperature_profile_file_path", "mutex"),
            ("constant_temperature", "mutex")
        ],
        type=Bool,
        json_type="boolean",
        unit="-"
    ),
    "constant_temperature" => (
        default=nothing,
        description="Constant temperature value",
        display_name="Constant temperature",
        required=false,
        conditionals=[
            ("temperature_profile_file_path", "mutex"),
            ("temperature_from_global_file", "mutex")
        ],
        type=Float64,
        json_type="number",
        unit="°C"
    ),
    "max_power_profile_file_path" => (
        default=nothing,
        description="Path to a profile file with maximum power values",
        display_name="Max. power profile file",
        required=false,
        conditionals=[
            ("constant_power", "mutex"),
        ],
        type=String,
        json_type="string",
        unit="-"
    ),
    "constant_power" => (
        default=nothing,
        description="Constant maximum power to supply",
        display_name="Constant max. power",
        required=false,
        conditionals=[
            ("max_power_profile_file_path", "mutex"),
        ],
        type=Float64,
        json_type="number",
        unit="W"
    ),
    "scale" => (
        default=1.0,
        description="Scaling factor for the max. power profile",
        display_name="Max. power scale",
        required=false,
        conditionals=[("max_power_profile_file_path", "is_not_nothing")],
        type=Float64,
        json_type="number",
        unit="-"
    ),
    "temperature_reduction_model" => (
        default="none",
        description="The temperature reduction model to use",
        display_name="Temp. reduction model",
        required=false,
        options=["none", "constant", "lmtd"],
        type=String,
        json_type="string",
        unit="-"
    ),
    "min_source_in_temperature" => (
        default=nothing,
        description="Minimum source input temperature",
        display_name="Min. source temp.",
        required=false,
        conditionals=[("temperature_reduction_model", "has_value", "lmtd")],
        type=Float64,
        json_type="number",
        unit="°C"
    ),
    "max_source_in_temperature" => (
        default=nothing,
        description="Maximum source input temperature",
        display_name="Max. source temp.",
        required=false,
        conditionals=[("temperature_reduction_model", "has_value", "lmtd")],
        type=Float64,
        json_type="number",
        unit="°C"
    ),
    "minimal_reduction" => (
        default=2.0,
        description="Minimal (or constant) reduction of temperature",
        display_name="Min. reduction",
        required=false,
        conditionals=[("temperature_reduction_model", "is_one_of", ("lmtd", "constant"))],
        type=Float64,
        json_type="number",
        unit="K"
    ),
)
#! format: on

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

    function GenericHeatSource(uac::String, config::Dict{String,Any}, sim_params::Dict{String,Any})
        return new(SSOT_parameter_constructor(GenericHeatSource, uac, config, sim_params)...)
    end
end

function component_parameters(x::Type{GenericHeatSource})::Dict{String,NamedTuple}
    return deepcopy(GENERIC_HEAT_SOURCE_PARAMETERS) # return a copy to prevent external modification
end

function extract_parameter(x::Type{GenericHeatSource}, config::Dict{String,Any}, param_name::String,
                           param_def::NamedTuple, sim_params::Dict{String,Any}, uac::String)
    if param_name == "constant_temperature" || param_name == "temperature_profile_file_path"
        constant_temperature,
        temperature_profile = get_parameter_profile_from_config(config,
                                                                sim_params,
                                                                "temperature",
                                                                "temperature_profile_file_path",
                                                                "temperature_from_global_file",
                                                                "constant_temperature",
                                                                uac)
        return param_name == "constant_temperature" ? constant_temperature : temperature_profile
    end

    return extract_parameter(Component, config, param_name, param_def, sim_params, uac)
end

function validate_config(x::Type{GenericHeatSource}, config::Dict{String,Any}, extracted::Dict{String,Any},
                         uac::String, sim_params::Dict{String,Any})
    validate_config(Component, extracted, uac, sim_params, component_parameters(GenericHeatSource))
end

function init_from_params(x::Type{GenericHeatSource}, uac::String, params::Dict{String,Any},
                          raw_params::Dict{String,Any}, sim_params::Dict{String,Any})::Tuple
    medium = Symbol(params["medium"])
    register_media([medium])

    max_power_profile = params["max_power_profile_file_path"] !== nothing ?
                        Profile(params["max_power_profile_file_path"], sim_params) :
                        nothing

    # return tuple in the order expected by new()
    return (uac,
            Controller(params["control_parameters"]),
            sf_flexible_source,
            medium,
            InterfaceMap(medium => nothing),
            InterfaceMap(medium => nothing),
            max_power_profile,
            params["temperature_profile_file_path"], # temperature_profile, might be from global weather data
            params["scale"],
            params["constant_power"],
            params["constant_temperature"],
            params["temperature_reduction_model"],
            params["min_source_in_temperature"],
            params["max_source_in_temperature"],
            nothing,                     # avg_source_in_temperature
            params["minimal_reduction"], # lmtd_min
            0.0,                         # max_energy
            nothing,                     # temperature_src_in
            nothing)                     # temperature_snk_out
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
        temperature_max, _ = check_temperatures_source(exchanges, unit.temperature_snk_out, unit.max_energy)
    end

    if sum(energy_demand; init=0.0) < 0.0
        add!(outface, abs.(energy_demand), temperature_min, temperature_max)
    end
end

function output_values(unit::GenericHeatSource)::Vector{String}
    if unit.temperature_profile === nothing && unit.constant_temperature === nothing
        return [string(unit.medium) * ":OUT",
                "Max_Energy"]
    else
        return [string(unit.medium) * ":OUT",
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
