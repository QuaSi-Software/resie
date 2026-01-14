#! format: off
const FLEXIBLE_SUPPLY_PARAMETERS = Dict(
    "medium" => (
        description="Medium of the supply (e.g. electricity, heat, gas, etc.)",
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
        unit="W"
    ),
)
#! format: on

"""
Implementation of a component modeling an abstract flexible supply of some medium.

This is particularly useful for testing, but can also be used to model any flexible
component or other equipment unit that processes energy in a given medium. The component
might still have a maximum power draw in a single time step, but can provide any fraction
of this to connected components.
"""
mutable struct FlexibleSupply <: Component
    uac::String
    controller::Controller
    sys_function::SystemFunction
    medium::Symbol

    input_interfaces::InterfaceMap
    output_interfaces::InterfaceMap

    max_power_profile::Union{Profile,Nothing}
    temperature_profile::Union{Profile,Nothing}
    scaling_factor::Float64

    max_energy::Float64
    temperature::Temperature
    constant_power::Union{Nothing,Float64}
    constant_temperature::Temperature

    function FlexibleSupply(uac::String, config::Dict{String,Any}, sim_params::Dict{String,Any})
        constructor_errored = false

        # extract all parameters using the parameter dictionary as the source of truth
        extracted_params = Dict{String,Any}()
        for (param_name, param_def) in FLEXIBLE_SUPPLY_PARAMETERS
            try
                extracted_params[param_name] = extract_parameter(FlexibleSupply, config, param_name,
                                                                 param_def, sim_params, uac)
            catch e
                @error "$(sprint(showerror, e))"
                constructor_errored = true
            end
        end

        try
            # extract control_parameters, which is essentially a subconfig
            extracted_params["control_parameters"] = extract_control_parameters(Component, config)

            # validate configuration, e.g. for interdependencies and allowed values
            validate_config(FlexibleSupply, config, extracted_params, uac, sim_params)
        catch e
            @error "$(sprint(showerror, e))"
            constructor_errored = true
        end

        # we delayed throwing the errors upwards so that many errors are caught at once
        if constructor_errored
            throw(InputError("Can't construct component $uac because of errors parsing " *
                             "the parameter config. Check the error log for more " *
                             "information on each error."))
        end

        # initialize and construct the object
        init_values = init_from_params(FlexibleSupply, uac, extracted_params, config, sim_params)
        return new(init_values...)
    end
end

function component_parameters(x::Type{FlexibleSupply})::Dict{String,NamedTuple}
    return deepcopy(FLEXIBLE_SUPPLY_PARAMETERS) # return a copy to prevent external modification
end

function extract_parameter(x::Type{FlexibleSupply}, config::Dict{String,Any}, param_name::String,
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

    return extract_parameter(Component, config, param_name, param_def, sim_params)
end

function validate_config(x::Type{FlexibleSupply}, config::Dict{String,Any}, extracted::Dict{String,Any},
                         uac::String, sim_params::Dict{String,Any})
    validate_config(Component, extracted, uac, sim_params, component_parameters(FlexibleSupply))
end

function init_from_params(x::Type{FlexibleSupply}, uac::String, params::Dict{String,Any},
                          raw_params::Dict{String,Any}, sim_params::Dict{String,Any})::Tuple
    medium = Symbol(params["medium"])
    register_media([medium])

    max_power_profile = params["max_power_profile_file_path"] !== nothing ?
                        Profile(params["max_power_profile_file_path"], sim_params) :
                        nothing

    # return tuple in the order expected by new()
    return (uac,                                     # uac
            Controller(params["control_parameters"]),
            sf_flexible_source,                      # sys_function
            medium,                                  # medium
            InterfaceMap(medium => nothing),         # input_interfaces
            InterfaceMap(medium => nothing),         # output_interfaces
            max_power_profile,                       # max_power_profile
            params["temperature_profile_file_path"], # temperature_profile, might be from global weather data
            params["scale"],                         # scaling_factor
            0.0,                                     # max_energy
            nothing,                                 # temperature
            params["constant_power"],                # constant_power
            params["constant_temperature"])          # constant_temperature
end

function initialise!(unit::FlexibleSupply, sim_params::Dict{String,Any})
    set_storage_transfer!(unit.output_interfaces[unit.medium],
                          load_storages(unit.controller, unit.medium))
end

function control(unit::FlexibleSupply,
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
        unit.temperature = unit.constant_temperature
    elseif unit.temperature_profile !== nothing
        unit.temperature = Profiles.value_at_time(unit.temperature_profile, sim_params)
    end
    set_max_energy!(unit.output_interfaces[unit.medium], unit.max_energy, nothing, unit.temperature)
end

function process(unit::FlexibleSupply, sim_params::Dict{String,Any})
    outface = unit.output_interfaces[unit.medium]
    exchanges = balance_on(outface, outface.target)

    # if we get the exchanges from a bus, the temperature check has already been performed
    if outface.target.sys_function == EnergySystems.sf_bus
        energy_demand = [e.balance + e.energy_potential for e in exchanges]
        temperature_min = [nothing for _ in exchanges]
        temperature_max = [unit.temperature for _ in exchanges]
    else # check temperature
        energy_demand,
        temperature_min,
        temperature_max, _ = check_temperatures_source(exchanges, unit.temperature, unit.max_energy)
    end

    if sum(energy_demand; init=0.0) < 0.0
        add!(outface, abs.(energy_demand), temperature_min, temperature_max)
    end
end

function output_values(unit::FlexibleSupply)::Vector{String}
    if unit.temperature_profile === nothing && unit.constant_temperature === nothing
        return [string(unit.medium) * ":OUT",
                "Max_Energy"]
    else
        return [string(unit.medium) * ":OUT",
                "Max_Energy",
                "Temperature"]
    end
end

function output_value(unit::FlexibleSupply, key::OutputKey)::Float64
    if key.value_key == "OUT"
        return calculate_energy_flow(unit.output_interfaces[key.medium])
    elseif key.value_key == "Max_Energy"
        return unit.max_energy
    elseif key.value_key == "Temperature"
        return unit.temperature
    end
    throw(KeyError(key.value_key))
end

export FlexibleSupply
