#! format: off
const STORAGE_COMPONENT_PARAMETERS = Dict(
    "medium" => (
        description="Medium of the storage (e.g. electricity, heat, gas, etc.)",
        display_name="Medium",
        required=true,
        type=String,
        json_type="string",
        unit="-"
    ),
    "capacity" => (
        default=nothing,
        description="Energy capacity of the storage",
        display_name="Capacity",
        required=true,
        validations=[
            ("self", "value_gte_num", 0.0)
        ],
        type=Float64,
        json_type="number",
        unit="Wh"
    ),
    "initial_load" => (
        default=0.0,
        description="Initial load relative to capacity",
        display_name="Initial load",
        required=false,
        validations=[
            ("self", "value_gte_num", 0.0),
            ("self", "value_lte_num", 1.0),
        ],
        type=Float64,
        json_type="number",
        unit="-"
    ),
)

const STORAGE_ECONOMIC_PARAMETERS = get_economic_standard_params("storage",
    Dict{String,Any}(
            "lifetime_years" => 20,
            "capex_specific" => nothing,
            "capex_price_change_rate_per_year" => 0.01,
            "maintenance_inspection_rate_per_year" => 0.02,
            "maintenance_inspection_price_change_rate_per_year" =>  0.005,
            "repair_rate_per_year" => 0.05,
            "repair_price_change_rate_per_year" =>  0.005,
            "operational_labour_hours_per_year" =>  0.5,
            "subsidy_rate_of_capex" =>  nothing,
            "subsidy_max" =>  nothing
    ),
    Dict{String,Any}(
            "capex_specific" => "€/Wh"
    )
)

const STORAGE_EMISSIONS_PARAMETERS = get_emissions_standard_params("storage",
    Dict{String,Any}(
        "lifetime_years" => 20,
        "embodied_emissions_specific" => "const:0.0",
        "embodied_emissions_change_rate_per_year" => 0.0
    ),
    Dict{String,Any}(
        "embodied_emissions_specific" => "g CO2/Wh"
    ),
)
#! format: on

"""
Implementation of a component modeling a generic storage of a chosen medium.

This is particularly useful for testing, but can also be used to model any storage or other
equipment unit that stores energy in a given medium.
"""
mutable struct Storage <: Component
    uac::String
    controller::Controller
    sys_function::SystemFunction

    input_interfaces::InterfaceMap
    output_interfaces::InterfaceMap

    medium::Symbol

    economic_parameters::Dict{String,Any}
    emissions_parameters::Dict{String,Any}

    capacity::Float64
    load::Float64
    load_end_of_last_timestep::Float64
    losses::Float64

    # indicates if the process step has already been performed in the current time step
    process_done::Bool
    # indicates if the load step has already been performed in the current time step
    load_done::Bool

    function Storage(uac::String, config::Dict{String,Any}, sim_params::Dict{String,Any})
        return new(SSOT_parameter_constructor(Storage, uac, config, sim_params)...)
    end
end

function component_parameters(x::Type{Storage})::Dict{String,Any}
    return deepcopy(STORAGE_COMPONENT_PARAMETERS)
end

function economic_parameters(x::Type{Storage})::Dict{String,Any}
    return deepcopy(STORAGE_ECONOMIC_PARAMETERS)
end

function emissions_parameters(x::Type{Storage})::Dict{String,Any}
    return deepcopy(STORAGE_EMISSIONS_PARAMETERS)
end

function extract_parameter(x::Type{Storage}, config::Dict{String,Any}, param_name::String, param_def::NamedTuple,
                           sim_params::Dict{String,Any}, uac::String)
    return extract_parameter(Component, config, param_name, param_def, sim_params, uac)
end

function validate_config(x::Type{Storage}, config::Dict{String,Any}, extracted::Dict{String,Any}, uac::String,
                         sim_params::Dict{String,Any}, param_type::String)
    if param_type == "economy"
        parameter = economic_parameters(Storage)
        uac = uac * " - economic_parameters"
    elseif param_type == "emissions"
        parameter = emissions_parameters(Storage)
        uac = uac * " - emissions_parameters"
    elseif param_type == "component"
        parameter = component_parameters(Storage)
    end
    validate_config(Component, extracted, uac, sim_params, parameter)
end

function init_from_params(x::Type{Storage}, uac::String, params::Dict{String,Any},
                          raw_params::Dict{String,Any}, sim_params::Dict{String,Any})::Tuple
    medium = Symbol(params["medium"])

    # return tuple in the order expected by new()
    return (uac,                             # uac
            Controller(params["control_parameters"]),
            sf_storage,                      # sys_function
            InterfaceMap(medium => nothing), # input_interfaces
            InterfaceMap(medium => nothing), # output_interfaces
            medium,                          # medium
            params["economic_parameters"],
            params["emissions_parameters"],
            params["capacity"],              # capacity
            params["initial_load"] * params["capacity"], # load
            0.0,                             # load_end_of_last_timestep
            0.0,                             # losses
            false,                           # process_done
            false)                           # load_done
end

function initialise!(unit::Storage, sim_params::Dict{String,Any})
    set_storage_transfer!(unit.input_interfaces[unit.medium],
                          unload_storages(unit.controller, unit.medium))
    set_storage_transfer!(unit.output_interfaces[unit.medium],
                          load_storages(unit.controller, unit.medium))

    unit.load_end_of_last_timestep = copy(unit.load)
end

function control(unit::Storage,
                 components::Grouping,
                 sim_params::Dict{String,Any})
    update(unit.controller)

    set_max_energy!(unit.input_interfaces[unit.medium], unit.capacity - unit.load)
    set_max_energy!(unit.output_interfaces[unit.medium], unit.load)
end

function process(unit::Storage, sim_params::Dict{String,Any})
    outface = unit.output_interfaces[unit.medium]
    exchanges = balance_on(outface, outface.target)
    energy_demand = balance(exchanges) + energy_potential(exchanges)

    if energy_demand >= 0.0 || unit.capacity <= sim_params["epsilon"]
        handle_component_update!(unit, "process", sim_params)
        set_max_energy!(unit.output_interfaces[unit.medium], 0.0)
        return # process is only concerned with moving energy to the target
    end

    if unit.load > abs(energy_demand)
        unit.load += energy_demand
        add!(outface, abs(energy_demand))
    else
        add!(outface, unit.load)
        unit.load = 0.0
    end
    handle_component_update!(unit, "process", sim_params)
end

function handle_component_update!(unit::Storage, step::String, sim_params::Dict{String,Any})
    if step == "process"
        unit.process_done = true
    elseif step == "load"
        unit.load_done = true
    end
    if unit.process_done && unit.load_done
        # update component
        unit.load_end_of_last_timestep = copy(unit.load)
        # reset 
        unit.process_done = false
        unit.load_done = false
    end
end

function load(unit::Storage, sim_params::Dict{String,Any})
    inface = unit.input_interfaces[unit.medium]
    exchanges = balance_on(inface, inface.source)
    energy_available = balance(exchanges) + energy_potential(exchanges)

    if energy_available <= 0.0 || unit.capacity <= sim_params["epsilon"]
        handle_component_update!(unit, "load", sim_params)
        set_max_energy!(unit.input_interfaces[unit.medium], 0.0)
        return # load is only concerned with receiving energy from the source
    end

    diff = unit.capacity - unit.load
    if diff > energy_available
        unit.load += energy_available
        sub!(inface, energy_available)
    else
        unit.load = unit.capacity
        sub!(inface, diff)
    end

    handle_component_update!(unit, "load", sim_params)
end

function get_reference_for_capex_and_embodied_emissions(unit::Storage)
    return unit.capacity # [Wh]
end

function output_values(unit::Storage)::Vector{String}
    return [string(unit.medium) * ":IN",
            string(unit.medium) * ":OUT",
            "Load",
            "Load%",
            "Capacity",
            "LossesGains"]
end

function output_value(unit::Storage, key::OutputKey)::Float64
    if key.value_key == "IN"
        return calculate_energy_flow(unit.input_interfaces[key.medium])
    elseif key.value_key == "OUT"
        return calculate_energy_flow(unit.output_interfaces[key.medium])
    elseif key.value_key == "Load"
        return unit.load
    elseif key.value_key == "Load%"
        return 100 * unit.load / unit.capacity
    elseif key.value_key == "Capacity"
        return unit.capacity
    elseif key.value_key == "LossesGains"
        return -unit.losses
    end
    throw(KeyError(key.value_key))
end

export Storage
