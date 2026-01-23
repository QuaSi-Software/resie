#! format: off
const BUFFER_TANK_PARAMETERS = Dict(
    "medium" => (
        default="m_h_w_ht1",
        description="Heat medium of the buffer tank",
        display_name="Medium",
        required=false,
        type=String,
        json_type="string",
        unit="-"
    ),
    "model_type" => (
        default="ideally_stratified",
        description="Operation model: 'ideally_stratified', 'balanced', or 'ideally_mixed'",
        display_name="Model type",
        required=false,
        type=String,
        json_type="string",
        options=["ideally_stratified", "balanced", "ideally_mixed"],
        unit="-"
    ),
    "consider_losses" => (
        default=false,
        description="Toggle to en-/disable the calculation of losses to the ambient",
        display_name="Consider losses to the ambient?",
        required=false,
        type=Bool,
        json_type="boolean",
        unit="-"
    ),
    "ambient_temperature_profile_file_path" => (
        default=nothing,
        description="Path to a temperature profile file",
        display_name="Ambient temp. profile",
        required=false,
        conditionals=[
            ("consider_losses", "is_true"),
            ("ambient_temperature_from_global_file", "mutex"),
            ("constant_ambient_temperature", "mutex")
        ],
        type=String,
        json_type="string",
        unit="-"
    ),
    "ambient_temperature_from_global_file" => (
        default=false,
        description="If true, take the ambient temperature profile from the global weather data file",
        display_name="Ambient temp. from global file",
        required=false,
        conditionals=[
            ("consider_losses", "is_true"),
            ("ambient_temperature_profile_file_path", "mutex"),
            ("constant_ambient_temperature", "mutex")
        ],
        type=Bool,
        json_type="boolean",
        unit="-"
    ),
    "constant_ambient_temperature" => (
        default=nothing,
        description="Constant ambient temperature value",
        display_name="Constant ambient temp.",
        required=false,
        conditionals=[
            ("consider_losses", "is_true"),
            ("ambient_temperature_profile_file_path", "mutex"),
            ("ambient_temperature_from_global_file", "mutex")
        ],
        type=Float64,
        json_type="number",
        unit="°C"
    ),
    "capacity" => (
        default=nothing,
        description="Capacity of the tank as energy value",
        display_name="Capacity (energy)",
        required=false,
        conditionals=[
            ("volume", "mutex"),
        ],
        validations=[
            ("self", "value_gt_num_or_nothing", 0.0),
        ],
        type=Float64,
        json_type="number",
        unit="Wh"
    ),
    "volume" => (
        default=nothing,
        description="Capacity of the tank as volume",
        display_name="Capacity (volume)",
        required=false,
        conditionals=[
            ("volume", "mutex"),
        ],
        validations=[
            ("self", "value_gt_num_or_nothing", 0.0),
        ],
        type=Float64,
        json_type="number",
        unit="m^3"
    ),
    "initial_load" => (
        default=0.0,
        description="Initial load as relative value",
        display_name="Initial load (relative)",
        required=false,
        validations=[
            ("self", "value_gte_num", 0.0),
            ("self", "value_lte_num", 1.0),
        ],
        type=Float64,
        json_type="number",
        unit="-"
    ),
    "rho_medium" => (
        default=1000.0,
        description="Density of the fluid",
        display_name="Density fluid",
        required=false,
        validations=[
            ("self", "value_gt_num", 0.0),
        ],
        type=Float64,
        json_type="number",
        unit="kg/m^3"
    ),
    "cp_medium" => (
        default=4.18,
        description="Specific heat capacity of the fluid",
        display_name="Spec. heat capacity fluid",
        required=false,
        validations=[
            ("self", "value_gt_num", 0.0),
        ],
        type=Float64,
        json_type="number",
        unit="kJ/kg*K"
    ),
    "high_temperature" => (
        default=75.0,
        description="Upper temperature of the buffer tank to the supply line",
        display_name="Upper temperature",
        required=false,
        type=Float64,
        json_type="number",
        unit="°C"
    ),
    "low_temperature" => (
        default=20.0,
        description="Lower temperature of the buffer tank from the return line",
        display_name="Lower temperature",
        required=false,
        validations=[
            ("self", "value_lt_rel", "high_temperature"),
        ],
        type=Float64,
        json_type="number",
        unit="°C"
    ),
    "max_load_rate" => (
        default=nothing,
        description="Maximum load rate of the tank relative to the capacity",
        display_name="Max. load rate",
        required=false,
        validations=[
            ("self", "value_gt_num_or_nothing", 0.0),
        ],
        type=Floathing,
        json_type="number",
        unit="1/h"
    ),
    "max_unload_rate" => (
        default=nothing,
        description="Maximum unload rate of the tank relative to the capacity",
        display_name="Max. unload rate",
        required=false,
        validations=[
            ("self", "value_gt_num_or_nothing", 0.0),
        ],
        type=Floathing,
        json_type="number",
        unit="1/h"
    ),
    "h_to_r" => (
        default=2.0,
        description="Ratio of height to radius of the cylinder",
        display_name="Height-radius ratio",
        required=false,
        validations=[
            ("self", "value_gt_num", 0.0),
        ],
        type=Float64,
        json_type="number",
        unit="-"
    ),
    "thermal_transmission_lid" => (
        default=1.2,
        description="Thermal transmission rate of the lid",
        display_name="Thermal transmission lid",
        required=false,
        conditionals=[
            ("consider_losses", "is_true"),
        ],
        validations=[
            ("self", "value_gte_num", 0.0),
        ],
        type=Float64,
        json_type="number",
        unit="W/m^2*K"
    ),
    "thermal_transmission_barrel" => (
        default=1.2,
        description="Thermal transmission rate of the barrel",
        display_name="Thermal transmission barrel",
        required=false,
        conditionals=[
            ("consider_losses", "is_true"),
        ],
        validations=[
            ("self", "value_gte_num", 0.0),
        ],
        type=Float64,
        json_type="number",
        unit="W/m^2*K"
    ),
    "thermal_transmission_bottom" => (
        default=1.2,
        description="Thermal transmission rate of the bottom",
        display_name="Thermal transmission bottom",
        required=false,
        conditionals=[
            ("consider_losses", "is_true"),
        ],
        validations=[
            ("self", "value_gte_num", 0.0),
        ],
        type=Float64,
        json_type="number",
        unit="W/m^2*K"
    ),
    "ground_temperature" => (
        default=12.0,
        description="Constant temperature of the ground",
        display_name="Ground temperature",
        required=false,
        conditionals=[
            ("consider_losses", "is_true"),
        ],
        type=Float64,
        json_type="number",
        unit="°C"
    ),
    "switch_point" => (
        default=0.15,
        description="Switch point for the balanced model as relative load",
        display_name="Switch point",
        required=false,
        conditionals=[
            ("model_type", "has_value", "balanced"),
        ],
        validations=[
            ("self", "value_gt_num", 0.0),
            ("self", "value_lt_num", 1.0),
        ],
        type=Float64,
        json_type="number",
        unit="-"
    ),
)
#! format: on

"""
Implementation of a buffer tank holding hot water for heating or DHW purposes.

This is a simplified model of a buffer tank with three different model types:
* ideally_stratified: An ideally stratified cylindrical tank providing thermal energy always 
                      at the tank's upper temperature. Losses, if activated, are only 
                      energy losses reducing the amount of available energy (no exergy losses). 
* ideally_mixed:      An ideally mixed cylindrical tank providing thermal energy at a temperature
                      between the tank's upper and lower temperature. Note that the upper 
                      temperature is only supplied at 100% load. Losses, if activated, are
                      considered as energy and exergy losses reducing the energy and the 
                      current supply temperature.
* balanced:           The balanced model of the cylindrical buffer tank is a mix of the 
                      ideally stratified and ideally mixed model. At a load higher than the 
                      user-defined switch point, the ideally stratified model is used. At a 
                      load less than the switch point, the model switches to the ideally mixed model, 
                      representing a more realistic temperature behaviour of the energy supply.
Note that all three models can only be loaded with energy at a temperature of at least the 
upper temperature of the tank.

"""
mutable struct BufferTank <: Component
    uac::String
    controller::Controller
    sys_function::SystemFunction

    input_interfaces::InterfaceMap
    output_interfaces::InterfaceMap
    medium::Symbol

    model_type::Symbol

    capacity::Floathing
    volume::Floathing
    rho_medium::Float64
    cp_medium::Float64
    high_temperature::Float64
    low_temperature::Float64
    max_load_rate::Floathing
    max_unload_rate::Floathing
    # maximum input energy per time step [Wh]
    max_input_energy::Floathing
    # maximum output energy per time step [Wh]
    max_output_energy::Floathing
    consider_losses::Bool

    ## for losses
    h_to_r::Float64
    # surface of the lid and the bottom of the cylinder [m^2]
    surface_lid_bottom::Float64
    # surface of the barrel of the cylinder [m^2]
    surface_barrel::Float64
    thermal_transmission_lid::Float64
    thermal_transmission_barrel::Float64
    thermal_transmission_bottom::Float64
    ambient_temperature_profile::Union{Profile,Nothing}
    ambient_temperature::Temperature
    ground_temperature::Temperature

    # for model type "balanced"
    switch_point::Float64

    # current_max_output_temperature at the beginning of the time step
    current_max_output_temperature::Float64
    initial_load::Float64
    # current load, set to inital_load at the beginning [Wh]
    load::Float64
    # stores the load of the previous time step without losses
    load_end_of_last_timestep::Float64
    # losses in current time step [Wh]
    losses::Float64
    # bool indicating if the process step has already been performed in the current time step
    process_done::Bool
    # bool indicating if the load step has already been performed in the current time step
    load_done::Bool

    function BufferTank(uac::String, config::Dict{String,Any}, sim_params::Dict{String,Any})
        new(SSOT_parameter_constructor(BufferTank, uac, config, sim_params)...)
    end
end

function component_parameters(x::Type{BufferTank})::Dict{String,NamedTuple}
    return deepcopy(BUFFER_TANK_PARAMETERS) # return a copy to prevent external modification
end

function extract_parameter(x::Type{BufferTank}, config::Dict{String,Any}, param_name::String, param_def::NamedTuple,
                           sim_params::Dict{String,Any}, uac::String)
    if param_name == "constant_ambient_temperature" || param_name == "ambient_temperature_profile_file_path"
        constant_temperature,
        temperature_profile = get_parameter_profile_from_config(config,
                                                                sim_params,
                                                                "ambient_emperature",
                                                                "ambient_emperature_profile_file_path",
                                                                "ambient_emperature_from_global_file",
                                                                "constant_ambient_temperature",
                                                                uac)
        return param_name == "constant_ambient_temperature" ? constant_temperature : temperature_profile
    end

    return extract_parameter(Component, config, param_name, param_def, sim_params, uac)
end

function validate_config(x::Type{BufferTank}, config::Dict{String,Any}, extracted::Dict{String,Any}, uac::String,
                         sim_params::Dict{String,Any})
    validate_config(Component, extracted, uac, sim_params, component_parameters(BufferTank))
end

function init_from_params(x::Type{BufferTank}, uac::String, params::Dict{String,Any},
                          raw_params::Dict{String,Any}, sim_params::Dict{String,Any})::Tuple
    medium = Symbol(params["medium"])
    register_media([medium])

    return (uac,
            Controller(params["control_parameters"]),
            sf_storage,
            InterfaceMap(medium => nothing),
            InterfaceMap(medium => nothing),
            medium,
            Symbol(params["model_type"]),
            params["capacity"],
            params["volume"],
            params["rho_medium"],
            params["cp_medium"],
            params["high_temperature"],
            params["low_temperature"],
            params["max_load_rate"],
            params["max_unload_rate"],
            nothing, # max_input_energy
            nothing, # max_output_energy
            params["consider_losses"],
            params["h_to_r"],
            0.0,     # surface_lid_bottom
            0.0,     # surface_barrel
            params["thermal_transmission_lid"],
            params["thermal_transmission_barrel"],
            params["thermal_transmission_bottom"],
            params["consider_losses"] ? params["ambient_temperature_profile"] : nothing,
            params["consider_losses"] ? params["constant_ambient_temperature"] : nothing,
            params["ground_temperature"],
            params["switch_point"],
            0.0,     # current_max_output_temperature
            params["initial_load"],
            0.0,     # load
            0.0,     # load_end_of_last_timestep,
            0.0,     # losses
            false,   # process_done
            false)   # load_done
end

function initialise!(unit::BufferTank, sim_params::Dict{String,Any})
    set_storage_transfer!(unit.input_interfaces[unit.medium],
                          unload_storages(unit.controller, unit.medium))
    set_storage_transfer!(unit.output_interfaces[unit.medium],
                          load_storages(unit.controller, unit.medium))

    # calculate volume and capacity
    if unit.capacity === nothing && unit.volume === nothing
        @error "For the buffer tank $(unit.uac), either a volume or a capacity has to be given, but none of them is given."
        throw(InputError())
    elseif unit.capacity === nothing && unit.volume !== nothing
        unit.capacity = unit.volume * unit.rho_medium / 3.6 * unit.cp_medium *
                        (unit.high_temperature - unit.low_temperature)             # [Wh]
    elseif unit.capacity !== nothing && unit.volume === nothing
        if unit.consider_losses
            unit.volume = unit.capacity /
                          (unit.rho_medium / 3.6 * unit.cp_medium * (unit.high_temperature - unit.low_temperature))  # [m^3]
        end
    else
        @error "For the buffer tank $(unit.uac), either a volume or a capacity has to be given, but both are given."
        throw(InputError())
    end

    # calculate maximum input and output energy
    if unit.max_load_rate === nothing
        unit.max_input_energy = Inf
    else
        unit.max_input_energy = unit.max_load_rate * unit.capacity * (sim_params["time_step_seconds"] / 60 / 60)  # [Wh per timestep]
    end
    if unit.max_unload_rate === nothing
        unit.max_output_energy = Inf
    else
        unit.max_output_energy = unit.max_unload_rate * unit.capacity * (sim_params["time_step_seconds"] / 60 / 60)  # [Wh per timestep]
    end

    # calculate surfaces of the buffer tank cylinder
    if unit.consider_losses
        radius = cbrt(unit.volume / (unit.h_to_r * pi))  # [m]
        height = radius * unit.h_to_r                    # [m]
        unit.surface_barrel = 2 * pi * radius * height   # [m^2]
        unit.surface_lid_bottom = radius^2 * pi          # [m^2]
    end

    # set initial state
    unit.load = unit.initial_load * unit.capacity
    unit.load_end_of_last_timestep = copy(unit.load)
end

function control(unit::BufferTank,
                 components::Grouping,
                 sim_params::Dict{String,Any})
    update(unit.controller)

    unit.current_max_output_temperature = temperature_at_load(unit)

    set_max_energy!(unit.input_interfaces[unit.medium], min(unit.capacity - unit.load, unit.max_input_energy),
                    unit.high_temperature, nothing)
    set_max_energy!(unit.output_interfaces[unit.medium], min(unit.load, unit.max_output_energy), nothing,
                    unit.current_max_output_temperature)

    if unit.ambient_temperature_profile !== nothing && unit.consider_losses
        unit.ambient_temperature = Profiles.value_at_time(unit.ambient_temperature_profile, sim_params)
    end
end

function temperature_at_load(unit::BufferTank)::Temperature
    if unit.model_type == :ideally_stratified
        # always high temperature is available
        return unit.high_temperature
    elseif unit.model_type == :balanced
        # When the storage is loaded above the switch_point, the high temperature is available. 
        # At loads below the switch_point, a linear course between high and low temperature can be supplied.
        partial_load = min(1.0, unit.load / (unit.capacity * unit.switch_point))
        return (unit.high_temperature - unit.low_temperature) * partial_load + unit.low_temperature
    elseif unit.model_type == :ideally_mixed
        # A linear course between high and low temperature is supplied, depending on the current load.
        return unit.low_temperature + (unit.load / unit.capacity) * (unit.high_temperature - unit.low_temperature)
    end
end

function calculate_losses!(unit::BufferTank, sim_params)
    if !unit.consider_losses
        unit.load_end_of_last_timestep = copy(unit.load)
        return
    end

    function calculate_energy_loss(unit, sim_params)
        if unit.load <= sim_params["epsilon"]
            # no gains or losses if storage is completely empty for ideally_stratified model
            return 0.0
        else
            barrel_surface = unit.surface_barrel * unit.load / unit.capacity
            temperature_difference = unit.high_temperature - unit.ambient_temperature

            return (unit.thermal_transmission_lid * unit.surface_lid_bottom * temperature_difference +
                    unit.thermal_transmission_barrel * barrel_surface * temperature_difference) *
                   sim_params["time_step_seconds"] / 60 / 60
        end
    end

    if unit.model_type == :ideally_stratified
        # losses only through layer with high temperature (energy losses)
        unit.losses = calculate_energy_loss(unit, sim_params)
    elseif unit.model_type == :balanced
        # losses above switch_point result in energy decrease, below switch_point in exergy losses
        if unit.load / unit.capacity > unit.switch_point
            unit.losses = calculate_energy_loss(unit, sim_params)
        else
            current_tank_temperature = unit.low_temperature +
                                       (unit.load / (unit.switch_point * unit.capacity)) *
                                       (unit.high_temperature - unit.low_temperature)
            temperature_difference_air = current_tank_temperature - unit.ambient_temperature
            barrel_surface = unit.surface_barrel * unit.switch_point

            unit.losses = (unit.thermal_transmission_lid * unit.surface_lid_bottom * temperature_difference_air +
                           unit.thermal_transmission_barrel * barrel_surface * temperature_difference_air) *
                          sim_params["time_step_seconds"] / 60 / 60
        end
    elseif unit.model_type == :ideally_mixed
        # losses are exergy losses through the whole storage
        current_tank_temperature = unit.low_temperature +
                                   (unit.load / unit.capacity) * (unit.high_temperature - unit.low_temperature)
        temperature_difference_air = current_tank_temperature - unit.ambient_temperature
        temperature_difference_ground = current_tank_temperature - unit.ground_temperature

        unit.losses = (unit.thermal_transmission_lid * unit.surface_lid_bottom * temperature_difference_air +
                       unit.thermal_transmission_bottom * unit.surface_lid_bottom * temperature_difference_ground +
                       unit.thermal_transmission_barrel * unit.surface_barrel * temperature_difference_air) *
                      sim_params["time_step_seconds"] / 60 / 60
    end

    # save load at the end of the current time step before applying the losses for the control modules
    unit.load_end_of_last_timestep = copy(unit.load)

    # update load of storage and limit losses to current load
    unit.losses = min(unit.losses, unit.load)
    unit.load = unit.load - unit.losses
end

function process(unit::BufferTank, sim_params::Dict{String,Any})
    outface = unit.output_interfaces[unit.medium]
    exchanges = balance_on(outface, outface.target)
    energy_demanded = balance(exchanges) + energy_potential(exchanges)

    # shortcut if there is no energy demanded
    if energy_demanded >= -sim_params["epsilon"]
        set_max_energy!(unit.output_interfaces[unit.medium], 0.0)
        handle_component_update!(unit, "process", sim_params)
        return
    end

    for exchange in exchanges
        demanded_on_interface = exchange.balance + exchange.energy_potential

        if demanded_on_interface >= -sim_params["epsilon"]
            continue
        end

        if (exchange.temperature_min !== nothing
            &&
            exchange.temperature_min > unit.current_max_output_temperature)
            # we can only supply energy at a temperature at or below the tank's current
            # output temperature
            continue
        end

        used_heat = min(abs(energy_demanded), abs(demanded_on_interface))

        if unit.load > used_heat
            unit.load -= used_heat
            add!(outface, used_heat, nothing, unit.current_max_output_temperature)
            energy_demanded += used_heat
        else
            add!(outface, unit.load, nothing, unit.current_max_output_temperature)
            energy_demanded += unit.load
            unit.load = 0.0
        end
    end

    handle_component_update!(unit, "process", sim_params)
end

function handle_component_update!(unit::BufferTank, step::String, sim_params::Dict{String,Any})
    if step == "process"
        unit.process_done = true
    elseif step == "load"
        unit.load_done = true
    end
    if unit.process_done && unit.load_done
        # update component
        calculate_losses!(unit, sim_params)
        # reset 
        unit.process_done = false
        unit.load_done = false
    end
end

function load(unit::BufferTank, sim_params::Dict{String,Any})
    inface = unit.input_interfaces[unit.medium]
    exchanges = balance_on(inface, inface.source)
    energy_available = balance(exchanges) + energy_potential(exchanges)

    # shortcut if there is no energy to be used
    if energy_available <= sim_params["epsilon"]
        handle_component_update!(unit, "load", sim_params)
        set_max_energy!(unit.input_interfaces[unit.medium], 0.0)
        return
    end

    for exchange in exchanges
        exchange_energy_available = exchange.balance + exchange.energy_potential

        if exchange_energy_available < sim_params["epsilon"]
            continue
        end

        if exchange.temperature_max !== nothing &&
           exchange.temperature_max < unit.high_temperature
            # we can only take in energy if it's at a higher/equal temperature than the
            # tank's upper limit for temperatures
            continue
        end

        used_heat = min(energy_available, exchange_energy_available)
        diff = unit.capacity - unit.load

        if diff > used_heat
            unit.load += used_heat
            sub!(inface, used_heat, unit.high_temperature, nothing)
            energy_available -= used_heat
        else
            unit.load = unit.capacity
            sub!(inface, diff, unit.high_temperature, nothing)
            energy_available -= diff
        end
    end

    handle_component_update!(unit, "load", sim_params)
end

function output_values(unit::BufferTank)::Vector{String}
    return [string(unit.medium) * ":IN",
            string(unit.medium) * ":OUT",
            "Load",
            "Load%",
            "Capacity",
            "LossesGains",
            "CurrentMaxOutTemp"]
end

function output_value(unit::BufferTank, key::OutputKey)::Float64
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
    elseif key.value_key == "CurrentMaxOutTemp"
        return unit.current_max_output_temperature
    end
    throw(KeyError(key.value_key))
end

export BufferTank
