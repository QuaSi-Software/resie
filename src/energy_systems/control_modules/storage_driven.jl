#! format: off
CONMOD_STORAGE_DRIVEN_PARAMS = Dict(
    "storage_uac" => (
        default=nothing,
        description="UAC of the storage component whose level controls this component.",
        display_name="Storage UAC",
        required=true,
        type=String,
        json_type="string",
        unit="-"
    ),
    "control_mode" => (
        default="charge",
        description="How the controlled component interacts with the monitored storage. " *
                    "Use 'charge' when the component fills the storage: it turns on below " *
                    "the lower threshold and off above the upper threshold. Use 'discharge' " *
                    "when the component consumes energy from the storage: it turns on above " *
                    "the upper threshold and off below the lower threshold.",
        display_name="Storage control mode",
        required=false,
        options=["charge", "discharge"],
        type=String,
        json_type="string",
        unit="-"
    ),
    "low_threshold" => (
        default=0.2,
        description="Lower relative storage-level boundary. In 'charge' mode, the component " *
                    "turns on below this value. In 'discharge' mode, it turns off at or below " *
                    "this value.",
        display_name="Lower threshold",
        required=false,
        validations=[
            ("self", "value_gte_num", 0.0),
            ("self", "value_lt_num", 1.0),
        ],
        type=Float64,
        json_type="number",
        unit="-"
    ),
    "high_threshold" => (
        default=0.95,
        description="Upper relative storage-level boundary. In 'charge' mode, the component " *
                    "turns off at or above this value. In 'discharge' mode, it turns on above " *
                    "this value.",
        display_name="Upper threshold",
        required=false,
        validations=[
            ("self", "value_gt_num", 0.0),
            ("self", "value_lte_num", 1.0),
            ("self", "value_gt_rel", "low_threshold"),
        ],
        type=Float64,
        json_type="number",
        unit="-"
    ),
    "min_run_time" => (
        default=1800,
        description="Minimum time for which the component remains on after starting. The " *
                    "normal turn-off threshold is ignored until this time has passed. A full " *
                    "storage in 'charge' mode or an empty storage in 'discharge' mode still " *
                    "stops the component immediately. This value should ideally be a multiple " *
                    "of the simulation time step.",
        display_name="Min. run time",
        required=false,
        validations=[
            ("self", "value_gte_num", 0.0),
        ],
        type=Float64,
        json_type="number",
        unit="s"
    ),
)
#! format: on

"""
Control module for enabling a component according to the level of one linked
storage component.

The module uses two storage-level boundaries and one control mode:

- `charge`: enable below `low_threshold` and disable at or above
  `high_threshold`. Use this when the controlled component fills the storage.
- `discharge`: enable above `high_threshold` and disable at or below
  `low_threshold`. Use this when the controlled component consumes energy
  from the storage.

The normal turn-off condition is applied only after `min_run_time`. A physical
storage limit overrides the minimum run time: full storage in `charge` mode and
empty storage in `discharge` mode.
"""
mutable struct CM_StorageDriven <: ControlModule
    name::String
    parameters::Dict{String,Any}
    state_machine::StateMachine

    function CM_StorageDriven(parameters::Dict{String,Any},
                              components::Grouping,
                              sim_params::Dict{String,Any},
                              unit_uac::String)
        default_parameters = Dict{String,Any}(
            "name" => "storage_driven",
            "control_mode" => "charge",
            "low_threshold" => 0.2,
            "high_threshold" => 0.95,
            "min_run_time" => 1800,
            "storage_uac" => nothing,
        )
        params = Base.merge(default_parameters, parameters)

        if !(params["storage_uac"] !== nothing
             && params["storage_uac"] in keys(components)
             && components[params["storage_uac"]] isa StorageComponent)
            @error "Required storage component `$(params["storage_uac"])` for control module storage_driven not given"
            throw(InputError())
        end
        params["storage"] = components[params["storage_uac"]]

        function should_turn_on(_)
            storage_level = params["storage"].load_end_of_last_timestep / params["storage"].capacity
            if params["control_mode"] == "charge"
                return storage_level < params["low_threshold"]
            else
                return storage_level > params["high_threshold"]
            end
        end
        turn_on_condition = should_turn_on

        function should_turn_off(state_machine)
            minimum_run_time_reached = state_machine.time_in_state * sim_params["time_step_seconds"] >=
                                       params["min_run_time"]

            if params["control_mode"] == "charge"
                storage_is_full = params["storage"].load_end_of_last_timestep >= params["storage"].capacity
                high_threshold_reached = params["storage"].load_end_of_last_timestep >=
                                         params["storage"].capacity * params["high_threshold"]

                return storage_is_full || (high_threshold_reached && minimum_run_time_reached)
            else
                storage_is_empty = params["storage"].load_end_of_last_timestep <= 0.0
                low_threshold_reached = params["storage"].load_end_of_last_timestep <=
                                        params["storage"].capacity * params["low_threshold"]

                return storage_is_empty || (low_threshold_reached && minimum_run_time_reached)
            end
        end
        turn_off_condition = should_turn_off

        state_machine = StateMachine(UInt(1),
                                     Dict{UInt,String}(
                                         1 => "Off",
                                         2 => "On",
                                     ),
                                     Dict{UInt,TruthTable}(
                                         1 => TruthTable(; conditions=[turn_on_condition],
                                                         table_data=Dict{Tuple,UInt}(
                                                                                     (false,) => 1,
                                                                                     (true,) => 2,
                                                    )),
                                         2 => TruthTable(; conditions=[turn_off_condition],
                                                         table_data=Dict{Tuple,UInt}(
                                                                                     (false,) => 2,
                                                                                     (true,) => 1,
                                                    )),
                                     ))

        return new("storage_driven", params, state_machine)
    end
end

# method for control module name on type-level
control_module_name(::Type{CM_StorageDriven})::String = "storage_driven"

# method for parameter definitions on type-level
control_module_parameters(x::Type{CM_StorageDriven})::Dict{String,NamedTuple} = CONMOD_STORAGE_DRIVEN_PARAMS

function has_method_for(mod::CM_StorageDriven, func::ControlModuleFunction)::Bool
    return func == cmf_upper_plr_limit
end

function update(mod::CM_StorageDriven)
    move_state(mod.state_machine)
end

function upper_plr_limit(mod::CM_StorageDriven, sim_params::Dict{String,Any})::Float64
    return mod.state_machine.state == 2 ? 1.0 : 0.0
end
