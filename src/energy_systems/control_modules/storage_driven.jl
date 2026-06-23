#! format: off
CONMOD_STORAGE_DRIVEN_PARAMS = Dict(
    "storage_uac" => (
        default=nothing,
        description="UAC of the storage component.",
        display_name="Storage UAC",
        required=true,
        type=String,
        json_type="string",
        unit="-"
    ),
    "low_threshold" => (
        default=0.2,
        description="Lower relative storage level at which the storage should begin to " *
                    "be filled.",
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
        description="Upper relative storage level at which the storage should stop " *
                    "being filled.",
        display_name="Lower threshold",
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
        description="Minimum run time of the hysteresis in the 'on' condition. This will " *
                    "override the upper threshold, but will in turn be overridden if the " *
                    "storage is completely full. Should ideally be a multiple of the " *
                    "simulation time step, as fractional time steps are not considered.",
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
Control module for running a component depending on the state of a linked storage component.
In particular it switches to a state of allowing operation of the component when the load
of the linked storage falls below the lower threshold. The module stays in this state until
the load has reached the upper threshold and the minimum run time has passed.
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
            "low_threshold" => 0.2,
            "high_threshold" => 0.95,
            "min_run_time" => 1800,
            "storage_uac" => nothing,
        )
        params = Base.merge(default_parameters, parameters)

        if !(params["storage_uac"] !== nothing
             && params["storage_uac"] in keys(components)
             && components[params["storage_uac"]] isa StorageComponent)
            @error "Required storage component for control module storage_driven not given"
        end
        params["storage"] = components[params["storage_uac"]]

        state_machine = StateMachine(UInt(1),           # state
                                     Dict{UInt,String}( # state_names
                                         1 => "Off",
                                         2 => "Load",
                                     ),
                                     Dict{UInt,TruthTable}( # transitions
                                         1 => TruthTable(;  # State: Off
                                                         conditions=[function (state_machine)
                                                                         return params["storage"].load_end_of_last_timestep <
                                                                                params["storage"].capacity *
                                                                                params["low_threshold"]
                                                                     end],
                                                         table_data=Dict{Tuple,UInt}(
                                                             (false,) => 1,
                                                             (true,) => 2,
                                                         )),
                                         2 => TruthTable(;  # State: Load
                                                         conditions=[function (state_machine)
                                                                         return params["storage"].load_end_of_last_timestep >=
                                                                                params["storage"].capacity *
                                                                                params["high_threshold"]
                                                                     end,
                                                                     function (state_machine)
                                                                         return state_machine.time_in_state *
                                                                                sim_params["time_step_seconds"] >=
                                                                                params["min_run_time"]
                                                                     end],
                                                         table_data=Dict{Tuple,UInt}(
                                                             (true, true) => 1,
                                                             (false, true) => 2,
                                                             (true, false) => 2,
                                                             (false, false) => 2,
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
