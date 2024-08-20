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
                              sim_params::Dict{String,Any})
        default_parameters = Dict{String,Any}(
            "name" => "storage_driven",
            "low_threshold" => 0.2,
            "high_treshold" => 0.95,
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
                                                                         return params["storage"].load <
                                                                                params["storage"].capacity *
                                                                                params["low_threshold"]
                                                                     end],
                                                         table_data=Dict{Tuple,UInt}(
                                                             (false,) => 1,
                                                             (true,) => 2,
                                                         )),
                                         2 => TruthTable(;  # State: Load
                                                         conditions=[function (state_machine)
                                                                         return params["storage"].load >=
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

function has_method_for(mod::CM_StorageDriven, func::ControlModuleFunction)::Bool
    return func == cmf_upper_plr_limit
end

function update(mod::CM_StorageDriven)
    move_state(mod.state_machine)
end

function upper_plr_limit(mod::CM_StorageDriven, sim_params::Dict{String,Any})::Float64
    return mod.state_machine.state == 2 ? 1.0 : 0.0
end
