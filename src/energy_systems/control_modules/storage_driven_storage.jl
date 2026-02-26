"""
Control module for setting a storage to be allowed to be charged depending on it's own 
charge state.
In particular it switches to a state of allowing charging when the load of the storage falls 
below the lower threshold. The module stays in this state until the load has reached the 
upper threshold and the minimum run time has passed.
"""

mutable struct CM_StorageDrivenStorage <: ControlModule
    name::String
    parameters::Dict{String,Any}
    state_machine::StateMachine

    function CM_StorageDrivenStorage(parameters::Dict{String,Any},
                              components::Grouping,
                              sim_params::Dict{String,Any},
                              unit_uac::String)
        default_parameters = Dict{String,Any}(
            "name" => "storage_driven_storage",
            "low_threshold" => 0.2,
            "high_threshold" => 0.95,
            "min_run_time" => 0,
            "charge_is_allowed" => true
        )
        params = Base.merge(default_parameters, parameters)

        params["storage"] = components[unit_uac]

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
                                                                     end,
                                                                     function (state_machine)
                                                                        return params["charge_is_allowed"]
                                                                     end],
                                                         table_data=Dict{Tuple,UInt}(
                                                             (false,true) => 1,
                                                             (true, true) => 2,
                                                             (false,false) => 1,
                                                             (true, false) => 1,
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
                                                                     end,
                                                                     function (state_machine)
                                                                        return params["charge_is_allowed"]
                                                                     end],
                                                         table_data=Dict{Tuple,UInt}(
                                                             (true, true, true) => 1,
                                                             (false,true, true) => 2,
                                                             (true, false,true) => 2,
                                                             (false,false,true) => 2,
                                                             (true, true, false) => 1,
                                                             (false,true, false) => 1,
                                                             (true, false,false) => 1,
                                                             (false,false,false) => 1,
                                                         )),
                                     ))

        return new(params["name"], params, state_machine)
    end
end

# method for control module name on type-level
control_module_name(::Type{CM_StorageDrivenStorage})::String = "storage_driven_storage"

function has_method_for(mod::CM_StorageDrivenStorage, func::ControlModuleFunction)::Bool
    return func == cmf_charge_is_allowed
end

function update(mod::CM_StorageDrivenStorage)
    move_state(mod.state_machine)
end

function charge_is_allowed(mod::CM_StorageDrivenStorage, sim_params::Dict{String,Any})::Bool
    return mod.state_machine.state == 2 ? true : false
end
