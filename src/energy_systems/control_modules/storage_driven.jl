
# strt_desc_storage_driven = """Storage-driven
# ------------------------
# Run a component depending on the state of a linked thermal buffer tank. In particular
# this strategy enables processing when the tank has fallen below a certain threshold and
# disables it when the tank has reached another threshold.
# """

mutable struct CM_StorageDriven <: ControlModule
    name::String
    parameters::Dict{String,Any}
    state_machine::StateMachine

    function CM_StorageDriven(
        parameters::Dict{String,Any},
        components::Grouping,
        sim_params::Dict{String,Any}
    )
        default_parameters=Dict{String,Any}(
            "name" => "storage_driven",
            "low_threshold" => 0.2,
            "high_treshold" => 0.95,
            "min_run_time" => 1800,
            "storage_uac" => nothing
        )
        params = Base.merge(default_parameters, parameters)

        if !(
            params["storage_uac"] !== nothing
            && params["storage_uac"] in keys(components)
            && components[params["storage_uac"]] isa StorageComponent
        )
            @error "Required storage component for control module storage_driven not given"
        end
        params["storage"] = components[params["storage_uac"]]

        state_machine = StateMachine(
            UInt(1), # state
            Dict{UInt,String}( # state_names
                1 => "Off",
                2 => "Load",
            ),
            Dict{UInt,TruthTable}( # transitions
                1 => TruthTable( # State: Off
                    conditions=[
                        function(state_machine)
                            return params["storage"].load <
                                params["storage"].capacity * params["low_threshold"]
                        end
                    ],
                    table_data=Dict{Tuple,UInt}(
                        (false,) => 1,
                        (true,) => 2,
                    )
                ),
                2 => TruthTable( # State: Load
                    conditions=[
                        function(state_machine)
                            return params["storage"].load >= params["storage"].capacity *
                                params["high_threshold"]
                        end,
                        function(state_machine)
                            return state_machine.time_in_state *
                                sim_params["time_step_seconds"] >= params["min_run_time"]
                        end
                    ],
                    table_data=Dict{Tuple,UInt}(
                        (true, true) => 1,
                        (false, true) => 2,
                        (true, false) => 2,
                        (false, false) => 2,
                    )
                ),
            )
        )

        return new("storage_driven", params, state_machine)
    end
end

function update(mod::CM_StorageDriven)
    move_state(mod.state_machine)
end

function upper_plr_limit(mod::CM_StorageDriven, sim_params::Dict{String,Any})
    return mod.state_machine.state == 2 ? 1.0 : 0.0
end