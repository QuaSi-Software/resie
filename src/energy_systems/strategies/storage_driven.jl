function strt_sm_storage_driven(parameters :: Dict{String, Any}) :: StateMachine
    return StateMachine(
        state=UInt(1),
        state_names=Dict{UInt, String}(
            1 => "Off",
            2 => "Load",
        ),
        time_in_state=UInt(0),
        transitions=Dict{UInt, TruthTable}(
            1 => TruthTable( # State: Off
                conditions=[
                    Condition(
                        "Buffer < X%",
                        Dict{String, Any}(
                            "percentage" => parameters["low_threshold"]
                        )
                    ),
                ],
                table_data=Dict{Tuple, UInt}(
                    (false,) => 1,
                    (true,) => 2,
                )
            ),

            2 => TruthTable( # State: Load
                conditions=[
                    Condition(
                        "Buffer >= X%",
                        Dict{String, Any}(
                            "percentage" => parameters["high_threshold"]
                        )
                    ),
                    Condition(
                        "Min run time",
                        Dict{String, Any}()
                    ),
                    Condition(
                        "Would overfill thermal buffer",
                        Dict{String, Any}()
                    ),
                ],
                table_data=Dict{Tuple, UInt}(
                    (false, false, false) => 2,
                    (false, true, false) => 2,
                    (true, false, false) => 2,
                    (true, true, false) => 1,
                    (false, false, true) => 1,
                    (false, true, true) => 1,
                    (true, false, true) => 1,
                    (true, true, true) => 1,
                )
            ),
        )
    )
end

STRT_SM_FUNCS["storage_driven"] = strt_sm_storage_driven
STRT_SM_PARAMS["storage_driven"] = Dict{String, Any}(
    "low_threshold" => 0.2,
    "high_threshold" => 0.95,
)