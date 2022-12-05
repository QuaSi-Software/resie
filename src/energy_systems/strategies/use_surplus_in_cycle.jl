function strt_sm_use_surplus_in_cycle(parameters :: Dict{String, Any}) :: StateMachine
    StateMachine(
        state=UInt(1),
        state_names=Dict{UInt, String}(
            1 => "Load",
            2 => "Produce",
        ),
        time_in_state=UInt(0),
        transitions=Dict{UInt, TruthTable}(
            1 => TruthTable( # State: Load
                conditions=[
                    Condition(
                        "Buffer >= X%",
                        Dict{String, Any}(
                            "percentage" => parameters["percentage"]
                        )
                    ),
                    Condition(
                        "HP is running",
                        Dict{String, Any}()
                    )
                ],
                table_data=Dict{Tuple, UInt}(
                    (false,false) => 2,
                    (false,true) => 2,
                    (true,false) => 1,
                    (true,true) => 2,
                )
            ),

            2 => TruthTable( # State: Produce
                conditions=[
                    Condition(
                        "Buffer >= X%",
                        Dict{String, Any}(
                            "percentage" => parameters["percentage"]
                        )
                    ),
                    Condition(
                        "HP is running",
                        Dict{String, Any}()
                    )
                ],
                table_data=Dict{Tuple, UInt}(
                        (false,false) => 2,
                        (false,true) => 2,
                        (true,false) => 1,
                        (true,true) => 2,
                    )
            ),
        )
    )
end

STRT_SM_FUNCS["use_surplus_in_cycle"] = strt_sm_use_surplus_in_cycle
STRT_SM_PARAMS["use_surplus_in_cycle"] = Dict{String, Any}(
    "percentage" => 0.95,
)