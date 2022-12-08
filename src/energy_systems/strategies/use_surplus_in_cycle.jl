function strt_sm_use_surplus_in_cycle(parameters :: Dict{String, Any}) :: StateMachine
    StateMachine(
        UInt(1), # state
        Dict{UInt, String}( # state_names
            1 => "Load",
            2 => "Produce",
        ),
        Dict{UInt, TruthTable}( # transitions
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