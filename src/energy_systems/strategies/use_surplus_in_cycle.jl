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

strt_desc_use_surplus_in_cycle = """Use surplus in cycle
------------------------
A special strategy used for ensuring that a long-term thermal storage system does not store
energy when a heat pump, which is fed by the LTTSS, is running at the same time. This would
make little sense as this cycle can only destroy energy and would never be implemented in
a real energy network.
"""

OP_STRATS["use_surplus_in_cycle"] = OperationalStrategyType(
    name="use_surplus_in_cycle",
    description=strt_desc_use_surplus_in_cycle,
    sm_constructor=strt_sm_use_surplus_in_cycle,
    conditions=[
        "Buffer >= X%",
        "HP is running",
        "Buffer >= X%",
    ],
    strategy_parameters=Dict{String, Any}(
        "percentage" => 0.95,
    ),
    required_systems=EnSysRequirements()
)
