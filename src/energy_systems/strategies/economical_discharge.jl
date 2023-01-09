function strt_sm_economical_discharge(parameters :: Dict{String, Any}) :: StateMachine
    return StateMachine(
        UInt(1), # state
        Dict{UInt, String}( # state_names
            1 => "Load",
            2 => "Discharge",
        ),
        Dict{UInt, TruthTable}( # transitions
            1 => TruthTable( # State: Load
                conditions=[
                    Condition(
                        "Little PV power",
                        Dict{String, Any}(
                            "threshold" => parameters["pv_threshold"]
                        )
                    ),
                    Condition(
                        "Sufficient charge",
                        Dict{String, Any}(
                            "threshold" => parameters["min_charge"]
                        )
                    )
                ],
                table_data=Dict{Tuple, UInt}(
                    (false,false) => 1,
                    (false,true) => 1,
                    (true,false) => 1,
                    (true,true) => 2,
                )
            ),

            2 => TruthTable( # State: Discharge
                conditions=[
                    Condition(
                        "Little PV power",
                        Dict{String, Any}(
                            "threshold" => parameters["pv_threshold"]
                        )
                    ),
                    Condition(
                        "Sufficient charge",
                        Dict{String, Any}(
                            "threshold" => parameters["discharge_limit"]
                        )
                    )
                ],
                table_data=Dict{Tuple, UInt}(
                        (false,false) => 1,
                        (false,true) => 1,
                        (true,false) => 1,
                        (true,true) => 2,
                    )
            ),
        )
    )
end

strt_desc_economical_discharge = """Economical discharge
------------------------
Discharges a battery when a linked PV plant produces little power and the battery has
sufficient charge to provide substantial amounts of energy.
"""

OP_STRATS["economical_discharge"] = OperationalStrategyType(
    name="economical_discharge",
    description=strt_desc_economical_discharge,
    sm_constructor=strt_sm_economical_discharge,
    conditions=[
        "Little PV power",
        "Sufficient charge"
    ],
    strategy_parameters=Dict{String, Any}(
        "pv_threshold" => 0.15,
        "min_charge" => 0.2,
        "discharge_limit" => 0.05,
    ),
    required_systems=EnSysRequirements()
)
