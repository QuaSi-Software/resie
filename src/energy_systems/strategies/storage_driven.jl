function strt_sm_storage_driven(parameters::Dict{String,Any})::StateMachine
    return StateMachine(
        UInt(1), # state
        Dict{UInt,String}( # state_names
            1 => "Off",
            2 => "Load",
        ),
        Dict{UInt,TruthTable}( # transitions
            1 => TruthTable( # State: Off
                conditions=[
                    Condition(
                        "Buffer < X%",
                        Dict{String,Any}(
                            "percentage" => parameters["low_threshold"]
                        )
                    ),
                ],
                table_data=Dict{Tuple,UInt}(
                    (false,) => 1,
                    (true,) => 2,
                )
            ), 2 => TruthTable( # State: Load
                conditions=[
                    Condition(
                        "Buffer >= X%",
                        Dict{String,Any}(
                            "percentage" => parameters["high_threshold"]
                        )
                    ),
                    Condition(
                        "Min run time",
                        Dict{String,Any}()
                    ),
                    Condition(
                        "Would overfill thermal buffer",
                        Dict{String,Any}()
                    ),
                ],
                table_data=Dict{Tuple,UInt}(
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

strt_desc_storage_driven = """Storage-driven
------------------------
Run an energy system depending on the state of a linked thermal buffer tank. In particular
this strategy enables production when the tank has fallen below a certain threshold and
disables it when the tank has reached another threshold.
"""

OP_STRATS["storage_driven"] = OperationalStrategyType(
    name="storage_driven",
    description=strt_desc_storage_driven,
    sm_constructor=strt_sm_storage_driven,
    conditions=[
        "Buffer < X%",
        "Buffer >= X%",
        "Min run time",
        "Would overfill thermal buffer"
    ],
    strategy_parameters=Dict{String,Any}(
        "low_threshold" => 0.2,
        "high_threshold" => 0.95,
        "name" => "storage_driven",
        "load_storages" => true,
        "unload_storages" => true,
        "m_el_in" => true,
        "m_el_out" => true,
        "m_gas_in" => true,
        "m_h2_out" => true,
        "m_o2_out" => true,
        "m_heat_out" => true,
        "m_heat_in" => true,
    ),
    required_systems=EnSysRequirements()
)
