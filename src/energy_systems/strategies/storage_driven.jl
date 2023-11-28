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
end

strt_desc_storage_driven = """Storage-driven
------------------------
Run a component depending on the state of a linked thermal buffer tank. In particular
this strategy enables processing when the tank has fallen below a certain threshold and
disables it when the tank has reached another threshold.
"""

OP_STRATS["storage_driven"] = OperationalStrategyType(
    name="storage_driven",
    description=strt_desc_storage_driven,
    sm_constructor=strt_sm_storage_driven,
    conditions=[
        "Buffer < X%",
        "Buffer >= X%",
        "Min run time"
    ],
    strategy_parameters=Dict{String,Any}(
        "low_threshold" => 0.2,
        "high_threshold" => 0.95,
        "name" => "storage_driven",
        "load_storages" => true,
        "unload_storages" => true,
        "operation_profile_path" => nothing,
        "consider_m_el_in" => true,
        "consider_m_el_out" => true,
        "consider_m_gas_in" => true,
        "consider_m_fuel_in" => true,
        "consider_m_h2_out" => true,
        "consider_m_o2_out" => true,
        "consider_m_heat_out" => true,
        "consider_m_heat_in" => true,
    ),
    required_components=EnSysRequirements()
)
