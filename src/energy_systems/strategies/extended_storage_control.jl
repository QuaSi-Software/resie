function strt_sm_extended_storage_control(parameters::Dict{String,Any})::StateMachine
    return StateMachine()
end

strt_desc_extended_storage_control = """Extended storage control
------------------------
This strategy offers the possibility to decide if a storage is allowed to either fill any connected
storage in the system framework or if a storage can be filled by any connected storage in the system framework.
While the connection matrix of each bus only defines the loading/unloading of directly connected
storages, the funcionality provided within this strategy is valid system-wide.
"""

OP_STRATS["extended_storage_control"] = OperationalStrategyType(
    name="extended_storage_control",
    description=strt_desc_extended_storage_control,
    sm_constructor=strt_sm_extended_storage_control,
    conditions=[],
    strategy_parameters=Dict{String,Any}(
        "name" => "extended_storage_control",
        "load_any_storage" => false,
        "unload_any_storage" => false
    ),
    required_systems=EnSysRequirements(
        "receiver" => (EnergySystem, nothing)
    )
)
