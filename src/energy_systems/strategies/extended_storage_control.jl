function strt_sm_extended_storage_control(parameters::Dict{String,Any})::StateMachine
    return StateMachine()
end

strt_desc_extended_storage_control = """Extended storage control
------------------------
This strategy offers the possibility to decide if a storage is allowed to either fill any connected
storage in the system framework or if a storage can be filled by any connected storage in the system framework.
While the connection matrix of each bus only defines the loading/unloading of directly connected
storages, the funcionality provided within this strategy is valid system-wide.
NOTE: If the properties of load_any_storage and unload_any_storage of different storages in a system don't 
      fit together, unexpected results can be retrieved. If a storage is allowed to load all other storages, 
      the restrictions of the other storages are ignored. If one of them is set to "unload_any_storage=false", 
      this can not be garanteed as is highly depends on the order of operation.
      To avoid this issue, it is recommended to only use the flag "load_any_storage" and leave "unload_any_storage"
      by it's default value. "load_any_storage" affects the produce() step of the storage and if set to true, 
      the storage potential is added to the energy that should be provided. This energy is then written to the balance
      of the interface, so the load() step of storages downstream of the first storage can not detect, from where the
      energy came from. The loading of only a selection of downstream storages from one upstream storage can not be 
      modeled reliably as this depends on the order of operation of the produce() and load() steps of all storages.
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
