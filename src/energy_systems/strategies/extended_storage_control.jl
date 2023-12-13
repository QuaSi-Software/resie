function strt_sm_extended_storage_control(cond_params::Dict{String,Any})::StateMachine
    return StateMachine()
end

strt_desc_extended_storage_control = """Extended storage control
------------------------
This strategy offers the possibility to decide if a storage is allowed to fill any connected storage 
in the energy system. While the connection matrix of each bus only defines the loading of directly
connected storages, the funcionality provided within this strategy is valid system-wide.
NOTE: This is only implemented in a downstream way with regards to the ouput energies. Controlling if a 
      storage can be loaded from other storages upstream is not possible reliably due to the architecture of The
      simulation model (balances in interfaces are writte from upstream components without the information
      of the source of this energy to decouple components from each other) and therefore not implemented. 
"""

OP_STRATS["extended_storage_control"] = OperationalStrategyType(
    name="extended_storage_control",
    description=strt_desc_extended_storage_control,
    sm_constructor=strt_sm_extended_storage_control,
    conditions=[],
    strategy_parameters=Dict{String,Any}(
        "name" => "extended_storage_control",
        "load_any_storage" => false
    ),
    required_components=EnSysRequirements(
        "receiver" => (Component, nothing)
    )
)
