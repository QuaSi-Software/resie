function strt_sm_demand_driven(parameters::Dict{String,Any})::StateMachine
    return StateMachine()
end

strt_desc_demand_driven = """Demand-driven
------------------------
Ensures that the processing of a component is done after the demand from components that receive
energy from the controlled component has been calculated. This also works if a bus acts as
intermediary between components. As this strategy does not implement any actual control, it is
mostly used for ensuring the components are calculated in the right order as well as
indicating to the process implementation that a given demand must be met (as opposed to
strategy "supply_driven").
"""

OP_STRATS["demand_driven"] = OperationalStrategyType(
    name="demand_driven",
    description=strt_desc_demand_driven,
    sm_constructor=strt_sm_demand_driven,
    conditions=[],
    strategy_parameters=Dict{String,Any}(
        "name" => "demand_driven",
        "load_storages" => true,
        "unload_storages" => true,
        "operation_profile_path" => nothing,
        "m_el_in" => true,
        "m_el_out" => true,
        "m_gas_in" => true,
        "m_h2_out" => true,
        "m_o2_out" => true,
        "m_heat_out" => true,
        "m_heat_in" => true
    ),
    required_components=EnSysRequirements(
        "receiver" => (Component, nothing)
    )
)
