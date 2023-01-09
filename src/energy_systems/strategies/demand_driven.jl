function strt_sm_demand_driven(parameters :: Dict{String, Any}) :: StateMachine
    return StateMachine()
end

strt_desc_demand_driven = """Demand-driven
------------------------
Ensures that the production of a system is done after the demand from systems that receive
energy from the controlled system has been calculated. This also works if a bus acts as
intermediary between systems. As this strategy does not implement any actual control, it is
mostly used for ensuring the systems are calculated in the right order as well as
indicating to the production implementation that a given demand must be met (as opposed to
strategy "supply_driven").
"""

OP_STRATS["demand_driven"] = OperationalStrategyType(
    name="demand_driven",
    description=strt_desc_demand_driven,
    sm_constructor=strt_sm_demand_driven,
    conditions=[],
    strategy_parameters=Dict{String, Any}(),
    required_systems=EnSysRequirements()
)
