function strt_sm_supply_driven(parameters :: Dict{String, Any}) :: StateMachine
    return StateMachine()
end

strt_desc_supply_driven = """Supply-driven
------------------------
Ensures that the production of a system is done after the production from systems that
provide energy to the controlled system has been calculated. This also works if a bus acts as
intermediary between systems. As this strategy does not implement any actual control, it is
mostly used for ensuring the systems are calculated in the right order as well as
indicating to the production implementation that a given supply must be used (as opposed to
strategy "demand_driven").
"""

OP_STRATS["supply_driven"] = OperationalStrategyType(
    name="supply_driven",
    description=strt_desc_supply_driven,
    sm_constructor=strt_sm_supply_driven,
    conditions=[],
    strategy_parameters=Dict{String, Any}(),
    required_systems=EnSysRequirements()
)

STRT_SM_FUNCS["supply_driven"] = strt_sm_supply_driven
STRT_SM_PARAMS["supply_driven"] = Dict{String, Any}()
