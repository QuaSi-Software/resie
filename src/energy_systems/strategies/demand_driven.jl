function strt_sm_demand_driven(cond_params::Dict{String,Any})::StateMachine
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
        "operation_profile_path" => nothing,
        "consider_m_el_in" => true,
        "consider_m_el_out" => true,
        "consider_m_gas_in" => true,
        "consider_m_fuel_in" => true,
        "consider_m_h2_out" => true,
        "consider_m_o2_out" => true,
        "consider_m_heat_out" => true,
        "consider_m_heat_ht_out" => true,
        "consider_m_heat_lt_out" => true,
        "consider_m_heat_in" => true
    ),
    required_components=EnSysRequirements(
        "receiver" => (Component, nothing)
    )
)
