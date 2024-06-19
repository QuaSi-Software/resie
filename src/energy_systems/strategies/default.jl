function strt_sm_default(cond_params::Dict{String,Any})::StateMachine
    return StateMachine()
end

strt_desc_default = """Default
------------------------
Default operational strategy that all components have if no other is specified. Supports
several mechanisms of controlling operation, including:
  * Setting a profile that defines when a component can operate. This only works components
    specifying how the profile is used.
  * Setting specific inputs or outputs to be assumed infinite when calculating the energies
    that can be utilised in a time step
  * Defining for any input or output if the energy supplied/requested over that interface is
    allowed to be used for un-/loading storage components
"""

OP_STRATS["default"] = OperationalStrategyType(
    name="default",
    description=strt_desc_default,
    sm_constructor=strt_sm_default,
    conditions=[],
    strategy_parameters=Dict{String,Any}(
        "name" => "default",
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
    required_components=EnSysRequirements()
)
