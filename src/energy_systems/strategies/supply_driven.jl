function strt_sm_supply_driven(parameters :: Dict{String, Any}) :: StateMachine
    return StateMachine()
end

STRT_SM_FUNCS["supply_driven"] = strt_sm_supply_driven
STRT_SM_PARAMS["supply_driven"] = Dict{String, Any}()
