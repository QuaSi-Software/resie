function strt_sm_demand_driven(parameters :: Dict{String, Any}) :: StateMachine
    return StateMachine()
end

STRT_SM_FUNCS["demand_driven"] = strt_sm_demand_driven
STRT_SM_PARAMS["demand_driven"] = Dict{String, Any}()
