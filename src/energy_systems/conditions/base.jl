function load_condition_prototypes()
    CONDITION_PROTOTYPES["Buffer < X%"] = ConditionPrototype(
        "Buffer < X%", # name
        Dict{String,Any}( # parameters
            "percentage" => 0.5
        ),
        EnSysRequirements( # required_components
            "buffer" => (BufferTank, nothing)
        ),
        function (condition, unit, simulation_parameters) # check_function
            return rel(condition, "buffer").load <
                   condition.parameters["percentage"] * rel(condition, "buffer").capacity
        end
    )

    CONDITION_PROTOTYPES["Buffer >= X%"] = ConditionPrototype(
        "Buffer >= X%", # name
        Dict{String,Any}( # parameters
            "percentage" => 0.5
        ),
        EnSysRequirements( # required_components
            "buffer" => (BufferTank, nothing)
        ),
        function (condition, unit, simulation_parameters) # check_function
            return rel(condition, "buffer").load >=
                   condition.parameters["percentage"] * rel(condition, "buffer").capacity
        end
    )

    CONDITION_PROTOTYPES["Min run time"] = ConditionPrototype(
        "Min run time", # name
        Dict{String,Any}(), # parameters
        EnSysRequirements(), # required_components
        function (condition, unit, simulation_parameters) # check_function
            return unit.controller.state_machine.time_in_state *
                   simulation_parameters["time_step_seconds"] >=
                   unit.min_run_time
        end
    )

    CONDITION_PROTOTYPES["Would overfill thermal buffer"] = ConditionPrototype(
        "Would overfill thermal buffer", # name
        Dict{String,Any}(), # parameters
        EnSysRequirements( # required_components
            "buffer" => (BufferTank, nothing)
        ),
        function (condition, unit, simulation_parameters) # check_function
            return rel(condition, "buffer").capacity - rel(condition, "buffer").load <
                   unit.power * unit.min_power_fraction
        end
    )

    CONDITION_PROTOTYPES["Little PV power"] = ConditionPrototype(
        "Little PV power", # name
        Dict{String,Any}( # parameters
            "threshold" => 1000
        ),
        EnSysRequirements( # required_components
            "pv_plant" => (PVPlant, nothing)
        ),
        function (condition, unit, simulation_parameters) # check_function
            outface = rel(condition, "pv_plant").output_interfaces[:m_e_ac_230v]
            return (outface.balance != 0.0 ?
                    outface.sum_abs_change :
                    outface.sum_abs_change * 0.5
            ) <
                   condition.parameters["threshold"] * rel(condition, "pv_plant").supply * 0.25
        end
    )

    CONDITION_PROTOTYPES["Sufficient charge"] = ConditionPrototype(
        "Sufficient charge", # name
        Dict{String,Any}( # parameters
            "threshold" => 0.2
        ),
        EnSysRequirements(), # required_components
        function (condition, unit, simulation_parameters) # check_function
            return unit.load >= condition.parameters["threshold"] * unit.capacity
        end
    )

    CONDITION_PROTOTYPES["HP is running"] = ConditionPrototype(
        "HP is running", # name
        Dict{String,Any}(), # parameters
        EnSysRequirements( # required_components
            "heat_pump" => (HeatPump, nothing)
        ),
        function (condition, unit, simulation_parameters) # check_function
            return rel(condition, "heat_pump").output_interfaces[:m_h_w_ht1].sum_abs_change >
                   simulation_parameters["epsilon"]
        end
    )
end
