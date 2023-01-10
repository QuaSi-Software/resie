function load_condition_prototypes()
    CONDITION_PROTOTYPES["Buffer < X%"] = ConditionPrototype(
        "Buffer < X%", # name
        Dict{String, Any}( # parameters
            "percentage" => 0.5
        ),
        EnSysRequirements( # required_systems
            "buffer" => (BufferTank, nothing)
        )
    )

    CONDITION_PROTOTYPES["Buffer >= X%"] = ConditionPrototype(
        "Buffer >= X%", # name
        Dict{String, Any}( # parameters
            "percentage" => 0.5
        ),
        EnSysRequirements( # required_systems
            "buffer" => (BufferTank, nothing)
        )
    )

    CONDITION_PROTOTYPES["Min run time"] = ConditionPrototype(
        "Min run time", # name
        Dict{String, Any}(), # parameters
        EnSysRequirements() # required_systems
    )

    CONDITION_PROTOTYPES["Would overfill thermal buffer"] = ConditionPrototype(
        "Would overfill thermal buffer", # name
        Dict{String, Any}(), # parameters
        EnSysRequirements( # required_systems
            "buffer" => (BufferTank, nothing)
        )
    )

    CONDITION_PROTOTYPES["Little PV power"] = ConditionPrototype(
        "Little PV power", # name
        Dict{String, Any}( # parameters
            "threshold" => 1000
        ),
        EnSysRequirements( # required_systems
            "pv_plant" => (PVPlant, nothing)
        )
    )

    CONDITION_PROTOTYPES["Much PV power"] = ConditionPrototype(
        "Much PV power", # name
        Dict{String, Any}( # parameters
            "threshold" => 1000
        ),
        EnSysRequirements( # required_systems
            "pv_plant" => (PVPlant, nothing)
        )
    )

    CONDITION_PROTOTYPES["Sufficient charge"] = ConditionPrototype(
        "Sufficient charge", # name
        Dict{String, Any}( # parameters
            "threshold" => 0.2
        ),
        EnSysRequirements() # required_systems
    )

    CONDITION_PROTOTYPES["Insufficient charge"] = ConditionPrototype(
        "Insufficient charge", # name
        Dict{String, Any}( # parameters
            "threshold" => 0.05
        ),
        EnSysRequirements() # required_systems
    )

    CONDITION_PROTOTYPES["HP is running"] = ConditionPrototype(
        "HP is running", # name
        Dict{String, Any}(), # parameters
        EnSysRequirements( # required_systems
            "heat_pump" => (HeatPump, nothing)
        )
    )
end
