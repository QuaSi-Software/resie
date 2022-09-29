Base.@kwdef mutable struct CHPP <: ControlledSystem
    controller :: StateMachine
    inputs :: Dict{MediumCategory, ControlledSystem}
    outputs :: Dict{MediumCategory, ControlledSystem}
    accepted_inputs :: Vector{MediumCategory}
    accepted_outputs :: Vector{MediumCategory}

    last_produced_e :: Float64 = 0.0
    last_produced_h :: Float64 = 0.0

    power :: Float64
    electricity_fraction :: Float64 = 0.4
    min_power_fraction :: Float64
    min_run_time :: UInt = 1800
end

function make_CHPP(strategy :: String, power :: Float64) :: CHPP
    if strategy == "Ensure storage"
        return CHPP(
            StateMachine( # CHPP.controller
                state=UInt(1),
                state_names=Dict{UInt, String}(
                    1 => "Off",
                    2 => "Load",
                ),
                time_in_state=UInt(0),
                transitions=Dict{UInt, TruthTable}(
                    1 => TruthTable( # State: Off
                        conditions=[
                            Condition(
                                "Buffer < X%",
                                Dict{String, Any}(
                                    "percentage" => 0.2
                                )
                            ),
                        ],
                        table_data=Dict{Tuple, UInt}(
                            (false,) => 1,
                            (true,) => 2,
                        )
                    ),

                    2 => TruthTable( # State: Load
                        conditions=[
                            Condition(
                                "Buffer >= X%",
                                Dict{String, Any}(
                                    "percentage" => 0.9
                                )
                            ),
                            Condition(
                                "Min run time",
                                Dict{String, Any}()
                            ),
                            Condition(
                                "Would overfill thermal buffer",
                                Dict{String, Any}()
                            ),
                        ],
                        table_data=Dict{Tuple, UInt}(
                            (false, false, false) => 2,
                            (false, true, false) => 2,
                            (true, false, false) => 2,
                            (true, true, false) => 1,
                            (false, false, true) => 1,
                            (false, true, true) => 1,
                            (true, false, true) => 1,
                            (true, true, true) => 1,
                        )
                    ),
                )
            ),
            Dict{MediumCategory, ControlledSystem}(), # CHPP.inputs
            Dict{MediumCategory, ControlledSystem}(), # CHPP.outputs
            [m_c_g_natgas], # CHPP.accepted_inputs
            [m_h_w_60c, m_e_ac_230v], # CHPP.accepted_outputs
            0.0, # CHPP.last_produced_e
            0.0, # CHPP.last_produced_h
            power, # CHPP.power
            0.4, # CHPP.electricity_fraction
            0.2, # CHPP.min_power_fraction
            1800 # CHPP.min_run_time
        )
    else
        return CHPP(controller=StateMachine(), power=power)
    end
end

function specific_values(unit :: CHPP, time :: Int) :: Vector{Tuple}
    return [
        ("Production E", "$(unit.last_produced_e)"),
        ("Production H", "$(unit.last_produced_h)")
    ]
end

export CHPP, make_CHPP, specific_values