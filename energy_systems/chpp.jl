Base.@kwdef mutable struct CHPP <: ControlledSystem
    controller :: StateMachine

    input_interfaces :: Dict{MediumCategory, SystemInterface}
    output_interfaces :: Dict{MediumCategory, SystemInterface}

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
            Dict{MediumCategory, SystemInterface}( # input_interfaces
                m_c_g_natgas => SystemInterface()
            ),
            Dict{MediumCategory, SystemInterface}( # output_interfaces
                m_h_w_60c => SystemInterface(),
                m_e_ac_230v => SystemInterface()
            ),
            0.0, # last_produced_e
            0.0, # last_produced_h
            power, # power
            0.4, # electricity_fraction
            0.2, # min_power_fraction
            1800 # min_run_time
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