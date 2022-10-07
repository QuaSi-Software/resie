const TIME_STEP = UInt(900)

include("energy_systems/base.jl")

using .EnergySystems

function print_system_state(systems :: Grouping, time :: Int)
    println("Time is ", time)
    for unit in each(systems)
        pprint(unit, time)
        print(" | ")
    end
    print("\n")
end

function reset_file(systems :: Grouping)
    open("./out.csv", "w") do file_handle
        write(file_handle, "Time")
        for unit in each(systems)
            for val in specific_values(unit, Int(0))
                write(file_handle, ";$(typeof(unit)) $(val[1])")
            end
        end
        write(file_handle, "\n")
    end
end

function write_to_file(systems :: Grouping, time :: Int)
    open("./out.csv", "a") do file_handle
        write(file_handle, "$time")
        for unit in each(systems)
            for val in specific_values(unit, time)
                write(file_handle, replace(
                    replace(";$(val[2])", "/" => ";"),
                    "." => ","
                ))
            end
        end
        write(file_handle, "\n")
    end
end

function run_simulation()
    systems = Grouping(
        "TST_01_HZG_01_GRI" => make_GridConnection(EnergySystems.m_c_g_natgas, true),
        "TST_01_ELT_01_GRI" => make_GridConnection(EnergySystems.m_e_ac_230v, true),
        "TST_01_ELT_01_GRO" => make_GridConnection(EnergySystems.m_e_ac_230v, false),
        "TST_01_HZG_01_BFT" => make_BufferTank(40000.0, 20000.0),
        "TST_01_HZG_01_CHP" => make_CHPP("Ensure storage", 12500.0),
        "TST_01_HZG_01_HTP" => make_HeatPump("Ensure storage", 20000.0, 3.0),
        "TST_01_ELT_01_PVP" => make_PVPlant(15000.0),
        "TST_01_ELT_01_BUS" => make_Bus(EnergySystems.m_e_ac_230v),
        "TST_01_HZG_01_BUS" => make_Bus(EnergySystems.m_h_w_60c),
        "TST_01_HZG_01_DEM" => make_Demand(EnergySystems.m_h_w_60c, 10000.0),
        "TST_01_ELT_01_DEM" => make_Demand(EnergySystems.m_e_ac_230v, 15000.0),
    )

    control_order = [
        "TST_01_HZG_01_GRI",
        "TST_01_ELT_01_GRI",
        "TST_01_ELT_01_GRO",
        "TST_01_HZG_01_BFT",
        "TST_01_HZG_01_CHP",
        "TST_01_HZG_01_HTP",
        "TST_01_ELT_01_PVP",
        "TST_01_ELT_01_BUS",
        "TST_01_HZG_01_BUS",
        "TST_01_HZG_01_DEM",
        "TST_01_ELT_01_DEM",
    ]

    production_order = [
        "TST_01_HZG_01_DEM",
        "TST_01_ELT_01_DEM",
        "TST_01_ELT_01_PVP",
        "TST_01_HZG_01_BUS",
        "TST_01_HZG_01_BFT",
        "TST_01_HZG_01_CHP",
        "TST_01_HZG_01_HTP",
        "TST_01_ELT_01_BUS",
        "TST_01_HZG_01_GRI",
        "TST_01_ELT_01_GRI",
        "TST_01_ELT_01_GRO",
    ]

    link_control_with(
        systems["TST_01_HZG_01_CHP"],
        Grouping("TST_01_HZG_01_BFT" => systems["TST_01_HZG_01_BFT"])
    )
    link_control_with(
        systems["TST_01_HZG_01_HTP"],
        Grouping("TST_01_HZG_01_BFT" => systems["TST_01_HZG_01_BFT"])
    )

    link_production_with(
        systems["TST_01_HZG_01_GRI"],
        Grouping("TST_01_HZG_01_CHP" => systems["TST_01_HZG_01_CHP"])
    )
    link_production_with(
        systems["TST_01_ELT_01_GRI"],
        Grouping("TST_01_ELT_01_BUS" => systems["TST_01_ELT_01_BUS"])
    )
    link_production_with(
        systems["TST_01_HZG_01_BFT"],
        Grouping("TST_01_HZG_01_BUS" => systems["TST_01_HZG_01_BUS"])
    )
    link_production_with(
        systems["TST_01_HZG_01_CHP"],
        Grouping(
            "TST_01_HZG_01_BUS" => systems["TST_01_HZG_01_BUS"],
            "TST_01_ELT_01_BUS" => systems["TST_01_ELT_01_BUS"]
        )
    )
    link_production_with(
        systems["TST_01_HZG_01_HTP"],
        Grouping("TST_01_HZG_01_BUS" => systems["TST_01_HZG_01_BUS"])
    )
    link_production_with(
        systems["TST_01_ELT_01_PVP"],
        Grouping("TST_01_ELT_01_BUS" => systems["TST_01_ELT_01_BUS"])
    )
    link_production_with(
        systems["TST_01_HZG_01_BUS"],
        Grouping(
            "TST_01_HZG_01_DEM" => systems["TST_01_HZG_01_DEM"],
            "TST_01_HZG_01_BFT" => systems["TST_01_HZG_01_BFT"]
        )
    )
    link_production_with(
        systems["TST_01_ELT_01_BUS"],
        Grouping(
            "TST_01_ELT_01_DEM" => systems["TST_01_ELT_01_DEM"],
            "TST_01_HZG_01_HTP" => systems["TST_01_HZG_01_HTP"],
            "TST_01_ELT_01_GRO" => systems["TST_01_ELT_01_GRO"]
        )
    )

    parameters = Dict{String, Any}(
        "time" => 0,
        "time_step_seconds" => TIME_STEP,
        "price_factor" => 0.5,
        "epsilon" => 1e-9
    )

    print_system_state(systems, parameters["time"])
    reset_file(systems)

    for i = 1:(96*7)
        # perform the simulation
        control(systems, control_order, parameters)
        produce(systems, production_order, parameters)

        warnings = check_balances(systems, parameters["epsilon"])
        if length(warnings) > 0
            print("Time is $(parameters["time"])\n")
            for (key, balance) in warnings
                print("Warning: Balance for system $key was not zero: $balance\n")
            end
        end

        # output
        write_to_file(systems, parameters["time"])

        # simulation update
        parameters["time"] += Int(TIME_STEP)
    end
end

run_simulation()
