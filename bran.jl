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
        "TST_01_HZG_01_GRD" => GridConnection(medium=EnergySystems.m_c_g_natgas),
        "TST_01_ELT_01_GRD" => GridConnection(medium=EnergySystems.m_e_ac_230v),
        "TST_01_HZG_01_BFT" => BufferTank(capacity=40000.0, load=20000.0),
        "TST_01_HZG_01_CHP" => make_CHPP("Ensure storage", 12500.0),
        "TST_01_HZG_01_HTP" => make_HeatPump("Ensure storage", 20000.0, 3.0),
        "TST_01_ELT_01_PVP" => PVPlant(amplitude=15000.0),
        "TST_01_ELT_01_BUS" => Bus(medium=EnergySystems.m_e_ac_230v),
        "TST_01_HZG_01_DEM" => Demand(medium=EnergySystems.m_h_w_60c, load=10000),
        "TST_01_ELT_01_DEM" => Demand(medium=EnergySystems.m_e_ac_230v, load=15000),
    )

    link_with(systems["TST_01_HZG_01_CHP"], Grouping("TST_01_HZG_01_BFT" => systems["TST_01_HZG_01_BFT"]))
    link_with(systems["TST_01_HZG_01_HTP"], Grouping("TST_01_HZG_01_BFT" => systems["TST_01_HZG_01_BFT"]))

    parameters = Dict{String, Any}(
        "time" => 0,
        "time_step_seconds" => TIME_STEP,
        "price_factor" => 0.5
    )

    print_system_state(systems, parameters["time"])
    reset_file(systems)

    for i = 1:(96*7)
        for unit in each(systems)
            control(unit, systems, parameters)
        end

        # production
        produce(systems, parameters)

        # output and simulation update
        # print_system_state(system, parameters["time"])
        write_to_file(systems, parameters["time"])
        parameters["time"] += Int(TIME_STEP)
    end
end

run_simulation()
