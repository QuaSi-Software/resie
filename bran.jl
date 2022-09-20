const TIME_STEP = UInt(900)

include("energy_systems/base.jl")

using .EnergySystems

function print_system_state(system :: Vector{ControlledSystem}, time :: Int)
    println("Time is ", time)
    for unit in system
        pprint(unit, time)
        print(" | ")
    end
    print("\n")
end

function reset_file(system :: Vector{ControlledSystem})
    open("./out.csv", "w") do file_handle
        write(file_handle, "Time")
        for unit in system
            for val in specific_values(unit, Int(0))
                write(file_handle, ";$(typeof(unit)) $(val[1])")
            end
        end
        write(file_handle, "\n")
    end
end

function write_to_file(system :: Vector{ControlledSystem}, time :: Int)
    open("./out.csv", "a") do file_handle
        write(file_handle, "$time")
        for unit in system
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
    buffer = BufferTank(capacity=40000.0, load=20000.0)
    system = [
        GridConnection(medium=EnergySystems.m_c_g_natgas),
        GridConnection(medium=EnergySystems.m_e_ac_230v),
        buffer,
        make_CHPP("Ensure storage", 12500.0, buffer),
        make_HeatPump("Ensure storage", 20000.0, 3.0, buffer),
        PVPlant(amplitude=15000.0),
        Bus(medium=EnergySystems.m_e_ac_230v),
        Demand(medium=EnergySystems.m_h_w_60c, load=10000),
        Demand(medium=EnergySystems.m_e_ac_230v, load=15000),
    ]
    parameters = Dict{String, Any}(
        "time" => 0,
        "time_step_seconds" => TIME_STEP,
        "price_factor" => 0.5
    )

    print_system_state(system, parameters["time"])
    reset_file(system)

    for i = 1:(96*7)
        for unit in system
            # control
            move_state(unit, system, parameters)
        end

        # production
        produce(system, parameters)

        # output and simulation update
        # print_system_state(system, parameters["time"])
        write_to_file(system, parameters["time"])
        parameters["time"] += Int(TIME_STEP)
    end
end

run_simulation()
