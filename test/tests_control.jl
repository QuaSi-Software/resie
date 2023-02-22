using Debugger
using Test
using Resie.EnergySystems

function setup_control_tests()
    systems_config = Dict{String, Any}(
        "TST_BT_01" => Dict{String, Any}(
            "type" => "BufferTank",
            "control_refs" => [],
            "production_refs" => [],
            "capacity" => 40000,
            "load" => 20000,
            "strategy" => Dict{String, Any}(
                "name" => "Default"
            )
        ),
        "TST_HP_01" => Dict{String, Any}(
            "type" => "HeatPump",
            "control_refs" => ["TST_BT_01"],
            "production_refs" => [],
            "strategy" => Dict{String, Any}(
                "name" => "storage_driven",
                "high_threshold" => 0.9,
                "low_threshold" => 0.2
            ),
            "power" => 20000,
            "fixed_cop" => 3.0
        ),
    )
    systems = Resie.load_systems(systems_config)

    simulation_params = Dict{String, Any}(
        "time" => 0,
        "time_step_seconds" => 900,
        "epsilon" => 1e-9
    )

    return systems, simulation_params
end

function test_move_state_default_strategy()
    systems, simulation_params = setup_control_tests()
    buffer_tank = systems["TST_BT_01"]

    @test buffer_tank.controller.state_machine.state == 1
    EnergySystems.move_state(buffer_tank, systems, simulation_params)
    @test buffer_tank.controller.state_machine.state == 1
end

function test_move_state_storage_driven_hp()
    systems, simulation_params = setup_control_tests()
    buffer_tank = systems["TST_BT_01"]
    heat_pump = systems["TST_HP_01"]

    buffer_tank.load = 20000
    heat_pump.controller.state_machine.state = 1
    EnergySystems.move_state(heat_pump, systems, simulation_params)
    @test heat_pump.controller.state_machine.state == 1
    buffer_tank.load = 0
    EnergySystems.move_state(heat_pump, systems, simulation_params)
    @test heat_pump.controller.state_machine.state == 2
end

@testset "control_tests" begin
    @testset "move_state_default_strategy" begin
        test_move_state_default_strategy()
    end

    @testset "move_state_storage_driven_hp" begin
        test_move_state_storage_driven_hp()
    end
end