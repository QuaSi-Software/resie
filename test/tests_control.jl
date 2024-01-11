using Debugger
using Test
using Resie
using Resie.EnergySystems

function setup_control_tests()
    components_config = Dict{String,Any}(
        "TST_BT_01" => Dict{String,Any}(
            "type" => "BufferTank",
            "control_refs" => [],
            "output_refs" => [],
            "capacity" => 40000,
            "load" => 20000,
            "strategy" => Dict{String,Any}(
                "name" => "Default"
            )
        ),
        "TST_HP_01" => Dict{String,Any}(
            "type" => "HeatPump",
            "control_refs" => ["TST_BT_01"],
            "output_refs" => [],
            "strategy" => Dict{String,Any}(
                "name" => "storage_driven",
                "high_threshold" => 0.9,
                "low_threshold" => 0.2
            ),
            "power_th" => 20000,
            "constant_cop" => 3.0
        ),
    )

    simulation_params = Dict{String,Any}(
        "time" => 0,
        "time_step_seconds" => 900,
        "epsilon" => 1e-9
    )

    components = Resie.load_components(components_config, simulation_params)


    return components, simulation_params
end

function test_move_state_default_strategy()
    components, simulation_params = setup_control_tests()
    buffer_tank = components["TST_BT_01"]

    @test buffer_tank.controller.state_machine.state == 1
    EnergySystems.move_state(buffer_tank, components, simulation_params)
    @test buffer_tank.controller.state_machine.state == 1
end

function test_move_state_storage_driven_hp()
    components, simulation_params = setup_control_tests()
    buffer_tank = components["TST_BT_01"]
    heat_pump = components["TST_HP_01"]

    buffer_tank.load = 20000
    heat_pump.controller.state_machine.state = 1
    EnergySystems.move_state(heat_pump, components, simulation_params)
    @test heat_pump.controller.state_machine.state == 1
    buffer_tank.load = 0
    EnergySystems.move_state(heat_pump, components, simulation_params)
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