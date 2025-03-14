function setup_control_tests()
    components_config = Dict{String,Any}(
        "TST_GRI_01" => Dict{String,Any}(
            "type" => "GridConnection",
            "output_refs" => ["TST_HP_01"],
            "is_source" => true,
            "medium" => "m_h_w_lt1",
        ),
        "TST_GRI_02" => Dict{String,Any}(
            "type" => "GridConnection",
            "output_refs" => ["TST_HP_01"],
            "is_source" => true,
            "medium" => "m_e_ac_230v",
        ),
        "TST_HP_01" => Dict{String,Any}(
            "type" => "HeatPump",
            "output_refs" => ["TST_BT_01"],
            "control_modules" => [Dict{String,Any}(
                                      "name" => "storage_driven",
                                      "high_threshold" => 0.9,
                                      "low_threshold" => 0.2,
                                      "storage_uac" => "TST_BT_01",
                                  )],
            "power_th" => 20000,
            "cop_function" => "const:3.0",
            "power_losses_factor" => 1.0,
            "heat_losses_factor" => 1.0,
        ),
        "TST_BT_01" => Dict{String,Any}(
            "type" => "BufferTank",
            "output_refs" => ["TST_DEM_01"],
            "medium" => "m_h_w_ht1",
            "model_type" => "ideally_stratified",
            "capacity" => 40000,
            "initial_load" => 0.5,
        ),
        "TST_DEM_01" => Dict{String,Any}(
            "type" => "Demand",
            "output_refs" => [],
            "medium" => "m_h_w_ht1",
            "constant_demand" => 5000,
            "constant_temperature" => 60,
        ),
    )

    simulation_params = get_default_sim_params()

    components = Resie.load_components(components_config, simulation_params)

    return components, simulation_params
end

function test_move_state_storage_driven_hp()
    components, simulation_params = setup_control_tests()
    buffer_tank = components["TST_BT_01"]
    heat_pump = components["TST_HP_01"]

    buffer_tank.load_end_of_last_timestep = 20000
    EnergySystems.update(heat_pump.controller)
    @test heat_pump.controller.modules[1].state_machine.state == 1
    buffer_tank.load_end_of_last_timestep = 0
    EnergySystems.update(heat_pump.controller)
    @test heat_pump.controller.modules[1].state_machine.state == 2
end

@testset "move_state_storage_driven_hp" begin
    test_move_state_storage_driven_hp()
end
