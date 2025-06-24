using Test
using Resie
using Resie.EnergySystems

include("../test_util.jl")

function get_energy_system_fuel_boiler_direct()
    components_config = Dict{String,Any}(
        "TST_GRI_01" => Dict{String,Any}(
            "type" => "GridConnection",
            "output_refs" => ["TST_FB_01"],
            "is_source" => true,
            "medium" => "m_c_g_natgas",
        ),
        "TST_FB_01" => Dict{String,Any}(
            "type" => "FuelBoiler",
            "m_fuel_in" => "m_c_g_natgas",
            "output_refs" => ["TST_DEM_01"],
            "control_modules" => [Dict{String,Any}(
                                      "name" => "profile_limited",
                                      "profile_path" => "./profiles/tests/operation_triangular.prf",
                                  )],
            "power_th" => 20000,
        ),
        "TST_DEM_01" => Dict{String,Any}(
            "type" => "BoundedSink",
            "output_refs" => [],
            "medium" => "m_h_w_ht1",
            "constant_power" => 20000,
            "constant_temperature" => 65,
        ),
    )
    return components_config
end

function test_boiler_direct_demand_follows_profile()
    components_config = get_energy_system_fuel_boiler_direct()
    simulation_params = get_default_sim_params()
    components = Resie.load_components(components_config, simulation_params)
    setup_mock_run!(components, simulation_params)
    grid = components["TST_GRI_01"]
    boiler = components["TST_FB_01"]
    sink = components["TST_DEM_01"]

    # first timestep, profile limits to 100%
    EnergySystems.reset(grid)
    EnergySystems.reset(boiler)
    EnergySystems.reset(sink)

    EnergySystems.control(boiler, components, simulation_params)
    EnergySystems.control(grid, components, simulation_params)
    EnergySystems.control(sink, components, simulation_params)

    EnergySystems.process(boiler, simulation_params)
    EnergySystems.process(sink, simulation_params)
    EnergySystems.process(grid, simulation_params)

    @test boiler.output_interfaces[boiler.m_heat_out].sum_abs_change * 0.5 ≈ 5000.0
    @test boiler.output_interfaces[boiler.m_heat_out].balance ≈ 0.0
    @test boiler.input_interfaces[boiler.m_fuel_in].sum_abs_change * 0.5 ≈ 5500.0

    # second timestep, profile limits to 50%
    simulation_params["time"] = 900
    simulation_params["current_date"] = DateTime(2024, 1, 1, 0, 15)
    EnergySystems.reset(grid)
    EnergySystems.reset(boiler)
    EnergySystems.reset(sink)

    EnergySystems.control(boiler, components, simulation_params)
    EnergySystems.control(grid, components, simulation_params)
    EnergySystems.control(sink, components, simulation_params)

    EnergySystems.process(boiler, simulation_params)
    EnergySystems.process(sink, simulation_params)
    EnergySystems.process(grid, simulation_params)

    @test boiler.output_interfaces[boiler.m_heat_out].sum_abs_change * 0.5 ≈ 2500.0
    @test boiler.output_interfaces[boiler.m_heat_out].balance ≈ 0.0
    @test boiler.input_interfaces[boiler.m_fuel_in].sum_abs_change * 0.5 ≈ 2750.0

    # third timestep, profile forbids operation
    simulation_params["time"] = 1800
    simulation_params["current_date"] = DateTime(2024, 1, 1, 0, 30)
    EnergySystems.reset(grid)
    EnergySystems.reset(boiler)
    EnergySystems.reset(sink)

    EnergySystems.control(boiler, components, simulation_params)
    EnergySystems.control(grid, components, simulation_params)
    EnergySystems.control(sink, components, simulation_params)

    EnergySystems.process(boiler, simulation_params)
    EnergySystems.process(sink, simulation_params)
    EnergySystems.process(grid, simulation_params)

    @test boiler.output_interfaces[boiler.m_heat_out].sum_abs_change * 0.5 ≈ 0.0
    @test boiler.output_interfaces[boiler.m_heat_out].balance ≈ 0.0
    @test boiler.input_interfaces[boiler.m_fuel_in].sum_abs_change * 0.5 ≈ 0.0
end

@testset "boiler_direct_demand_follows_profile" begin
    test_boiler_direct_demand_follows_profile()
end

function get_energy_system_heat_pump_cascade()
    components_config = Dict{String,Any}(
        "TST_SRC_01" => Dict{String,Any}(
            "type" => "BoundedSupply",
            "output_refs" => ["TST_HP_01"],
            "is_source" => true,
            "medium" => "m_h_w_lt1",
            "constant_temperature" => 10,
            "constant_power" => 40000,
        ),
        "TST_GRI_01" => Dict{String,Any}(
            "type" => "GridConnection",
            "output_refs" => ["TST_HP_01"],
            "is_source" => true,
            "medium" => "m_e_ac_230v",
        ),
        "TST_GRI_02" => Dict{String,Any}(
            "type" => "GridConnection",
            "output_refs" => ["TST_HP_02"],
            "is_source" => true,
            "medium" => "m_e_ac_230v",
        ),
        "TST_HP_01" => Dict{String,Any}(
            "type" => "HeatPump",
            "model_type" => "simplified",
            "m_heat_in" => "m_h_w_lt1",
            "m_heat_out" => "m_h_w_lt2",
            "output_refs" => ["TST_HP_02"],
            "power_th" => 4000,
            "cop_function" => "const:2.0",
            "output_temperature" => 50,
            "power_losses_factor" => 1.0,
            "heat_losses_factor" => 1.0,
            "fudge_factor" => 1.0,
        ),
        "TST_HP_02" => Dict{String,Any}(
            "type" => "HeatPump",
            "model_type" => "simplified",
            "m_heat_in" => "m_h_w_lt2",
            "output_refs" => ["TST_GRO_01"],
            "control_modules" => [Dict{String,Any}(
                                      "name" => "profile_limited",
                                      "profile_path" => "./profiles/tests/operation_triangular.prf",
                                  )],
            "power_th" => 8000,
            "cop_function" => "const:2.0",
            "power_losses_factor" => 1.0,
            "heat_losses_factor" => 1.0,
            "fudge_factor" => 1.0,
        ),
        "TST_GRO_01" => Dict{String,Any}(
            "type" => "GridConnection",
            "output_refs" => [],
            "medium" => "m_h_w_ht1",
            "is_source" => false,
            "output_temperature" => 90,
        ),
    )
    return components_config
end

function test_heat_pump_cascade_follows_profile()
    components_config = get_energy_system_heat_pump_cascade()
    simulation_params = get_default_sim_params()
    components = Resie.load_components(components_config, simulation_params)
    setup_mock_run!(components, simulation_params)
    source = components["TST_SRC_01"]
    grid_1 = components["TST_GRI_01"]
    grid_2 = components["TST_GRI_02"]
    hp_1 = components["TST_HP_01"]
    hp_2 = components["TST_HP_02"]
    sink = components["TST_GRO_01"]

    # first timestep, profile limits hp2 to 100%
    EnergySystems.reset(source)
    EnergySystems.reset(grid_1)
    EnergySystems.reset(grid_2)
    EnergySystems.reset(hp_1)
    EnergySystems.reset(hp_2)
    EnergySystems.reset(sink)

    EnergySystems.control(hp_1, components, simulation_params)
    EnergySystems.control(hp_2, components, simulation_params)
    EnergySystems.control(source, components, simulation_params)
    EnergySystems.control(grid_1, components, simulation_params)
    EnergySystems.control(grid_2, components, simulation_params)
    EnergySystems.control(sink, components, simulation_params)

    EnergySystems.potential(hp_1, simulation_params)
    EnergySystems.potential(hp_2, simulation_params)
    EnergySystems.process(hp_2, simulation_params)
    EnergySystems.process(hp_1, simulation_params)
    EnergySystems.process(sink, simulation_params)
    EnergySystems.process(source, simulation_params)
    EnergySystems.process(grid_1, simulation_params)
    EnergySystems.process(grid_2, simulation_params)

    @test hp_1.output_interfaces[hp_1.m_heat_out].sum_abs_change * 0.5 ≈ 1000.0
    @test hp_1.output_interfaces[hp_1.m_heat_out].balance ≈ 0.0
    @test hp_1.input_interfaces[hp_1.m_heat_in].sum_abs_change * 0.5 ≈ 500.0
    @test hp_1.input_interfaces[hp_1.m_el_in].sum_abs_change * 0.5 ≈ 500.0

    @test hp_2.output_interfaces[hp_2.m_heat_out].sum_abs_change * 0.5 ≈ 2000.0
    @test hp_2.output_interfaces[hp_2.m_heat_out].balance ≈ 0.0
    @test hp_2.input_interfaces[hp_2.m_heat_in].sum_abs_change * 0.5 ≈ 1000.0
    @test hp_2.input_interfaces[hp_2.m_el_in].sum_abs_change * 0.5 ≈ 1000.0

    # second timestep, profile limits to 50%
    simulation_params["time"] = 900
    simulation_params["current_date"] = DateTime(2024, 1, 1, 0, 15)
    EnergySystems.reset(source)
    EnergySystems.reset(grid_1)
    EnergySystems.reset(grid_2)
    EnergySystems.reset(hp_1)
    EnergySystems.reset(hp_2)
    EnergySystems.reset(sink)

    EnergySystems.control(hp_1, components, simulation_params)
    EnergySystems.control(hp_2, components, simulation_params)
    EnergySystems.control(source, components, simulation_params)
    EnergySystems.control(grid_1, components, simulation_params)
    EnergySystems.control(grid_2, components, simulation_params)
    EnergySystems.control(sink, components, simulation_params)

    EnergySystems.potential(hp_1, simulation_params)
    EnergySystems.potential(hp_2, simulation_params)
    EnergySystems.process(hp_2, simulation_params)
    EnergySystems.process(hp_1, simulation_params)
    EnergySystems.process(sink, simulation_params)
    EnergySystems.process(source, simulation_params)
    EnergySystems.process(grid_1, simulation_params)
    EnergySystems.process(grid_2, simulation_params)

    @test hp_1.output_interfaces[hp_1.m_heat_out].sum_abs_change * 0.5 ≈ 500.0
    @test hp_1.output_interfaces[hp_1.m_heat_out].balance ≈ 0.0
    @test hp_1.input_interfaces[hp_1.m_heat_in].sum_abs_change * 0.5 ≈ 250.0
    @test hp_1.input_interfaces[hp_1.m_el_in].sum_abs_change * 0.5 ≈ 250.0

    @test hp_2.output_interfaces[hp_2.m_heat_out].sum_abs_change * 0.5 ≈ 1000.0
    @test hp_2.output_interfaces[hp_2.m_heat_out].balance ≈ 0.0
    @test hp_2.input_interfaces[hp_2.m_heat_in].sum_abs_change * 0.5 ≈ 500.0
    @test hp_2.input_interfaces[hp_2.m_el_in].sum_abs_change * 0.5 ≈ 500.0

    # third timestep, profile forbids operation
    simulation_params["time"] = 1800
    simulation_params["current_date"] = DateTime(2024, 1, 1, 0, 30)
    EnergySystems.reset(source)
    EnergySystems.reset(grid_1)
    EnergySystems.reset(grid_2)
    EnergySystems.reset(hp_1)
    EnergySystems.reset(hp_2)
    EnergySystems.reset(sink)

    EnergySystems.control(hp_1, components, simulation_params)
    EnergySystems.control(hp_2, components, simulation_params)
    EnergySystems.control(source, components, simulation_params)
    EnergySystems.control(grid_1, components, simulation_params)
    EnergySystems.control(grid_2, components, simulation_params)
    EnergySystems.control(sink, components, simulation_params)

    EnergySystems.potential(hp_1, simulation_params)
    EnergySystems.potential(hp_2, simulation_params)
    EnergySystems.process(hp_2, simulation_params)
    EnergySystems.process(hp_1, simulation_params)
    EnergySystems.process(sink, simulation_params)
    EnergySystems.process(source, simulation_params)
    EnergySystems.process(grid_1, simulation_params)
    EnergySystems.process(grid_2, simulation_params)

    @test hp_1.output_interfaces[hp_1.m_heat_out].sum_abs_change * 0.5 ≈ 0.0
    @test hp_1.output_interfaces[hp_1.m_heat_out].balance ≈ 0.0
    @test hp_1.input_interfaces[hp_1.m_heat_in].sum_abs_change * 0.5 ≈ 0.0
    @test hp_1.input_interfaces[hp_1.m_el_in].sum_abs_change * 0.5 ≈ 0.0

    @test hp_2.output_interfaces[hp_2.m_heat_out].sum_abs_change * 0.5 ≈ 0.0
    @test hp_2.output_interfaces[hp_2.m_heat_out].balance ≈ 0.0
    @test hp_2.input_interfaces[hp_2.m_heat_in].sum_abs_change * 0.5 ≈ 0.0
    @test hp_2.input_interfaces[hp_2.m_el_in].sum_abs_change * 0.5 ≈ 0.0
end

@testset "heat_pump_cascade_follows_profile" begin
    test_heat_pump_cascade_follows_profile()
end
