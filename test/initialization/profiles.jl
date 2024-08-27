using Debugger
using Test
using Resie
using Resie.EnergySystems
using Resie.Profiles

include("../test_util.jl")

function energy_system()::Dict{String,Any}
    return Dict{String,Any}(
        "TST_GRI_01" => Dict{String,Any}(
            "type" => "GridConnection",
            "medium" => "m_h_w_ht1",
            "output_refs" => ["TST_DEM_01"],
            "is_source" => true,
        ),
        "TST_DEM_01" => Dict{String,Any}(
            "type" => "Demand",
            "medium" => "m_h_w_ht1",
            "output_refs" => [],
            "energy_profile_file_path" => "./profiles/tests/heating_demand_short.prf",
            "temperature_profile_file_path" => "./profiles/tests/temperature_short.prf",
            "scale" => 1,
        ),
    )
end

function test_profile_aggregation_two()
    components_config = energy_system()

    simulation_parameters = get_default_sim_params()
    simulation_parameters["time_step_seconds"] = 1800

    components = Resie.load_components(components_config, simulation_parameters)
    demand = components["TST_DEM_01"]

    @test demand.energy_profile.time_step == 1800
    @test demand.energy_profile.data_type == "extensive"
    @test demand.energy_profile.data == Dict(DateTime("2024-01-01T00:00:00") => 22.0,
                                             DateTime("2024-01-01T00:30:00") => 60.0,
                                             DateTime("2024-01-01T01:00:00") => 25.5,
                                             DateTime("2024-01-01T01:30:00") => 16.0,
                                             DateTime("2024-01-01T02:00:00") => 43.0)
    @test sum(values(demand.energy_profile.data)) == 166.5

    @test demand.temperature_profile.time_step == 1800
    @test demand.temperature_profile.data_type == "intensive"
    @test demand.temperature_profile.data == Dict(DateTime("2024-01-01T00:00:00") => 64.0,
                                                  DateTime("2024-01-01T00:30:00") => 69.0,
                                                  DateTime("2024-01-01T01:00:00") => 72.0,
                                                  DateTime("2024-01-01T01:30:00") => 58.5,
                                                  DateTime("2024-01-01T02:00:00") => 68.0)
end

function test_profile_aggregation_four()
    components_config = energy_system()

    simulation_parameters = get_default_sim_params()
    simulation_parameters["time_step_seconds"] = 3600

    components = Resie.load_components(components_config, simulation_parameters)
    demand = components["TST_DEM_01"]

    @test demand.energy_profile.time_step == 3600
    @test demand.energy_profile.data_type == "extensive"
    @test demand.energy_profile.data == Dict(DateTime("2024-01-01T00:00:00") => 82.0,
                                             DateTime("2024-01-01T01:00:00") => 41.5,
                                             DateTime("2024-01-01T02:00:00") => 43.0)
    @test sum(values(demand.energy_profile.data)) == 166.5

    @test demand.temperature_profile.time_step == 3600
    @test demand.temperature_profile.data_type == "intensive"
    @test demand.temperature_profile.data == Dict(DateTime("2024-01-01T00:00:00") => 66.5,
                                                  DateTime("2024-01-01T01:00:00") => 65.25,
                                                  DateTime("2024-01-01T02:00:00") => 68.0)
end

function test_profile_segmentation_half()
    components_config = energy_system()

    simulation_parameters = get_default_sim_params()
    simulation_parameters["time_step_seconds"] = 450

    components = Resie.load_components(components_config, simulation_parameters)
    demand = components["TST_DEM_01"]

    @test demand.energy_profile.time_step == 450
    @test demand.energy_profile.data_type == "extensive"
    expected = Dict(DateTime("2024-01-01T00:00:00") => 5,
                    DateTime("2024-01-01T00:07:30") => 5,
                    DateTime("2024-01-01T00:15:00") => 6,
                    DateTime("2024-01-01T00:22:30") => 6,
                    DateTime("2024-01-01T00:30:00") => 30,
                    DateTime("2024-01-01T00:37:30") => 30,
                    DateTime("2024-01-01T00:45:00") => 0,
                    DateTime("2024-01-01T00:52:30") => 0,
                    DateTime("2024-01-01T01:00:00") => 1.25,
                    DateTime("2024-01-01T01:07:30") => 1.25,
                    DateTime("2024-01-01T01:15:00") => 11.5,
                    DateTime("2024-01-01T01:22:30") => 11.5,
                    DateTime("2024-01-01T01:30:00") => 2.5,
                    DateTime("2024-01-01T01:37:30") => 2.5,
                    DateTime("2024-01-01T01:45:00") => 5.5,
                    DateTime("2024-01-01T01:52:30") => 5.5,
                    DateTime("2024-01-01T02:00:00") => 0.5,
                    DateTime("2024-01-01T02:07:30") => 0.5,
                    DateTime("2024-01-01T02:15:00") => 21,
                    DateTime("2024-01-01T02:22:30") => 21)
    @test all([demand.energy_profile.data[k] - expected[k] for k in keys(expected)] .<
              simulation_parameters["epsilon"]) == true
    @test all([demand.energy_profile.data[k] - expected[k] for k in keys(demand.energy_profile.data)] .<
              simulation_parameters["epsilon"]) == true

    @test sum(values(demand.energy_profile.data)) == 166.5

    @test demand.temperature_profile.time_step == 450
    @test demand.temperature_profile.data_type == "intensive"
    expected = Dict(DateTime("2024-01-01T00:00:00") => 65,
                    DateTime("2024-01-01T00:07:30") => 64,
                    DateTime("2024-01-01T00:15:00") => 63,
                    DateTime("2024-01-01T00:22:30") => 66.5,
                    DateTime("2024-01-01T00:30:00") => 70,
                    DateTime("2024-01-01T00:37:30") => 69,
                    DateTime("2024-01-01T00:45:00") => 68,
                    DateTime("2024-01-01T00:52:30") => 70,
                    DateTime("2024-01-01T01:00:00") => 72,
                    DateTime("2024-01-01T01:07:30") => 72,
                    DateTime("2024-01-01T01:15:00") => 72,
                    DateTime("2024-01-01T01:22:30") => 64,
                    DateTime("2024-01-01T01:30:00") => 56,
                    DateTime("2024-01-01T01:37:30") => 58.5,
                    DateTime("2024-01-01T01:45:00") => 61,
                    DateTime("2024-01-01T01:52:30") => 63.5,
                    DateTime("2024-01-01T02:00:00") => 66,
                    DateTime("2024-01-01T02:07:30") => 68,
                    DateTime("2024-01-01T02:15:00") => 70)
    @test all([demand.temperature_profile.data[k] - expected[k] for k in keys(expected)] .<
              simulation_parameters["epsilon"]) == true
    @test all([demand.temperature_profile.data[k] - expected[k] for k in keys(demand.temperature_profile.data)] .<
              simulation_parameters["epsilon"]) == true
end

function test_profile_segmentation_third()
    components_config = energy_system()

    simulation_parameters = get_default_sim_params()
    simulation_parameters["time_step_seconds"] = 300

    components = Resie.load_components(components_config, simulation_parameters)
    demand = components["TST_DEM_01"]

    @test demand.energy_profile.time_step == 300
    @test demand.energy_profile.data_type == "extensive"
    expected = Dict(DateTime("2024-01-01T00:00:00") => 10 / 3,
                    DateTime("2024-01-01T00:05:00") => 10 / 3,
                    DateTime("2024-01-01T00:10:00") => 10 / 3,
                    DateTime("2024-01-01T00:15:00") => 12 / 3,
                    DateTime("2024-01-01T00:20:00") => 12 / 3,
                    DateTime("2024-01-01T00:25:00") => 12 / 3,
                    DateTime("2024-01-01T00:30:00") => 60 / 3,
                    DateTime("2024-01-01T00:35:00") => 60 / 3,
                    DateTime("2024-01-01T00:40:00") => 60 / 3,
                    DateTime("2024-01-01T00:45:00") => 0 / 3,
                    DateTime("2024-01-01T00:50:00") => 0 / 3,
                    DateTime("2024-01-01T00:55:00") => 0 / 3,
                    DateTime("2024-01-01T01:00:00") => 2.5 / 3,
                    DateTime("2024-01-01T01:05:00") => 2.5 / 3,
                    DateTime("2024-01-01T01:10:00") => 2.5 / 3,
                    DateTime("2024-01-01T01:15:00") => 23 / 3,
                    DateTime("2024-01-01T01:20:00") => 23 / 3,
                    DateTime("2024-01-01T01:25:00") => 23 / 3,
                    DateTime("2024-01-01T01:30:00") => 5 / 3,
                    DateTime("2024-01-01T01:35:00") => 5 / 3,
                    DateTime("2024-01-01T01:40:00") => 5 / 3,
                    DateTime("2024-01-01T01:45:00") => 11 / 3,
                    DateTime("2024-01-01T01:50:00") => 11 / 3,
                    DateTime("2024-01-01T01:55:00") => 11 / 3,
                    DateTime("2024-01-01T02:00:00") => 1 / 3,
                    DateTime("2024-01-01T02:05:00") => 1 / 3,
                    DateTime("2024-01-01T02:10:00") => 1 / 3,
                    DateTime("2024-01-01T02:15:00") => 42 / 3,
                    DateTime("2024-01-01T02:20:00") => 42 / 3,
                    DateTime("2024-01-01T02:25:00") => 42 / 3)
    @test all([demand.energy_profile.data[k] - expected[k] for k in keys(expected)] .<
              simulation_parameters["epsilon"]) == true
    @test all([demand.energy_profile.data[k] - expected[k] for k in keys(demand.energy_profile.data)] .<
              simulation_parameters["epsilon"]) == true
    @test sum(values(demand.energy_profile.data)) ≈ 166.5

    @test demand.temperature_profile.time_step == 300
    @test demand.temperature_profile.data_type == "intensive"
    expected = Dict(DateTime("2024-01-01T00:00:00") => 65,
                    DateTime("2024-01-01T00:05:00") => 64.3333333333333,
                    DateTime("2024-01-01T00:10:00") => 63.6666666666666,
                    DateTime("2024-01-01T00:15:00") => 63,
                    DateTime("2024-01-01T00:20:00") => 65.3333333333333,
                    DateTime("2024-01-01T00:25:00") => 67.6666666666666,
                    DateTime("2024-01-01T00:30:00") => 70,
                    DateTime("2024-01-01T00:35:00") => 69.3333333333333,
                    DateTime("2024-01-01T00:40:00") => 68.6666666666666,
                    DateTime("2024-01-01T00:45:00") => 68,
                    DateTime("2024-01-01T00:50:00") => 69.3333333333333,
                    DateTime("2024-01-01T00:55:00") => 70.6666666666666,
                    DateTime("2024-01-01T01:00:00") => 72,
                    DateTime("2024-01-01T01:05:00") => 72,
                    DateTime("2024-01-01T01:10:00") => 72,
                    DateTime("2024-01-01T01:15:00") => 72,
                    DateTime("2024-01-01T01:20:00") => 66.6666666666666,
                    DateTime("2024-01-01T01:25:00") => 61.3333333333333,
                    DateTime("2024-01-01T01:30:00") => 56,
                    DateTime("2024-01-01T01:35:00") => 57.6666666666666,
                    DateTime("2024-01-01T01:40:00") => 59.3333333333333,
                    DateTime("2024-01-01T01:45:00") => 61,
                    DateTime("2024-01-01T01:50:00") => 62.6666666666666,
                    DateTime("2024-01-01T01:55:00") => 64.3333333333333,
                    DateTime("2024-01-01T02:00:00") => 66,
                    DateTime("2024-01-01T02:05:00") => 67.3333333333333,
                    DateTime("2024-01-01T02:10:00") => 68.6666666666666,
                    DateTime("2024-01-01T02:15:00") => 70)
    @test all([demand.temperature_profile.data[k] - expected[k] for k in keys(expected)] .<
              simulation_parameters["epsilon"]) == true
    @test all([demand.temperature_profile.data[k] - expected[k] for k in keys(demand.temperature_profile.data)] .<
              simulation_parameters["epsilon"]) == true
end

@testset "profiles" begin
    test_profile_aggregation_two()
    test_profile_aggregation_four()
    test_profile_segmentation_half()
    test_profile_segmentation_third()
end
