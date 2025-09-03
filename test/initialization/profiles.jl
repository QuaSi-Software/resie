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

    simulation_parameters = Dict{String,Any}(
        "time" => 0,
        "current_date" => DateTime("01.01.2024 00:00", "dd.mm.yyy HH:MM"),
        "time_step_seconds" => 1800,
        "number_of_time_steps" => 6,
        "start_date" => DateTime("01.01.2024 00:00", "dd.mm.yyy HH:MM"),
        "end_date" => DateTime("01.01.2024 02:00", "dd.mm.yyy HH:MM"),
        "force_profiles_to_repeat" => false,
        "epsilon" => 1e-9,
        "latitude" => nothing,
        "longitude" => nothing,
        "watt_to_wh" => function (w)
            return w * 0.5
        end,
        "wh_to_watts" => function (w)
            return w * 2.0
        end,
    )

    components = Resie.load_components(components_config, simulation_parameters)
    demand = components["TST_DEM_01"]

    @test demand.energy_profile.time_step == 1800
    @test demand.energy_profile.data_type == :extensive
    @test demand.energy_profile.data == Dict(DateTime("2024-01-01T00:00:00") => 22.0,
                                             DateTime("2024-01-01T00:30:00") => 60.0,
                                             DateTime("2024-01-01T01:00:00") => 25.5,
                                             DateTime("2024-01-01T01:30:00") => 16.0,
                                             DateTime("2024-01-01T02:00:00") => 43.0)
    @test sum(values(demand.energy_profile.data)) == 166.5

    @test demand.temperature_profile.time_step == 1800
    @test demand.temperature_profile.data_type == :intensive
    @test demand.temperature_profile.data == Dict(DateTime("2024-01-01T00:00:00") => 64.0,
                                                  DateTime("2024-01-01T00:30:00") => 69.0,
                                                  DateTime("2024-01-01T01:00:00") => 72.0,
                                                  DateTime("2024-01-01T01:30:00") => 58.5,
                                                  DateTime("2024-01-01T02:00:00") => 68.0)
end

function test_profile_aggregation_four()
    components_config = energy_system()

    simulation_parameters = Dict{String,Any}(
        "time" => 0,
        "current_date" => DateTime("01.01.2024 00:00", "dd.mm.yyy HH:MM"),
        "time_step_seconds" => 3600,
        "number_of_time_steps" => 3,
        "start_date" => DateTime("01.01.2024 00:00", "dd.mm.yyy HH:MM"),
        "end_date" => DateTime("01.01.2024 02:00", "dd.mm.yyy HH:MM"),
        "force_profiles_to_repeat" => false,
        "epsilon" => 1e-9,
        "latitude" => nothing,
        "longitude" => nothing,
        "watt_to_wh" => function (w)
            return w * 1.0
        end,
        "wh_to_watts" => function (w)
            return w * 1.0
        end,
    )

    components = Resie.load_components(components_config, simulation_parameters)
    demand = components["TST_DEM_01"]

    @test demand.energy_profile.time_step == 3600
    @test demand.energy_profile.data_type == :extensive
    @test demand.energy_profile.data == Dict(DateTime("2024-01-01T00:00:00") => 82.0,
                                             DateTime("2024-01-01T01:00:00") => 41.5,
                                             DateTime("2024-01-01T02:00:00") => 127.0)
    @test sum(values(demand.energy_profile.data)) == 166.5 + 2 * 42

    @test demand.temperature_profile.time_step == 3600
    @test demand.temperature_profile.data_type == :intensive
    @test demand.temperature_profile.data == Dict(DateTime("2024-01-01T00:00:00") => 66.5,
                                                  DateTime("2024-01-01T01:00:00") => 65.25,
                                                  DateTime("2024-01-01T02:00:00") => 69.0)
end

function test_profile_segmentation_half()
    components_config = energy_system()

    simulation_parameters = Dict{String,Any}(
        "time" => 0,
        "current_date" => DateTime("01.01.2024 00:00", "dd.mm.yyy HH:MM"),
        "time_step_seconds" => 450,
        "number_of_time_steps" => 18,
        "start_date" => DateTime("01.01.2024 00:00", "dd.mm.yyy HH:MM"),
        "end_date" => DateTime("01.01.2024 02:15", "dd.mm.yyy HH:MM"),
        "force_profiles_to_repeat" => false,
        "epsilon" => 1e-9,
        "latitude" => nothing,
        "longitude" => nothing,
        "watt_to_wh" => function (w)
            return w * 0.125
        end,
        "wh_to_watts" => function (w)
            return w * 8.0
        end,
    )

    components = Resie.load_components(components_config, simulation_parameters)
    demand = components["TST_DEM_01"]

    @test demand.energy_profile.time_step == 450
    @test demand.energy_profile.data_type == :extensive
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
                    DateTime("2024-01-01T02:15:00") => 21)
    @test all([demand.energy_profile.data[k] - expected[k] for k in keys(expected)] .<
              simulation_parameters["epsilon"]) == true
    @test all([demand.energy_profile.data[k] - expected[k] for k in keys(demand.energy_profile.data)] .<
              simulation_parameters["epsilon"]) == true

    @test sum(values(demand.energy_profile.data)) == 166.5 - 21

    @test demand.temperature_profile.time_step == 450
    @test demand.temperature_profile.data_type == :intensive
    expected = Dict(DateTime("2024-01-01T00:00:00") => 65.0,
                    DateTime("2024-01-01T00:07:30") => 65.0,
                    DateTime("2024-01-01T00:15:00") => 63.0,
                    DateTime("2024-01-01T00:22:30") => 63.0,
                    DateTime("2024-01-01T00:30:00") => 70.0,
                    DateTime("2024-01-01T00:37:30") => 70.0,
                    DateTime("2024-01-01T00:45:00") => 68.0,
                    DateTime("2024-01-01T00:52:30") => 68.0,
                    DateTime("2024-01-01T01:00:00") => 72.0,
                    DateTime("2024-01-01T01:07:30") => 72.0,
                    DateTime("2024-01-01T01:15:00") => 72.0,
                    DateTime("2024-01-01T01:22:30") => 72.0,
                    DateTime("2024-01-01T01:30:00") => 56.0,
                    DateTime("2024-01-01T01:37:30") => 56.0,
                    DateTime("2024-01-01T01:45:00") => 61.0,
                    DateTime("2024-01-01T01:52:30") => 61.0,
                    DateTime("2024-01-01T02:00:00") => 66.0,
                    DateTime("2024-01-01T02:07:30") => 66.0,
                    DateTime("2024-01-01T02:15:00") => 70.0)
    @test all([demand.temperature_profile.data[k] - expected[k] for k in keys(expected)] .<
              simulation_parameters["epsilon"]) == true
    @test all([demand.temperature_profile.data[k] - expected[k] for k in keys(demand.temperature_profile.data)] .<
              simulation_parameters["epsilon"]) == true
end

function test_profile_segmentation_third()
    components_config = energy_system()

    simulation_parameters = Dict{String,Any}(
        "time" => 0,
        "current_date" => DateTime("01.01.2024 00:00", "dd.mm.yyy HH:MM"),
        "time_step_seconds" => 300,
        "number_of_time_steps" => 27,
        "start_date" => DateTime("01.01.2024 00:00", "dd.mm.yyy HH:MM"),
        "end_date" => DateTime("01.01.2024 02:15", "dd.mm.yyy HH:MM"),
        "force_profiles_to_repeat" => false,
        "epsilon" => 1e-9,
        "latitude" => nothing,
        "longitude" => nothing,
        "watt_to_wh" => function (w)
            return w * 0.08333333333333
        end,
        "wh_to_watts" => function (w)
            return w * 12.0
        end,
    )

    components = Resie.load_components(components_config, simulation_parameters)
    demand = components["TST_DEM_01"]

    @test demand.energy_profile.time_step == 300
    @test demand.energy_profile.data_type == :extensive
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
                    DateTime("2024-01-01T02:15:00") => 42 / 3)
    @test all([demand.energy_profile.data[k] - expected[k] for k in keys(expected)] .<
              simulation_parameters["epsilon"]) == true
    @test all([demand.energy_profile.data[k] - expected[k] for k in keys(demand.energy_profile.data)] .<
              simulation_parameters["epsilon"]) == true
    @test sum(values(demand.energy_profile.data)) ≈ 166.5 - 42 / 3 * 2

    @test demand.temperature_profile.time_step == 300
    @test demand.temperature_profile.data_type == :intensive
    expected = Dict(DateTime("2024-01-01T00:00:00") => 65.0,
                    DateTime("2024-01-01T00:05:00") => 65.0,
                    DateTime("2024-01-01T00:10:00") => 65.0,
                    DateTime("2024-01-01T00:15:00") => 63.0,
                    DateTime("2024-01-01T00:20:00") => 63.0,
                    DateTime("2024-01-01T00:25:00") => 63.0,
                    DateTime("2024-01-01T00:30:00") => 70.0,
                    DateTime("2024-01-01T00:35:00") => 70.0,
                    DateTime("2024-01-01T00:40:00") => 70.0,
                    DateTime("2024-01-01T00:45:00") => 68.0,
                    DateTime("2024-01-01T00:50:00") => 68.0,
                    DateTime("2024-01-01T00:55:00") => 68.0,
                    DateTime("2024-01-01T01:00:00") => 72.0,
                    DateTime("2024-01-01T01:05:00") => 72.0,
                    DateTime("2024-01-01T01:10:00") => 72.0,
                    DateTime("2024-01-01T01:15:00") => 72.0,
                    DateTime("2024-01-01T01:20:00") => 72.0,
                    DateTime("2024-01-01T01:25:00") => 72.0,
                    DateTime("2024-01-01T01:30:00") => 56.0,
                    DateTime("2024-01-01T01:35:00") => 56.0,
                    DateTime("2024-01-01T01:40:00") => 56.0,
                    DateTime("2024-01-01T01:45:00") => 61.0,
                    DateTime("2024-01-01T01:50:00") => 61.0,
                    DateTime("2024-01-01T01:55:00") => 61.0,
                    DateTime("2024-01-01T02:00:00") => 66.0,
                    DateTime("2024-01-01T02:05:00") => 66.0,
                    DateTime("2024-01-01T02:10:00") => 66.0,
                    DateTime("2024-01-01T02:15:00") => 70.0)
    @test all([demand.temperature_profile.data[k] - expected[k] for k in keys(expected)] .<
              simulation_parameters["epsilon"]) == true
    @test all([demand.temperature_profile.data[k] - expected[k] for k in keys(demand.temperature_profile.data)] .<
              simulation_parameters["epsilon"]) == true
end

@testset "profiles_stepwise" begin
    test_profile_aggregation_two()
    test_profile_aggregation_four()
    test_profile_segmentation_half()
    test_profile_segmentation_third()
end

function energy_system_linear_classic()::Dict{String,Any}
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
            "energy_profile_file_path" => "./profiles/tests/heating_demand_short_linear.prf",
            "temperature_profile_file_path" => "./profiles/tests/temperature_short_linear.prf",
            "scale" => 1,
        ),
    )
end

function test_profile_segmentation_half_linear_classic()
    components_config = energy_system_linear_classic()

    simulation_parameters = Dict{String,Any}(
        "time" => 0,
        "current_date" => DateTime("01.01.2024 00:00", "dd.mm.yyy HH:MM"),
        "time_step_seconds" => 450,
        "number_of_time_steps" => 18,
        "start_date" => DateTime("01.01.2024 00:00", "dd.mm.yyy HH:MM"),
        "end_date" => DateTime("01.01.2024 02:15", "dd.mm.yyy HH:MM"),
        "force_profiles_to_repeat" => false,
        "epsilon" => 1e-9,
        "latitude" => nothing,
        "longitude" => nothing,
        "watt_to_wh" => function (w)
            return w * 0.125
        end,
        "wh_to_watts" => function (w)
            return w * 8.0
        end,
    )

    components = Resie.load_components(components_config, simulation_parameters)
    demand = components["TST_DEM_01"]

    @test demand.energy_profile.time_step == 450
    @test demand.energy_profile.data_type == :extensive
    expected = Dict(DateTime("2024-01-01T00:00:00") => 5.0,
                    DateTime("2024-01-01T00:07:30") => 5.5,
                    DateTime("2024-01-01T00:15:00") => 6.0,
                    DateTime("2024-01-01T00:22:30") => 18.0,
                    DateTime("2024-01-01T00:30:00") => 30.0,
                    DateTime("2024-01-01T00:37:30") => 15.0,
                    DateTime("2024-01-01T00:45:00") => 0.0,
                    DateTime("2024-01-01T00:52:30") => 0.625,
                    DateTime("2024-01-01T01:00:00") => 1.25,
                    DateTime("2024-01-01T01:07:30") => 6.375,
                    DateTime("2024-01-01T01:15:00") => 11.5,
                    DateTime("2024-01-01T01:22:30") => 7.0,
                    DateTime("2024-01-01T01:30:00") => 2.5,
                    DateTime("2024-01-01T01:37:30") => 4.0,
                    DateTime("2024-01-01T01:45:00") => 5.5,
                    DateTime("2024-01-01T01:52:30") => 3.0,
                    DateTime("2024-01-01T02:00:00") => 0.5,
                    DateTime("2024-01-01T02:07:30") => 10.75,
                    DateTime("2024-01-01T02:15:00") => 21.0)
    @test all([demand.energy_profile.data[k] - expected[k] for k in keys(expected)] .<
              simulation_parameters["epsilon"]) == true
    @test all([demand.energy_profile.data[k] - expected[k] for k in keys(demand.energy_profile.data)] .<
              simulation_parameters["epsilon"]) == true

    @test sum(values(demand.energy_profile.data)) == 153.5

    @test demand.temperature_profile.time_step == 450
    @test demand.temperature_profile.data_type == :intensive
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

function test_profile_segmentation_third_linear_classic()
    components_config = energy_system_linear_classic()

    simulation_parameters = Dict{String,Any}(
        "time" => 0,
        "current_date" => DateTime("01.01.2024 00:00", "dd.mm.yyy HH:MM"),
        "time_step_seconds" => 300,
        "number_of_time_steps" => 27,
        "start_date" => DateTime("01.01.2024 00:00", "dd.mm.yyy HH:MM"),
        "end_date" => DateTime("01.01.2024 02:15", "dd.mm.yyy HH:MM"),
        "force_profiles_to_repeat" => false,
        "epsilon" => 1e-9,
        "latitude" => nothing,
        "longitude" => nothing,
        "watt_to_wh" => function (w)
            return w * 0.0833333333333
        end,
        "wh_to_watts" => function (w)
            return w * 12.0
        end,
    )

    components = Resie.load_components(components_config, simulation_parameters)
    demand = components["TST_DEM_01"]

    @test demand.energy_profile.time_step == 300
    @test demand.energy_profile.data_type == :extensive
    expected = Dict(DateTime("2024-01-01T00:00:00") => 3.333333333333333,
                    DateTime("2024-01-01T00:05:00") => 3.5555555555555554,
                    DateTime("2024-01-01T00:10:00") => 3.7777777777777777,
                    DateTime("2024-01-01T00:15:00") => 4.0,
                    DateTime("2024-01-01T00:20:00") => 9.333333333333332,
                    DateTime("2024-01-01T00:25:00") => 14.666666666666666,
                    DateTime("2024-01-01T00:30:00") => 20.0,
                    DateTime("2024-01-01T00:35:00") => 13.333333333333332,
                    DateTime("2024-01-01T00:40:00") => 6.666666666666666,
                    DateTime("2024-01-01T00:45:00") => 0.0,
                    DateTime("2024-01-01T00:50:00") => 0.27777777777777773,
                    DateTime("2024-01-01T00:55:00") => 0.5555555555555555,
                    DateTime("2024-01-01T01:00:00") => 0.8333333333333333,
                    DateTime("2024-01-01T01:05:00") => 3.1111111111111107,
                    DateTime("2024-01-01T01:10:00") => 5.388888888888888,
                    DateTime("2024-01-01T01:15:00") => 7.666666666666666,
                    DateTime("2024-01-01T01:20:00") => 5.666666666666666,
                    DateTime("2024-01-01T01:25:00") => 3.666666666666666,
                    DateTime("2024-01-01T01:30:00") => 1.6666666666666665,
                    DateTime("2024-01-01T01:35:00") => 2.333333333333333,
                    DateTime("2024-01-01T01:40:00") => 2.9999999999999996,
                    DateTime("2024-01-01T01:45:00") => 3.6666666666666665,
                    DateTime("2024-01-01T01:50:00") => 2.5555555555555554,
                    DateTime("2024-01-01T01:55:00") => 1.4444444444444442,
                    DateTime("2024-01-01T02:00:00") => 0.3333333333333333,
                    DateTime("2024-01-01T02:05:00") => 4.888888888888888,
                    DateTime("2024-01-01T02:10:00") => 9.444444444444443,
                    DateTime("2024-01-01T02:15:00") => 14.0)
    @test all([demand.energy_profile.data[k] - expected[k] for k in keys(expected)] .<
              simulation_parameters["epsilon"]) == true
    @test all([demand.energy_profile.data[k] - expected[k] for k in keys(demand.energy_profile.data)] .<
              simulation_parameters["epsilon"]) == true
    @test sum(values(demand.energy_profile.data)) ≈ 149.16666666666666

    @test demand.temperature_profile.time_step == 300
    @test demand.temperature_profile.data_type == :intensive
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

@testset "profiles_linear_classic" begin
    test_profile_segmentation_half_linear_classic()
    test_profile_segmentation_third_linear_classic()
end

function energy_system_linear_time_preserving()::Dict{String,Any}
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
            "energy_profile_file_path" => "./profiles/tests/heating_demand_short_linear_time_preserving.prf",
            "temperature_profile_file_path" => "./profiles/tests/temperature_short_linear_time_preserving.prf",
            "scale" => 1,
        ),
    )
end

function test_profile_segmentation_half_linear_time_preserving()
    components_config = energy_system_linear_time_preserving()

    simulation_parameters = Dict{String,Any}(
        "time" => 0,
        "current_date" => DateTime("01.01.2024 00:00", "dd.mm.yyy HH:MM"),
        "time_step_seconds" => 450,
        "number_of_time_steps" => 18,
        "start_date" => DateTime("01.01.2024 00:00", "dd.mm.yyy HH:MM"),
        "end_date" => DateTime("01.01.2024 02:15", "dd.mm.yyy HH:MM"),
        "force_profiles_to_repeat" => false,
        "epsilon" => 1e-9,
        "latitude" => nothing,
        "longitude" => nothing,
        "watt_to_wh" => function (w)
            return w * 0.125
        end,
        "wh_to_watts" => function (w)
            return w * 8.0
        end,
    )

    components = Resie.load_components(components_config, simulation_parameters)
    demand = components["TST_DEM_01"]

    @test demand.energy_profile.time_step == 450
    @test demand.energy_profile.data_type == :extensive
    expected = Dict(DateTime("2024-01-01T00:00:00") => 5.0,
                    DateTime("2024-01-01T00:07:30") => 5.25,
                    DateTime("2024-01-01T00:15:00") => 5.75,
                    DateTime("2024-01-01T00:22:30") => 12.0,
                    DateTime("2024-01-01T00:30:00") => 24.0,
                    DateTime("2024-01-01T00:37:30") => 22.5,
                    DateTime("2024-01-01T00:45:00") => 7.5,
                    DateTime("2024-01-01T00:52:30") => 0.3125,
                    DateTime("2024-01-01T01:00:00") => 0.9375,
                    DateTime("2024-01-01T01:07:30") => 3.8125,
                    DateTime("2024-01-01T01:15:00") => 8.9375,
                    DateTime("2024-01-01T01:22:30") => 9.25,
                    DateTime("2024-01-01T01:30:00") => 4.75,
                    DateTime("2024-01-01T01:37:30") => 3.25,
                    DateTime("2024-01-01T01:45:00") => 4.75,
                    DateTime("2024-01-01T01:52:30") => 4.25,
                    DateTime("2024-01-01T02:00:00") => 1.75,
                    DateTime("2024-01-01T02:07:30") => 5.625,
                    DateTime("2024-01-01T02:15:00") => 10.75)
    @test all([demand.energy_profile.data[k] - expected[k] for k in keys(expected)] .<
              simulation_parameters["epsilon"]) == true
    @test all([demand.energy_profile.data[k] - expected[k] for k in keys(demand.energy_profile.data)] .<
              simulation_parameters["epsilon"]) == true

    @test sum(values(demand.energy_profile.data)) == 140.375

    @test demand.temperature_profile.time_step == 450
    @test demand.temperature_profile.data_type == :intensive
    expected = Dict(DateTime("2024-01-01T00:00:00") => 65.0,
                    DateTime("2024-01-01T00:07:30") => 64.5,
                    DateTime("2024-01-01T00:15:00") => 63.5,
                    DateTime("2024-01-01T00:22:30") => 64.75,
                    DateTime("2024-01-01T00:30:00") => 68.25,
                    DateTime("2024-01-01T00:37:30") => 69.5,
                    DateTime("2024-01-01T00:45:00") => 68.5,
                    DateTime("2024-01-01T00:52:30") => 69.0,
                    DateTime("2024-01-01T01:00:00") => 71.0,
                    DateTime("2024-01-01T01:07:30") => 72.0,
                    DateTime("2024-01-01T01:15:00") => 72.0,
                    DateTime("2024-01-01T01:22:30") => 68.0,
                    DateTime("2024-01-01T01:30:00") => 60.0,
                    DateTime("2024-01-01T01:37:30") => 57.25,
                    DateTime("2024-01-01T01:45:00") => 59.75,
                    DateTime("2024-01-01T01:52:30") => 62.25,
                    DateTime("2024-01-01T02:00:00") => 64.75,
                    DateTime("2024-01-01T02:07:30") => 67.0,
                    DateTime("2024-01-01T02:15:00") => 68.0)
    @test all([demand.temperature_profile.data[k] - expected[k] for k in keys(expected)] .<
              simulation_parameters["epsilon"]) == true
    @test all([demand.temperature_profile.data[k] - expected[k] for k in keys(demand.temperature_profile.data)] .<
              simulation_parameters["epsilon"]) == true
end

function test_profile_segmentation_third_linear_time_preserving()
    components_config = energy_system_linear_time_preserving()

    simulation_parameters = Dict{String,Any}(
        "time" => 0,
        "current_date" => DateTime("01.01.2024 00:00", "dd.mm.yyy HH:MM"),
        "time_step_seconds" => 300,
        "number_of_time_steps" => 27,
        "start_date" => DateTime("01.01.2024 00:00", "dd.mm.yyy HH:MM"),
        "end_date" => DateTime("01.01.2024 02:15", "dd.mm.yyy HH:MM"),
        "force_profiles_to_repeat" => false,
        "epsilon" => 1e-9,
        "latitude" => nothing,
        "longitude" => nothing,
        "watt_to_wh" => function (w)
            return w * 0.0833333333333
        end,
        "wh_to_watts" => function (w)
            return w * 12.0
        end,
    )

    components = Resie.load_components(components_config, simulation_parameters)
    demand = components["TST_DEM_01"]

    @test demand.energy_profile.time_step == 300
    @test demand.energy_profile.data_type == :extensive
    expected = Dict(DateTime("2024-01-01T00:00:00") => 3.333333333333333,
                    DateTime("2024-01-01T00:05:00") => 3.333333333333333,
                    DateTime("2024-01-01T00:10:00") => 3.5555555555555554,
                    DateTime("2024-01-01T00:15:00") => 3.7777777777777777,
                    DateTime("2024-01-01T00:20:00") => 4.0,
                    DateTime("2024-01-01T00:25:00") => 9.333333333333332,
                    DateTime("2024-01-01T00:30:00") => 14.666666666666666,
                    DateTime("2024-01-01T00:35:00") => 20.0,
                    DateTime("2024-01-01T00:40:00") => 13.333333333333332,
                    DateTime("2024-01-01T00:45:00") => 6.666666666666666,
                    DateTime("2024-01-01T00:50:00") => 0.0,
                    DateTime("2024-01-01T00:55:00") => 0.2777777777777777,
                    DateTime("2024-01-01T01:00:00") => 0.5555555555555555,
                    DateTime("2024-01-01T01:05:00") => 0.8333333333333333,
                    DateTime("2024-01-01T01:10:00") => 3.1111111111111107,
                    DateTime("2024-01-01T01:15:00") => 5.388888888888888,
                    DateTime("2024-01-01T01:20:00") => 7.666666666666666,
                    DateTime("2024-01-01T01:25:00") => 5.666666666666666,
                    DateTime("2024-01-01T01:30:00") => 3.666666666666666,
                    DateTime("2024-01-01T01:35:00") => 1.6666666666666665,
                    DateTime("2024-01-01T01:40:00") => 2.333333333333333,
                    DateTime("2024-01-01T01:45:00") => 2.9999999999999996,
                    DateTime("2024-01-01T01:50:00") => 3.6666666666666665,
                    DateTime("2024-01-01T01:55:00") => 2.5555555555555554,
                    DateTime("2024-01-01T02:00:00") => 1.4444444444444442,
                    DateTime("2024-01-01T02:05:00") => 0.3333333333333333,
                    DateTime("2024-01-01T02:10:00") => 4.888888888888888,
                    DateTime("2024-01-01T02:15:00") => 7.166666666666667)
    @test all([demand.energy_profile.data[k] - expected[k] for k in keys(expected)] .<
              simulation_parameters["epsilon"]) == true
    @test all([demand.energy_profile.data[k] - expected[k] for k in keys(demand.energy_profile.data)] .<
              simulation_parameters["epsilon"]) == true
    @test sum(values(demand.energy_profile.data)) ≈ 136.22222222222222

    @test demand.temperature_profile.time_step == 300
    @test demand.temperature_profile.data_type == :intensive
    expected = Dict(DateTime("2024-01-01T00:00:00") => 65.0,
                    DateTime("2024-01-01T00:05:00") => 65.0,
                    DateTime("2024-01-01T00:10:00") => 64.33333333333333,
                    DateTime("2024-01-01T00:15:00") => 63.666666666666664,
                    DateTime("2024-01-01T00:20:00") => 63.0,
                    DateTime("2024-01-01T00:25:00") => 65.33333333333333,
                    DateTime("2024-01-01T00:30:00") => 67.66666666666666,
                    DateTime("2024-01-01T00:35:00") => 70.0,
                    DateTime("2024-01-01T00:40:00") => 69.33333333333333,
                    DateTime("2024-01-01T00:45:00") => 68.66666666666666,
                    DateTime("2024-01-01T00:50:00") => 68.0,
                    DateTime("2024-01-01T00:55:00") => 69.33333333333333,
                    DateTime("2024-01-01T01:00:00") => 70.66666666666666,
                    DateTime("2024-01-01T01:05:00") => 72.0,
                    DateTime("2024-01-01T01:10:00") => 72.0,
                    DateTime("2024-01-01T01:15:00") => 72.0,
                    DateTime("2024-01-01T01:20:00") => 72.0,
                    DateTime("2024-01-01T01:25:00") => 66.66666666666666,
                    DateTime("2024-01-01T01:30:00") => 61.33333333333333,
                    DateTime("2024-01-01T01:35:00") => 56.0,
                    DateTime("2024-01-01T01:40:00") => 57.66666666666666,
                    DateTime("2024-01-01T01:45:00") => 59.33333333333333,
                    DateTime("2024-01-01T01:50:00") => 61.0,
                    DateTime("2024-01-01T01:55:00") => 62.666666666666664,
                    DateTime("2024-01-01T02:00:00") => 64.33333333333333,
                    DateTime("2024-01-01T02:05:00") => 66.0,
                    DateTime("2024-01-01T02:10:00") => 67.33333333333333,
                    DateTime("2024-01-01T02:15:00") => 68.0)
    @test all([demand.temperature_profile.data[k] - expected[k] for k in keys(expected)] .<
              simulation_parameters["epsilon"]) == true
    @test all([demand.temperature_profile.data[k] - expected[k] for k in keys(demand.temperature_profile.data)] .<
              simulation_parameters["epsilon"]) == true
end

@testset "profiles_linear_time_preserving" begin
    test_profile_segmentation_half_linear_time_preserving()
    test_profile_segmentation_third_linear_time_preserving()
end
