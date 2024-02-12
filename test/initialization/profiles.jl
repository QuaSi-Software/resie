using Debugger
using Test
using Resie
using Resie.EnergySystems
using Resie.Profiles

function energy_system()::Dict{String,Any}
    return Dict{String,Any}(
        "TST_GRI_01" => Dict{String,Any}(
            "type" => "GridConnection",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "output_refs" => ["TST_DEM_01"],
            "is_source" => true,
        ),
        "TST_DEM_01" => Dict{String,Any}(
            "type" => "Demand",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "output_refs" => [],
            "energy_profile_file_path" => "./profiles/tests/heating_demand_short.prf",
            "temperature_profile_file_path" => "./profiles/tests/temperature_short.prf",
            "scale" => 1
        ),
    )
end


function test_profile_aggregation_two()
    components_config = energy_system()

    simulation_parameters = Dict{String,Any}(
        "time_step_seconds" => 1800,
        "time" => 0,
        "epsilon" => 1e-9
    )

    components = Resie.load_components(components_config, simulation_parameters)
    demand = components["TST_DEM_01"]

    @test demand.energy_profile.time_step == 1800
    @test demand.energy_profile.is_power == false
    @test demand.energy_profile.data == [22, 60, 25.5, 16, 43]
    @test sum(demand.energy_profile.data) == 166.5

    @test demand.temperature_profile.time_step == 1800
    @test demand.temperature_profile.is_power == true
    @test demand.temperature_profile.data == [64, 69, 72, 58.5, 68]
end

function test_profile_aggregation_four()
    components_config = energy_system()

    simulation_parameters = Dict{String,Any}(
        "time_step_seconds" => 3600,
        "time" => 0,
        "epsilon" => 1e-9
    )

    components = Resie.load_components(components_config, simulation_parameters)
    demand = components["TST_DEM_01"]

    @test demand.energy_profile.time_step == 3600
    @test demand.energy_profile.is_power == false
    @test demand.energy_profile.data == [82, 41.5, 43]
    @test sum(demand.energy_profile.data) == 166.5

    @test demand.temperature_profile.time_step == 3600
    @test demand.temperature_profile.is_power == true
    @test demand.temperature_profile.data == [66.5, 65.25, 68]
end

function test_profile_segmentation_half()
    components_config = energy_system()

    simulation_parameters = Dict{String,Any}(
        "time_step_seconds" => 450,
        "time" => 0,
        "epsilon" => 1e-9
    )

    components = Resie.load_components(components_config, simulation_parameters)
    demand = components["TST_DEM_01"]

    @test demand.energy_profile.time_step == 450
    @test demand.energy_profile.is_power == false
    @test demand.energy_profile.data == [5,5,6,6,30,30,0,0,1.25,1.25,11.5,11.5,2.5,2.5,5.5,5.5,0.5,0.5,21,21]
    @test sum(demand.energy_profile.data) == 166.5

    @test demand.temperature_profile.time_step == 450
    @test demand.temperature_profile.is_power == true
    @test demand.temperature_profile.data == [65,64,63,66.5,70,69,68,70,72,72,72,64,56,58.5,61,63.5,66,68,70,70]
end

function test_profile_segmentation_third()
    components_config = energy_system()

    simulation_parameters = Dict{String,Any}(
        "time_step_seconds" => 300,
        "time" => 0,
        "epsilon" => 1e-9
    )

    components = Resie.load_components(components_config, simulation_parameters)
    demand = components["TST_DEM_01"]

    @test demand.energy_profile.time_step == 300
    @test demand.energy_profile.is_power == false
    @test all(demand.energy_profile.data .- [10,10,10,12,12,12,60,60,60,0,0,0,2.5,2.5,2.5,23,23,23,5,5,5,11,11,11,1,1,1,42,42,42]./3 .< simulation_parameters["epsilon"]) == true
    @test sum(demand.energy_profile.data) == 166.5

    @test demand.temperature_profile.time_step == 300
    @test demand.temperature_profile.is_power == true
    @test all(demand.temperature_profile.data .- [  65,
                                                64.33333333333334,
                                                63.66666666666667,
                                                63,
                                                65.33333333333334,
                                                67.66666666666667,
                                                70,
                                                69.33333333333334,
                                                68.66666666666667,
                                                68,
                                                69.33333333333334,
                                                70.66666666666667,
                                                72,
                                                72,
                                                72,
                                                72,
                                                66.66666666666667,
                                                61.33333333333334,
                                                56,
                                                57.66666666666667,
                                                59.33333333333334,
                                                61,
                                                62.66666666666667,
                                                64.33333333333334,
                                                66,
                                                67.33333333333334,
                                                68.66666666666667,
                                                70,
                                                70,
                                                70] .< simulation_parameters["epsilon"]) == true
end

@testset "profiles" begin
    test_profile_aggregation_two()
    test_profile_aggregation_four()
    test_profile_segmentation_half()
    test_profile_segmentation_third()
end