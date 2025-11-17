using Debugger
using Test
using Resie
using Resie.EnergySystems
using Resie.Profiles

include("../../test_util.jl")

function test_one_to_one_grid()
    components_config = Dict{String,Any}(
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
            "constant_demand" => 4000,
            "constant_temperature" => 55,
        ),
    )

    simulation_parameters = get_default_sim_params()

    components = Resie.load_components(components_config, simulation_parameters)
    demand = components["TST_DEM_01"]
    grid = components["TST_GRI_01"]

    EnergySystems.reset(demand)
    EnergySystems.reset(grid)

    @test demand.demand == 0.0
    @test demand.temperature === nothing

    EnergySystems.control(demand, components, simulation_parameters)
    EnergySystems.control(grid, components, simulation_parameters)

    @test demand.demand == 1000.0
    @test demand.temperature == 55.0

    EnergySystems.process(demand, simulation_parameters)

    @test demand.input_interfaces[demand.medium].balance == -1000.0
    @test demand.input_interfaces[demand.medium].max_energy.temperature_min == [55.0]

    EnergySystems.process(grid, simulation_parameters)

    @test demand.input_interfaces[demand.medium].balance == 0.0
    @test demand.input_interfaces[demand.medium].max_energy.temperature_min == [55.0]
end

@testset "one_to_one_grid" begin
    test_one_to_one_grid()
end

function test_one_to_one_flexible_source()
    components_config = Dict{String,Any}(
        "TST_SRC_01" => Dict{String,Any}(
            "type" => "FlexibleSupply",
            "medium" => "m_h_w_ht1",
            "output_refs" => ["TST_DEM_01"],
            "constant_temperature" => 50,
            "constant_power" => 4000,
        ),
        "TST_DEM_01" => Dict{String,Any}(
            "type" => "Demand",
            "medium" => "m_h_w_ht1",
            "output_refs" => [],
            "constant_demand" => 4000,
            "constant_temperature" => 55,
        ),
    )

    simulation_parameters = get_default_sim_params()

    components = Resie.load_components(components_config, simulation_parameters)
    demand = components["TST_DEM_01"]
    source = components["TST_SRC_01"]

    EnergySystems.reset(demand)
    EnergySystems.reset(source)

    @test demand.demand == 0.0
    @test demand.temperature === nothing

    EnergySystems.control(demand, components, simulation_parameters)
    EnergySystems.control(source, components, simulation_parameters)
    EnergySystems.process(demand, simulation_parameters)

    @test demand.input_interfaces[demand.medium].balance == -1000.0
    @test demand.input_interfaces[demand.medium].max_energy.temperature_min == [55.0]

    EnergySystems.process(source, simulation_parameters)

    @test demand.input_interfaces[demand.medium].balance == -1000.0
    @test demand.input_interfaces[demand.medium].max_energy.temperature_min == [55.0]
end

@testset "one_to_one_flexible_source" begin
    test_one_to_one_flexible_source()
end
