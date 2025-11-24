using Debugger
using Test
using Resie
using Resie.EnergySystems

include("../test_util.jl")

function energy_system()::Dict{String,Any}
    return Dict{String,Any}(
        "TST_GRI_01" => Dict{String,Any}(
            "type" => "GridConnection",
            "medium" => "m_e_ac_230v",
            "output_refs" => ["TST_BUS_01"],
            "is_source" => true,
        ),
        "TST_PVP_01" => Dict{String,Any}(
            "type" => "PVPlant",
            "output_refs" => ["TST_BUS_01"],
            "energy_profile_file_path" => "./profiles/tests/source_power_pv.prf",
            "scale" => 20000,
        ),
        "TST_BUS_01" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_e_ac_230v",
        ),
        "TST_BAT_01" => Dict{String,Any}(
            "type" => "Battery",
            "output_refs" => ["TST_BUS_01"],
            "capacity" => 10000,
            "load" => 5000,
        ),
        "TST_DEM_01" => Dict{String,Any}(
            "type" => "Demand",
            "medium" => "m_e_ac_230v",
            "output_refs" => [],
            "energy_profile_file_path" => "./profiles/tests/demand_electricity.prf",
            "scale" => 1,
            "constant_demand" => 1000,
        ),
    )
end

function test_load_no_connections()
    components_config = energy_system()

    simulation_parameters = get_default_sim_params()

    exception_occured = false
    try
        components = Resie.load_components(components_config, simulation_parameters)
        bus = components["TST_BUS_01"]
        @test length(bus.connectivity.input_order) == 0
        @test length(bus.connectivity.output_order) == 0
        @test bus.connectivity.energy_flow === nothing
    catch e
        exception_occured = true
        @test e isa MethodError
        @test String(nameof(e.f)) == "set_storage_transfer!"
    end

    @test exception_occured
end

@testset "load_no_connections" begin
    test_load_no_connections()
end

function test_load_given_lists_empty()
    components_config = energy_system()
    components_config["TST_BUS_01"]["connections"] = Dict{String,Any}(
        "input_order" => [],
        "output_order" => [],
    )

    simulation_parameters = get_default_sim_params()
    exception_occured = false
    try
        components = Resie.load_components(components_config, simulation_parameters)
        bus = components["TST_BUS_01"]
        @test length(bus.connectivity.input_order) == 0
        @test length(bus.connectivity.output_order) == 0
        @test bus.connectivity.energy_flow === nothing
    catch e
        exception_occured = true
        @test e isa MethodError
        @test String(nameof(e.f)) == "set_storage_transfer!"
    end

    @test exception_occured
end

@testset "load_given_lists_empty" begin
    test_load_given_lists_empty()
end

function test_fully_specified()
    components_config = energy_system()
    components_config["TST_BUS_01"]["connections"] = Dict{String,Any}(
        "input_order" => ["TST_PVP_01",
                          "TST_BAT_01",
                          "TST_GRI_01"],
        "output_order" => ["TST_DEM_01",
                           "TST_BAT_01"],
        "energy_flow" => [[1, 1],
                          [1, 0],
                          [1, 0]],
    )

    simulation_parameters = get_default_sim_params()

    components = Resie.load_components(components_config, simulation_parameters)
    bus = components["TST_BUS_01"]
    @test length(bus.connectivity.input_order) == 3
    bus.connectivity.input_order[1] == "TST_PVP_01"
    bus.connectivity.input_order[2] == "TST_BAT_01"
    bus.connectivity.input_order[3] == "TST_GRI_01"
    @test length(bus.connectivity.output_order) == 2
    bus.connectivity.output_order[1] == "TST_DEM_01"
    bus.connectivity.output_order[2] == "TST_BAT_01"
    @test length(bus.connectivity.energy_flow) == 3
    @test bus.connectivity.energy_flow[1][1] != 0
    @test bus.connectivity.energy_flow[1][2] != 0
    @test bus.connectivity.energy_flow[2][1] != 0
    @test bus.connectivity.energy_flow[2][2] == 0
    @test bus.connectivity.energy_flow[3][1] != 0
    @test bus.connectivity.energy_flow[3][2] == 0
end

@testset "fully_specified" begin
    test_fully_specified()
end
