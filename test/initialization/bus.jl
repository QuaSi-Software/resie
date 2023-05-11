using Debugger
using Test
using Resie
using Resie.EnergySystems

function system_topology()::Dict{String,Any}
    return Dict{String,Any}(
        "TST_GRI_01" => Dict{String,Any}(
            "type" => "GridConnection",
            "medium" => "m_e_ac_230v",
            "control_refs" => [],
            "output_refs" => ["TST_BUS_01"],
            "is_source" => true,
        ),
        "TST_PVP_01" => Dict{String,Any}(
            "type" => "PVPlant",
            "control_refs" => [],
            "output_refs" => ["TST_BUS_01"],
            "energy_profile_file_path" => "./profiles/tests/source_power_pv.prf",
            "scale" => 20000
        ),
        "TST_BUS_01" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_e_ac_230v",
            "control_refs" => [],
            "output_refs" => ["TST_DEM_01", "TST_BAT_01"],
        ),
        "TST_BAT_01" => Dict{String,Any}(
            "type" => "Battery",
            "control_refs" => ["TST_PVP_01"],
            "output_refs" => ["TST_BUS_01"],
            "strategy" => Dict{String,Any}(
                "name" => "economical_discharge",
                "pv_threshold" => 0.15,
                "min_charge" => 0.2,
                "discharge_limit" => 0.05
            ),
            "capacity" => 10000,
            "load" => 5000
        ),
        "TST_DEM_01" => Dict{String,Any}(
            "type" => "Demand",
            "medium" => "m_e_ac_230v",
            "control_refs" => [],
            "output_refs" => [],
            "energy_profile_file_path" => "./profiles/tests/demand_electricity.prf",
            "scale" => 1,
            "static_load" => 1000,
        ),
    )
end

function test_load_no_connection_matrix()
    systems_config = system_topology()
    systems = Resie.load_systems(systems_config)
    bus = systems["TST_BUS_01"]
    @test length(bus.connectivity.input_order) == 0
    @test length(bus.connectivity.output_order) == 2
    @test bus.connectivity.output_order[1] == "TST_DEM_01"
    @test bus.connectivity.output_order[2] == "TST_BAT_01"
    @test bus.connectivity.storage_loading === nothing
end

@testset "load_no_connection_matrix" begin
    test_load_no_connection_matrix()
end

function test_load_given_lists_empty()
    systems_config = system_topology()
    systems_config["TST_BUS_01"]["connection_matrix"] = Dict{String,Any}(
        "input_order" => [],
        "output_order" => [],
    )
    systems = Resie.load_systems(systems_config)
    bus = systems["TST_BUS_01"]
    @test length(bus.connectivity.input_order) == 0
    @test length(bus.connectivity.output_order) == 0
    @test bus.connectivity.storage_loading === nothing
end

@testset "load_given_lists_empty" begin
    test_load_given_lists_empty()
end

function test_fully_specified()
    systems_config = system_topology()
    systems_config["TST_BUS_01"]["connection_matrix"] = Dict{String,Any}(
        "input_order" => [
            "TST_PVP_01",
            "TST_GRI_01",
        ],
        "output_order" => [
            "TST_DEM_01",
            "TST_BAT_01"
        ],
        "storage_loading" => [
            [1, 1],
            [1, 0],
            [1, 0]
        ]
    )
    systems = Resie.load_systems(systems_config)
    bus = systems["TST_BUS_01"]
    @test length(bus.connectivity.input_order) == 2
    bus.connectivity.input_order[1] == "TST_PVP_01"
    bus.connectivity.input_order[2] == "TST_GRI_01"
    @test length(bus.connectivity.output_order) == 2
    bus.connectivity.output_order[1] == "TST_DEM_01"
    bus.connectivity.output_order[2] == "TST_BAT_01"
    @test length(bus.connectivity.storage_loading) == 3
    @test bus.connectivity.storage_loading[1][1]
    @test bus.connectivity.storage_loading[1][2]
    @test bus.connectivity.storage_loading[2][1]
    @test !bus.connectivity.storage_loading[2][2]
    @test bus.connectivity.storage_loading[3][1]
    @test !bus.connectivity.storage_loading[3][2]
end

@testset "fully_specified" begin
    test_fully_specified()
end