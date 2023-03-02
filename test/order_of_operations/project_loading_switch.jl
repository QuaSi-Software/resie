using Debugger
using Test
using Resie
using Resie.EnergySystems

function test_ooo_storage_loading_switch()
    systems_config = Dict{String,Any}(
        "TST_GRI_01" => Dict{String,Any}(
            "type" => "GridConnection",
            "medium" => "m_c_g_natgas",
            "control_refs" => [],
            "production_refs" => ["TST_GBO_01"],
            "is_source" => true,
        ),
        "TST_GRI_02" => Dict{String,Any}(
            "type" => "GridConnection",
            "medium" => "m_c_g_natgas",
            "control_refs" => [],
            "production_refs" => ["TST_GBO_02"],
            "is_source" => true,
        ),
        "TST_GBO_01" => Dict{String,Any}(
            "type" => "GasBoiler",
            "control_refs" => ["TST_BFT_01"],
            "production_refs" => [
                "TST_BUS_01"
            ],
            "strategy" => Dict{String,Any}(
                "name" => "demand_driven"
            ),
            "power" => 40000
        ),
        "TST_GBO_02" => Dict{String,Any}(
            "type" => "GasBoiler",
            "control_refs" => ["TST_BFT_01"],
            "production_refs" => [
                "TST_BUS_01"
            ],
            "strategy" => Dict{String,Any}(
                "name" => "demand_driven"
            ),
            "power" => 40000
        ),
        "TST_BUS_01" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "production_refs" => ["TST_DEM_01", "TST_BFT_01"],
            "input_priorities" => ["TST_GBO_01", "TST_BFT_01", "TST_GBO_02"]
        ),
        "TST_BFT_01" => Dict{String,Any}(
            "type" => "BufferTank",
            "control_refs" => [],
            "production_refs" => [
                "TST_BUS_01"
            ],
            "capacity" => 40000,
            "load" => 0
        ),
        "TST_DEM_01" => Dict{String,Any}(
            "type" => "Demand",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "production_refs" => [],
            "energy_profile_file_path" => "./profiles/tests/demand_heating_energy.prf",
            "temperature_profile_file_path" => "./profiles/tests/demand_heating_temperature.prf",
            "scale" => 1000
        ),
    )

    expected = [
        ("TST_DEM_01", EnergySystems.s_reset),
        ("TST_BUS_01", EnergySystems.s_reset),
        ("TST_GBO_02", EnergySystems.s_reset),
        ("TST_GBO_01", EnergySystems.s_reset),
        ("TST_BFT_01", EnergySystems.s_reset),
        ("TST_GRI_02", EnergySystems.s_reset),
        ("TST_GRI_01", EnergySystems.s_reset),
        ("TST_DEM_01", EnergySystems.s_control),
        ("TST_BUS_01", EnergySystems.s_control),
        ("TST_GBO_01", EnergySystems.s_control),
        ("TST_BFT_01", EnergySystems.s_control),
        ("TST_GBO_02", EnergySystems.s_control),
        ("TST_GRI_02", EnergySystems.s_control),
        ("TST_GRI_01", EnergySystems.s_control),
        ("TST_DEM_01", EnergySystems.s_produce),
        ("TST_BUS_01", EnergySystems.s_produce),
        ("TST_GBO_01", EnergySystems.s_produce),
        ("TST_BFT_01", EnergySystems.s_produce),
        ("TST_GBO_02", EnergySystems.s_produce),
        ("TST_BFT_01", EnergySystems.s_load),
        ("TST_GRI_02", EnergySystems.s_produce),
        ("TST_GRI_01", EnergySystems.s_produce),
        ("TST_BUS_01", EnergySystems.s_distribute),
    ]

    systems = Resie.load_systems(systems_config)
    ooo = Resie.order_of_operations(systems)
    @test all(ooo .== expected)
end

@testset "ooo_storage_loading_switch" begin
    test_ooo_storage_loading_switch()
end