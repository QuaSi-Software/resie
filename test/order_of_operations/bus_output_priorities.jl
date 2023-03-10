using Debugger
using Test
using Resie
using Resie.EnergySystems

function test_ooo_bus_output_priorities()
    systems_config = Dict{String,Any}(
        "TST_GRI_01" => Dict{String,Any}(
            "type" => "GridConnection",
            "medium" => "m_c_g_natgas",
            "control_refs" => [],
            "production_refs" => ["TST_GBO_01"],
            "is_source" => true,
        ),
        "TST_GBO_01" => Dict{String,Any}(
            "type" => "GasBoiler",
            "control_refs" => ["TST_BUS_01"],
            "production_refs" => ["TST_BUS_01"],
            "strategy" => Dict{String,Any}(
                "name" => "demand_driven",
            ),
            "power" => 10000
        ),
        "TST_BUS_01" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "production_refs" => ["TST_BUS_02", "TST_BUS_03"],
            "connection_matrix" => Dict{String, Any}(
                "input_order" => [
                    "TST_GBO_01",
                ],
                "output_order" => [
                    "TST_BUS_02",
                    "TST_BUS_03"
                ],
            )
        ),
        "TST_BUS_02" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "production_refs" => ["TST_DEM_01"],
            "connection_matrix" => Dict{String, Any}(
                "input_order" => [
                    "TST_BUS_01",
                ],
                "output_order" => [
                    "TST_DEM_01",
                ],
            )
        ),
        "TST_BUS_03" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "production_refs" => ["TST_DEM_02"],
            "connection_matrix" => Dict{String, Any}(
                "input_order" => [
                    "TST_BUS_01",
                ],
                "output_order" => [
                    "TST_DEM_02",
                ],
            )
        ),
        "TST_DEM_01" => Dict{String,Any}(
            "type" => "Demand",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "production_refs" => [],
            "static_load" => 1000,
            "static_temperature" => 60,
            "scale" => 1
        ),
        "TST_DEM_02" => Dict{String,Any}(
            "type" => "Demand",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "production_refs" => [],
            "static_load" => 1000,
            "static_temperature" => 60,
            "scale" => 1
        ),
    )

    expected = [
        ("TST_DEM_02", EnergySystems.s_reset),
        ("TST_DEM_01", EnergySystems.s_reset),
        ("TST_BUS_02", EnergySystems.s_reset),
        ("TST_BUS_01", EnergySystems.s_reset),
        ("TST_BUS_03", EnergySystems.s_reset),
        ("TST_GBO_01", EnergySystems.s_reset),
        ("TST_GRI_01", EnergySystems.s_reset),
        ("TST_DEM_02", EnergySystems.s_control),
        ("TST_DEM_01", EnergySystems.s_control),
        ("TST_BUS_02", EnergySystems.s_control),
        ("TST_BUS_01", EnergySystems.s_control),
        ("TST_BUS_03", EnergySystems.s_control),
        ("TST_GBO_01", EnergySystems.s_control),
        ("TST_GRI_01", EnergySystems.s_control),
        ("TST_DEM_02", EnergySystems.s_produce),
        ("TST_DEM_01", EnergySystems.s_produce),
        ("TST_BUS_02", EnergySystems.s_produce),
        ("TST_BUS_01", EnergySystems.s_produce),
        ("TST_BUS_03", EnergySystems.s_produce),
        ("TST_GBO_01", EnergySystems.s_produce),
        ("TST_GRI_01", EnergySystems.s_produce),
        ("TST_BUS_02", EnergySystems.s_distribute),
        ("TST_BUS_03", EnergySystems.s_distribute),
        ("TST_BUS_01", EnergySystems.s_distribute),
    ]

    systems = Resie.load_systems(systems_config)
    ooo = Resie.order_of_operations(systems)
    @test all(ooo .== expected)
end


@testset "ooo_bus_output_priorities" begin
    test_ooo_bus_output_priorities()
end