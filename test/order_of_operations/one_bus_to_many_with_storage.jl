
using Debugger
using Test
using Resie
using Resie.EnergySystems

include("../test_util.jl")

function test_ooo_one_bus_to_many_with_storage()
    components_config = Dict{String,Any}(
        "TST_SRC_01" => Dict{String,Any}(
            "type" => "BoundedSupply",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "output_refs" => ["TST_BUS_01"],
            "is_source" => true,
            "static_power" => 400,
        ),
        "TST_BUS_01" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "output_refs" => [
                "TST_BUS_02",
                "TST_BUS_02",
            ],
            "connection_matrix" => Dict{String, Any}(
                "input_order" => [
                    "TST_SRC_01",
                ],
                "output_order" => [
                    "TST_BUS_02",
                    "TST_BUS_03",
                ],
            )
        ),
        "TST_BUS_02" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "output_refs" => [
                "TST_DEM_01",
                "TST_TES_01",
            ],
            "connection_matrix" => Dict{String, Any}(
                "input_order" => [
                    "TST_BUS_01",
                    "TST_TES_01",
                ],
                "output_order" => [
                    "TST_DEM_01",
                    "TST_TES_01",
                ],
            )
        ),
        "TST_BUS_03" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "output_refs" => [
                "TST_DEM_02",
                "TST_TES_02",
            ],
            "connection_matrix" => Dict{String, Any}(
                "input_order" => [
                    "TST_BUS_01",
                    "TST_TES_02",
                ],
                "output_order" => [
                    "TST_DEM_02",
                    "TST_TES_02",
                ],
            )
        ),
        "TST_TES_01" => Dict{String,Any}(
            "type" => "Storage",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "output_refs" => ["TST_BUS_02"],
            "capacity" => 1000,
            "load" => 500,
        ),
        "TST_TES_02" => Dict{String,Any}(
            "type" => "Storage",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "output_refs" => ["TST_BUS_03"],
            "capacity" => 1000,
            "load" => 500,
        ),
        "TST_DEM_01" => Dict{String,Any}(
            "type" => "Demand",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "output_refs" => [],
            "static_load" => 100.0,
            "static_temperature" => 55.0,
        ),
        "TST_DEM_02" => Dict{String,Any}(
            "type" => "Demand",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "output_refs" => [],
            "static_load" => 100.0,
            "static_temperature" => 55.0,
        ),
    )

    expected = [
        ("TST_DEM_02", EnergySystems.s_reset),
        ("TST_DEM_01", EnergySystems.s_reset),
        ("TST_BUS_02", EnergySystems.s_reset),
        ("TST_BUS_01", EnergySystems.s_reset),
        ("TST_BUS_03", EnergySystems.s_reset),
        ("TST_TES_02", EnergySystems.s_reset),
        ("TST_TES_01", EnergySystems.s_reset),
        ("TST_SRC_01", EnergySystems.s_reset),
        ("TST_DEM_02", EnergySystems.s_control),
        ("TST_DEM_01", EnergySystems.s_control),
        ("TST_BUS_02", EnergySystems.s_control),
        ("TST_BUS_01", EnergySystems.s_control),
        ("TST_BUS_03", EnergySystems.s_control),
        ("TST_TES_02", EnergySystems.s_control),
        ("TST_TES_01", EnergySystems.s_control),
        ("TST_SRC_01", EnergySystems.s_control),
        ("TST_DEM_02", EnergySystems.s_process),
        ("TST_DEM_01", EnergySystems.s_process),
        ("TST_BUS_02", EnergySystems.s_process),
        ("TST_BUS_01", EnergySystems.s_process),
        ("TST_BUS_03", EnergySystems.s_process),
        ("TST_TES_01", EnergySystems.s_process),
        ("TST_TES_02", EnergySystems.s_process),
        ("TST_TES_01", EnergySystems.s_load),
        ("TST_TES_02", EnergySystems.s_load),
        ("TST_SRC_01", EnergySystems.s_process),
        ("TST_BUS_02", EnergySystems.s_distribute),
        ("TST_BUS_03", EnergySystems.s_distribute),
        ("TST_BUS_01", EnergySystems.s_distribute),
    ]

    components = Resie.load_components(components_config)
    ooo = Resie.calculate_order_of_operations(components)
    @test pwc_steps_astr(expected, ooo) == ""
end

@testset "ooo_one_bus_to_many_with_storage" begin
    test_ooo_one_bus_to_many_with_storage()
end