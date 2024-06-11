using Debugger
using Test
using Resie
using Resie.EnergySystems

include("../test_util.jl")

function test_ooo_bus_output_priorities()
    components_config = Dict{String,Any}(
        "TST_GRI_01" => Dict{String,Any}(
            "type" => "GridConnection",
            "medium" => "m_c_g_natgas",
            "control_refs" => [],
            "output_refs" => ["TST_GBO_01"],
            "is_source" => true,
        ),
        "TST_GBO_01" => Dict{String,Any}(
            "type" => "FuelBoiler",
            "m_fuel_in" => "m_c_g_natgas",
            "control_refs" => ["TST_BUS_01"],
            "output_refs" => ["TST_BUS_01"],
            "strategy" => Dict{String,Any}(
                "name" => "demand_driven",
            ),
            "power_th" => 10000,
            "efficiency_fuel_in" => "const:1.0",
        ),
        "TST_BUS_01" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "connections" => Dict{String, Any}(
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
            "connections" => Dict{String, Any}(
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
            "connections" => Dict{String, Any}(
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
            "output_refs" => [],
            "constant_demand" => 1000,
            "constant_temperature" => 60,
            "scale" => 1
        ),
        "TST_DEM_02" => Dict{String,Any}(
            "type" => "Demand",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "output_refs" => [],
            "constant_demand" => 1000,
            "constant_temperature" => 60,
            "scale" => 1
        ),
    )

    expected = [
        ("TST_DEM_02", EnergySystems.s_reset),
        ("TST_DEM_01", EnergySystems.s_reset),
        ("TST_BUS_02", EnergySystems.s_reset),
        ("TST_BUS_01", EnergySystems.s_reset),
        ("Proxy-TST_BUS_01|TST_BUS_03|TST_BUS_02", EnergySystems.s_reset),
        ("TST_BUS_03", EnergySystems.s_reset),
        ("TST_GBO_01", EnergySystems.s_reset),
        ("TST_GRI_01", EnergySystems.s_reset),
        ("TST_DEM_02", EnergySystems.s_control),
        ("TST_DEM_01", EnergySystems.s_control),
        ("TST_BUS_02", EnergySystems.s_control),
        ("TST_BUS_01", EnergySystems.s_control),
        ("Proxy-TST_BUS_01|TST_BUS_03|TST_BUS_02", EnergySystems.s_control),
        ("TST_BUS_03", EnergySystems.s_control),
        ("TST_GBO_01", EnergySystems.s_control),
        ("TST_GRI_01", EnergySystems.s_control),
        ("TST_DEM_02", EnergySystems.s_process),
        ("TST_DEM_01", EnergySystems.s_process),
        ("TST_BUS_02", EnergySystems.s_process),
        ("TST_BUS_01", EnergySystems.s_process),
        ("Proxy-TST_BUS_01|TST_BUS_03|TST_BUS_02", EnergySystems.s_process),
        ("TST_BUS_03", EnergySystems.s_process),
        ("TST_GBO_01", EnergySystems.s_process),
        ("TST_GRI_01", EnergySystems.s_process),
        ("Proxy-TST_BUS_01|TST_BUS_03|TST_BUS_02", EnergySystems.s_distribute),
        ("TST_BUS_02", EnergySystems.s_distribute),
        ("TST_BUS_01", EnergySystems.s_distribute),
        ("TST_BUS_03", EnergySystems.s_distribute),
    ]

    simulation_parameters = Dict{String,Any}(
        "time_step_seconds" => 900,
        "time" => 0,
        "epsilon" => 1e-9
    )

    components = Resie.load_components(components_config, simulation_parameters)
    ooo = Resie.calculate_order_of_operations(components)
    @test pwc_steps_astr(expected, ooo) == ""
end


@testset "ooo_bus_output_priorities" begin
    test_ooo_bus_output_priorities()
end