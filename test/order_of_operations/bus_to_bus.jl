
using Debugger
using Test
using Resie
using Resie.EnergySystems

include("../test_util.jl")

function test_ooo_bus_to_bus()
    components_config = Dict{String,Any}(
        "TST_GRI_01" => Dict{String,Any}(
            "type" => "GridConnection",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "output_refs" => ["TST_BUS_01"],
            "is_source" => true,
        ),
        "TST_BUS_01" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "connections" => Dict{String, Any}(
                "input_order" => [
                    "TST_BFT_01",
                    "TST_GRI_01"
                ],
                "output_order" => [
                    "TST_BUS_02",
                    "TST_BFT_01"
                ],
            )
        ),
        "TST_BFT_01" => Dict{String,Any}(
            "type" => "BufferTank",
            "control_refs" => [],
            "output_refs" => [
                "TST_BUS_01"
            ],
            "capacity" => 40000,
            "load" => 20000
        ),
        "TST_BUS_02" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "connections" => Dict{String, Any}(
                "input_order" => [
                    "TST_BFT_02",
                    "TST_BUS_01"
                ],
                "output_order" => [
                    "TST_DEM_01",
                    "TST_BFT_02"
                ],
            )
        ),
        "TST_BFT_02" => Dict{String,Any}(
            "type" => "BufferTank",
            "control_refs" => [],
            "output_refs" => [
                "TST_BUS_02"
            ],
            "capacity" => 20000,
            "load" => 10000
        ),
        "TST_DEM_01" => Dict{String,Any}(
            "type" => "Demand",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "output_refs" => [],
            "energy_profile_file_path" => "./profiles/tests/demand_heating_energy.prf",
            "temperature_profile_file_path" => "./profiles/tests/demand_heating_temperature.prf",
            "scale" => 1000
        ),
    )

    expected = [
        ("TST_DEM_01", EnergySystems.s_reset),
        ("TST_BUS_02", EnergySystems.s_reset),
        ("TST_BUS_01", EnergySystems.s_reset),
        ("Proxy-TST_BUS_01|TST_BUS_02", EnergySystems.s_reset),
        ("TST_BFT_02", EnergySystems.s_reset),
        ("TST_BFT_01", EnergySystems.s_reset),
        ("TST_GRI_01", EnergySystems.s_reset),
        ("TST_DEM_01", EnergySystems.s_control),
        ("TST_BUS_02", EnergySystems.s_control),
        ("TST_BUS_01", EnergySystems.s_control),
        ("Proxy-TST_BUS_01|TST_BUS_02", EnergySystems.s_control),
        ("TST_BFT_02", EnergySystems.s_control),
        ("TST_BFT_01", EnergySystems.s_control),
        ("TST_GRI_01", EnergySystems.s_control),
        ("TST_DEM_01", EnergySystems.s_process),
        ("TST_BUS_02", EnergySystems.s_process),
        ("Proxy-TST_BUS_01|TST_BUS_02", EnergySystems.s_process),
        ("TST_BFT_02", EnergySystems.s_process),
        ("TST_BFT_01", EnergySystems.s_process),
        ("TST_BUS_01", EnergySystems.s_process),
        ("TST_BFT_02", EnergySystems.s_load),
        ("TST_BFT_01", EnergySystems.s_load),
        ("TST_GRI_01", EnergySystems.s_process),
        ("Proxy-TST_BUS_01|TST_BUS_02", EnergySystems.s_distribute),
        ("TST_BUS_02", EnergySystems.s_distribute),
        ("TST_BUS_01", EnergySystems.s_distribute),
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

@testset "ooo_bus_to_bus" begin
    test_ooo_bus_to_bus()
end