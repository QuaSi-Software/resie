using Debugger
using Test
using Resie
using Resie.EnergySystems

include("../test_util.jl")

function test_ooo_storage_loading_switch()
    components_config = Dict{String,Any}(
        "TST_GRI_01" => Dict{String,Any}(
            "type" => "GridConnection",
            "medium" => "m_c_g_natgas",
            "control_refs" => [],
            "output_refs" => ["TST_GBO_01"],
            "is_source" => true,
        ),
        "TST_GRI_02" => Dict{String,Any}(
            "type" => "GridConnection",
            "medium" => "m_c_g_natgas",
            "control_refs" => [],
            "output_refs" => ["TST_GBO_02"],
            "is_source" => true,
        ),
        "TST_GBO_01" => Dict{String,Any}(
            "type" => "FuelBoiler",
            "m_fuel_in" => "m_c_g_natgas",
            "control_refs" => ["TST_BFT_01"],
            "output_refs" => [
                "TST_BUS_01"
            ],
            "strategy" => Dict{String,Any}(
                "name" => "demand_driven"
            ),
            "power_th" => 40000
        ),
        "TST_GBO_02" => Dict{String,Any}(
            "type" => "FuelBoiler",
            "m_fuel_in" => "m_c_g_natgas",
            "control_refs" => ["TST_BFT_01"],
            "output_refs" => [
                "TST_BUS_01"
            ],
            "strategy" => Dict{String,Any}(
                "name" => "demand_driven"
            ),
            "power_th" => 40000
        ),
        "TST_BUS_01" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "connections" => Dict{String, Any}(
                "input_order" => [
                    "TST_GBO_01",
                    "TST_BFT_01",
                    "TST_GBO_02"
                ],
                "output_order" => [
                    "TST_DEM_01",
                    "TST_BFT_01"
                ],
                "energy_flow" => [
                    [1, 1],
                    [1, 0],
                    [1, 0]
                ]
            )
        ),
        "TST_BFT_01" => Dict{String,Any}(
            "type" => "BufferTank",
            "control_refs" => [],
            "output_refs" => [
                "TST_BUS_01"
            ],
            "capacity" => 40000,
            "load" => 0
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
        ("TST_BUS_01", EnergySystems.s_reset),
        ("TST_GBO_02", EnergySystems.s_reset),
        ("TST_GBO_01", EnergySystems.s_reset),
        ("TST_BFT_01", EnergySystems.s_reset),
        ("TST_GRI_02", EnergySystems.s_reset),
        ("TST_GRI_01", EnergySystems.s_reset),
        ("TST_DEM_01", EnergySystems.s_control),
        ("TST_BUS_01", EnergySystems.s_control),
        ("TST_GBO_02", EnergySystems.s_control),
        ("TST_GBO_01", EnergySystems.s_control),
        ("TST_BFT_01", EnergySystems.s_control),
        ("TST_GRI_02", EnergySystems.s_control),
        ("TST_GRI_01", EnergySystems.s_control),
        ("TST_DEM_01", EnergySystems.s_process),
        ("TST_BUS_01", EnergySystems.s_process),
        ("TST_GBO_01", EnergySystems.s_potential),
        ("TST_GBO_02", EnergySystems.s_potential),
        ("TST_GBO_01", EnergySystems.s_process),
        ("TST_BFT_01", EnergySystems.s_process),
        ("TST_GBO_02", EnergySystems.s_process),
        ("TST_BFT_01", EnergySystems.s_load),
        ("TST_GRI_02", EnergySystems.s_process),
        ("TST_GRI_01", EnergySystems.s_process),
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

@testset "ooo_storage_loading_switch" begin
    test_ooo_storage_loading_switch()
end