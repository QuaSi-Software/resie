using Debugger
using Test
using Resie
using Resie.EnergySystems

include("../test_util.jl")

function test_base_order()
    components_config = Dict{String,Any}(
        "TST_01_ELT_01_PVP" => Dict{String,Any}(
            "type" => "PVPlant",
            "output_refs" => [
                "TST_01_ELT_01_BUS"
            ],
            "energy_profile_file_path" => "./profiles/tests/source_power_pv.prf",
            "scale" => 20000
        ),
        "TST_01_HZG_01_DEM" => Dict{String,Any}(
            "type" => "Demand",
            "medium" => "m_h_w_ht1",
            "output_refs" => [],
            "energy_profile_file_path" => "./profiles/tests/demand_electricity.prf",
            "temperature_profile_file_path" => "./profiles/tests/demand_heating_temperature.prf",
            "scale" => 10000.0
        ),
        "TST_01_ELT_01_DEM" => Dict{String,Any}(
            "type" => "Demand",
            "medium" => "m_e_ac_230v",
            "output_refs" => [],
            "energy_profile_file_path" => "./profiles/tests/demand_heating_energy.prf",
            "scale" => 15000
        ),
        "TST_01_HZG_01_BUS" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_h_w_ht1",
            "connections" => Dict{String,Any}(
                "input_order" => [
                    "TST_01_HZG_01_CHP",
                    "TST_01_HZG_01_HTP",
                    "TST_01_HZG_01_BFT"
                ],
                "output_order" => [
                    "TST_01_HZG_01_DEM",
                    "TST_01_HZG_01_BFT"
                ]
            )
        ),
        "TST_01_ELT_01_BUS" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_e_ac_230v",
            "connections" => Dict{String,Any}(
                "input_order" => [
                    "TST_01_ELT_01_PVP",
                    "TST_01_HZG_01_CHP",
                    "TST_01_ELT_01_BAT",
                    "TST_01_ELT_01_GRI"
                ],
                "output_order" => [
                    "TST_01_ELT_01_DEM",
                    "TST_01_HZG_01_HTP",
                    "TST_01_ELT_01_BAT",
                    "TST_01_ELT_01_GRO"
                ]
            )
        ),
        "TST_01_HZG_01_CHP" => Dict{String,Any}(
            "type" => "CHPP",
            "output_refs" => [
                "TST_01_HZG_01_BUS",
                "TST_01_ELT_01_BUS"
            ],
            "control_modules" => [
                Dict{String,Any}(
                    "name" => "storage_driven",
                    "high_threshold" => 0.9,
                    "low_threshold" => 0.2,
                    "storage_uac" => "TST_01_HZG_01_BFT"
                )
            ],
            "power_el" => 12500
        ),
        "TST_01_HZG_01_HTP" => Dict{String,Any}(
            "type" => "HeatPump",
            "output_refs" => [
                "TST_01_HZG_01_BUS"
            ],
            "control_modules" => [
                Dict{String,Any}(
                    "name" => "storage_driven",
                    "high_threshold" => 0.5,
                    "low_threshold" => 0.1,
                    "storage_uac" => "TST_01_HZG_01_BFT"
                )
            ],
            "power_th" => 20000,
            "constant_cop" => 3.0
        ),
        "TST_01_HZG_01_BFT" => Dict{String,Any}(
            "type" => "BufferTank",
            "output_refs" => [
                "TST_01_HZG_01_BUS"
            ],
            "capacity" => 40000,
            "load" => 20000
        ),
        "TST_01_ELT_01_BAT" => Dict{String,Any}(
            "type" => "Battery",
            "output_refs" => [
                "TST_01_ELT_01_BUS"
            ],
            "control_modules" => [
                Dict{String,Any}(
                    "name" => "economical_discharge",
                    "pv_threshold" => 750.0,
                    "min_charge" => 0.2,
                    "discharge_limit" => 0.05,
                    "battery_uac" => "TST_01_ELT_01_BAT",
                    "pv_plant_uac" => "TST_01_ELT_01_PVP"
                )
            ],
            "capacity" => 10000,
            "load" => 5000
        ),
        "TST_01_HZG_01_GRI" => Dict{String,Any}(
            "type" => "GridConnection",
            "medium" => "m_c_g_natgas",
            "output_refs" => [
                "TST_01_HZG_01_CHP"
            ],
            "is_source" => true
        ),
        "TST_01_HZG_02_SRC" => Dict{String,Any}(
            "type" => "BoundedSupply",
            "medium" => "m_h_w_lt1",
            "output_refs" => [
                "TST_01_HZG_01_HTP"
            ],
            "max_power_profile_file_path" => "./profiles/tests/source_heat_max_power.prf",
            "temperature_profile_file_path" => "./profiles/tests/source_heat_temperature.prf",
            "scale" => 25000
        ),
        "TST_01_ELT_01_GRI" => Dict{String,Any}(
            "type" => "GridConnection",
            "medium" => "m_e_ac_230v",
            "output_refs" => [
                "TST_01_ELT_01_BUS"
            ],
            "is_source" => true
        ),
        "TST_01_ELT_01_GRO" => Dict{String,Any}(
            "type" => "GridConnection",
            "medium" => "m_e_ac_230v",
            "output_refs" => [],
            "is_source" => false
        )
    )

    expected = [
        [1300, ("TST_01_ELT_01_PVP", EnergySystems.s_reset)],
        [1299, ("TST_01_ELT_01_DEM", EnergySystems.s_reset)],
        [1298, ("TST_01_HZG_01_DEM", EnergySystems.s_reset)],
        [1297, ("TST_01_ELT_01_BUS", EnergySystems.s_reset)],
        [1296, ("TST_01_HZG_01_BUS", EnergySystems.s_reset)],
        [1295, ("TST_01_HZG_01_CHP", EnergySystems.s_reset)],
        [1294, ("TST_01_HZG_01_HTP", EnergySystems.s_reset)],
        [1293, ("TST_01_HZG_01_BFT", EnergySystems.s_reset)],
        [1292, ("TST_01_ELT_01_BAT", EnergySystems.s_reset)],
        [1291, ("TST_01_ELT_01_GRI", EnergySystems.s_reset)],
        [1290, ("TST_01_HZG_01_GRI", EnergySystems.s_reset)],
        [1289, ("TST_01_HZG_02_SRC", EnergySystems.s_reset)],
        [1288, ("TST_01_ELT_01_GRO", EnergySystems.s_reset)],
        [1287, ("TST_01_ELT_01_PVP", EnergySystems.s_control)],
        [1286, ("TST_01_ELT_01_DEM", EnergySystems.s_control)],
        [1285, ("TST_01_HZG_01_DEM", EnergySystems.s_control)],
        [1284, ("TST_01_ELT_01_BUS", EnergySystems.s_control)],
        [1283, ("TST_01_HZG_01_BUS", EnergySystems.s_control)],
        [1282, ("TST_01_HZG_01_CHP", EnergySystems.s_control)],
        [1281, ("TST_01_HZG_01_HTP", EnergySystems.s_control)],
        [1280, ("TST_01_HZG_01_BFT", EnergySystems.s_control)],
        [1279, ("TST_01_ELT_01_BAT", EnergySystems.s_control)],
        [1278, ("TST_01_ELT_01_GRI", EnergySystems.s_control)],
        [1277, ("TST_01_HZG_01_GRI", EnergySystems.s_control)],
        [1276, ("TST_01_HZG_02_SRC", EnergySystems.s_control)],
        [1275, ("TST_01_ELT_01_GRO", EnergySystems.s_control)],
        [1274, ("TST_01_ELT_01_PVP", EnergySystems.s_process)],
        [1273, ("TST_01_ELT_01_DEM", EnergySystems.s_process)],
        [1272, ("TST_01_HZG_01_DEM", EnergySystems.s_process)],
        [1271, ("TST_01_ELT_01_BUS", EnergySystems.s_process)],
        [1270, ("TST_01_HZG_01_BUS", EnergySystems.s_process)],
        [1269, ("TST_01_HZG_01_CHP", EnergySystems.s_process)],
        [1268, ("TST_01_HZG_01_HTP", EnergySystems.s_process)],
        [1267, ("TST_01_HZG_01_BFT", EnergySystems.s_process)],
        [1266, ("TST_01_ELT_01_BAT", EnergySystems.s_process)],
        [1265, ("TST_01_HZG_01_BFT", EnergySystems.s_load)],
        [1264, ("TST_01_ELT_01_BAT", EnergySystems.s_load)],
        [1263, ("TST_01_ELT_01_GRI", EnergySystems.s_process)],
        [1262, ("TST_01_HZG_01_GRI", EnergySystems.s_process)],
        [1261, ("TST_01_HZG_02_SRC", EnergySystems.s_process)],
        [1260, ("TST_01_ELT_01_GRO", EnergySystems.s_process)],
        [1259, ("TST_01_ELT_01_BUS", EnergySystems.s_distribute)],
        [1258, ("TST_01_HZG_01_BUS", EnergySystems.s_distribute)],
    ]

    simulation_parameters = Dict{String,Any}(
        "time_step_seconds" => 900,
        "time" => 0,
        "epsilon" => 1e-9
    )
    
    components = Resie.load_components(components_config, simulation_parameters)
    by_function = Resie.categorize_by_function(components)
    steps = Resie.base_order(by_function)
    @test pwc_steps_astr(expected, steps) == ""
end


@testset "base_order" begin
    test_base_order()
end