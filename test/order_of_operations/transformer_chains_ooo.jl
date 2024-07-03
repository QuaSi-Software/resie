using Debugger
using Test
using Resie
using Resie.EnergySystems

include("../test_util.jl")

function test_base_order()
    components_config = Dict{String,Any}(
        "TST_01_ELT_01_PVP" => Dict{String,Any}(
            "type" => "PVPlant",
            "control_refs" => [],
            "output_refs" => [
                "TST_01_ELT_01_BUS"
            ],
            "energy_profile_file_path" => "./profiles/tests/source_power_pv.prf",
            "scale" => 20000
        ),
        "TST_01_HZG_01_DEM" => Dict{String,Any}(
            "type" => "Demand",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "output_refs" => [],
            "energy_profile_file_path" => "./profiles/tests/demand_electricity.prf",
            "temperature_profile_file_path" => "./profiles/tests/demand_heating_temperature.prf",
            "scale" => 10000.0
        ),
        "TST_01_ELT_01_DEM" => Dict{String,Any}(
            "type" => "Demand",
            "medium" => "m_e_ac_230v",
            "control_refs" => [],
            "output_refs" => [],
            "energy_profile_file_path" => "./profiles/tests/demand_heating_energy.prf",
            "scale" => 15000
        ),
        "TST_01_HZG_01_BUS" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
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
            "control_refs" => [],
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
            "control_refs" => ["TST_01_HZG_01_BFT"],
            "output_refs" => [
                "TST_01_HZG_01_BUS",
                "TST_01_ELT_01_BUS"
            ],
            "strategy" => Dict{String,Any}(
                "name" => "storage_driven",
                "high_threshold" => 0.9,
                "low_threshold" => 0.2
            ),
            "power_gas" => 12500
        ),
        "TST_01_HZG_01_HTP" => Dict{String,Any}(
            "type" => "HeatPump",
            "control_refs" => ["TST_01_HZG_01_BFT"],
            "output_refs" => [
                "TST_01_HZG_01_BUS"
            ],
            "strategy" => Dict{String,Any}(
                "name" => "storage_driven",
                "high_threshold" => 0.5,
                "low_threshold" => 0.1
            ),
            "power_th" => 20000,
            "constant_cop" => 3.0
        ),
        "TST_01_HZG_01_BFT" => Dict{String,Any}(
            "type" => "BufferTank",
            "control_refs" => [],
            "output_refs" => [
                "TST_01_HZG_01_BUS"
            ],
            "capacity" => 40000,
            "load" => 20000
        ),
        "TST_01_ELT_01_BAT" => Dict{String,Any}(
            "type" => "Battery",
            "control_refs" => ["TST_01_ELT_01_PVP"],
            "output_refs" => [
                "TST_01_ELT_01_BUS"
            ],
            "strategy" => Dict{String,Any}(
                "name" => "economical_discharge",
                "pv_threshold" => 0.15,
                "min_charge" => 0.2,
                "discharge_limit" => 0.05
            ),
            "capacity" => 10000,
            "load" => 5000
        ),
        "TST_01_HZG_01_GRI" => Dict{String,Any}(
            "type" => "GridConnection",
            "medium" => "m_c_g_natgas",
            "control_refs" => [],
            "output_refs" => [
                "TST_01_HZG_01_CHP"
            ],
            "is_source" => true
        ),
        "TST_01_HZG_02_SRC" => Dict{String,Any}(
            "type" => "BoundedSupply",
            "medium" => "m_h_w_lt1",
            "control_refs" => [],
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
            "control_refs" => [],
            "output_refs" => [
                "TST_01_ELT_01_BUS"
            ],
            "is_source" => true
        ),
        "TST_01_ELT_01_GRO" => Dict{String,Any}(
            "type" => "GridConnection",
            "medium" => "m_e_ac_230v",
            "control_refs" => [],
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
        [1269, ("TST_01_HZG_01_CHP", EnergySystems.s_potential)],
        [1268, ("TST_01_HZG_01_HTP", EnergySystems.s_potential)],
        [1267, ("TST_01_HZG_01_CHP", EnergySystems.s_process)],
        [1266, ("TST_01_HZG_01_HTP", EnergySystems.s_process)],
        [1265, ("TST_01_HZG_01_BFT", EnergySystems.s_process)],
        [1264, ("TST_01_ELT_01_BAT", EnergySystems.s_process)],
        [1263, ("TST_01_HZG_01_BFT", EnergySystems.s_load)],
        [1262, ("TST_01_ELT_01_BAT", EnergySystems.s_load)],
        [1261, ("TST_01_ELT_01_GRI", EnergySystems.s_process)],
        [1260, ("TST_01_HZG_01_GRI", EnergySystems.s_process)],
        [1259, ("TST_01_HZG_02_SRC", EnergySystems.s_process)],
        [1258, ("TST_01_ELT_01_GRO", EnergySystems.s_process)],
        [1257, ("TST_01_ELT_01_BUS", EnergySystems.s_distribute)],
        [1256, ("TST_01_HZG_01_BUS", EnergySystems.s_distribute)],
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

function test_ooo_middle_bus()
    components_config = Dict{String,Any}(
        "TST_DEM_01" => Dict{String,Any}(
            "type" => "Demand",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "output_refs" => [],
            "energy_profile_file_path" => "./profiles/tests/demand_heating_energy.prf",
            "constant_temperature" => 95,
            "scale" => 400
        ),
        "TST_DEM_02" => Dict{String,Any}(
            "type" => "Demand",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "output_refs" => [],
            "energy_profile_file_path" => "./profiles/tests/demand_heating_energy.prf",
            "constant_temperature" => 95,
            "scale" => 400
        ),
        "TST_SRC_01" => Dict{String,Any}(
            "type" => "BoundedSupply",
            "medium" => "m_h_w_lt1",
            "control_refs" => [],
            "output_refs" => [
                "TST_HP_01"
            ],
            "max_power_profile_file_path" => "./profiles/tests/demand_electricity.prf",
            "constant_temperature" => 5,
            "scale" => 1000
        ),
        "TST_SRC_02" => Dict{String,Any}(
            "type" => "BoundedSupply",
            "medium" => "m_h_w_lt1",
            "control_refs" => [],
            "output_refs" => [
                "TST_HP_03"
            ],
            "max_power_profile_file_path" => "./profiles/tests/demand_electricity.prf",
            "constant_temperature" => 5,
            "scale" => 1000
        ),
        "TST_GRI_EL" => Dict{String,Any}(
            "type" => "GridConnection",
            "medium" => "m_e_ac_230v",
            "control_refs" => [],
            "output_refs" => [
                "TST_BUS_EL"
            ],
            "is_source" => true
        ),
        "TST_HP_01" => Dict{String,Any}(
            "type" => "HeatPump",
            "m_heat_in" => "m_h_w_lt1",
            "m_heat_out" => "m_h_w_ht1",
            "control_refs" => [],
            "output_refs" => [
                "TST_HP_02"
            ],
            "strategy" => Dict{String,Any}(
                "name" => "demand_driven"
            ),
            "power_th" => 1200,
            "output_temperature" => 10,
            "min_power_fraction" => 0.0
        ),
        "TST_HP_02" => Dict{String,Any}(
            "type" => "HeatPump",
            "m_heat_in" => "m_h_w_ht1",
            "m_heat_out" => "m_h_w_ht1",
            "control_refs" => [],
            "output_refs" => [
                "TST_BUS_TH"
            ],
            "strategy" => Dict{String,Any}(
                "name" => "demand_driven"
            ),
            "power_th" => 1200,
            "output_temperature" => 25,
            "min_power_fraction" => 0.0
        ),
        "TST_HP_03" => Dict{String,Any}(
            "type" => "HeatPump",
            "m_heat_in" => "m_h_w_lt1",
            "m_heat_out" => "m_h_w_ht1",
            "control_refs" => [],
            "output_refs" => [
                "TST_HP_04"
            ],
            "strategy" => Dict{String,Any}(
                "name" => "demand_driven"
            ),
            "power_th" => 1200,
            "output_temperature" => 20,
            "min_power_fraction" => 0.0
        ),
        "TST_HP_04" => Dict{String,Any}(
            "type" => "HeatPump",
            "m_heat_in" => "m_h_w_ht1",
            "m_heat_out" => "m_h_w_ht1",
            "control_refs" => [],
            "output_refs" => [
                "TST_BUS_TH"
            ],
            "strategy" => Dict{String,Any}(
                "name" => "demand_driven"
            ),
            "power_th" => 1200,
            "output_temperature" => 30,
            "min_power_fraction" => 0.0
        ),
        "TST_HP_05" => Dict{String,Any}(
            "type" => "HeatPump",
            "m_heat_in" => "m_h_w_ht1",
            "m_heat_out" => "m_h_w_ht1",
            "control_refs" => [],
            "output_refs" => [
                "TST_HP_06"
            ],
            "strategy" => Dict{String,Any}(
                "name" => "demand_driven"
            ),
            "power_th" => 1200,
            "output_temperature" => 60,
            "min_power_fraction" => 0.0
        ),
        "TST_HP_06" => Dict{String,Any}(
            "type" => "HeatPump",
            "m_heat_in" => "m_h_w_ht1",
            "m_heat_out" => "m_h_w_ht1",
            "control_refs" => [],
            "output_refs" => [
                "TST_DEM_01"
            ],
            "strategy" => Dict{String,Any}(
                "name" => "demand_driven"
            ),
            "power_th" => 1200,
            "min_power_fraction" => 0.0
        ),
        "TST_HP_07" => Dict{String,Any}(
            "type" => "HeatPump",
            "m_heat_in" => "m_h_w_ht1",
            "m_heat_out" => "m_h_w_ht1",
            "control_refs" => [],
            "output_refs" => [
                "TST_HP_08"
            ],
            "strategy" => Dict{String,Any}(
                "name" => "demand_driven"
            ),
            "power_th" => 1200,
            "output_temperature" => 70,
            "min_power_fraction" => 0.0
        ),
        "TST_HP_08" => Dict{String,Any}(
            "type" => "HeatPump",
            "m_heat_in" => "m_h_w_ht1",
            "m_heat_out" => "m_h_w_ht1",
            "control_refs" => [],
            "output_refs" => [
                "TST_DEM_02"
            ],
            "strategy" => Dict{String,Any}(
                "name" => "demand_driven"
            ),
            "power_th" => 1200,
            "min_power_fraction" => 0.0
        ),
        "TST_BUS_EL" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_e_ac_230v",
            "control_refs" => [],
            "connections" => Dict{String,Any}(
                "input_order" => [
                    "TST_GRI_EL"
                ],
                "output_order" => [
                    "TST_HP_01",
                    "TST_HP_02",
                    "TST_HP_03",
                    "TST_HP_04",
                    "TST_HP_05",
                    "TST_HP_06",
                    "TST_HP_07",
                    "TST_HP_08"
                ],
                "energy_flow" => [
                    [1, 1, 1, 1, 1, 1, 1, 1]
                ]
            )
        ),
        "TST_BUS_TH" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "connections" => Dict{String,Any}(
                "input_order" => [
                    "TST_HP_02",
                    "TST_HP_04"
                ],
                "output_order" => [
                    "TST_HP_07",
                    "TST_HP_05"
                ],
                "energy_flow" => [
                    [1, 1],
                    [1, 1]
                ]
            )
        )
    )

    expected = [
        "TST_DEM_01 s_reset",
        "TST_DEM_02 s_reset",
        "TST_BUS_TH s_reset",
        "TST_BUS_EL s_reset",
        "TST_HP_01 s_reset",
        "TST_HP_06 s_reset",
        "TST_HP_03 s_reset",
        "TST_HP_05 s_reset",
        "TST_HP_04 s_reset",
        "TST_HP_02 s_reset",
        "TST_HP_08 s_reset",
        "TST_HP_07 s_reset",
        "TST_SRC_01 s_reset",
        "TST_SRC_02 s_reset",
        "TST_GRI_EL s_reset",
        "TST_DEM_01 s_control",
        "TST_DEM_02 s_control",
        "TST_BUS_TH s_control",
        "TST_BUS_EL s_control",
        "TST_HP_01 s_control",
        "TST_HP_06 s_control",
        "TST_HP_03 s_control",
        "TST_HP_05 s_control",
        "TST_HP_04 s_control",
        "TST_HP_02 s_control",
        "TST_HP_08 s_control",
        "TST_HP_07 s_control",
        "TST_SRC_01 s_control",
        "TST_SRC_02 s_control",
        "TST_GRI_EL s_control",
        "TST_DEM_01 s_process",
        "TST_DEM_02 s_process",
        "TST_BUS_TH s_process",
        "TST_BUS_EL s_process",
        "TST_HP_01 s_potential",
        "TST_HP_02 s_potential",
        "TST_HP_03 s_potential",
        "TST_HP_04 s_potential",
        "TST_HP_08 s_potential",
        "TST_HP_07 s_potential",
        "TST_HP_06 s_potential",
        "TST_HP_05 s_potential",
        "TST_HP_02 s_process",
        "TST_HP_01 s_process",
        "TST_HP_04 s_process",
        "TST_HP_03 s_process",
        "TST_HP_07 s_process",
        "TST_HP_08 s_process",
        "TST_HP_05 s_process",
        "TST_HP_06 s_process",
        "TST_SRC_01 s_process",
        "TST_SRC_02 s_process",
        "TST_GRI_EL s_process",
        "TST_BUS_TH s_distribute",
        "TST_BUS_EL s_distribute"
    ]

    simulation_parameters = Dict{String,Any}(
        "time_step_seconds" => 900,
        "time" => 0,
        "epsilon" => 1e-9
    )
    
    components = Resie.load_components(components_config, simulation_parameters)
    steps = Resie.calculate_order_of_operations(components)
    ooo = []
    for step in steps
        push!(ooo, "$(step[1]) $(step[2])")
    end
    @test pwc_ooo_astr(expected, ooo) == ""
end

function test_ooo_middle_bus_different_order()
    components_config = Dict{String,Any}(
        "TST_DEM_01" => Dict{String,Any}(
            "type" => "Demand",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "output_refs" => [],
            "energy_profile_file_path" => "./profiles/tests/demand_heating_energy.prf",
            "constant_temperature" => 95,
            "scale" => 400
        ),
        "TST_DEM_02" => Dict{String,Any}(
            "type" => "Demand",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "output_refs" => [],
            "energy_profile_file_path" => "./profiles/tests/demand_heating_energy.prf",
            "constant_temperature" => 95,
            "scale" => 400
        ),
        "TST_SRC_01" => Dict{String,Any}(
            "type" => "BoundedSupply",
            "medium" => "m_h_w_lt1",
            "control_refs" => [],
            "output_refs" => [
                "TST_HP_01"
            ],
            "max_power_profile_file_path" => "./profiles/tests/demand_electricity.prf",
            "constant_temperature" => 5,
            "scale" => 1000
        ),
        "TST_SRC_02" => Dict{String,Any}(
            "type" => "BoundedSupply",
            "medium" => "m_h_w_lt1",
            "control_refs" => [],
            "output_refs" => [
                "TST_HP_03"
            ],
            "max_power_profile_file_path" => "./profiles/tests/demand_electricity.prf",
            "constant_temperature" => 5,
            "scale" => 1000
        ),
        "TST_GRI_EL" => Dict{String,Any}(
            "type" => "GridConnection",
            "medium" => "m_e_ac_230v",
            "control_refs" => [],
            "output_refs" => [
                "TST_BUS_EL"
            ],
            "is_source" => true
        ),
        "TST_HP_01" => Dict{String,Any}(
            "type" => "HeatPump",
            "m_heat_in" => "m_h_w_lt1",
            "m_heat_out" => "m_h_w_ht1",
            "control_refs" => [],
            "output_refs" => [
                "TST_HP_02"
            ],
            "strategy" => Dict{String,Any}(
                "name" => "demand_driven"
            ),
            "power_th" => 1200,
            "output_temperature" => 10,
            "min_power_fraction" => 0.0
        ),
        "TST_HP_02" => Dict{String,Any}(
            "type" => "HeatPump",
            "m_heat_in" => "m_h_w_ht1",
            "m_heat_out" => "m_h_w_ht1",
            "control_refs" => [],
            "output_refs" => [
                "TST_BUS_TH"
            ],
            "strategy" => Dict{String,Any}(
                "name" => "demand_driven"
            ),
            "power_th" => 1200,
            "output_temperature" => 25,
            "min_power_fraction" => 0.0
        ),
        "TST_HP_03" => Dict{String,Any}(
            "type" => "HeatPump",
            "m_heat_in" => "m_h_w_lt1",
            "m_heat_out" => "m_h_w_ht1",
            "control_refs" => [],
            "output_refs" => [
                "TST_HP_04"
            ],
            "strategy" => Dict{String,Any}(
                "name" => "demand_driven"
            ),
            "power_th" => 1200,
            "output_temperature" => 20,
            "min_power_fraction" => 0.0
        ),
        "TST_HP_04" => Dict{String,Any}(
            "type" => "HeatPump",
            "m_heat_in" => "m_h_w_ht1",
            "m_heat_out" => "m_h_w_ht1",
            "control_refs" => [],
            "output_refs" => [
                "TST_BUS_TH"
            ],
            "strategy" => Dict{String,Any}(
                "name" => "demand_driven"
            ),
            "power_th" => 1200,
            "output_temperature" => 30,
            "min_power_fraction" => 0.0
        ),
        "TST_HP_05" => Dict{String,Any}(
            "type" => "HeatPump",
            "m_heat_in" => "m_h_w_ht1",
            "m_heat_out" => "m_h_w_ht1",
            "control_refs" => [],
            "output_refs" => [
                "TST_HP_06"
            ],
            "strategy" => Dict{String,Any}(
                "name" => "demand_driven"
            ),
            "power_th" => 1200,
            "output_temperature" => 60,
            "min_power_fraction" => 0.0
        ),
        "TST_HP_06" => Dict{String,Any}(
            "type" => "HeatPump",
            "m_heat_in" => "m_h_w_ht1",
            "m_heat_out" => "m_h_w_ht1",
            "control_refs" => [],
            "output_refs" => [
                "TST_DEM_01"
            ],
            "strategy" => Dict{String,Any}(
                "name" => "demand_driven"
            ),
            "power_th" => 1200,
            "min_power_fraction" => 0.0
        ),
        "TST_HP_07" => Dict{String,Any}(
            "type" => "HeatPump",
            "m_heat_in" => "m_h_w_ht1",
            "m_heat_out" => "m_h_w_ht1",
            "control_refs" => [],
            "output_refs" => [
                "TST_HP_08"
            ],
            "strategy" => Dict{String,Any}(
                "name" => "demand_driven"
            ),
            "power_th" => 1200,
            "output_temperature" => 70,
            "min_power_fraction" => 0.0
        ),
        "TST_HP_08" => Dict{String,Any}(
            "type" => "HeatPump",
            "m_heat_in" => "m_h_w_ht1",
            "m_heat_out" => "m_h_w_ht1",
            "control_refs" => [],
            "output_refs" => [
                "TST_DEM_02"
            ],
            "strategy" => Dict{String,Any}(
                "name" => "demand_driven"
            ),
            "power_th" => 1200,
            "min_power_fraction" => 0.0
        ),
        "TST_BUS_EL" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_e_ac_230v",
            "control_refs" => [],
            "connections" => Dict{String,Any}(
                "input_order" => [
                    "TST_GRI_EL"
                ],
                "output_order" => [
                    "TST_HP_01",
                    "TST_HP_02",
                    "TST_HP_03",
                    "TST_HP_04",
                    "TST_HP_05",
                    "TST_HP_06",
                    "TST_HP_07",
                    "TST_HP_08"
                ],
                "energy_flow" => [
                    [1, 1, 1, 1, 1, 1, 1, 1]
                ]
            )
        ),
        "TST_BUS_TH" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "connections" => Dict{String,Any}(
                "input_order" => [
                    "TST_HP_04",
                    "TST_HP_02"
                ],
                "output_order" => [
                    "TST_HP_05",
                    "TST_HP_07"
                ],
                "energy_flow" => [
                    [1, 1],
                    [1, 1]
                ]
            )
        )
    )

    expected = [
        "TST_DEM_01 s_reset",
        "TST_DEM_02 s_reset",
        "TST_BUS_TH s_reset",
        "TST_BUS_EL s_reset",
        "TST_HP_01 s_reset",
        "TST_HP_06 s_reset",
        "TST_HP_03 s_reset",
        "TST_HP_05 s_reset",
        "TST_HP_04 s_reset",
        "TST_HP_02 s_reset",
        "TST_HP_08 s_reset",
        "TST_HP_07 s_reset",
        "TST_SRC_01 s_reset",
        "TST_SRC_02 s_reset",
        "TST_GRI_EL s_reset",
        "TST_DEM_01 s_control",
        "TST_DEM_02 s_control",
        "TST_BUS_TH s_control",
        "TST_BUS_EL s_control",
        "TST_HP_01 s_control",
        "TST_HP_06 s_control",
        "TST_HP_03 s_control",
        "TST_HP_05 s_control",
        "TST_HP_04 s_control",
        "TST_HP_02 s_control",
        "TST_HP_08 s_control",
        "TST_HP_07 s_control",
        "TST_SRC_01 s_control",
        "TST_SRC_02 s_control",
        "TST_GRI_EL s_control",
        "TST_DEM_01 s_process",
        "TST_DEM_02 s_process",
        "TST_BUS_TH s_process",
        "TST_BUS_EL s_process",
        "TST_HP_03 s_potential",
        "TST_HP_04 s_potential",
        "TST_HP_01 s_potential",
        "TST_HP_02 s_potential",
        "TST_HP_06 s_potential",
        "TST_HP_05 s_potential",
        "TST_HP_08 s_potential",
        "TST_HP_07 s_potential",
        "TST_HP_04 s_process",
        "TST_HP_03 s_process",
        "TST_HP_02 s_process",
        "TST_HP_01 s_process",
        "TST_HP_05 s_process",
        "TST_HP_06 s_process",
        "TST_HP_07 s_process",
        "TST_HP_08 s_process",
        "TST_SRC_01 s_process",
        "TST_SRC_02 s_process",
        "TST_GRI_EL s_process",
        "TST_BUS_TH s_distribute",
        "TST_BUS_EL s_distribute"
    ]

    simulation_parameters = Dict{String,Any}(
        "time_step_seconds" => 900,
        "time" => 0,
        "epsilon" => 1e-9
    )
    
    components = Resie.load_components(components_config, simulation_parameters)
    steps = Resie.calculate_order_of_operations(components)
    ooo = []
    for step in steps
        push!(ooo, "$(step[1]) $(step[2])")
    end
    @test pwc_ooo_astr(expected, ooo) == ""
end

function test_ooo_middle_transformer()
    components_config = Dict{String,Any}(
        "TST_DEM_01" => Dict{String,Any}(
            "type" => "Demand",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "output_refs" => [],
            "energy_profile_file_path" => "./profiles/examples/general/dem_heat_nrg_var_hi-amp.prf",
            "constant_temperature" => 85,
            "scale" => 500
        ),
        "TST_DEM_02" => Dict{String,Any}(
            "type" => "Demand",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "output_refs" => [],
            "energy_profile_file_path" => "./profiles/examples/general/dem_heat_nrg_var_hi-amp.prf",
            "constant_temperature" => 90,
            "scale" => 500
        ),
        "TST_GRI_O2" => Dict{String,Any}(
            "type" => "GridConnection",
            "medium" => "m_c_g_o2",
            "control_refs" => [],
            "input_refs" => [
                "TST_01_ELY_01"
            ],
            "is_source" => false
        ),
        "TST_SRC_01" => Dict{String,Any}(
            "type" => "BoundedSupply",
            "medium" => "m_h_w_lt1",
            "control_refs" => [],
            "output_refs" => [
                "TST_HP_01"
            ],
            "constant_power" => 500000,
            "constant_temperature" => 20
        ),
        "TST_GRI_01" => Dict{String,Any}(
            "type" => "GridConnection",
            "medium" => "m_e_ac_230v",
            "control_refs" => [],
            "output_refs" => [
                "TST_BUS_EL"
            ],
            "is_source" => true
        ),
        "TST_HP_01" => Dict{String,Any}(
            "type" => "HeatPump",
            "control_refs" => [],
            "output_refs" => [
                "TST_BUS_TH"
            ],
            "strategy" => Dict{String,Any}(
                "name" => "demand_driven"
            ),
            "power_th" => 1200,
            "min_power_fraction" => 0.0
        ),
        "TST_HP_02" => Dict{String,Any}(
            "type" => "HeatPump",
            "control_refs" => [],
            "output_refs" => [
                "TST_BUS_TH"
            ],
            "strategy" => Dict{String,Any}(
                "name" => "demand_driven"
            ),
            "power_th" => 9000,
            "min_power_fraction" => 0.0
        ),
        "TST_01_ELY_01" => Dict{String,Any}(
            "type" => "Electrolyser",
            "control_refs" => [],
            "output_refs" => [
                "TST_HP_02",
                "TST_BUS_H2",
                "TST_GRI_O2"
            ],
            "strategy" => Dict{String,Any}(
                "name" => "demand_driven"
            ),
            "power_el" => 40000,
            "min_power_fraction" => 0.0,
            "output_temperature" => 45
        ),
        "ESS_FBH2_01" => Dict{String,Any}(
            "type" => "FuelBoiler",
            "m_fuel_in" => "m_c_g_h2",
            "control_refs" => [],
            "output_refs" => [
                "TST_DEM_02"
            ],
            "strategy" => Dict{String,Any}(
                "name" => "demand_driven"
            ),
            "min_power_fraction" => 0,
            "power_th" => 1500
        ),
        "TST_BUS_EL" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_e_ac_230v",
            "control_refs" => [],
            "connections" => Dict{String,Any}(
                "input_order" => [
                    "TST_GRI_01"
                ],
                "output_order" => [
                    "TST_HP_01",
                    "TST_HP_02",
                    "TST_01_ELY_01"
                ],
                "energy_flow" => [
                    [1, 1, 1]
                ]
            )
        ),
        "TST_BUS_TH" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "connections" => Dict{String,Any}(
                "input_order" => [
                    "TST_HP_02",
                    "TST_HP_01"
                ],
                "output_order" => [
                    "TST_DEM_01"
                ],
                "energy_flow" => [
                    [1, 1, 1],
                    [1, 1, 1]
                ]
            )
        ),
        "TST_BUS_H2" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_c_g_h2",
            "control_refs" => [],
            "connections" => Dict{String,Any}(
                "input_order" => [
                    "TST_01_ELY_01"
                ],
                "output_order" => [
                    "ESS_FBH2_01"
                ],
                "energy_flow" => [
                    [1]
                ]
            )
        )
    )
    expected = [
        "TST_DEM_01 s_reset",
        "TST_DEM_02 s_reset",
        "TST_BUS_TH s_reset",
        "TST_BUS_H2 s_reset",
        "TST_BUS_EL s_reset",
        "TST_HP_01 s_reset",
        "TST_01_ELY_01 s_reset",
        "TST_HP_02 s_reset",
        "ESS_FBH2_01 s_reset",
        "TST_SRC_01 s_reset",
        "TST_GRI_01 s_reset",
        "TST_GRI_O2 s_reset",
        "TST_DEM_01 s_control",
        "TST_DEM_02 s_control",
        "TST_BUS_TH s_control",
        "TST_BUS_H2 s_control",
        "TST_BUS_EL s_control",
        "TST_HP_01 s_control",
        "TST_01_ELY_01 s_control",
        "TST_HP_02 s_control",
        "ESS_FBH2_01 s_control",
        "TST_SRC_01 s_control",
        "TST_GRI_01 s_control",
        "TST_GRI_O2 s_control",
        "TST_DEM_01 s_process",
        "TST_DEM_02 s_process",
        "TST_BUS_TH s_process",
        "TST_BUS_H2 s_process",
        "TST_BUS_EL s_process",
        "ESS_FBH2_01 s_potential",
        "TST_HP_01 s_potential",
        "TST_HP_02 s_potential",
        "TST_01_ELY_01 s_process",
        "ESS_FBH2_01 s_process",
        "TST_HP_02 s_process",
        "TST_HP_01 s_process",
        "TST_SRC_01 s_process",
        "TST_GRI_01 s_process",
        "TST_GRI_O2 s_process",
        "TST_BUS_TH s_distribute",
        "TST_BUS_H2 s_distribute",
        "TST_BUS_EL s_distribute"
    ]

    simulation_parameters = Dict{String,Any}(
        "time_step_seconds" => 900,
        "time" => 0,
        "epsilon" => 1e-9
    )
    
    components = Resie.load_components(components_config, simulation_parameters)
    steps = Resie.calculate_order_of_operations(components)
    ooo = []
    for step in steps
        push!(ooo, "$(step[1]) $(step[2])")
    end
    @test pwc_ooo_astr(expected, ooo) == ""
end

function test_ooo_parallels()
    components_config = Dict{String,Any}(
        "TST_DEM_01" => Dict{String,Any}(
            "type" => "Demand",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "output_refs" => [],
            "energy_profile_file_path" => "./profiles/examples/general/dem_heat_nrg_var_hi-amp.prf",
            "constant_temperature" => 85,
            "scale" => 1000
        ),
        "TST_SRC_01" => Dict{String,Any}(
            "type" => "BoundedSupply",
            "medium" => "m_h_w_lt1",
            "control_refs" => [],
            "output_refs" => [
                "TST_BUS_TH1"
            ],
            "constant_power" => 500000,
            "constant_temperature" => 5
        ),
        "TST_GRI_01" => Dict{String,Any}(
            "type" => "GridConnection",
            "medium" => "m_e_ac_230v",
            "control_refs" => [],
            "output_refs" => [
                "TST_BUS_EL"
            ],
            "is_source" => true
        ),
        "TST_HP_01" => Dict{String,Any}(
            "type" => "HeatPump",
            "m_heat_in" => "m_h_w_lt1",
            "m_heat_out" => "m_h_w_ht1",
            "control_refs" => [],
            "output_refs" => [
                "TST_HP_02"
            ],
            "strategy" => Dict{String,Any}(
                "name" => "demand_driven"
            ),
            "power_th" => 500,
            "output_temperature" => 30,
            "min_power_fraction" => 0.0
        ),
        "TST_HP_02" => Dict{String,Any}(
            "type" => "HeatPump",
            "m_heat_in" => "m_h_w_ht1",
            "m_heat_out" => "m_h_w_ht1",
            "control_refs" => [],
            "output_refs" => [
                "TST_BUS_TH2"
            ],
            "strategy" => Dict{String,Any}(
                "name" => "demand_driven"
            ),
            "power_th" => 500,
            "min_power_fraction" => 0.0
        ),
        "TST_HP_03" => Dict{String,Any}(
            "type" => "HeatPump",
            "m_heat_in" => "m_h_w_lt1",
            "m_heat_out" => "m_h_w_ht1",
            "control_refs" => [],
            "output_refs" => [
                "TST_HP_04"
            ],
            "strategy" => Dict{String,Any}(
                "name" => "demand_driven"
            ),
            "power_th" => 500,
            "output_temperature" => 60,
            "min_power_fraction" => 0.0
        ),
        "TST_HP_04" => Dict{String,Any}(
            "type" => "HeatPump",
            "m_heat_in" => "m_h_w_ht1",
            "m_heat_out" => "m_h_w_ht1",
            "control_refs" => [],
            "output_refs" => [
                "TST_BUS_TH2"
            ],
            "strategy" => Dict{String,Any}(
                "name" => "demand_driven"
            ),
            "power_th" => 600,
            "min_power_fraction" => 0.0
        ),
        "TST_HP_05" => Dict{String,Any}(
            "type" => "HeatPump",
            "m_heat_in" => "m_h_w_lt1",
            "m_heat_out" => "m_h_w_ht1",
            "control_refs" => [],
            "output_refs" => [
                "TST_BUS_TH2"
            ],
            "strategy" => Dict{String,Any}(
                "name" => "demand_driven"
            ),
            "power_th" => 1500,
            "min_power_fraction" => 0.0
        ),
        "TST_BUS_EL" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_e_ac_230v",
            "control_refs" => [],
            "connections" => Dict{String,Any}(
                "input_order" => [
                    "TST_GRI_01"
                ],
                "output_order" => [
                    "TST_HP_01",
                    "TST_HP_02",
                    "TST_HP_03",
                    "TST_HP_04",
                    "TST_HP_05"
                ],
                "energy_flow" => [
                    [1, 1, 1, 1, 1]
                ]
            )
        ),
        "TST_BUS_TH1" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_h_w_lt1",
            "control_refs" => [],
            "connections" => Dict{String,Any}(
                "input_order" => [
                    "TST_SRC_01"
                ],
                "output_order" => [
                    "TST_HP_01",
                    "TST_HP_03",
                    "TST_HP_05"
                ],
                "energy_flow" => [
                    [1, 1, 1]
                ]
            )
        ),
        "TST_BUS_TH2" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "connections" => Dict{String,Any}(
                "input_order" => [
                    "TST_HP_02",
                    "TST_HP_04",
                    "TST_HP_05"
                ],
                "output_order" => [
                    "TST_DEM_01"
                ],
                "energy_flow" => [
                    [1],
                    [1],
                    [1]
                ]
            )
        )
    )
    expected = [
        "TST_DEM_01 s_reset",
        "TST_BUS_TH2 s_reset",
        "TST_BUS_EL s_reset",
        "TST_BUS_TH1 s_reset",
        "TST_HP_01 s_reset",
        "TST_HP_03 s_reset",
        "TST_HP_05 s_reset",
        "TST_HP_04 s_reset",
        "TST_HP_02 s_reset",
        "TST_SRC_01 s_reset",
        "TST_GRI_01 s_reset",
        "TST_DEM_01 s_control",
        "TST_BUS_TH2 s_control",
        "TST_BUS_EL s_control",
        "TST_BUS_TH1 s_control",
        "TST_HP_01 s_control",
        "TST_HP_03 s_control",
        "TST_HP_05 s_control",
        "TST_HP_04 s_control",
        "TST_HP_02 s_control",
        "TST_SRC_01 s_control",
        "TST_GRI_01 s_control",
        "TST_DEM_01 s_process",
        "TST_BUS_TH2 s_process",
        "TST_BUS_EL s_process",
        "TST_BUS_TH1 s_process",
        "TST_HP_02 s_potential",
        "TST_HP_01 s_potential",
        "TST_HP_02 s_potential",
        "TST_HP_04 s_potential",
        "TST_HP_03 s_potential",
        "TST_HP_04 s_potential",
        "TST_HP_05 s_potential",
        "TST_HP_01 s_process",
        "TST_HP_02 s_process",
        "TST_HP_03 s_process",
        "TST_HP_04 s_process",
        "TST_HP_05 s_process",
        "TST_SRC_01 s_process",
        "TST_GRI_01 s_process",
        "TST_BUS_TH2 s_distribute",
        "TST_BUS_EL s_distribute",
        "TST_BUS_TH1 s_distribute"
    ]

    simulation_parameters = Dict{String,Any}(
        "time_step_seconds" => 900,
        "time" => 0,
        "epsilon" => 1e-9
    )
    
    components = Resie.load_components(components_config, simulation_parameters)
    steps = Resie.calculate_order_of_operations(components)
    ooo = []
    for step in steps
        push!(ooo, "$(step[1]) $(step[2])")
    end
    @test pwc_ooo_astr(expected, ooo) == ""
end

function test_ooo_parallels_different_order()
    components_config = Dict{String,Any}(
        "TST_DEM_01" => Dict{String,Any}(
            "type" => "Demand",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "output_refs" => [],
            "energy_profile_file_path" => "./profiles/examples/general/dem_heat_nrg_var_hi-amp.prf",
            "constant_temperature" => 85,
            "scale" => 1000
        ),
        "TST_SRC_01" => Dict{String,Any}(
            "type" => "BoundedSupply",
            "medium" => "m_h_w_lt1",
            "control_refs" => [],
            "output_refs" => [
                "TST_BUS_TH1"
            ],
            "constant_power" => 500000,
            "constant_temperature" => 5
        ),
        "TST_GRI_01" => Dict{String,Any}(
            "type" => "GridConnection",
            "medium" => "m_e_ac_230v",
            "control_refs" => [],
            "output_refs" => [
                "TST_BUS_EL"
            ],
            "is_source" => true
        ),
        "TST_HP_01" => Dict{String,Any}(
            "type" => "HeatPump",
            "m_heat_in" => "m_h_w_lt1",
            "m_heat_out" => "m_h_w_ht1",
            "control_refs" => [],
            "output_refs" => [
                "TST_HP_02"
            ],
            "strategy" => Dict{String,Any}(
                "name" => "demand_driven"
            ),
            "power_th" => 500,
            "output_temperature" => 30,
            "min_power_fraction" => 0.0
        ),
        "TST_HP_02" => Dict{String,Any}(
            "type" => "HeatPump",
            "m_heat_in" => "m_h_w_ht1",
            "m_heat_out" => "m_h_w_ht1",
            "control_refs" => [],
            "output_refs" => [
                "TST_BUS_TH2"
            ],
            "strategy" => Dict{String,Any}(
                "name" => "demand_driven"
            ),
            "power_th" => 500,
            "min_power_fraction" => 0.0
        ),
        "TST_HP_03" => Dict{String,Any}(
            "type" => "HeatPump",
            "m_heat_in" => "m_h_w_lt1",
            "m_heat_out" => "m_h_w_ht1",
            "control_refs" => [],
            "output_refs" => [
                "TST_HP_04"
            ],
            "strategy" => Dict{String,Any}(
                "name" => "demand_driven"
            ),
            "power_th" => 500,
            "output_temperature" => 60,
            "min_power_fraction" => 0.0
        ),
        "TST_HP_04" => Dict{String,Any}(
            "type" => "HeatPump",
            "m_heat_in" => "m_h_w_ht1",
            "m_heat_out" => "m_h_w_ht1",
            "control_refs" => [],
            "output_refs" => [
                "TST_BUS_TH2"
            ],
            "strategy" => Dict{String,Any}(
                "name" => "demand_driven"
            ),
            "power_th" => 600,
            "min_power_fraction" => 0.0
        ),
        "TST_HP_05" => Dict{String,Any}(
            "type" => "HeatPump",
            "m_heat_in" => "m_h_w_lt1",
            "m_heat_out" => "m_h_w_ht1",
            "control_refs" => [],
            "output_refs" => [
                "TST_BUS_TH2"
            ],
            "strategy" => Dict{String,Any}(
                "name" => "demand_driven"
            ),
            "power_th" => 1500,
            "min_power_fraction" => 0.0
        ),
        "TST_BUS_EL" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_e_ac_230v",
            "control_refs" => [],
            "connections" => Dict{String,Any}(
                "input_order" => [
                    "TST_GRI_01"
                ],
                "output_order" => [
                    "TST_HP_01",
                    "TST_HP_02",
                    "TST_HP_03",
                    "TST_HP_04",
                    "TST_HP_05"
                ],
                "energy_flow" => [
                    [1, 1, 1, 1, 1]
                ]
            )
        ),
        "TST_BUS_TH1" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_h_w_lt1",
            "control_refs" => [],
            "connections" => Dict{String,Any}(
                "input_order" => [
                    "TST_SRC_01"
                ],
                "output_order" => [
                    "TST_HP_03",
                    "TST_HP_05",
                    "TST_HP_01"
                ],
                "energy_flow" => [
                    [1, 1, 1]
                ]
            )
        ),
        "TST_BUS_TH2" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "connections" => Dict{String,Any}(
                "input_order" => [
                    "TST_HP_04",
                    "TST_HP_05",
                    "TST_HP_02"
                ],
                "output_order" => [
                    "TST_DEM_01"
                ],
                "energy_flow" => [
                    [1],
                    [1],
                    [1]
                ]
            )
        )
    )
    expected = [
        "TST_DEM_01 s_reset",
        "TST_BUS_TH2 s_reset",
        "TST_BUS_EL s_reset",
        "TST_BUS_TH1 s_reset",
        "TST_HP_01 s_reset",
        "TST_HP_03 s_reset",
        "TST_HP_05 s_reset",
        "TST_HP_04 s_reset",
        "TST_HP_02 s_reset",
        "TST_SRC_01 s_reset",
        "TST_GRI_01 s_reset",
        "TST_DEM_01 s_control",
        "TST_BUS_TH2 s_control",
        "TST_BUS_EL s_control",
        "TST_BUS_TH1 s_control",
        "TST_HP_01 s_control",
        "TST_HP_03 s_control",
        "TST_HP_05 s_control",
        "TST_HP_04 s_control",
        "TST_HP_02 s_control",
        "TST_SRC_01 s_control",
        "TST_GRI_01 s_control",
        "TST_DEM_01 s_process",
        "TST_BUS_TH2 s_process",
        "TST_BUS_EL s_process",
        "TST_BUS_TH1 s_process",
        "TST_HP_04 s_potential",
        "TST_HP_03 s_potential",
        "TST_HP_04 s_potential",
        "TST_HP_05 s_potential",
        "TST_HP_02 s_potential",
        "TST_HP_01 s_potential",
        "TST_HP_03 s_process",
        "TST_HP_04 s_process",
        "TST_HP_05 s_process",
        "TST_HP_01 s_process",
        "TST_HP_02 s_process",
        "TST_SRC_01 s_process",
        "TST_GRI_01 s_process",
        "TST_BUS_TH2 s_distribute",
        "TST_BUS_EL s_distribute",
        "TST_BUS_TH1 s_distribute"
    ]

    simulation_parameters = Dict{String,Any}(
        "time_step_seconds" => 900,
        "time" => 0,
        "epsilon" => 1e-9
    )
    
    components = Resie.load_components(components_config, simulation_parameters)
    steps = Resie.calculate_order_of_operations(components)
    ooo = []
    for step in steps
        push!(ooo, "$(step[1]) $(step[2])")
    end
    @test pwc_ooo_astr(expected, ooo) == ""
end

function test_ooo_parallels_in_chain()
    components_config = Dict{String,Any}(
        "TST_DEM_01" => Dict{String,Any}(
            "type" => "Demand",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "output_refs" => [],
            "energy_profile_file_path" => "./profiles/tests/demand_heating_energy.prf",
            "constant_temperature" => 70,
            "scale" => 500
        ),
        "TST_DEM_02" => Dict{String,Any}(
            "type" => "Demand",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "output_refs" => [],
            "energy_profile_file_path" => "./profiles/tests/demand_heating_energy.prf",
            "constant_temperature" => 90,
            "scale" => 500
        ),
        "TST_DEM_03" => Dict{String,Any}(
            "type" => "Demand",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "output_refs" => [],
            "energy_profile_file_path" => "./profiles/tests/demand_heating_energy.prf",
            "constant_temperature" => 95,
            "scale" => 500
        ),
        "TST_GRI_H2" => Dict{String,Any}(
            "type" => "GridConnection",
            "medium" => "m_c_g_h2",
            "control_refs" => [],
            "input_refs" => [
                "TST_01_ELY_01"
            ],
            "is_source" => false
        ),
        "TST_GRI_O2" => Dict{String,Any}(
            "type" => "GridConnection",
            "medium" => "m_c_g_o2",
            "control_refs" => [],
            "input_refs" => [
                "TST_01_ELY_01"
            ],
            "is_source" => false
        ),
        "TST_SRC_01" => Dict{String,Any}(
            "type" => "BoundedSupply",
            "medium" => "m_h_w_lt1",
            "control_refs" => [],
            "output_refs" => [
                "TST_BUS_00"
            ],
            "max_power_profile_file_path" => "./profiles/tests/demand_electricity.prf",
            "temperature_profile_file_path" => "./profiles/tests/demand_heating_temperature.prf",
            "scale" => 1000
        ),
        "TST_GRI_00" => Dict{String,Any}(
            "type" => "GridConnection",
            "medium" => "m_e_ac_230v",
            "control_refs" => [],
            "output_refs" => [
                "TST_01_ELY_01"
            ],
            "is_source" => true
        ),
        "TST_GRI_01" => Dict{String,Any}(
            "type" => "GridConnection",
            "medium" => "m_e_ac_230v",
            "control_refs" => [],
            "output_refs" => [
                "TST_HP_01"
            ],
            "is_source" => true
        ),
        "TST_SRC_1b" => Dict{String,Any}(
            "type" => "BoundedSupply",
            "medium" => "m_e_ac_230v",
            "control_refs" => [],
            "output_refs" => [
                "TST_HP_01b"
            ],
            "is_source" => true,
            "constant_power" => 400
        ),
        "TST_GRI_02" => Dict{String,Any}(
            "type" => "GridConnection",
            "medium" => "m_e_ac_230v",
            "control_refs" => [],
            "output_refs" => [
                "TST_HP_02"
            ],
            "is_source" => true
        ),
        "TST_GRI_03" => Dict{String,Any}(
            "type" => "GridConnection",
            "medium" => "m_e_ac_230v",
            "control_refs" => [],
            "output_refs" => [
                "TST_HP_03"
            ],
            "is_source" => true
        ),
        "TST_GRI_04" => Dict{String,Any}(
            "type" => "GridConnection",
            "medium" => "m_e_ac_230v",
            "control_refs" => [],
            "output_refs" => [
                "TST_HP_04"
            ],
            "is_source" => true
        ),
        "TST_HP_01" => Dict{String,Any}(
            "type" => "HeatPump",
            "control_refs" => [],
            "output_refs" => [
                "TST_HP_01b"
            ],
            "strategy" => Dict{String,Any}(
                "name" => "demand_driven"
            ),
            "power_th" => 1200,
            "output_temperature" => 60,
            "min_power_fraction" => 0.0
        ),
        "TST_HP_01b" => Dict{String,Any}(
            "type" => "HeatPump",
            "control_refs" => [],
            "m_heat_in" => "m_h_w_ht1",
            "output_refs" => [
                "TST_BUS_01"
            ],
            "strategy" => Dict{String,Any}(
                "name" => "demand_driven"
            ),
            "power_th" => 9000,
            "min_power_fraction" => 0.0,
            "input_temperature" => 60
        ),
        "TST_HP_02" => Dict{String,Any}(
            "type" => "HeatPump",
            "control_refs" => [],
            "output_refs" => [
                "TST_BUS_01"
            ],
            "strategy" => Dict{String,Any}(
                "name" => "demand_driven"
            ),
            "power_th" => 9000,
            "min_power_fraction" => 0.0
        ),
        "TST_HP_03" => Dict{String,Any}(
            "type" => "HeatPump",
            "control_refs" => [],
            "m_heat_in" => "m_h_w_ht1",
            "output_refs" => [
                "TST_DEM_02"
            ],
            "strategy" => Dict{String,Any}(
                "name" => "demand_driven"
            ),
            "power_th" => 9000,
            "min_power_fraction" => 0.0,
            "input_temperature" => 70
        ),
        "TST_HP_04" => Dict{String,Any}(
            "type" => "HeatPump",
            "control_refs" => [],
            "m_heat_in" => "m_h_w_ht1",
            "output_refs" => [
                "TST_DEM_03"
            ],
            "strategy" => Dict{String,Any}(
                "name" => "demand_driven"
            ),
            "power_th" => 9000,
            "min_power_fraction" => 0.0,
            "input_temperature" => 70
        ),

        "TST_01_ELY_01" => Dict{String,Any}(
            "type" => "Electrolyser",
            "control_refs" => [],
            "output_refs" => [
                "TST_BUS_00",
                "TST_GRI_H2",
                "TST_GRI_O2"
            ],
            "strategy" => Dict{String,Any}(
                "name" => "demand_driven"
            ),
            "power_el" => 40000,
            "min_power_fraction" => 0.0,
            "output_temperature" => 45
        ),
        "TST_BUS_00" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_h_w_lt1",
            "control_refs" => [],
            "connections" => Dict{String,Any}(
                "input_order" => [
                    "TST_01_ELY_01",
                    "TST_SRC_01"
                ],
                "output_order" => [
                    "TST_HP_01",
                    "TST_HP_02"
                ],
                "energy_flow" => [
                    [1, 1],
                    [1, 1]
                ]
            )
        ),
        "TST_BUS_01" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "connections" => Dict{String,Any}(
                "input_order" => [
                    "TST_HP_01b",
                    "TST_HP_02"
                ],
                "output_order" => [
                    "TST_HP_03",
                    "TST_HP_04",
                    "TST_DEM_01"
                ],
                "energy_flow" => [
                    [1, 1, 1],
                    [1, 1, 1]
                ]
            )
        )
    )

    expected = [
        "TST_DEM_01 s_reset",
        "TST_DEM_03 s_reset",
        "TST_DEM_02 s_reset",
        "TST_BUS_00 s_reset",
        "TST_BUS_01 s_reset",
        "TST_HP_01 s_reset",
        "TST_HP_03 s_reset",
        "TST_01_ELY_01 s_reset",
        "TST_HP_04 s_reset",
        "TST_HP_02 s_reset",
        "TST_HP_01b s_reset",
        "TST_GRI_00 s_reset",
        "TST_GRI_03 s_reset",
        "TST_GRI_02 s_reset",
        "TST_SRC_01 s_reset",
        "TST_SRC_1b s_reset",
        "TST_GRI_01 s_reset",
        "TST_GRI_04 s_reset",
        "TST_GRI_O2 s_reset",
        "TST_GRI_H2 s_reset",
        "TST_DEM_01 s_control",
        "TST_DEM_03 s_control",
        "TST_DEM_02 s_control",
        "TST_BUS_00 s_control",
        "TST_BUS_01 s_control",
        "TST_HP_01 s_control",
        "TST_HP_03 s_control",
        "TST_01_ELY_01 s_control",
        "TST_HP_04 s_control",
        "TST_HP_02 s_control",
        "TST_HP_01b s_control",
        "TST_GRI_00 s_control",
        "TST_GRI_03 s_control",
        "TST_GRI_02 s_control",
        "TST_SRC_01 s_control",
        "TST_SRC_1b s_control",
        "TST_GRI_01 s_control",
        "TST_GRI_04 s_control",
        "TST_GRI_O2 s_control",
        "TST_GRI_H2 s_control",
        "TST_DEM_01 s_process",
        "TST_DEM_03 s_process",
        "TST_DEM_02 s_process",
        "TST_BUS_00 s_process",
        "TST_BUS_01 s_process",
        "TST_01_ELY_01 s_potential",
        "TST_HP_03 s_potential",
        "TST_HP_04 s_potential",
        "TST_HP_01b s_potential",
        "TST_HP_01 s_potential",
        "TST_HP_01b s_potential",
        "TST_HP_02 s_potential",
        "TST_01_ELY_01 s_process",
        "TST_HP_01 s_process",
        "TST_HP_01b s_process",
        "TST_HP_02 s_process",
        "TST_HP_03 s_process",
        "TST_HP_04 s_process",
        "TST_GRI_00 s_process",
        "TST_GRI_03 s_process",
        "TST_GRI_02 s_process",
        "TST_SRC_01 s_process",
        "TST_SRC_1b s_process",
        "TST_GRI_01 s_process",
        "TST_GRI_04 s_process",
        "TST_GRI_O2 s_process",
        "TST_GRI_H2 s_process",
        "TST_BUS_00 s_distribute",
        "TST_BUS_01 s_distribute"
    ]

    simulation_parameters = Dict{String,Any}(
        "time_step_seconds" => 900,
        "time" => 0,
        "epsilon" => 1e-9
    )
    
    components = Resie.load_components(components_config, simulation_parameters)
    steps = Resie.calculate_order_of_operations(components)
    ooo = []
    for step in steps
        push!(ooo, "$(step[1]) $(step[2])")
    end
    @test pwc_ooo_astr(expected, ooo) == ""
end

function test_ooo_parallels_in_a_row()
    components_config = Dict{String,Any}(
        "TST_DEM_01" => Dict{String,Any}(
            "type" => "Demand",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "output_refs" => [],
            "energy_profile_file_path" => "./profiles/examples/general/dem_heat_nrg_var_hi-amp.prf",
            "constant_temperature" => 85,
            "scale" => 1000
        ),
        "TST_SRC_01" => Dict{String,Any}(
            "type" => "BoundedSupply",
            "medium" => "m_h_w_lt1",
            "control_refs" => [],
            "output_refs" => [
                "TST_BUS_TH1"
            ],
            "constant_power" => 500000,
            "constant_temperature" => 5
        ),
        "TST_GRI_01" => Dict{String,Any}(
            "type" => "GridConnection",
            "medium" => "m_e_ac_230v",
            "control_refs" => [],
            "output_refs" => [
                "TST_BUS_EL"
            ],
            "is_source" => true
        ),
        "TST_HP_01" => Dict{String,Any}(
            "type" => "HeatPump",
            "m_heat_in" => "m_h_w_lt1",
            "m_heat_out" => "m_h_w_ht1",
            "control_refs" => [],
            "output_refs" => [
                "TST_HP_02"
            ],
            "strategy" => Dict{String,Any}(
                "name" => "demand_driven"
            ),
            "power_th" => 500,
            "output_temperature" => 30,
            "min_power_fraction" => 0.0
        ),
        "TST_HP_02" => Dict{String,Any}(
            "type" => "HeatPump",
            "m_heat_in" => "m_h_w_ht1",
            "m_heat_out" => "m_h_w_ht1",
            "control_refs" => [],
            "output_refs" => [
                "TST_BUS_TH2"
            ],
            "strategy" => Dict{String,Any}(
                "name" => "demand_driven"
            ),
            "power_th" => 500,
            "output_temperature" => 40,
            "min_power_fraction" => 0.0
        ),
        "TST_HP_03" => Dict{String,Any}(
            "type" => "HeatPump",
            "m_heat_in" => "m_h_w_lt1",
            "m_heat_out" => "m_h_w_ht1",
            "control_refs" => [],
            "output_refs" => [
                "TST_HP_04"
            ],
            "strategy" => Dict{String,Any}(
                "name" => "demand_driven"
            ),
            "power_th" => 500,
            "output_temperature" => 35,
            "min_power_fraction" => 0.0
        ),
        "TST_HP_04" => Dict{String,Any}(
            "type" => "HeatPump",
            "m_heat_in" => "m_h_w_ht1",
            "m_heat_out" => "m_h_w_ht1",
            "control_refs" => [],
            "output_refs" => [
                "TST_BUS_TH2"
            ],
            "strategy" => Dict{String,Any}(
                "name" => "demand_driven"
            ),
            "power_th" => 600,
            "output_temperature" => 45,
            "min_power_fraction" => 0.0
        ),
        "TST_HP_05" => Dict{String,Any}(
            "type" => "HeatPump",
            "m_heat_in" => "m_h_w_ht1",
            "m_heat_out" => "m_h_w_ht1",
            "control_refs" => [],
            "output_refs" => [
                "TST_HP_06"
            ],
            "strategy" => Dict{String,Any}(
                "name" => "demand_driven"
            ),
            "power_th" => 600,
            "output_temperature" => 65,
            "min_power_fraction" => 0.0
        ),
        "TST_HP_06" => Dict{String,Any}(
            "type" => "HeatPump",
            "m_heat_in" => "m_h_w_ht1",
            "m_heat_out" => "m_h_w_ht1",
            "control_refs" => [],
            "output_refs" => [
                "TST_BUS_TH3"
            ],
            "strategy" => Dict{String,Any}(
                "name" => "demand_driven"
            ),
            "power_th" => 600,
            "min_power_fraction" => 0.0
        ),
        "TST_HP_07" => Dict{String,Any}(
            "type" => "HeatPump",
            "m_heat_in" => "m_h_w_ht1",
            "m_heat_out" => "m_h_w_ht1",
            "control_refs" => [],
            "output_refs" => [
                "TST_HP_08"
            ],
            "strategy" => Dict{String,Any}(
                "name" => "demand_driven"
            ),
            "power_th" => 600,
            "output_temperature" => 75,
            "min_power_fraction" => 0.0
        ),
        "TST_HP_08" => Dict{String,Any}(
            "type" => "HeatPump",
            "m_heat_in" => "m_h_w_ht1",
            "m_heat_out" => "m_h_w_ht1",
            "control_refs" => [],
            "output_refs" => [
                "TST_BUS_TH3"
            ],
            "strategy" => Dict{String,Any}(
                "name" => "demand_driven"
            ),
            "power_th" => 600,
            "min_power_fraction" => 0.0
        ),
        "TST_BUS_EL" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_e_ac_230v",
            "control_refs" => [],
            "connections" => Dict{String,Any}(
                "input_order" => [
                    "TST_GRI_01"
                ],
                "output_order" => [
                    "TST_HP_01",
                    "TST_HP_02",
                    "TST_HP_03",
                    "TST_HP_04",
                    "TST_HP_05",
                    "TST_HP_06",
                    "TST_HP_07",
                    "TST_HP_08",
                ],
                "energy_flow" => [
                    [1, 1, 1, 1, 1, 1, 1, 1]
                ]
            )
        ),
        "TST_BUS_TH1" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_h_w_lt1",
            "control_refs" => [],
            "connections" => Dict{String,Any}(
                "input_order" => [
                    "TST_SRC_01"
                ],
                "output_order" => [
                    "TST_HP_01",
                    "TST_HP_03"
                ],
                "energy_flow" => [
                    [1, 1]
                ]
            )
        ),
        "TST_BUS_TH2" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "connections" => Dict{String,Any}(
                "input_order" => [
                    "TST_HP_02",
                    "TST_HP_04"
                ],
                "output_order" => [
                    "TST_HP_05",
                    "TST_HP_07"
                ],
                "energy_flow" => [
                    [1, 1],
                    [1, 1]
                ]
            )
        ),
        "TST_BUS_TH3" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "connections" => Dict{String,Any}(
                "input_order" => [
                    "TST_HP_06",
                    "TST_HP_08"
                ],
                "output_order" => [
                    "TST_DEM_01"
                ],
                "energy_flow" => [
                    [1],
                    [1]
                ]
            )
        )
    )

    expected = [
        "TST_DEM_01 s_reset",
        "TST_BUS_TH2 s_reset",
        "TST_BUS_TH3 s_reset",
        "TST_BUS_EL s_reset",
        "TST_BUS_TH1 s_reset",
        "TST_HP_01 s_reset",
        "TST_HP_06 s_reset",
        "TST_HP_03 s_reset",
        "TST_HP_05 s_reset",
        "TST_HP_04 s_reset",
        "TST_HP_02 s_reset",
        "TST_HP_08 s_reset",
        "TST_HP_07 s_reset",
        "TST_SRC_01 s_reset",
        "TST_GRI_01 s_reset",
        "TST_DEM_01 s_control",
        "TST_BUS_TH2 s_control",
        "TST_BUS_TH3 s_control",
        "TST_BUS_EL s_control",
        "TST_BUS_TH1 s_control",
        "TST_HP_01 s_control",
        "TST_HP_06 s_control",
        "TST_HP_03 s_control",
        "TST_HP_05 s_control",
        "TST_HP_04 s_control",
        "TST_HP_02 s_control",
        "TST_HP_08 s_control",
        "TST_HP_07 s_control",
        "TST_SRC_01 s_control",
        "TST_GRI_01 s_control",
        "TST_DEM_01 s_process",
        "TST_BUS_TH2 s_process",
        "TST_BUS_TH3 s_process",
        "TST_BUS_EL s_process",
        "TST_BUS_TH1 s_process",
        "TST_HP_06 s_potential",
        "TST_HP_05 s_potential",
        "TST_HP_06 s_potential",
        "TST_HP_08 s_potential",
        "TST_HP_07 s_potential",
        "TST_HP_01 s_potential",
        "TST_HP_02 s_potential",
        "TST_HP_01 s_potential",
        "TST_HP_03 s_potential",
        "TST_HP_04 s_potential",
        "TST_HP_01 s_process",
        "TST_HP_02 s_process",
        "TST_HP_03 s_process",
        "TST_HP_04 s_process",
        "TST_HP_05 s_process",
        "TST_HP_06 s_process",
        "TST_HP_07 s_process",
        "TST_HP_08 s_process",
        "TST_SRC_01 s_process",
        "TST_GRI_01 s_process",
        "TST_BUS_TH2 s_distribute",
        "TST_BUS_TH3 s_distribute",
        "TST_BUS_EL s_distribute",
        "TST_BUS_TH1 s_distribute"
    ]

    simulation_parameters = Dict{String,Any}(
        "time_step_seconds" => 900,
        "time" => 0,
        "epsilon" => 1e-9
    )
    
    components = Resie.load_components(components_config, simulation_parameters)
    steps = Resie.calculate_order_of_operations(components)
    ooo = []
    for step in steps
        push!(ooo, "$(step[1]) $(step[2])")
    end
    @test pwc_ooo_astr(expected, ooo) == ""
end


function test_ooo_connected_middle_busses()
    components_config = Dict{String,Any}(
        "TST_DEM_01" => Dict{String,Any}(
            "type" => "Demand",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "output_refs" => [],
            "energy_profile_file_path" => "./profiles/tests/demand_heating_energy.prf",
            "constant_temperature" => 95,
            "scale" => 400
        ),
        "TST_DEM_02" => Dict{String,Any}(
            "type" => "Demand",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "output_refs" => [],
            "energy_profile_file_path" => "./profiles/tests/demand_heating_energy.prf",
            "constant_temperature" => 95,
            "scale" => 400
        ),
        "TST_DEM_03" => Dict{String,Any}(
            "type" => "Demand",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "output_refs" => [],
            "energy_profile_file_path" => "./profiles/tests/demand_heating_energy.prf",
            "constant_temperature" => 115,
            "scale" => 400
        ),
        "TST_SRC_01" => Dict{String,Any}(
            "type" => "BoundedSupply",
            "medium" => "m_h_w_lt1",
            "control_refs" => [],
            "output_refs" => [
                "TST_HP_01"
            ],
            "max_power_profile_file_path" => "./profiles/tests/demand_electricity.prf",
            "constant_temperature" => 5,
            "scale" => 1000
        ),
        "TST_SRC_02" => Dict{String,Any}(
            "type" => "BoundedSupply",
            "medium" => "m_h_w_lt1",
            "control_refs" => [],
            "output_refs" => [
                "TST_HP_07"
            ],
            "max_power_profile_file_path" => "./profiles/tests/demand_electricity.prf",
            "constant_temperature" => 5,
            "scale" => 1000
        ),
        "TST_GRI_EL" => Dict{String,Any}(
            "type" => "GridConnection",
            "medium" => "m_e_ac_230v",
            "control_refs" => [],
            "output_refs" => [
                "TST_BUS_EL"
            ],
            "is_source" => true
        ),
        "TST_HP_01" => Dict{String,Any}(
            "type" => "HeatPump",
            "m_heat_in" => "m_h_w_lt1",
            "m_heat_out" => "m_h_w_ht1",
            "control_refs" => [],
            "output_refs" => [
                "TST_HP_02"
            ],
            "strategy" => Dict{String,Any}(
                "name" => "demand_driven"
            ),
            "power_th" => 1200,
            "output_temperature" => 10,
            "min_power_fraction" => 0.0
        ),
        "TST_HP_02" => Dict{String,Any}(
            "type" => "HeatPump",
            "m_heat_in" => "m_h_w_ht1",
            "m_heat_out" => "m_h_w_ht1",
            "control_refs" => [],
            "output_refs" => [
                "TST_BUS_TH_01"
            ],
            "strategy" => Dict{String,Any}(
                "name" => "demand_driven"
            ),
            "power_th" => 1200,
            "output_temperature" => 25,
            "min_power_fraction" => 0.0
        ),
        "TST_HP_03" => Dict{String,Any}(
            "type" => "HeatPump",
            "m_heat_in" => "m_h_w_ht1",
            "m_heat_out" => "m_h_w_ht1",
            "control_refs" => [],
            "output_refs" => [
                "TST_HP_04"
            ],
            "strategy" => Dict{String,Any}(
                "name" => "demand_driven"
            ),
            "power_th" => 1200,
            "output_temperature" => 60,
            "min_power_fraction" => 0.0
        ),
        "TST_HP_04" => Dict{String,Any}(
            "type" => "HeatPump",
            "m_heat_in" => "m_h_w_ht1",
            "m_heat_out" => "m_h_w_ht1",
            "control_refs" => [],
            "output_refs" => [
                "TST_DEM_01"
            ],
            "strategy" => Dict{String,Any}(
                "name" => "demand_driven"
            ),
            "power_th" => 1200,
            "min_power_fraction" => 0.0
        ),
        "TST_HP_05" => Dict{String,Any}(
            "type" => "HeatPump",
            "m_heat_in" => "m_h_w_ht1",
            "m_heat_out" => "m_h_w_ht1",
            "control_refs" => [],
            "output_refs" => [
                "TST_HP_06"
            ],
            "strategy" => Dict{String,Any}(
                "name" => "demand_driven"
            ),
            "power_th" => 1200,
            "output_temperature" => 45,
            "min_power_fraction" => 0.0
        ),
        "TST_HP_06" => Dict{String,Any}(
            "type" => "HeatPump",
            "m_heat_in" => "m_h_w_ht1",
            "m_heat_out" => "m_h_w_ht1",
            "control_refs" => [],
            "output_refs" => [
                "TST_BUS_TH_02"
            ],
            "strategy" => Dict{String,Any}(
                "name" => "demand_driven"
            ),
            "power_th" => 1200,
            "output_temperature" => 55,
            "min_power_fraction" => 0.0
        ),
        "TST_HP_07" => Dict{String,Any}(
            "type" => "HeatPump",
            "m_heat_in" => "m_h_w_lt1",
            "m_heat_out" => "m_h_w_ht1",
            "control_refs" => [],
            "output_refs" => [
                "TST_HP_08"
            ],
            "strategy" => Dict{String,Any}(
                "name" => "demand_driven"
            ),
            "power_th" => 1200,
            "output_temperature" => 35,
            "min_power_fraction" => 0.0
        ),
        "TST_HP_08" => Dict{String,Any}(
            "type" => "HeatPump",
            "m_heat_in" => "m_h_w_ht1",
            "m_heat_out" => "m_h_w_ht1",
            "control_refs" => [],
            "output_refs" => [
                "TST_BUS_TH_02"
            ],
            "strategy" => Dict{String,Any}(
                "name" => "demand_driven"
            ),
            "power_th" => 1200,
            "output_temperature" => 60,
            "min_power_fraction" => 0.0
        ),
        "TST_HP_09" => Dict{String,Any}(
            "type" => "HeatPump",
            "m_heat_in" => "m_h_w_ht1",
            "m_heat_out" => "m_h_w_ht1",
            "control_refs" => [],
            "output_refs" => [
                "TST_HP_10"
            ],
            "strategy" => Dict{String,Any}(
                "name" => "demand_driven"
            ),
            "power_th" => 1200,
            "output_temperature" => 80,
            "min_power_fraction" => 0.0
        ),
        "TST_HP_10" => Dict{String,Any}(
            "type" => "HeatPump",
            "m_heat_in" => "m_h_w_ht1",
            "m_heat_out" => "m_h_w_ht1",
            "control_refs" => [],
            "output_refs" => [
                "TST_DEM_02"
            ],
            "strategy" => Dict{String,Any}(
                "name" => "demand_driven"
            ),
            "power_th" => 1200,
            "min_power_fraction" => 0.0
        ),
        "TST_HP_11" => Dict{String,Any}(
            "type" => "HeatPump",
            "m_heat_in" => "m_h_w_ht1",
            "m_heat_out" => "m_h_w_ht1",
            "control_refs" => [],
            "output_refs" => [
                "TST_HP_12"
            ],
            "strategy" => Dict{String,Any}(
                "name" => "demand_driven"
            ),
            "power_th" => 1200,
            "output_temperature" => 90,
            "min_power_fraction" => 0.0
        ),
        "TST_HP_12" => Dict{String,Any}(
            "type" => "HeatPump",
            "m_heat_in" => "m_h_w_ht1",
            "m_heat_out" => "m_h_w_ht1",
            "control_refs" => [],
            "output_refs" => [
                "TST_DEM_03"
            ],
            "strategy" => Dict{String,Any}(
                "name" => "demand_driven"
            ),
            "power_th" => 1200,
            "min_power_fraction" => 0.0
        ),
        "TST_BUS_EL" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_e_ac_230v",
            "control_refs" => [],
            "connections" => Dict{String,Any}(
                "input_order" => [
                    "TST_GRI_EL"
                ],
                "output_order" => [
                    "TST_HP_01",
                    "TST_HP_02",
                    "TST_HP_03",
                    "TST_HP_04",
                    "TST_HP_05",
                    "TST_HP_06",
                    "TST_HP_07",
                    "TST_HP_08",
                    "TST_HP_09",
                    "TST_HP_10",
                    "TST_HP_11",
                    "TST_HP_12"
                ],
                "energy_flow" => [
                    [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
                ]
            )
        ),
        "TST_BUS_TH_01" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "connections" => Dict{String,Any}(
                "input_order" => [
                    "TST_HP_02"
                ],
                "output_order" => [
                    "TST_HP_03",
                    "TST_HP_05"
                ],
                "energy_flow" => [
                    [1, 1]
                ]
            )
        ),
        "TST_BUS_TH_02" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "connections" => Dict{String,Any}(
                "input_order" => [
                    "TST_HP_06",
                    "TST_HP_08"
                ],
                "output_order" => [
                    "TST_HP_09",
                    "TST_HP_11"
                ],
                "energy_flow" => [
                    [1, 1],
                    [1, 1]
                ]
            )
        )
    )

    expected = [
        "TST_DEM_01 s_reset",
        "TST_DEM_03 s_reset",
        "TST_DEM_02 s_reset",
        "TST_BUS_TH_01 s_reset",
        "TST_BUS_EL s_reset",
        "TST_BUS_TH_02 s_reset",
        "TST_HP_01 s_reset",
        "TST_HP_09 s_reset",
        "TST_HP_12 s_reset",
        "TST_HP_06 s_reset",
        "TST_HP_03 s_reset",
        "TST_HP_05 s_reset",
        "TST_HP_10 s_reset",
        "TST_HP_11 s_reset",
        "TST_HP_04 s_reset",
        "TST_HP_02 s_reset",
        "TST_HP_08 s_reset",
        "TST_HP_07 s_reset",
        "TST_SRC_01 s_reset",
        "TST_SRC_02 s_reset",
        "TST_GRI_EL s_reset",
        "TST_DEM_01 s_control",
        "TST_DEM_03 s_control",
        "TST_DEM_02 s_control",
        "TST_BUS_TH_01 s_control",
        "TST_BUS_EL s_control",
        "TST_BUS_TH_02 s_control",
        "TST_HP_01 s_control",
        "TST_HP_09 s_control",
        "TST_HP_12 s_control",
        "TST_HP_06 s_control",
        "TST_HP_03 s_control",
        "TST_HP_05 s_control",
        "TST_HP_10 s_control",
        "TST_HP_11 s_control",
        "TST_HP_04 s_control",
        "TST_HP_02 s_control",
        "TST_HP_08 s_control",
        "TST_HP_07 s_control",
        "TST_SRC_01 s_control",
        "TST_SRC_02 s_control",
        "TST_GRI_EL s_control",
        "TST_DEM_01 s_process",
        "TST_DEM_03 s_process",
        "TST_DEM_02 s_process",
        "TST_BUS_TH_01 s_process",
        "TST_BUS_EL s_process",
        "TST_BUS_TH_02 s_process",
        "TST_HP_01 s_potential",
        "TST_HP_02 s_potential",
        "TST_HP_04 s_potential",
        "TST_HP_03 s_potential",
        "TST_HP_10 s_potential",
        "TST_HP_09 s_potential",
        "TST_HP_12 s_potential",
        "TST_HP_11 s_potential",
        "TST_HP_07 s_potential",
        "TST_HP_08 s_potential",
        "TST_HP_06 s_potential",
        "TST_HP_05 s_potential",
        "TST_HP_02 s_process",
        "TST_HP_01 s_process",
        "TST_HP_03 s_process",
        "TST_HP_04 s_process",
        "TST_HP_05 s_process",
        "TST_HP_06 s_process",
        "TST_HP_08 s_process",
        "TST_HP_07 s_process",
        "TST_HP_09 s_process",
        "TST_HP_10 s_process",
        "TST_HP_11 s_process",
        "TST_HP_12 s_process",
        "TST_SRC_01 s_process",
        "TST_SRC_02 s_process",
        "TST_GRI_EL s_process",
        "TST_BUS_TH_01 s_distribute",
        "TST_BUS_EL s_distribute",
        "TST_BUS_TH_02 s_distribute"
    ]

    simulation_parameters = Dict{String,Any}(
        "time_step_seconds" => 900,
        "time" => 0,
        "epsilon" => 1e-9
    )
    
    components = Resie.load_components(components_config, simulation_parameters)
    steps = Resie.calculate_order_of_operations(components)
    ooo = []
    for step in steps
        push!(ooo, "$(step[1]) $(step[2])")
    end
    @test pwc_ooo_astr(expected, ooo) == ""
end


function test_ooo_connected_middle_transformer()
    components_config = Dict{String,Any}(
        "TST_DEM_01" => Dict{String,Any}(
            "type" => "Demand",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "output_refs" => [],
            "energy_profile_file_path" => "./profiles/examples/general/dem_heat_nrg_var_hi-amp.prf",
            "constant_temperature" => 85,
            "scale" => 500
        ),
        "TST_DEM_02" => Dict{String,Any}(
            "type" => "Demand",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "output_refs" => [],
            "energy_profile_file_path" => "./profiles/examples/general/dem_heat_nrg_var_hi-amp.prf",
            "constant_temperature" => 90,
            "scale" => 500
        ),
        "TST_GRI_O2" => Dict{String,Any}(
            "type" => "GridConnection",
            "medium" => "m_c_g_o2",
            "control_refs" => [],
            "input_refs" => [
                "TST_01_ELY_01"
            ],
            "is_source" => false
        ),
        "TST_GRI_01" => Dict{String,Any}(
            "type" => "GridConnection",
            "medium" => "m_e_ac_230v",
            "control_refs" => [],
            "output_refs" => [
                "TST_BUS_EL"
            ],
            "is_source" => true
        ),
        "TST_GRI_02" => Dict{String,Any}(
            "type" => "GridConnection",
            "medium" => "m_e_ac_230v",
            "control_refs" => [],
            "input_refs" => [
                "TST_BUS_EL"
            ],
            "is_source" => false
        ),
        "TST_GRI_NG" => Dict{String,Any}(
            "type" => "GridConnection",
            "medium" => "m_c_g_h2",
            "control_refs" => [],
            "output_refs" => [
                "TST_CHP_01"
            ],
            "is_source" => true
        ),
        "TST_HP_01" => Dict{String,Any}(
            "type" => "HeatPump",
            "m_heat_in" => "m_h_w_lt1",
            "m_heat_out" => "m_h_w_ht1",
            "control_refs" => [],
            "output_refs" => [
                "TST_HP_02"
            ],
            "strategy" => Dict{String,Any}(
                "name" => "demand_driven"
            ),
            "power_th" => 1200,
            "output_temperature" => 50,
            "min_power_fraction" => 0.0
        ),
        "TST_HP_02" => Dict{String,Any}(
            "type" => "HeatPump",
            "m_heat_in" => "m_h_w_ht1",
            "m_heat_out" => "m_h_w_ht1",
            "control_refs" => [],
            "output_refs" => [
                "TST_BUS_TH"
            ],
            "strategy" => Dict{String,Any}(
                "name" => "demand_driven"
            ),
            "power_th" => 9000,
            "min_power_fraction" => 0.0
        ),
        "TST_01_ELY_01" => Dict{String,Any}(
            "type" => "Electrolyser",
            "control_refs" => [],
            "output_refs" => [
                "TST_HP_01",
                "TST_BUS_H2",
                "TST_GRI_O2"
            ],
            "strategy" => Dict{String,Any}(
                "name" => "demand_driven"
            ),
            "power_el" => 40000,
            "min_power_fraction" => 0.0,
            "output_temperature" => 45
        ),
        "ESS_FBH2_01" => Dict{String,Any}(
            "type" => "FuelBoiler",
            "m_fuel_in" => "m_c_g_h2",
            "control_refs" => [],
            "output_refs" => [
                "TST_DEM_02"
            ],
            "strategy" => Dict{String,Any}(
                "name" => "demand_driven"
            ),
            "min_power_fraction" => 0,
            "power_th" => 1500
        ),
        "TST_CHP_01" => Dict{String,Any}(
            "type" => "CHPP",
            "m_gas_in" => "m_c_g_h2",
            "control_refs" => [],
            "output_refs" => [
                "TST_BUS_TH",
                "TST_BUS_EL"
            ],
            "strategy" => Dict{String,Any}(
                "name" => "demand_driven",
                "load_storages" => false
            ),
            "power_gas" => 250000,
            "min_power_fraction" => 0.0
        ),
        "TST_BUS_EL" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_e_ac_230v",
            "control_refs" => [],
            "connections" => Dict{String,Any}(
                "input_order" => [
                    "TST_CHP_01",
                    "TST_GRI_01"
                ],
                "output_order" => [
                    "TST_HP_01",
                    "TST_HP_02",
                    "TST_01_ELY_01",
                    "TST_GRI_02"
                ],
                "energy_flow" => [
                    [1, 1, 1, 1],
                    [1, 1, 1, 0]
                ]
            )
        ),
        "TST_BUS_TH" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "connections" => Dict{String,Any}(
                "input_order" => [
                    "TST_HP_02",
                    "TST_CHP_01"
                ],
                "output_order" => [
                    "TST_DEM_01"
                ],
                "energy_flow" => [
                    [1],
                    [1]
                ]
            )
        ),
        "TST_BUS_H2" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_c_g_h2",
            "control_refs" => [],
            "connections" => Dict{String,Any}(
                "input_order" => [
                    "TST_01_ELY_01"
                ],
                "output_order" => [
                    "ESS_FBH2_01"
                ],
                "energy_flow" => [
                    [1]
                ]
            )
        )
    )

    expected = [
        "TST_DEM_01 s_reset",
        "TST_DEM_02 s_reset",
        "TST_BUS_TH s_reset",
        "TST_BUS_H2 s_reset",
        "TST_BUS_EL s_reset",
        "TST_HP_01 s_reset",
        "TST_CHP_01 s_reset",
        "TST_01_ELY_01 s_reset",
        "TST_HP_02 s_reset",
        "ESS_FBH2_01 s_reset",
        "TST_GRI_NG s_reset",
        "TST_GRI_01 s_reset",
        "TST_GRI_02 s_reset",
        "TST_GRI_O2 s_reset",
        "TST_DEM_01 s_control",
        "TST_DEM_02 s_control",
        "TST_BUS_TH s_control",
        "TST_BUS_H2 s_control",
        "TST_BUS_EL s_control",
        "TST_HP_01 s_control",
        "TST_CHP_01 s_control",
        "TST_01_ELY_01 s_control",
        "TST_HP_02 s_control",
        "ESS_FBH2_01 s_control",
        "TST_GRI_NG s_control",
        "TST_GRI_01 s_control",
        "TST_GRI_02 s_control",
        "TST_GRI_O2 s_control",
        "TST_DEM_01 s_process",
        "TST_DEM_02 s_process",
        "TST_BUS_TH s_process",
        "TST_BUS_H2 s_process",
        "TST_BUS_EL s_process",
        "TST_CHP_01 s_potential",
        "TST_HP_02 s_potential",
        "TST_HP_01 s_potential",
        "ESS_FBH2_01 s_potential",
        "TST_01_ELY_01 s_process",
        "TST_HP_01 s_process",
        "TST_HP_02 s_process",
        "TST_CHP_01 s_process",
        "ESS_FBH2_01 s_process",
        "TST_GRI_NG s_process",
        "TST_GRI_01 s_process",
        "TST_GRI_02 s_process",
        "TST_GRI_O2 s_process",
        "TST_BUS_TH s_distribute",
        "TST_BUS_H2 s_distribute",
        "TST_BUS_EL s_distribute"
    ]

    simulation_parameters = Dict{String,Any}(
        "time_step_seconds" => 900,
        "time" => 0,
        "epsilon" => 1e-9
    )
    
    components = Resie.load_components(components_config, simulation_parameters)
    steps = Resie.calculate_order_of_operations(components)
    ooo = []
    for step in steps
        push!(ooo, "$(step[1]) $(step[2])")
    end
    @test pwc_ooo_astr(expected, ooo) == ""
end

function test_ooo_connected_middle_transformer_variant()
    components_config = Dict{String,Any}(
        "TST_DEM_01" => Dict{String,Any}(
            "type" => "Demand",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "output_refs" => [],
            "energy_profile_file_path" => "./profiles/examples/general/dem_heat_nrg_var_hi-amp.prf",
            "constant_temperature" => 85,
            "scale" => 500
        ),
        "TST_DEM_02" => Dict{String,Any}(
            "type" => "Demand",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "output_refs" => [],
            "energy_profile_file_path" => "./profiles/examples/general/dem_heat_nrg_var_hi-amp.prf",
            "constant_temperature" => 90,
            "scale" => 500
        ),
        "TST_GRI_O2" => Dict{String,Any}(
            "type" => "GridConnection",
            "medium" => "m_c_g_o2",
            "control_refs" => [],
            "input_refs" => [
                "TST_01_ELY_01"
            ],
            "is_source" => false
        ),
        "TST_GRI_01" => Dict{String,Any}(
            "type" => "GridConnection",
            "medium" => "m_e_ac_230v",
            "control_refs" => [],
            "output_refs" => [
                "TST_BUS_EL"
            ],
            "is_source" => true
        ),
        "TST_GRI_02" => Dict{String,Any}(
            "type" => "GridConnection",
            "medium" => "m_e_ac_230v",
            "control_refs" => [],
            "input_refs" => [
                "TST_BUS_EL"
            ],
            "is_source" => false
        ),
        "TST_SRC_01" => Dict{String,Any}(
        "type" => "BoundedSupply",
        "medium" => "m_h_w_lt1",
        "control_refs"  => [],
        "output_refs" => ["TST_HP_03"],
        "constant_temperature" => 20,
        "constant_power" => 100000
        ),
        "TST_GRI_NG" => Dict{String,Any}(
            "type" => "GridConnection",
            "medium" => "m_c_g_h2",
            "control_refs" => [],
            "output_refs" => [
                "TST_CHP_01"
            ],
            "is_source" => true
        ),
        "TST_HP_01" => Dict{String,Any}(
            "type" => "HeatPump",
            "m_heat_in" => "m_h_w_lt1",
            "m_heat_out" => "m_h_w_ht1",
            "control_refs" => [],
            "output_refs" => [
                "TST_HP_02"
            ],
            "strategy" => Dict{String,Any}(
                "name" => "demand_driven"
            ),
            "power_th" => 1200,
            "output_temperature" => 50,
            "min_power_fraction" => 0.0
        ),
        "TST_HP_02" => Dict{String,Any}(
            "type" => "HeatPump",
            "m_heat_in" => "m_h_w_ht1",
            "m_heat_out" => "m_h_w_ht1",
            "control_refs" => [],
            "output_refs" => [
                "TST_BUS_TH"
            ],
            "strategy" => Dict{String,Any}(
                "name" => "demand_driven"
            ),
            "power_th" => 9000,
            "min_power_fraction" => 0.0
        ),
        "TST_HP_03" => Dict{String,Any}(
            "type" => "HeatPump",
            "m_heat_in" => "m_h_w_lt1",
            "m_heat_out" => "m_h_w_ht1",
            "control_refs" => [],
            "output_refs" => [
                "TST_BUS_TH2"
            ],
            "strategy" => Dict{String,Any}(
                "name" => "demand_driven"
            ),
            "power_th" => 9000,
            "min_power_fraction" => 0.0
        ),
        "TST_01_ELY_01" => Dict{String,Any}(
            "type" => "Electrolyser",
            "control_refs" => [],
            "output_refs" => [
                "TST_HP_01",
                "TST_BUS_H2",
                "TST_GRI_O2"
            ],
            "strategy" => Dict{String,Any}(
                "name" => "demand_driven"
            ),
            "power_el" => 40000,
            "min_power_fraction" => 0.0,
            "output_temperature" => 45
        ),
        "ESS_FBH2_01" => Dict{String,Any}(
            "type" => "FuelBoiler",
            "m_fuel_in" => "m_c_g_h2",
            "control_refs" => [],
            "output_refs" => [
                "TST_BUS_TH2"
            ],
            "strategy" => Dict{String,Any}(
                "name" => "demand_driven"
            ),
            "min_power_fraction" => 0,
            "power_th" => 1500
        ),
        "TST_CHP_01" => Dict{String,Any}(
            "type" => "CHPP",
            "m_gas_in" => "m_c_g_h2",
            "control_refs" => [],
            "output_refs" => [
                "TST_BUS_TH",
                "TST_BUS_EL"
            ],
            "strategy" => Dict{String,Any}(
                "name" => "demand_driven",
                "load_storages" => false
            ),
            "power_gas" => 250000,
            "min_power_fraction" => 0.0
        ),
        "TST_BUS_EL" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_e_ac_230v",
            "control_refs" => [],
            "connections" => Dict{String,Any}(
                "input_order" => [
                    "TST_CHP_01",
                    "TST_GRI_01"
                ],
                "output_order" => [
                    "TST_HP_01",
                    "TST_HP_02",
                    "TST_HP_03",
                    "TST_01_ELY_01",
                    "TST_GRI_02"
                ],
                "energy_flow" => [
                    [1, 1, 1, 1, 1],
                    [1, 1, 1, 1, 0]
                ]
            )
        ),
        "TST_BUS_TH" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "connections" => Dict{String,Any}(
                "input_order" => [
                    "TST_HP_02",
                    "TST_CHP_01"
                ],
                "output_order" => [
                    "TST_DEM_01"
                ],
                "energy_flow" => [
                    [1],
                    [1]
                ]
            )
        ),
        "TST_BUS_TH2" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "connections" => Dict{String,Any}(
                "input_order" => [
                    "ESS_FBH2_01",
                    "TST_HP_03"
                ],
                "output_order" => [
                    "TST_DEM_02"
                ],
                "energy_flow" => [
                    [1],
                    [1]
                ]
            )
        ),
        "TST_BUS_H2" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_c_g_h2",
            "control_refs" => [],
            "connections" => Dict{String,Any}(
                "input_order" => [
                    "TST_01_ELY_01"
                ],
                "output_order" => [
                    "ESS_FBH2_01"
                ],
                "energy_flow" => [
                    [1]
                ]
            )
        )
    )
    expected = [
        "TST_DEM_01 s_reset",
        "TST_DEM_02 s_reset",
        "TST_BUS_TH2 s_reset",
        "TST_BUS_TH s_reset",
        "TST_BUS_H2 s_reset",
        "TST_BUS_EL s_reset",
        "TST_HP_01 s_reset",
        "TST_HP_03 s_reset",
        "TST_CHP_01 s_reset",
        "TST_01_ELY_01 s_reset",
        "TST_HP_02 s_reset",
        "ESS_FBH2_01 s_reset",
        "TST_SRC_01 s_reset",
        "TST_GRI_NG s_reset",
        "TST_GRI_01 s_reset",
        "TST_GRI_02 s_reset",
        "TST_GRI_O2 s_reset",
        "TST_DEM_01 s_control",
        "TST_DEM_02 s_control",
        "TST_BUS_TH2 s_control",
        "TST_BUS_TH s_control",
        "TST_BUS_H2 s_control",
        "TST_BUS_EL s_control",
        "TST_HP_01 s_control",
        "TST_HP_03 s_control",
        "TST_CHP_01 s_control",
        "TST_01_ELY_01 s_control",
        "TST_HP_02 s_control",
        "ESS_FBH2_01 s_control",
        "TST_SRC_01 s_control",
        "TST_GRI_NG s_control",
        "TST_GRI_01 s_control",
        "TST_GRI_02 s_control",
        "TST_GRI_O2 s_control",
        "TST_DEM_01 s_process",
        "TST_DEM_02 s_process",
        "TST_BUS_TH2 s_process",
        "TST_BUS_TH s_process",
        "TST_BUS_H2 s_process",
        "TST_BUS_EL s_process",
        "TST_CHP_01 s_potential",
        "TST_HP_02 s_potential",
        "TST_HP_01 s_potential",
        "TST_HP_03 s_potential",
        "ESS_FBH2_01 s_potential",
        "TST_01_ELY_01 s_process",
        "TST_HP_01 s_process",
        "TST_HP_02 s_process",
        "TST_CHP_01 s_process",
        "ESS_FBH2_01 s_process",
        "TST_HP_03 s_process",
        "TST_SRC_01 s_process",
        "TST_GRI_NG s_process",
        "TST_GRI_01 s_process",
        "TST_GRI_02 s_process",
        "TST_GRI_O2 s_process",
        "TST_BUS_TH2 s_distribute",
        "TST_BUS_TH s_distribute",
        "TST_BUS_H2 s_distribute",
        "TST_BUS_EL s_distribute"
    ]

    simulation_parameters = Dict{String,Any}(
        "time_step_seconds" => 900,
        "time" => 0,
        "epsilon" => 1e-9
    )
    
    components = Resie.load_components(components_config, simulation_parameters)
    steps = Resie.calculate_order_of_operations(components)
    ooo = []
    for step in steps
        push!(ooo, "$(step[1]) $(step[2])")
    end
    @test pwc_ooo_astr(expected, ooo) == ""
end

function test_ooo_circle_grid_input_denied()
    components_config = Dict{String,Any}(
        "TST_DEM_01" => Dict{String,Any}(
            "type" => "Demand",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "output_refs" => [],
            "energy_profile_file_path" => "./profiles/examples/general/dem_heat_nrg_var_hi-amp.prf",
            "constant_temperature" => 85,
            "scale" => 1000
        ),
        "TST_SRC_01" => Dict{String,Any}(
            "type" => "BoundedSupply",
            "medium" => "m_h_w_lt1",
            "control_refs" => [],
            "output_refs" => [
                "TST_HP_01"
            ],
            "constant_power" => 500000,
            "constant_temperature" => 5
        ),
        "TST_GRI_01" => Dict{String,Any}(
            "type" => "GridConnection",
            "medium" => "m_e_ac_230v",
            "control_refs" => [],
            "output_refs" => [
                "TST_BUS_EL"
            ],
            "is_source" => true
        ),
        "TST_GRI_02" => Dict{String,Any}(
            "type" => "GridConnection",
            "medium" => "m_e_ac_230v",
            "control_refs" => [],
            "input_refs" => [
                "TST_BUS_EL"
            ],
            "is_source" => false
        ),
        "TST_GRI_natgas" => Dict{String,Any}(
            "type" => "GridConnection",
            "medium" => "m_c_natgas",
            "control_refs" => [],
            "output_refs" => [
                "TST_CHP_01"
            ],
            "is_source" => true
        ),
        "TST_HP_01" => Dict{String,Any}(
            "type" => "HeatPump",
            "m_heat_in" => "m_h_w_lt1",
            "m_heat_out" => "m_h_w_ht1",
            "control_refs" => [],
            "output_refs" => [
                "TST_BUS_TH"
            ],
            "strategy" => Dict{String,Any}(
                "name" => "demand_driven"
            ),
            "power_th" => 500,
            "min_power_fraction" => 0.0
        ),
        "TST_CHP_01" => Dict{String,Any}(
            "type" => "CHPP",
            "m_gas_in" => "m_c_natgas",
            "m_heat_out" => "m_h_w_ht1",
            "control_refs" => [],
            "output_refs" => [
                "TST_BUS_TH",
                "TST_BUS_EL"
            ],
            "strategy" => Dict{String,Any}(
                "name" => "demand_driven"
            ),
            "power_gas" => 250000,
            "min_power_fraction" => 0.0
        ),
        "TST_BUS_EL" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_e_ac_230v",
            "control_refs" => [],
            "connections" => Dict{String,Any}(
                "input_order" => [
                    "TST_CHP_01",
                    "TST_GRI_01"
                ],
                "output_order" => [
                    "TST_HP_01",
                    "TST_GRI_02"
                ],
                "energy_flow" => [
                    [1, 0],
                    [1, 0]
                ]
            )
        ),
        "TST_BUS_TH" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "connections" => Dict{String,Any}(
                "input_order" => [
                    "TST_HP_01",
                    "TST_CHP_01"
                ],
                "output_order" => [
                    "TST_DEM_01"
                ],
                "energy_flow" => [
                    [1],
                    [1]
                ]
            )
        )
    )

    expected = [
        "TST_DEM_01 s_reset",
        "TST_BUS_TH s_reset",
        "TST_BUS_EL s_reset",
        "TST_HP_01 s_reset",
        "TST_CHP_01 s_reset",
        "TST_SRC_01 s_reset",
        "TST_GRI_natgas s_reset",
        "TST_GRI_01 s_reset",
        "TST_GRI_02 s_reset",
        "TST_DEM_01 s_control",
        "TST_BUS_TH s_control",
        "TST_BUS_EL s_control",
        "TST_HP_01 s_control",
        "TST_CHP_01 s_control",
        "TST_SRC_01 s_control",
        "TST_GRI_natgas s_control",
        "TST_GRI_01 s_control",
        "TST_GRI_02 s_control",
        "TST_DEM_01 s_process",
        "TST_BUS_TH s_process",
        "TST_BUS_EL s_process",
        "TST_HP_01 s_potential",
        "TST_CHP_01 s_potential",
        "TST_HP_01 s_process",
        "TST_CHP_01 s_process",
        "TST_SRC_01 s_process",
        "TST_GRI_natgas s_process",
        "TST_GRI_01 s_process",
        "TST_GRI_02 s_process",
        "TST_BUS_TH s_distribute",
        "TST_BUS_EL s_distribute"
    ]

    simulation_parameters = Dict{String,Any}(
        "time_step_seconds" => 900,
        "time" => 0,
        "epsilon" => 1e-9
    )
    
    components = Resie.load_components(components_config, simulation_parameters)
    steps = Resie.calculate_order_of_operations(components)
    ooo = []
    for step in steps
        push!(ooo, "$(step[1]) $(step[2])")
    end
    @test pwc_ooo_astr(expected, ooo) == ""
end

function test_ooo_circle_grid_input_allowed()
    components_config = Dict{String,Any}(
        "TST_DEM_01" => Dict{String,Any}(
            "type" => "Demand",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "output_refs" => [],
            "energy_profile_file_path" => "./profiles/examples/general/dem_heat_nrg_var_hi-amp.prf",
            "constant_temperature" => 85,
            "scale" => 1000
        ),
        "TST_SRC_01" => Dict{String,Any}(
            "type" => "BoundedSupply",
            "medium" => "m_h_w_lt1",
            "control_refs" => [],
            "output_refs" => [
                "TST_HP_01"
            ],
            "constant_power" => 500000,
            "constant_temperature" => 5
        ),
        "TST_GRI_01" => Dict{String,Any}(
            "type" => "GridConnection",
            "medium" => "m_e_ac_230v",
            "control_refs" => [],
            "output_refs" => [
                "TST_BUS_EL"
            ],
            "is_source" => true
        ),
        "TST_GRI_02" => Dict{String,Any}(
            "type" => "GridConnection",
            "medium" => "m_e_ac_230v",
            "control_refs" => [],
            "input_refs" => [
                "TST_BUS_EL"
            ],
            "is_source" => false
        ),
        "TST_GRI_natgas" => Dict{String,Any}(
            "type" => "GridConnection",
            "medium" => "m_c_natgas",
            "control_refs" => [],
            "output_refs" => [
                "TST_CHP_01"
            ],
            "is_source" => true
        ),
        "TST_HP_01" => Dict{String,Any}(
            "type" => "HeatPump",
            "m_heat_in" => "m_h_w_lt1",
            "m_heat_out" => "m_h_w_ht1",
            "control_refs" => [],
            "output_refs" => [
                "TST_BUS_TH"
            ],
            "strategy" => Dict{String,Any}(
                "name" => "demand_driven"
            ),
            "power_th" => 500,
            "min_power_fraction" => 0.0
        ),
        "TST_CHP_01" => Dict{String,Any}(
            "type" => "CHPP",
            "m_gas_in" => "m_c_natgas",
            "m_heat_out" => "m_h_w_ht1",
            "control_refs" => [],
            "output_refs" => [
                "TST_BUS_TH",
                "TST_BUS_EL"
            ],
            "strategy" => Dict{String,Any}(
                "name" => "demand_driven"
            ),
            "power_gas" => 250000,
            "min_power_fraction" => 0.0
        ),
        "TST_BUS_EL" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_e_ac_230v",
            "control_refs" => [],
            "connections" => Dict{String,Any}(
                "input_order" => [
                    "TST_CHP_01",
                    "TST_GRI_01"
                ],
                "output_order" => [
                    "TST_HP_01",
                    "TST_GRI_02"
                ],
                "energy_flow" => [
                    [1, 1],
                    [1, 0]
                ]
            )
        ),
        "TST_BUS_TH" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "connections" => Dict{String,Any}(
                "input_order" => [
                    "TST_HP_01",
                    "TST_CHP_01"
                ],
                "output_order" => [
                    "TST_DEM_01"
                ],
                "energy_flow" => [
                    [1],
                    [1]
                ]
            )
        )
    )

    expected = [
        "TST_DEM_01 s_reset",
        "TST_BUS_TH s_reset",
        "TST_BUS_EL s_reset",
        "TST_HP_01 s_reset",
        "TST_CHP_01 s_reset",
        "TST_SRC_01 s_reset",
        "TST_GRI_natgas s_reset",
        "TST_GRI_01 s_reset",
        "TST_GRI_02 s_reset",
        "TST_DEM_01 s_control",
        "TST_BUS_TH s_control",
        "TST_BUS_EL s_control",
        "TST_HP_01 s_control",
        "TST_CHP_01 s_control",
        "TST_SRC_01 s_control",
        "TST_GRI_natgas s_control",
        "TST_GRI_01 s_control",
        "TST_GRI_02 s_control",
        "TST_DEM_01 s_process",
        "TST_BUS_TH s_process",
        "TST_BUS_EL s_process",
        "TST_HP_01 s_potential",
        "TST_CHP_01 s_potential",
        "TST_HP_01 s_process",
        "TST_CHP_01 s_process",
        "TST_SRC_01 s_process",
        "TST_GRI_natgas s_process",
        "TST_GRI_01 s_process",
        "TST_GRI_02 s_process",
        "TST_BUS_TH s_distribute",
        "TST_BUS_EL s_distribute"
    ]

    simulation_parameters = Dict{String,Any}(
        "time_step_seconds" => 900,
        "time" => 0,
        "epsilon" => 1e-9
    )
    
    components = Resie.load_components(components_config, simulation_parameters)
    steps = Resie.calculate_order_of_operations(components)
    ooo = []
    for step in steps
        push!(ooo, "$(step[1]) $(step[2])")
    end
    @test pwc_ooo_astr(expected, ooo) == ""
end



function test_ooo_circle_middle_transformer_input()
    components_config = Dict{String,Any}(
        "TST_DEM_01" => Dict{String,Any}(
            "type" => "Demand",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "output_refs" => [],
            "energy_profile_file_path" => "./profiles/examples/general/dem_heat_nrg_var_hi-amp.prf",
            "constant_temperature" => 85,
            "scale" => 1000
        ),
        "TST_SRC_01" => Dict{String,Any}(
            "type" => "BoundedSupply",
            "medium" => "m_h_w_lt1",
            "control_refs" => [],
            "output_refs" => [
                "TST_HP_01"
            ],
            "constant_power" => 500000,
            "constant_temperature" => 5
        ),
        "TST_GRI_01" => Dict{String,Any}(
            "type" => "GridConnection",
            "medium" => "m_e_ac_230v",
            "control_refs" => [],
            "output_refs" => [
                "TST_BUS_EL"
            ],
            "is_source" => true
        ),
        "TST_GRI_02" => Dict{String,Any}(
            "type" => "GridConnection",
            "medium" => "m_e_ac_230v",
            "control_refs" => [],
            "input_refs" => [
                "TST_BUS_EL"
            ],
            "is_source" => false
        ),
        "TST_GRI_natgas" => Dict{String,Any}(
            "type" => "GridConnection",
            "medium" => "m_c_natgas",
            "control_refs" => [],
            "output_refs" => [
                "TST_CHP_01"
            ],
            "is_source" => true
        ),
        "TST_HP_01" => Dict{String,Any}(
            "type" => "HeatPump",
            "m_heat_in" => "m_h_w_lt1",
            "m_heat_out" => "m_h_w_ht1",
            "control_refs" => [],
            "output_refs" => [
                "TST_HP_02"
            ],
            "strategy" => Dict{String,Any}(
                "name" => "demand_driven"
            ),
            "power_th" => 500,
            "output_temperature" => 65,
            "min_power_fraction" => 0.0
        ),
        "TST_HP_02" => Dict{String,Any}(
            "type" => "HeatPump",
            "m_heat_in" => "m_h_w_ht1",
            "m_heat_out" => "m_h_w_ht1",
            "control_refs" => [],
            "output_refs" => [
                "TST_BUS_TH"
            ],
            "strategy" => Dict{String,Any}(
                "name" => "demand_driven"
            ),
            "power_th" => 500,
            "min_power_fraction" => 0.0
        ),

        "TST_CHP_01" => Dict{String,Any}(
            "type" => "CHPP",
            "m_gas_in" => "m_c_natgas",
            "m_heat_out" => "m_h_w_ht1",
            "control_refs" => [],
            "output_refs" => [
                "TST_BUS_TH",
                "TST_BUS_EL"
            ],
            "strategy" => Dict{String,Any}(
                "name" => "demand_driven"
            ),
            "power_gas" => 250000,
            "min_power_fraction" => 0.0
        ),
        "TST_BUS_EL" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_e_ac_230v",
            "control_refs" => [],
            "connections" => Dict{String,Any}(
                "input_order" => [
                    "TST_CHP_01",
                    "TST_GRI_01"
                ],
                "output_order" => [
                    "TST_HP_01",
                    "TST_HP_02",
                    "TST_GRI_02"
                ],
                "energy_flow" => [
                    [1, 1, 0],
                    [1, 1, 0]
                ]
            )
        ),
        "TST_BUS_TH" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "connections" => Dict{String,Any}(
                "input_order" => [
                    "TST_HP_02",
                    "TST_CHP_01"
                ],
                "output_order" => [
                    "TST_DEM_01"
                ],
                "energy_flow" => [
                    [1],
                    [1]
                ]
            )
        )
    )

    expected = [
        "TST_DEM_01 s_reset",
        "TST_BUS_TH s_reset",
        "TST_BUS_EL s_reset",
        "TST_HP_02 s_reset",
        "TST_HP_01 s_reset",
        "TST_CHP_01 s_reset",
        "TST_SRC_01 s_reset",
        "TST_GRI_natgas s_reset",
        "TST_GRI_01 s_reset",
        "TST_GRI_02 s_reset",
        "TST_DEM_01 s_control",
        "TST_BUS_TH s_control",
        "TST_BUS_EL s_control",
        "TST_HP_02 s_control",
        "TST_HP_01 s_control",
        "TST_CHP_01 s_control",
        "TST_SRC_01 s_control",
        "TST_GRI_natgas s_control",
        "TST_GRI_01 s_control",
        "TST_GRI_02 s_control",
        "TST_DEM_01 s_process",
        "TST_BUS_TH s_process",
        "TST_BUS_EL s_process",
        "TST_HP_01 s_potential",
        "TST_CHP_01 s_potential",
        "TST_HP_02 s_process",
        "TST_HP_01 s_process",
        "TST_CHP_01 s_process",
        "TST_SRC_01 s_process",
        "TST_GRI_natgas s_process",
        "TST_GRI_01 s_process",
        "TST_GRI_02 s_process",
        "TST_BUS_TH s_distribute",
        "TST_BUS_EL s_distribute"
    ]

    simulation_parameters = Dict{String,Any}(
        "time_step_seconds" => 900,
        "time" => 0,
        "epsilon" => 1e-9
    )
    
    components = Resie.load_components(components_config, simulation_parameters)
    steps = Resie.calculate_order_of_operations(components)
    ooo = []
    for step in steps
        push!(ooo, "$(step[1]) $(step[2])")
    end
    @test pwc_ooo_astr(expected, ooo) == ""
end

function test_ooo_circle_variant()
    components_config = Dict{String,Any}(
        "TST_DEM_01" => Dict{String,Any}(
            "type" => "Demand",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "output_refs" => [],
            "energy_profile_file_path" => "./profiles/examples/general/dem_heat_nrg_var_hi-amp.prf",
            "constant_temperature" => 85,
            "scale" => 1000
        ),
        "TST_SRC_01" => Dict{String,Any}(
            "type" => "BoundedSupply",
            "medium" => "m_h_w_lt1",
            "control_refs" => [],
            "output_refs" => [
                "TST_HP_01"
            ],
            "constant_power" => 500000,
            "constant_temperature" => 5
        ),
        "TST_GRI_01" => Dict{String,Any}(
            "type" => "GridConnection",
            "medium" => "m_e_ac_230v",
            "control_refs" => [],
            "output_refs" => [
                "TST_BUS_EL"
            ],
            "is_source" => true
        ),
        "TST_GRI_02" => Dict{String,Any}(
            "type" => "GridConnection",
            "medium" => "m_e_ac_230v",
            "control_refs" => [],
            "input_refs" => [
                "TST_BUS_EL"
            ],
            "is_source" => false
        ),
        "TST_GRI_natgas_1" => Dict{String,Any}(
            "type" => "GridConnection",
            "medium" => "m_c_natgas",
            "control_refs" => [],
            "output_refs" => [
                "TST_CHP_01"
            ],
            "is_source" => true
        ),
        "TST_GRI_natgas_2" => Dict{String,Any}(
            "type" => "GridConnection",
            "medium" => "m_c_natgas",
            "control_refs" => [],
            "output_refs" => [
                "TST_CHP_02"
            ],
            "is_source" => true
        ),
        "TST_HP_01" => Dict{String,Any}(
            "type" => "HeatPump",
            "m_heat_in" => "m_h_w_lt1",
            "m_heat_out" => "m_h_w_ht1",
            "control_refs" => [],
            "output_refs" => [
                "TST_BUS_TH"
            ],
            "strategy" => Dict{String,Any}(
                "name" => "demand_driven"
            ),
            "power_th" => 500,
            "min_power_fraction" => 0.0
        ),
        "TST_CHP_01" => Dict{String,Any}(
            "type" => "CHPP",
            "m_gas_in" => "m_c_natgas",
            "m_heat_out" => "m_h_w_ht1",
            "control_refs" => [],
            "output_refs" => [
                "TST_BUS_TH",
                "TST_BUS_EL"
            ],
            "strategy" => Dict{String,Any}(
                "name" => "demand_driven"
            ),
            "power_gas" => 500,
            "min_power_fraction" => 0.0
        ),
        "TST_CHP_02" => Dict{String,Any}(
            "type" => "CHPP",
            "m_gas_in" => "m_c_natgas",
            "m_heat_out" => "m_h_w_ht1",
            "control_refs" => [],
            "output_refs" => [
                "TST_BUS_TH",
                "TST_BUS_EL"
            ],
            "strategy" => Dict{String,Any}(
                "name" => "demand_driven"
            ),
            "power_gas" => 120000,
            "min_power_fraction" => 0.0
        ),
        "TST_BUS_EL" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_e_ac_230v",
            "control_refs" => [],
            "connections" => Dict{String,Any}(
                "input_order" => [
                    "TST_CHP_01",
                    "TST_CHP_02",
                    "TST_GRI_01"
                ],
                "output_order" => [
                    "TST_HP_01",
                    "TST_GRI_02"
                ],
                "energy_flow" => [
                    [1, 1],
                    [1, 0],
                    [1, 0]
                ]
            )
        ),
        "TST_BUS_TH" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "connections" => Dict{String,Any}(
                "input_order" => [
                    "TST_HP_01",
                    "TST_CHP_01",
                    "TST_CHP_02"

                ],
                "output_order" => [
                    "TST_DEM_01"
                ],
                "energy_flow" => [
                    [1],
                    [1],
                    [1]
                ]
            )
        )
    )

    expected = [
        "TST_DEM_01 s_reset",
        "TST_BUS_TH s_reset",
        "TST_BUS_EL s_reset",
        "TST_HP_01 s_reset",
        "TST_CHP_01 s_reset",
        "TST_CHP_02 s_reset",
        "TST_GRI_natgas_2 s_reset",
        "TST_SRC_01 s_reset",
        "TST_GRI_01 s_reset",
        "TST_GRI_natgas_1 s_reset",
        "TST_GRI_02 s_reset",
        "TST_DEM_01 s_control",
        "TST_BUS_TH s_control",
        "TST_BUS_EL s_control",
        "TST_HP_01 s_control",
        "TST_CHP_01 s_control",
        "TST_CHP_02 s_control",
        "TST_GRI_natgas_2 s_control",
        "TST_SRC_01 s_control",
        "TST_GRI_01 s_control",
        "TST_GRI_natgas_1 s_control",
        "TST_GRI_02 s_control",
        "TST_DEM_01 s_process",
        "TST_BUS_TH s_process",
        "TST_BUS_EL s_process",
        "TST_HP_01 s_potential",
        "TST_CHP_02 s_potential",
        "TST_CHP_01 s_potential",
        "TST_HP_01 s_process",
        "TST_CHP_01 s_process",
        "TST_CHP_02 s_process",
        "TST_GRI_natgas_2 s_process",
        "TST_SRC_01 s_process",
        "TST_GRI_01 s_process",
        "TST_GRI_natgas_1 s_process",
        "TST_GRI_02 s_process",
        "TST_BUS_TH s_distribute",
        "TST_BUS_EL s_distribute"
    ]

    simulation_parameters = Dict{String,Any}(
        "time_step_seconds" => 900,
        "time" => 0,
        "epsilon" => 1e-9
    )
    
    components = Resie.load_components(components_config, simulation_parameters)
    steps = Resie.calculate_order_of_operations(components)
    ooo = []
    for step in steps
        push!(ooo, "$(step[1]) $(step[2])")
    end
    @test pwc_ooo_astr(expected, ooo) == ""
end

function test_ooo_circle_middle_transformer_interconnections()
    components_config = Dict{String,Any}(
        "TST_DEM_01" => Dict{String,Any}(
            "type" => "Demand",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "output_refs" => [],
            "energy_profile_file_path" => "./profiles/examples/general/dem_heat_nrg_var_hi-amp.prf",
            "constant_temperature" => 85,
            "scale" => 500
        ),
        "TST_DEM_02" => Dict{String,Any}(
            "type" => "Demand",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "output_refs" => [],
            "energy_profile_file_path" => "./profiles/examples/general/dem_heat_nrg_var_hi-amp.prf",
            "constant_temperature" => 50,
            "scale" => 500
        ),
        "TST_SRC_01" => Dict{String,Any}(
            "type" => "BoundedSupply",
            "medium" => "m_h_w_lt1",
            "control_refs" => [],
            "output_refs" => [
                "TST_HP_01"
            ],
            "constant_power" => 500000,
            "constant_temperature" => 20
        ),
        "TST_GRI_natgas" => Dict{String,Any}(
            "type" => "GridConnection",
            "medium" => "m_c_natgas",
            "control_refs" => [],
            "output_refs" => [
                "TST_CHP_01"
            ],
            "is_source" => true
        ),
        "TST_GRI_01" => Dict{String,Any}(
            "type" => "GridConnection",
            "medium" => "m_e_ac_230v",
            "control_refs" => [],
            "output_refs" => [
                "TST_BUS_EL"
            ],
            "is_source" => true
        ),
        "TST_GRI_02" => Dict{String,Any}(
            "type" => "GridConnection",
            "medium" => "m_e_ac_230v",
            "control_refs" => [],
            "input_refs" => [
                "TST_BUS_EL"
            ],
            "is_source" => false
        ),
        "TST_HP_01" => Dict{String,Any}(
            "type" => "HeatPump",
            "control_refs" => [],
            "output_refs" => [
                "TST_BUS_TH"
            ],
            "strategy" => Dict{String,Any}(
                "name" => "demand_driven"
            ),
            "power_th" => 1200,
            "min_power_fraction" => 0.0
        ),
        "TST_HP_02" => Dict{String,Any}(
            "type" => "HeatPump",
            "m_heat_in" => "m_h_w_ht1",
            "m_heat_out" => "m_h_w_ht1",
            "control_refs" => [],
            "output_refs" => [
                "TST_BUS_TH"
            ],
            "strategy" => Dict{String,Any}(
                "name" => "demand_driven"
            ),
            "power_th" => 9000,
            "input_temperature" => 30,
            "min_power_fraction" => 0.0
        ),
        "TST_HP_03" => Dict{String,Any}(
            "type" => "HeatPump",
            "m_heat_in" => "m_h_w_ht1",
            "m_heat_out" => "m_h_w_ht1",
            "control_refs" => [],
            "output_refs" => [
                "TST_DEM_01"
            ],
            "strategy" => Dict{String,Any}(
                "name" => "demand_driven"
            ),
            "power_th" => 9000,
            "input_temperature" => 50,
            "min_power_fraction" => 0.0
        ),
        "TST_CHP_01" => Dict{String,Any}(
            "type" => "CHPP",
            "m_gas_in" => "m_c_natgas",
            "m_heat_out" => "m_h_w_ht1",
            "control_refs" => [],
            "output_refs" => [
                "TST_HP_02",
                "TST_BUS_EL"
            ],
            "strategy" => Dict{String,Any}(
                "name" => "demand_driven"
            ),
            "power_gas" => 250000,
            "min_power_fraction" => 0.0
        ),
        "TST_BUS_EL" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_e_ac_230v",
            "control_refs" => [],
            "connections" => Dict{String,Any}(
                "input_order" => [
                    "TST_CHP_01",
                    "TST_GRI_01"
                ],
                "output_order" => [
                    "TST_HP_02",
                    "TST_HP_01",
                    "TST_HP_03",
                    "TST_GRI_02"
                ],
                "energy_flow" => [
                    [0, 1, 1, 0],
                    [1, 1, 1, 1]
                ]
            )
        ),
        "TST_BUS_TH" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "connections" => Dict{String,Any}(
                "input_order" => [
                    "TST_HP_01",
                    "TST_HP_02",
                    "TST_BFT_TH_01"
                ],
                "output_order" => [
                    "TST_HP_03",
                    "TST_DEM_02",
                    "TST_BFT_TH_01"
                ],
                "energy_flow" => [
                    [1, 1, 1],
                    [1, 1, 1],
                    [1, 1, 1]
                ]
            )
        ),
        "TST_BFT_TH_01" => Dict{String,Any}(
            "type" => "BufferTank",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "output_refs" => [
                "TST_BUS_TH"
            ],
            "capacity" => 150000,
            "load" => 70000,
            "high_temperature" => 50.0,
            "low_temperature" => 20.0,
            "use_adaptive_temperature" => false
        )
    )

    expected = [
        "TST_DEM_01 s_reset",
        "TST_DEM_02 s_reset",
        "TST_BUS_TH s_reset",
        "TST_BUS_EL s_reset",
        "TST_HP_01 s_reset",
        "TST_HP_03 s_reset",
        "TST_CHP_01 s_reset",
        "TST_HP_02 s_reset",
        "TST_BFT_TH_01 s_reset",
        "TST_GRI_natgas s_reset",
        "TST_SRC_01 s_reset",
        "TST_GRI_01 s_reset",
        "TST_GRI_02 s_reset",
        "TST_DEM_01 s_control",
        "TST_DEM_02 s_control",
        "TST_BUS_TH s_control",
        "TST_BUS_EL s_control",
        "TST_HP_01 s_control",
        "TST_HP_03 s_control",
        "TST_CHP_01 s_control",
        "TST_HP_02 s_control",
        "TST_BFT_TH_01 s_control",
        "TST_GRI_natgas s_control",
        "TST_SRC_01 s_control",
        "TST_GRI_01 s_control",
        "TST_GRI_02 s_control",
        "TST_DEM_01 s_process",
        "TST_DEM_02 s_process",
        "TST_BUS_TH s_process",
        "TST_BUS_EL s_process",
        "TST_HP_01 s_potential",
        "TST_HP_03 s_potential",
        "TST_HP_02 s_potential",
        "TST_CHP_01 s_process",
        "TST_HP_03 s_process",
        "TST_HP_01 s_process",
        "TST_HP_02 s_process",
        "TST_BFT_TH_01 s_process",
        "TST_BFT_TH_01 s_load",
        "TST_GRI_natgas s_process",
        "TST_SRC_01 s_process",
        "TST_GRI_01 s_process",
        "TST_GRI_02 s_process",
        "TST_BUS_TH s_distribute",
        "TST_BUS_EL s_distribute"
    ]

    simulation_parameters = Dict{String,Any}(
        "time_step_seconds" => 900,
        "time" => 0,
        "epsilon" => 1e-9
    )
    
    components = Resie.load_components(components_config, simulation_parameters)
    steps = Resie.calculate_order_of_operations(components)
    ooo = []
    for step in steps
        push!(ooo, "$(step[1]) $(step[2])")
    end
    @test pwc_ooo_astr(expected, ooo) == ""
end
@testset "transformer_chains_ooo" begin
    test_ooo_middle_bus()
    test_ooo_middle_bus_different_order()
    test_ooo_middle_transformer()
    test_ooo_parallels()
    test_ooo_parallels_different_order()
    test_ooo_parallels_in_chain()
    test_ooo_parallels_in_a_row()
    test_ooo_connected_middle_busses()
    test_ooo_connected_middle_transformer()
    test_ooo_connected_middle_transformer_variant()
    test_ooo_circle_grid_input_allowed()
    test_ooo_circle_grid_input_denied()
    test_ooo_circle_middle_transformer_input()
    test_ooo_circle_variant()
    test_ooo_circle_middle_transformer_interconnections()
end