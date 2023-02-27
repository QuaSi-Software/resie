using Debugger
using Test
using Resie
using Resie.EnergySystems

function test_load_from_dict()
    systems_config = Dict{String, Any}(
        "TST_BT_01" => Dict{String, Any}(
            "type" => "BufferTank",
            "control_refs" => [],
            "production_refs" => [],
            "capacity" => 40000,
            "load" => 20000,
            "strategy" => Dict{String, Any}(
                "name" => "Default"
            )
        ),
        "TST_HP_01" => Dict{String, Any}(
            "type" => "HeatPump",
            "control_refs" => ["TST_BT_01"],
            "production_refs" => [],
            "strategy" => Dict{String, Any}(
                "name" => "storage_driven",
                "high_threshold" => 0.5,
                "low_threshold" => 0.1
            ),
            "power" => 20000,
            "fixed_cop" => 3.0
        ),
    )

    systems = Resie.load_systems(systems_config)
    @test length(keys(systems)) == 2
    @test typeof(systems["TST_BT_01"]) == Resie.EnergySystems.BufferTank
    @test systems["TST_BT_01"].sys_function == Resie.EnergySystems.sf_storage
    @test systems["TST_HP_01"].power == 20000
end

function test_ooo_for_heat_pumps_wrong()
    systems_config = Dict{String, Any}(
        "TST_DEM_01" => Dict{String, Any}(
            "type" => "Demand",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "production_refs" => [],
            "energy_profile_file_path" => "./profiles/tests/demand_heating_energy.prf",
            "temperature_profile_file_path" => "./profiles/tests/demand_heating_temperature.prf",
            "scale" => 1500
        ),
        "TST_SRC_01" => Dict{String, Any}(
            "type" => "DispatchableSupply",
            "medium" => "m_h_w_lt1",
            "control_refs" => [],
            "production_refs" => ["TST_HP_01"],
            "max_power_profile_file_path" => "./profiles/tests/demand_heating_energy.prf",
            "temperature_profile_file_path" => "./profiles/tests/demand_heating_temperature.prf",
            "scale" => 6000
        ),
        "TST_GRI_01" => Dict{String, Any}(
            "type" => "GridConnection",
            "medium" => "m_e_ac_230v",
            "control_refs" => [],
            "production_refs" => ["TST_HP_01"],
            "is_source" => true,
        ),
        "TST_HP_01" => Dict{String, Any}(
            "type" => "HeatPump",
            "control_refs" => ["TST_DEM_01"],
            "production_refs" => ["TST_DEM_01"],
            "strategy" => Dict{String, Any}(
                "name" => "demand_driven",
            ),
            "power" => 12000,
            "fixed_cop" => 3.0
        ),
    )

    expected = [
        ("TST_DEM_01", EnergySystems.s_reset),
        ("TST_HP_01", EnergySystems.s_reset),
        ("TST_SRC_01", EnergySystems.s_reset),
        ("TST_GRI_01", EnergySystems.s_reset),
        ("TST_DEM_01", EnergySystems.s_control),
        ("TST_HP_01", EnergySystems.s_control),
        ("TST_SRC_01", EnergySystems.s_control),
        ("TST_GRI_01", EnergySystems.s_control),
        ("TST_DEM_01", EnergySystems.s_produce),
        ("TST_HP_01", EnergySystems.s_produce),
        ("TST_SRC_01", EnergySystems.s_produce),
        ("TST_GRI_01", EnergySystems.s_produce),
    ]

    systems = Resie.load_systems(systems_config)
    ooo = Resie.order_of_operations(systems)
    @test all(ooo .== expected)
end

function test_ooo_bus_to_bus()
    systems_config = Dict{String, Any}(
        "TST_GRI_01" => Dict{String, Any}(
            "type" => "GridConnection",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "production_refs" => ["TST_BUS_01"],
            "is_source" => true,
        ),
        "TST_BUS_01" => Dict{String, Any}(
            "type" => "Bus",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "production_refs" => ["TST_BUS_02", "TST_BFT_01"],
            "input_priorities" => ["TST_BFT_01", "TST_GRI_01"]
        ),
        "TST_BFT_01" => Dict{String, Any}(
            "type" => "BufferTank",
            "control_refs" => [],
            "production_refs" => [
                "TST_BUS_01"
            ],
            "capacity" => 40000,
            "load" => 20000
        ),
        "TST_BUS_02" => Dict{String, Any}(
            "type" => "Bus",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "production_refs" => ["TST_DEM_01", "TST_BFT_02"],
            "input_priorities" => ["TST_BFT_02", "TST_BUS_01"]
        ),
        "TST_BFT_02" => Dict{String, Any}(
            "type" => "BufferTank",
            "control_refs" => [],
            "production_refs" => [
                "TST_BUS_02"
            ],
            "capacity" => 20000,
            "load" => 10000
        ),
        "TST_DEM_01" => Dict{String, Any}(
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
        ("TST_BUS_02", EnergySystems.s_reset),
        ("TST_BUS_01", EnergySystems.s_reset),
        ("TST_BFT_02", EnergySystems.s_reset),
        ("TST_BFT_01", EnergySystems.s_reset),
        ("TST_GRI_01", EnergySystems.s_reset),
        ("TST_DEM_01", EnergySystems.s_control),
        ("TST_BUS_02", EnergySystems.s_control),
        ("TST_BFT_02", EnergySystems.s_control),
        ("TST_BUS_01", EnergySystems.s_control),
        ("TST_BFT_01", EnergySystems.s_control),
        ("TST_GRI_01", EnergySystems.s_control),
        ("TST_DEM_01", EnergySystems.s_produce),
        ("TST_BUS_02", EnergySystems.s_produce),
        ("TST_BFT_02", EnergySystems.s_produce),
        ("TST_BUS_01", EnergySystems.s_produce),
        ("TST_BFT_01", EnergySystems.s_produce),
        ("TST_BFT_02", EnergySystems.s_load),
        ("TST_BFT_01", EnergySystems.s_load),
        ("TST_GRI_01", EnergySystems.s_produce),
        ("TST_BUS_02", EnergySystems.s_distribute),
        ("TST_BUS_01", EnergySystems.s_distribute),
    ]

    systems = Resie.load_systems(systems_config)
    ooo = Resie.order_of_operations(systems)
    @test all(ooo .== expected)
end

function test_ooo_storage_loading_switch()
    systems_config = Dict{String, Any}(
        "TST_GRI_01" => Dict{String, Any}(
            "type" => "GridConnection",
            "medium" => "m_c_g_natgas",
            "control_refs" => [],
            "production_refs" => ["TST_GBO_01"],
            "is_source" => true,
        ),
        "TST_GRI_02" => Dict{String, Any}(
            "type" => "GridConnection",
            "medium" => "m_c_g_natgas",
            "control_refs" => [],
            "production_refs" => ["TST_GBO_02"],
            "is_source" => true,
        ),
        "TST_GBO_01" => Dict{String, Any}(
            "type" => "GasBoiler",
            "control_refs" => ["TST_BFT_01"],
            "production_refs" => [
                "TST_BUS_01"
            ],
            "strategy" => Dict{String, Any}(
                "name" => "demand_driven"
            ),
            "power" => 40000
        ),
        "TST_GBO_02" => Dict{String, Any}(
            "type" => "GasBoiler",
            "control_refs" => ["TST_BFT_01"],
            "production_refs" => [
                "TST_BUS_01"
            ],
            "strategy" => Dict{String, Any}(
                "name" => "demand_driven"
            ),
            "power" => 40000
        ),
        "TST_BUS_01" => Dict{String, Any}(
            "type" => "Bus",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "production_refs" => ["TST_DEM_01", "TST_BFT_01"],
            "input_priorities" => ["TST_GBO_01", "TST_BFT_01", "TST_GBO_02"]
        ),
        "TST_BFT_01" => Dict{String, Any}(
            "type" => "BufferTank",
            "control_refs" => [],
            "production_refs" => [
                "TST_BUS_01"
            ],
            "capacity" => 40000,
            "load" => 0
        ),
        "TST_DEM_01" => Dict{String, Any}(
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

function test_ooo_bus_output_priorities()
    systems_config = Dict{String, Any}(
        "TST_GRI_01" => Dict{String, Any}(
            "type" => "GridConnection",
            "medium" => "m_c_g_natgas",
            "control_refs" => [],
            "production_refs" => ["TST_GBO_01"],
            "is_source" => true,
        ),
        "TST_GBO_01" => Dict{String, Any}(
            "type" => "GasBoiler",
            "control_refs" => ["TST_BUS_01"],
            "production_refs" => ["TST_BUS_01"],
            "strategy" => Dict{String, Any}(
                "name" => "demand_driven",
            ),
            "power" => 10000
        ),
        "TST_BUS_01" => Dict{String, Any}(
            "type" => "Bus",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "production_refs" => ["TST_BUS_02", "TST_BUS_03"],
            "input_priorities" => ["TST_GBO_01"]
        ),
        "TST_BUS_02" => Dict{String, Any}(
            "type" => "Bus",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "production_refs" => ["TST_DEM_01"],
            "input_priorities" => ["TST_BUS_01"]
        ),
        "TST_BUS_03" => Dict{String, Any}(
            "type" => "Bus",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "production_refs" => ["TST_DEM_02"],
            "input_priorities" => ["TST_BUS_01"]
        ),
        "TST_DEM_01" => Dict{String, Any}(
            "type" => "Demand",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "production_refs" => [],
            "static_load" => 1000,
            "static_temperature" => 60,
            "scale" => 1
        ),
        "TST_DEM_02" => Dict{String, Any}(
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

@testset "project_loading_tests" begin
    @testset "load_from_dict" begin
        test_load_from_dict()
    end

    @testset "ooo_for_heat_pumps_wrong" begin
        test_ooo_for_heat_pumps_wrong()
    end

    @testset "ooo_bus_to_bus" begin
        test_ooo_bus_to_bus()
    end

    @testset "ooo_storage_loading_switch" begin
        test_ooo_storage_loading_switch()
    end

    @testset "ooo_bus_output_priorities" begin
        test_ooo_bus_output_priorities()
    end
end