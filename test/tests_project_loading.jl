using Resie
using Resie.EnergySystems

@testset "project_loading_tests" begin

    @testset "load_from_dict" begin
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
                "cop" => 3.0
            ),
        )

        systems = Resie.load_systems(systems_config)
        @test length(keys(systems)) == 2
        @test typeof(systems["TST_BT_01"]) == Resie.EnergySystems.BufferTank
        @test systems["TST_BT_01"].sys_function == Resie.EnergySystems.sf_storage
        @test systems["TST_HP_01"].power == 20000
    end

    @testset "test_ooo_for_heat_pumps_wrong" begin
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
                "cop" => 3.0
            ),
        )

        expected = [
            ("TST_DEM_01", EnergySystems.s_reset),
            ("TST_HP_01", EnergySystems.s_reset),
            ("TST_SRC_01", EnergySystems.s_reset),
            ("TST_GRI_01", EnergySystems.s_reset),
            ("TST_DEM_01", EnergySystems.s_control),
            ("TST_DEM_01", EnergySystems.s_produce),
            ("TST_HP_01", EnergySystems.s_control),
            ("TST_HP_01", EnergySystems.s_produce),
            ("TST_SRC_01", EnergySystems.s_control),
            ("TST_SRC_01", EnergySystems.s_produce),
            ("TST_GRI_01", EnergySystems.s_control),
            ("TST_GRI_01", EnergySystems.s_produce),
        ]

        systems = Resie.load_systems(systems_config)
        ooo = Resie.order_of_operations(systems)
        @test all(ooo .== expected)
    end
end