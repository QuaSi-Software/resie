using Bran

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

        systems = Bran.load_systems(systems_config)
        @test length(keys(systems)) == 2
        @test typeof(systems["TST_BT_01"]) == Bran.EnergySystems.BufferTank
        @test systems["TST_BT_01"].sys_function == Bran.EnergySystems.sf_storage
        @test systems["TST_HP_01"].power == 20000
    end
end