using Debugger
using Test
using Resie
using Resie.EnergySystems

function test_load_from_dict()
    systems_config = Dict{String,Any}(
        "TST_BT_01" => Dict{String,Any}(
            "type" => "BufferTank",
            "control_refs" => [],
            "production_refs" => [],
            "capacity" => 40000,
            "load" => 20000,
            "strategy" => Dict{String,Any}(
                "name" => "Default"
            )
        ),
        "TST_HP_01" => Dict{String,Any}(
            "type" => "HeatPump",
            "control_refs" => ["TST_BT_01"],
            "production_refs" => [],
            "strategy" => Dict{String,Any}(
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

@testset "project_loading_tests" begin
    @testset "load_from_dict" begin
        test_load_from_dict()
    end

    include("order_of_operations/bus_output_priorities.jl")
    include("order_of_operations/bus_to_bus.jl")
    include("order_of_operations/heat_pumps_wrong.jl")
    include("order_of_operations/project_loading_switch.jl")
end