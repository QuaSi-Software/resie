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

function test_load_custom_medium_categories()
    systems_config = Dict{String,Any}(
        "TST_ELY_01" => Dict{String,Any}(
            "type" => "Electrolyser",
            "control_refs" => [],
            "production_refs" => ["TST_HP_01"],
            "power" => 1000,
            "m_el_in" => "m_e_dc_1000v",
            "m_heat_out" => "m_h_w_55c",
            "m_h2_out" => "m_c_g_h2-pure",
            "m_o2_out" => "m_c_g_o2-impure",
            "strategy" => Dict{String,Any}(
                "name" => "Default"
            )
        ),
        "TST_HP_01" => Dict{String,Any}(
            "type" => "HeatPump",
            "control_refs" => [],
            "production_refs" => [],
            "strategy" => Dict{String,Any}(
                "name" => "demand_driven",
            ),
            "m_el_in" => "m_e_dc_1000v",
            "m_heat_in" => "m_h_w_55c",
            "m_heat_out" => "m_h_w_85c",
            "power" => 12000,
            "fixed_cop" => 3.0
        )
    )
    systems = Resie.load_systems(systems_config)
    electrolyser = systems["TST_ELY_01"]
    heat_pump = systems["TST_HP_01"]

    @test electrolyser.m_el_in == Symbol("m_e_dc_1000v")
    @test electrolyser.m_heat_out == Symbol("m_h_w_55c")
    @test electrolyser.m_h2_out == Symbol("m_c_g_h2-pure")
    @test electrolyser.m_o2_out == Symbol("m_c_g_o2-impure")

    @test heat_pump.m_el_in == Symbol("m_e_dc_1000v")
    @test heat_pump.m_heat_in == Symbol("m_h_w_55c")
    @test heat_pump.m_heat_out == Symbol("m_h_w_85c")

    @test Symbol("m_e_dc_1000v") in EnergySystems.medium_categories
    @test Symbol("m_h_w_55c") in EnergySystems.medium_categories
    @test Symbol("m_h_w_85c") in EnergySystems.medium_categories
    @test Symbol("m_c_g_h2-pure") in EnergySystems.medium_categories
    @test Symbol("m_c_g_o2-impure") in EnergySystems.medium_categories

    # check if iterface of user-defined medium is built up correctly
    @test electrolyser.output_interfaces[electrolyser.m_heat_out].target.uac == "TST_HP_01"
end

@testset "project_loading_tests" begin
    @testset "load_from_dict" begin
        test_load_from_dict()
    end

    @testset "load_custom_medium_categories" begin
        test_load_custom_medium_categories()
    end

    include("initialization/bus.jl")

    include("order_of_operations/base_order.jl")
    include("order_of_operations/load_order_of_operation.jl")
    include("order_of_operations/bus_output_priorities.jl")
    include("order_of_operations/bus_to_bus.jl")
    include("order_of_operations/heat_pumps_wrong.jl")
    include("order_of_operations/reorderings.jl")
    include("order_of_operations/storage_loading_switch.jl")
    include("order_of_operations/transformer_chains.jl")
end