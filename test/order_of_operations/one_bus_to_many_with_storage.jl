
using Debugger
using Test
using Resie
using Resie.EnergySystems

include("../test_util.jl")

function test_ooo_one_bus_to_many_with_storage()
    components_config = Dict{String,Any}(
        "TST_SRC_01" => Dict{String,Any}(
            "type" => "BoundedSupply",
            "medium" => "m_h_w_ht1",
            "output_refs" => ["TST_BUS_01"],
            "is_source" => true,
            "constant_power" => 400,
        ),
        "TST_BUS_01" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_h_w_ht1",
            "connections" => Dict{String,Any}(
                "input_order" => ["TST_SRC_01"],
                "output_order" => ["TST_BUS_02",
                                   "TST_BUS_03"],
            ),
        ),
        "TST_BUS_02" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_h_w_ht1",
            "connections" => Dict{String,Any}(
                "input_order" => ["TST_BUS_01",
                                  "TST_TES_01"],
                "output_order" => ["TST_DEM_01",
                                   "TST_TES_01"],
            ),
        ),
        "TST_BUS_03" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_h_w_ht1",
            "connections" => Dict{String,Any}(
                "input_order" => ["TST_BUS_01",
                                  "TST_TES_02"],
                "output_order" => ["TST_DEM_02",
                                   "TST_TES_02"],
            ),
        ),
        "TST_TES_01" => Dict{String,Any}(
            "type" => "Storage",
            "medium" => "m_h_w_ht1",
            "output_refs" => ["TST_BUS_02"],
            "capacity" => 1000,
            "load" => 500,
        ),
        "TST_TES_02" => Dict{String,Any}(
            "type" => "Storage",
            "medium" => "m_h_w_ht1",
            "output_refs" => ["TST_BUS_03"],
            "capacity" => 1000,
            "load" => 500,
        ),
        "TST_DEM_01" => Dict{String,Any}(
            "type" => "Demand",
            "medium" => "m_h_w_ht1",
            "output_refs" => [],
            "constant_demand" => 100.0,
            "constant_temperature" => 55.0,
        ),
        "TST_DEM_02" => Dict{String,Any}(
            "type" => "Demand",
            "medium" => "m_h_w_ht1",
            "output_refs" => [],
            "constant_demand" => 100.0,
            "constant_temperature" => 55.0,
        ),
    )

    expected = [("TST_DEM_02", EnergySystems.s_reset),
                ("TST_DEM_01", EnergySystems.s_reset),
                ("TST_BUS_02", EnergySystems.s_reset),
                ("TST_BUS_01", EnergySystems.s_reset),
                ("Proxy-TST_BUS_01|TST_BUS_03|TST_BUS_02", EnergySystems.s_reset),
                ("TST_BUS_03", EnergySystems.s_reset),
                ("TST_TES_02", EnergySystems.s_reset),
                ("TST_TES_01", EnergySystems.s_reset),
                ("TST_SRC_01", EnergySystems.s_reset),
                ("TST_DEM_02", EnergySystems.s_control),
                ("TST_DEM_01", EnergySystems.s_control),
                ("TST_BUS_02", EnergySystems.s_control),
                ("TST_BUS_01", EnergySystems.s_control),
                ("Proxy-TST_BUS_01|TST_BUS_03|TST_BUS_02", EnergySystems.s_control),
                ("TST_BUS_03", EnergySystems.s_control),
                ("TST_TES_02", EnergySystems.s_control),
                ("TST_TES_01", EnergySystems.s_control),
                ("TST_SRC_01", EnergySystems.s_control),
                ("TST_DEM_02", EnergySystems.s_process),
                ("TST_DEM_01", EnergySystems.s_process),
                ("TST_BUS_02", EnergySystems.s_process),
                ("TST_BUS_01", EnergySystems.s_process),
                ("Proxy-TST_BUS_01|TST_BUS_03|TST_BUS_02", EnergySystems.s_process),
                ("TST_BUS_03", EnergySystems.s_process),
                ("TST_SRC_01", EnergySystems.s_process),
                ("TST_TES_02", EnergySystems.s_process),
                ("TST_TES_01", EnergySystems.s_process),
                ("TST_TES_01", EnergySystems.s_load),
                ("TST_TES_02", EnergySystems.s_load),
                ("Proxy-TST_BUS_01|TST_BUS_03|TST_BUS_02", EnergySystems.s_distribute),
                ("TST_BUS_02", EnergySystems.s_distribute),
                ("TST_BUS_01", EnergySystems.s_distribute),
                ("TST_BUS_03", EnergySystems.s_distribute)]

    simulation_parameters = get_default_sim_params()

    components = Resie.load_components(components_config, simulation_parameters)
    ooo = Resie.calculate_order_of_operations(components)
    @test pwc_steps_astr(expected, ooo) == ""
end

@testset "ooo_one_bus_to_many_with_storage" begin
    test_ooo_one_bus_to_many_with_storage()
end
