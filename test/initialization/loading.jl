using Debugger
using Test
using Resie
using Resie.EnergySystems

function test_load_from_dict()
    components_config = Dict{String,Any}(
        "TST_GRI_01" => Dict{String,Any}(
            "type" => "GridConnection",
            "output_refs" => ["TST_HP_01"],
            "is_source" => true,
            "medium" => "m_h_w_lt1",
        ),
        "TST_GRI_02" => Dict{String,Any}(
            "type" => "GridConnection",
            "output_refs" => ["TST_HP_01"],
            "is_source" => true,
            "medium" => "m_e_ac_230v",
        ),
        "TST_HP_01" => Dict{String,Any}(
            "type" => "HeatPump",
            "output_refs" => ["TST_BT_01"],
            "control_modules" => [Dict{String,Any}(
                                      "name" => "storage_driven",
                                      "high_threshold" => 0.5,
                                      "low_threshold" => 0.1,
                                      "storage_uac" => "TST_BT_01",
                                  )],
            "power_th" => 20000,
            "constant_cop" => 3.0,
        ),
        "TST_BT_01" => Dict{String,Any}(
            "type" => "BufferTank",
            "output_refs" => ["TST_DEM_01"],
            "medium" => "m_h_w_ht1",
            "capacity" => 40000,
            "load" => 20000,
        ),
        "TST_DEM_01" => Dict{String,Any}(
            "type" => "Demand",
            "output_refs" => [],
            "medium" => "m_h_w_ht1",
            "constant_demand" => 5000,
            "constant_temperature" => 60,
        ),
    )
    simulation_params = Dict{String,Any}(
        "time" => 0,
        "time_step_seconds" => 900,
        "epsilon" => 1e-9,
    )
    components = Resie.load_components(components_config, simulation_params)
    @test length(keys(components)) == 5
    @test typeof(components["TST_BT_01"]) == Resie.EnergySystems.BufferTank
    @test components["TST_BT_01"].sys_function == Resie.EnergySystems.sf_storage
    @test components["TST_HP_01"].power_th == 20000
end

@testset "load_from_dict" begin
    test_load_from_dict()
end

function test_load_custom_medium_categories()
    components_config = Dict{String,Any}(
        "TST_GRI_01" => Dict{String,Any}(
            "type" => "GridConnection",
            "output_refs" => ["TST_ELY_01"],
            "is_source" => true,
            "medium" => "m_e_dc_1000v",
        ),
        "TST_GRI_02" => Dict{String,Any}(
            "type" => "GridConnection",
            "output_refs" => ["TST_HP_01"],
            "is_source" => true,
            "medium" => "m_e_dc_1000v",
        ),
        "TST_ELY_01" => Dict{String,Any}(
            "type" => "Electrolyser",
            "output_refs" => ["TST_HP_01",
                              "TST_GRO_02",
                              "TST_GRO_03"],
            "power_el" => 1000,
            "m_el_in" => "m_e_dc_1000v",
            "m_heat_ht_out" => "m_h_w_55c",
            "m_h2_out" => "m_c_g_h2-pure",
            "m_o2_out" => "m_c_g_o2-impure",
            "heat_lt_is_usable" => false,
        ),
        "TST_HP_01" => Dict{String,Any}(
            "type" => "HeatPump",
            "output_refs" => ["TST_GRO_01"],
            "m_el_in" => "m_e_dc_1000v",
            "m_heat_in" => "m_h_w_55c",
            "m_heat_out" => "m_h_w_85c",
            "power_th" => 12000,
            "constant_cop" => 3.0,
        ),
        "TST_GRO_01" => Dict{String,Any}(
            "type" => "GridConnection",
            "output_refs" => [],
            "is_source" => false,
            "medium" => "m_h_w_85c",
        ),
        "TST_GRO_02" => Dict{String,Any}(
            "type" => "GridConnection",
            "output_refs" => [],
            "is_source" => false,
            "medium" => "m_c_g_h2-pure",
        ),
        "TST_GRO_03" => Dict{String,Any}(
            "type" => "GridConnection",
            "output_refs" => [],
            "is_source" => false,
            "medium" => "m_c_g_o2-impure",
        ),
    )
    simulation_params = Dict{String,Any}(
        "time" => 0,
        "time_step_seconds" => 900,
        "epsilon" => 1e-9,
    )
    components = Resie.load_components(components_config, simulation_params)
    electrolyser = components["TST_ELY_01"]
    heat_pump = components["TST_HP_01"]

    @test electrolyser.m_el_in == Symbol("m_e_dc_1000v")
    @test electrolyser.m_heat_ht_out == Symbol("m_h_w_55c")
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
    @test electrolyser.output_interfaces[electrolyser.m_heat_ht_out].target.uac == "TST_HP_01"
end

@testset "load_custom_medium_categories" begin
    test_load_custom_medium_categories()
end
