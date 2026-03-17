using Debugger
using Test
using Resie
using Resie.EnergySystems

include("../test_util.jl")

function test_load_from_dict()
    components_config = Dict{String,Any}(
        "TST_GRI_01" => Dict{String,Any}(
            "type" => "GridInput",
            "output_refs" => ["TST_HP_01"],
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
            "cop_function" => "const:3.0",
            "power_losses_factor" => 1.0,
            "heat_losses_factor" => 1.0,
        ),
        "TST_BT_01" => Dict{String,Any}(
            "type" => "BufferTank",
            "output_refs" => ["TST_DEM_01"],
            "medium" => "m_h_w_ht1",
            "model_type" => "ideally_stratified",
            "capacity" => 40000,
            "initial_load" => 0.5,
        ),
        "TST_DEM_01" => Dict{String,Any}(
            "type" => "Demand",
            "output_refs" => [],
            "medium" => "m_h_w_ht1",
            "constant_demand" => 5000,
            "constant_temperature" => 60,
        ),
    )
    simulation_params = get_default_sim_params()
    components = Resie.load_components(components_config, simulation_params)
    @test length(keys(components)) == 5
    @test typeof(components["TST_BT_01"]) == Resie.EnergySystems.BufferTank
    @test components["TST_BT_01"].sys_function == Resie.EnergySystems.sf_storage
    @test components["TST_HP_01"].design_power_th == 20000
end

@testset "load_from_dict" begin
    test_load_from_dict()
end
