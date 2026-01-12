using Debugger
using Test
using Resie
using Resie.EnergySystems

include("../test_util.jl")

function test_mutex_parameters()
    components_config = Dict{String,Any}(
        "TST_GRI_01" => Dict{String,Any}(
            "type" => "GridInput",
            "output_refs" => ["TST_DEM_01"],
            "medium" => "m_h_w_lt1",
            "temperature_from_global_file" => true,
            "constant_temperature" => 50.0,
        ),
        "TST_DEM_01" => Dict{String,Any}(
            "type" => "Demand",
            "output_refs" => [],
            "medium" => "m_h_w_lt1",
            "constant_demand" => 500.0,
        ),
    )
    simulation_params = get_default_sim_params()

    construction_errored = false
    msg = ""
    try
        _ = Resie.load_components(components_config, simulation_params)
    catch e
        construction_errored = true
        msg = sprint(showerror, e)
    end
    @assert construction_errored
    @assert occursin("Can't construct component TST_GRI_01", msg)
end

@testset "mutex_parameters" begin
    test_mutex_parameters()
end

function test_parameter_with_options()
    components_config = Dict{String,Any}(
        "TST_GRI_01" => Dict{String,Any}(
            "type" => "GridInput",
            "output_refs" => ["TST_HP_01"],
            "medium" => "m_h_w_lt1",
            "constant_temperature" => 20.0,
        ),
        "TST_GRI_02" => Dict{String,Any}(
            "type" => "GridInput",
            "output_refs" => ["TST_HP_01"],
            "medium" => "m_e_ac_230v",
        ),
        "TST_HP_01" => Dict{String,Any}(
            "type" => "HeatPump",
            "output_refs" => ["TST_DEM_01"],
            "model_type" => "not_valid",
            "power_th" => 10000,
        ),
        "TST_DEM_01" => Dict{String,Any}(
            "type" => "Demand",
            "output_refs" => [],
            "medium" => "m_h_w_ht1",
            "constant_temperature" => 50.0,
            "constant_demand" => 2000.0,
        ),
    )
    simulation_params = get_default_sim_params()

    construction_errored = false
    msg = ""
    try
        _ = Resie.load_components(components_config, simulation_params)
    catch e
        construction_errored = true
        msg = sprint(showerror, e)
    end
    @assert construction_errored
    @assert occursin("Can't construct component TST_HP_01", msg)
end

@testset "parameter_with_options" begin
    test_parameter_with_options()
end
