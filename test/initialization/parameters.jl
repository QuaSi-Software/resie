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
            "temperature_from_global_file" => "temp key",
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

function test_validations_fail()
    components_config = Dict{String,Any}(
        "TST_GRI_01" => Dict{String,Any}(
            "type" => "FlexibleSupply",
            "output_refs" => ["TST_BUS_01"],
            "medium" => "m_h_w_ht1",
            "constant_temperature" => 70.0,
        ),
        "TST_BUS_01" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_h_w_ht1",
            "connections" => Dict{String,Any}(
                "input_order" => ["TST_STO_01",
                                  "TST_GRI_01"],
                "output_order" => ["TST_DEM_01",
                                   "TST_STO_01"],
                "energy_flow" => [[1, 0],
                                  [1, 1]],
            ),
        ),
        "TST_STO_01" => Dict{String,Any}(
            "type" => "Storage",
            "output_refs" => ["TST_BUS_01"],
            "medium" => "m_h_w_ht1",
            "capacity" => 10000,
            "initial_load" => 1.1,
        ),
        "TST_DEM_01" => Dict{String,Any}(
            "type" => "Demand",
            "output_refs" => [],
            "medium" => "m_h_w_ht1",
            "constant_temperature" => 50.0,
            "constant_demand" => 1000.0,
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
    @assert occursin("Can't construct component TST_STO_01", msg)
end

@testset "validations_fail" begin
    test_validations_fail()
end

function test_required_with_conditional()
    components_config = Dict{String,Any}(
        "TST_GRI_01" => Dict{String,Any}(
            "type" => "FlexibleSupply",
            "output_refs" => ["TST_BUS_01"],
            "medium" => "m_e_ac_230v",
        ),
        "TST_BUS_01" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_e_ac_230v",
            "connections" => Dict{String,Any}(
                "input_order" => ["TST_BAT_01",
                                  "TST_GRI_01"],
                "output_order" => ["TST_DEM_01",
                                   "TST_BAT_01"],
                "energy_flow" => [[1, 0],
                                  [1, 1]],
            ),
        ),
        "TST_BAT_01" => Dict{String,Any}(
            "type" => "Battery",
            "output_refs" => ["TST_BUS_01"],
            "capacity" => 10000,
            "initial_load" => 0.5,
            "model_type" => "simplified",
        ),
        "TST_DEM_01" => Dict{String,Any}(
            "type" => "Demand",
            "output_refs" => [],
            "medium" => "m_e_ac_230v",
            "constant_demand" => 1000.0,
        ),
    )
    simulation_params = get_default_sim_params()

    # first attempt should fail because required param charge_efficiency is not given for
    # simplified model
    construction_errored = false
    msg = ""
    try
        _ = Resie.load_components(components_config, simulation_params)
    catch e
        construction_errored = true
        msg = sprint(showerror, e)
    end
    @assert construction_errored
    @assert occursin("Can't construct component TST_BAT_01", msg)

    # second attempt should succeed despite required parameter V_n_bat not being given
    # because the conditional on model type != simplified removes that requirement
    components_config["TST_BAT_01"]["charge_efficiency"] = 0.99
    components_config["TST_BAT_01"]["discharge_efficiency"] = 0.99
    construction_errored = false
    try
        _ = Resie.load_components(components_config, simulation_params)
    catch e
        showerror(stdout, e)
        print(stacktrace(catch_backtrace()))
        construction_errored = true
    end
    @assert !construction_errored
end

@testset "required_with_conditional" begin
    test_required_with_conditional()
end
