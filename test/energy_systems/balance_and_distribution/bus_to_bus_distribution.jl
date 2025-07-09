using Debugger
using Test
using Resie
using Resie.EnergySystems
using Resie.Profiles

include("../../test_util.jl")

function test_short_chain_distribution()
    components_config = Dict{String,Any}(
        "TST_SRC_01" => Dict{String,Any}(
            "type" => "BoundedSupply",
            "medium" => "m_e_ac_230v",
            "output_refs" => ["TST_BUS_01"],
            "constant_power" => 8000,
        ),
        "TST_DEM_01" => Dict{String,Any}(
            "type" => "Demand",
            "medium" => "m_e_ac_230v",
            "output_refs" => [],
            "constant_demand" => 4000,
        ),
        "TST_DEM_02" => Dict{String,Any}(
            "type" => "Demand",
            "medium" => "m_e_ac_230v",
            "output_refs" => [],
            "constant_demand" => 2000,
        ),
        "TST_BUS_01" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_e_ac_230v",
            "connections" => Dict{String,Any}(
                "input_order" => ["TST_SRC_01"],
                "output_order" => ["TST_BUS_02"],
            ),
        ),
        "TST_BUS_02" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_e_ac_230v",
            "connections" => Dict{String,Any}(
                "input_order" => ["TST_BUS_01"],
                "output_order" => ["TST_DEM_01",
                                   "TST_DEM_02"],
            ),
        ),
    )

    simulation_parameters = get_default_sim_params()

    components = Resie.load_components(components_config, simulation_parameters)
    proxy = components["Proxy-TST_BUS_01|TST_BUS_02"]

    @test length(proxy.input_interfaces) == 1
    @test length(proxy.output_interfaces) == 2
    @test proxy.output_interfaces[2].target.uac == "TST_DEM_02"

    EnergySystems.reset(components["TST_BUS_01"])
    EnergySystems.reset(components["TST_BUS_02"])
    EnergySystems.reset(proxy)
    EnergySystems.reset(components["TST_SRC_01"])
    EnergySystems.reset(components["TST_DEM_01"])
    EnergySystems.reset(components["TST_DEM_02"])

    EnergySystems.control(components["TST_BUS_01"], components, simulation_parameters)
    EnergySystems.control(components["TST_BUS_02"], components, simulation_parameters)
    EnergySystems.control(proxy, components, simulation_parameters)
    EnergySystems.control(components["TST_SRC_01"], components, simulation_parameters)
    EnergySystems.control(components["TST_DEM_01"], components, simulation_parameters)
    EnergySystems.control(components["TST_DEM_02"], components, simulation_parameters)

    EnergySystems.process(components["TST_DEM_01"], simulation_parameters)
    EnergySystems.process(components["TST_DEM_02"], simulation_parameters)
    EnergySystems.process(components["TST_SRC_01"], simulation_parameters)
    EnergySystems.process(components["TST_BUS_01"], simulation_parameters)
    EnergySystems.process(components["TST_BUS_02"], simulation_parameters)
    EnergySystems.process(proxy, simulation_parameters)

    EnergySystems.distribute!(proxy)
    EnergySystems.distribute!(components["TST_BUS_01"])
    EnergySystems.distribute!(components["TST_BUS_02"])

    bus_1 = components["TST_BUS_01"]
    @test bus_1.output_interfaces[1].sum_abs_change == 1500.0 * 2

    bus_2 = components["TST_BUS_02"]
    @test bus_2.input_interfaces[1].sum_abs_change == 1500.0 * 2
end

@testset "short_chain_distribution" begin
    test_short_chain_distribution()
end

function test_long_chain_distribution()
    components_config = Dict{String,Any}(
        "TST_SRC_01" => Dict{String,Any}(
            "type" => "BoundedSupply",
            "medium" => "m_e_ac_230v",
            "output_refs" => ["TST_BUS_01"],
            "constant_power" => 8000,
        ),
        "TST_DEM_01" => Dict{String,Any}(
            "type" => "Demand",
            "medium" => "m_e_ac_230v",
            "output_refs" => [],
            "constant_demand" => 4000,
        ),
        "TST_DEM_02" => Dict{String,Any}(
            "type" => "Demand",
            "medium" => "m_e_ac_230v",
            "output_refs" => [],
            "constant_demand" => 2000,
        ),
        "TST_BUS_01" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_e_ac_230v",
            "connections" => Dict{String,Any}(
                "input_order" => ["TST_SRC_01"],
                "output_order" => ["TST_BUS_02"],
            ),
        ),
        "TST_BUS_02" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_e_ac_230v",
            "connections" => Dict{String,Any}(
                "input_order" => ["TST_BUS_01"],
                "output_order" => ["TST_BUS_03",
                                   "TST_BUS_04"],
            ),
        ),
        "TST_BUS_03" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_e_ac_230v",
            "connections" => Dict{String,Any}(
                "input_order" => ["TST_BUS_02"],
                "output_order" => ["TST_DEM_01"],
            ),
        ),
        "TST_BUS_04" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_e_ac_230v",
            "connections" => Dict{String,Any}(
                "input_order" => ["TST_BUS_02"],
                "output_order" => ["TST_DEM_02"],
            ),
        ),
    )

    simulation_parameters = get_default_sim_params()

    components = Resie.load_components(components_config, simulation_parameters)
    proxy = components["Proxy-TST_BUS_01|TST_BUS_02|TST_BUS_04|TST_BUS_03"]

    @test length(proxy.input_interfaces) == 1
    @test length(proxy.output_interfaces) == 2
    @test proxy.output_interfaces[2].target.uac == "TST_DEM_02"

    EnergySystems.reset(components["TST_BUS_01"])
    EnergySystems.reset(components["TST_BUS_02"])
    EnergySystems.reset(components["TST_BUS_03"])
    EnergySystems.reset(components["TST_BUS_04"])
    EnergySystems.reset(proxy)
    EnergySystems.reset(components["TST_SRC_01"])
    EnergySystems.reset(components["TST_DEM_01"])
    EnergySystems.reset(components["TST_DEM_02"])

    EnergySystems.control(components["TST_BUS_01"], components, simulation_parameters)
    EnergySystems.control(components["TST_BUS_02"], components, simulation_parameters)
    EnergySystems.control(components["TST_BUS_03"], components, simulation_parameters)
    EnergySystems.control(components["TST_BUS_04"], components, simulation_parameters)
    EnergySystems.control(proxy, components, simulation_parameters)
    EnergySystems.control(components["TST_SRC_01"], components, simulation_parameters)
    EnergySystems.control(components["TST_DEM_01"], components, simulation_parameters)
    EnergySystems.control(components["TST_DEM_02"], components, simulation_parameters)

    EnergySystems.process(components["TST_DEM_01"], simulation_parameters)
    EnergySystems.process(components["TST_DEM_02"], simulation_parameters)
    EnergySystems.process(components["TST_SRC_01"], simulation_parameters)
    EnergySystems.process(components["TST_BUS_01"], simulation_parameters)
    EnergySystems.process(components["TST_BUS_02"], simulation_parameters)
    EnergySystems.process(components["TST_BUS_03"], simulation_parameters)
    EnergySystems.process(components["TST_BUS_04"], simulation_parameters)
    EnergySystems.process(proxy, simulation_parameters)

    EnergySystems.distribute!(proxy)
    EnergySystems.distribute!(components["TST_BUS_01"])
    EnergySystems.distribute!(components["TST_BUS_02"])
    EnergySystems.distribute!(components["TST_BUS_03"])
    EnergySystems.distribute!(components["TST_BUS_04"])

    bus_1 = components["TST_BUS_01"]
    @test bus_1.output_interfaces[1].sum_abs_change == 1500.0 * 2

    bus_2 = components["TST_BUS_02"]
    @test bus_2.output_interfaces[1].sum_abs_change == 1000.0 * 2
    @test bus_2.output_interfaces[2].sum_abs_change == 500.0 * 2

    bus_3 = components["TST_BUS_03"]
    @test bus_3.input_interfaces[1].sum_abs_change == 1000.0 * 2

    bus_4 = components["TST_BUS_04"]
    @test bus_4.input_interfaces[1].sum_abs_change == 500.0 * 2
end

@testset "long_chain_distribution" begin
    test_long_chain_distribution()
end
