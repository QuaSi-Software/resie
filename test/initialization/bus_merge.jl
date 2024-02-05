using Test
using Resie
using Resie.EnergySystems

function energy_system_simple()::Dict{String,Any}
    return Dict{String,Any}(
        "TST_SRC_01" => Dict{String,Any}(
            "type" => "BoundedSupply",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "output_refs" => ["TST_BUS_01"],
            "constant_supply" => 1500,
            "constant_temperature" => 65
        ),
        "TST_SRC_02" => Dict{String,Any}(
            "type" => "BoundedSupply",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "output_refs" => ["TST_BUS_02"],
            "constant_supply" => 500,
            "constant_temperature" => 40
        ),
        "TST_BUS_01" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "connections" => Dict{String, Any}(
                "input_order" => [
                    "TST_SRC_01"
                ],
                "output_order" => [
                    "TST_DEM_01",
                    "TST_BUS_02"
                ],
                "energy_flow" => [
                    [1,1],
                ]
            )
        ),
        "TST_BUS_02" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "connections" => Dict{String, Any}(
                "input_order" => [
                    "TST_BUS_01",
                    "TST_SRC_02",
                ],
                "output_order" => [
                    "TST_DEM_02",
                ],
                "energy_flow" => [
                    [1],
                    [1],
                ]
            )
        ),
        "TST_DEM_01" => Dict{String,Any}(
            "type" => "Demand",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "output_refs" => [],
            "constant_demand" => 1000,
            "constant_temperature" => 65
        ),
        "TST_DEM_02" => Dict{String,Any}(
            "type" => "Demand",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "output_refs" => [],
            "constant_demand" => 1000,
            "constant_temperature" => 40
        ),
    )
end

function test_shallow_copy()
    components_config = energy_system_simple()

    simulation_parameters = Dict{String,Any}(
        "time_step_seconds" => 900,
        "time" => 0,
        "epsilon" => 1e-9
    )

    components = Resie.load_components(components_config, simulation_parameters)
    bus = components["TST_BUS_01"]
    new_bus = Bus(bus, true)

    bus.uac = "new_uac"
    @test new_bus.uac == "TST_BUS_01"

    bus.controller.parameter["foo"] = "bar"
    @test new_bus.controller.parameter["foo"] == "bar"
end

@testset "shallow_copy" begin
    test_shallow_copy()
end

function test_deep_copy()
    components_config = energy_system_simple()

    simulation_parameters = Dict{String,Any}(
        "time_step_seconds" => 900,
        "time" => 0,
        "epsilon" => 1e-9
    )

    components = Resie.load_components(components_config, simulation_parameters)
    bus = components["TST_BUS_01"]
    new_bus = EnergySystems.deepcopy(bus)

    bus.uac = "new_uac"
    @test new_bus.uac == "TST_BUS_01"

    bus.controller.parameter["foo"] = "bar"
    @test !haskey(new_bus.controller.parameter, "foo")

    bus.input_interfaces[1].source.uac = "new_uac"
    @test new_bus.input_interfaces[1].source.uac == "new_uac"
end

@testset "deep_copy" begin
    test_deep_copy()
end
