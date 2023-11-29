using Debugger
using Test
using Resie
using Resie.EnergySystems
using Resie.Profiles

EnergySystems.set_timestep(900)

function test_many_to_one()
    components_config = Dict{String,Any}(
        "TST_SRC_01" => Dict{String,Any}(
            "type" => "BoundedSupply",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "output_refs" => ["TST_BUS_TH_01"],
            "static_power" => 4000,
            "static_temperature" => 55,
        ),
        "TST_DEM_01" => Dict{String,Any}(
            "type" => "Demand",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "output_refs" => [],
            "static_load" => 500,
            "static_temperature" => 55,
        ),
        "TST_DEM_02" => Dict{String,Any}(
            "type" => "Demand",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "output_refs" => [],
            "static_load" => 500,
            "static_temperature" => 55,
        ),
        "TST_BUS_TH_01" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "output_refs" => [
                "TST_DEM_01",
                "TST_DEM_02",
            ],
            "connection_matrix" => Dict{String,Any}(
                "input_order" => [
                    "TST_SRC_01",
                ],
                "output_order" => [
                    "TST_DEM_01",
                    "TST_DEM_02",
                ],
                "storage_loading" => [
                    [1, 1],
                ],
            ),
        ),
    )
    components = Resie.load_components(components_config)
    source = components["TST_SRC_01"]
    demand_1 = components["TST_DEM_01"]
    demand_2 = components["TST_DEM_02"]
    bus = components["TST_BUS_TH_01"]

    simulation_parameters = Dict{String,Any}(
        "time_step_seconds" => 900,
        "time" => 0,
    )

    EnergySystems.reset(demand_1)
    EnergySystems.reset(demand_2)
    EnergySystems.reset(source)
    EnergySystems.reset(bus)

    @test demand_1.load == 0.0
    @test demand_1.temperature === nothing

    @test demand_2.load == 0.0
    @test demand_2.temperature === nothing

    EnergySystems.control(demand_1, components, simulation_parameters)
    EnergySystems.control(demand_2, components, simulation_parameters)
    EnergySystems.control(source, components, simulation_parameters)
    EnergySystems.control(bus, components, simulation_parameters)

    @test demand_1.load == 500.0
    @test demand_1.temperature == 55.0

    @test demand_2.load == 500.0
    @test demand_2.temperature == 55.0

    EnergySystems.process(demand_1, simulation_parameters)

    @test demand_1.input_interfaces[demand_1.medium].balance == -500.0
    @test demand_1.input_interfaces[demand_1.medium].temperature == 55.0

    EnergySystems.process(demand_2, simulation_parameters)

    @test demand_2.input_interfaces[demand_2.medium].balance == -500.0
    @test demand_2.input_interfaces[demand_2.medium].temperature == 55.0

    EnergySystems.process(bus, simulation_parameters)
    EnergySystems.process(source, simulation_parameters)

    @test source.output_interfaces[source.medium].balance == 1000.0
    @test source.output_interfaces[source.medium].temperature == 55.0

    blnc = EnergySystems.balance(bus)
    @test blnc == 0.0

    EnergySystems.distribute!(bus)

    @test demand_1.input_interfaces[demand_1.medium].balance == 0.0
    @test demand_1.input_interfaces[demand_1.medium].temperature === 55.0

    @test demand_2.input_interfaces[demand_2.medium].balance == 0.0
    @test demand_2.input_interfaces[demand_2.medium].temperature === 55.0

    @test source.output_interfaces[source.medium].balance == 0.0
    @test source.output_interfaces[source.medium].temperature == 55.0

end

@testset "many_to_one" begin
    test_many_to_one()
end