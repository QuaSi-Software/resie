using Debugger
using Test
using Resie
using Resie.EnergySystems
using Resie.Profiles

include("../../test_util.jl")

function test_many_to_one()
    components_config = Dict{String,Any}(
        "TST_SRC_01" => Dict{String,Any}(
            "type" => "BoundedSupply",
            "medium" => "m_h_w_ht1",
            "output_refs" => ["TST_BUS_TH_01"],
            "constant_power" => 3200,
            "constant_temperature" => 55,
        ),
        "TST_SRC_02" => Dict{String,Any}(
            "type" => "BoundedSupply",
            "medium" => "m_h_w_ht1",
            "output_refs" => ["TST_BUS_TH_01"],
            "constant_power" => 2000,
            "constant_temperature" => 55,
        ),
        "TST_DEM_01" => Dict{String,Any}(
            "type" => "Demand",
            "medium" => "m_h_w_ht1",
            "output_refs" => [],
            "constant_demand" => 4000,
            "constant_temperature" => 55,
        ),
        "TST_BUS_TH_01" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_h_w_ht1",
            "connections" => Dict{String,Any}(
                "input_order" => ["TST_SRC_01",
                                  "TST_SRC_02"],
                "output_order" => ["TST_DEM_01"],
                "energy_flow" => [[1],
                                  [1]],
            ),
        ),
    )

    simulation_parameters = get_default_sim_params()

    components = Resie.load_components(components_config, simulation_parameters)
    source_1 = components["TST_SRC_01"]
    source_2 = components["TST_SRC_02"]
    demand = components["TST_DEM_01"]
    bus = components["TST_BUS_TH_01"]

    # first timestep, everything works out
    EnergySystems.reset(demand)
    EnergySystems.reset(source_1)
    EnergySystems.reset(source_2)
    EnergySystems.reset(bus)

    @test demand.demand == 0.0
    @test demand.temperature === nothing

    EnergySystems.control(demand, components, simulation_parameters)
    EnergySystems.control(source_1, components, simulation_parameters)
    EnergySystems.control(source_2, components, simulation_parameters)
    EnergySystems.control(bus, components, simulation_parameters)

    @test demand.demand == 1000.0
    @test demand.temperature == 55.0

    EnergySystems.process(demand, simulation_parameters)

    @test demand.input_interfaces[demand.medium].balance == -1000.0
    @test demand.input_interfaces[demand.medium].temperature_min == 55.0

    EnergySystems.process(bus, simulation_parameters)
    EnergySystems.process(source_1, simulation_parameters)

    @test source_1.output_interfaces[source_1.medium].balance == 800.0
    @test source_1.output_interfaces[source_1.medium].temperature_max == 55.0

    EnergySystems.process(source_2, simulation_parameters)

    @test source_2.output_interfaces[source_2.medium].balance == 200.0
    @test source_2.output_interfaces[source_2.medium].temperature_max == 55.0

    blnc = EnergySystems.balance(bus)
    @test blnc == 0.0

    EnergySystems.distribute!(bus)

    @test demand.input_interfaces[demand.medium].balance == 0.0
    @test demand.input_interfaces[demand.medium].temperature_min === 55.0

    @test source_1.output_interfaces[source_1.medium].balance == 0.0
    @test source_1.output_interfaces[source_1.medium].temperature_max == 55.0

    @test source_2.output_interfaces[source_2.medium].balance == 0.0
    @test source_2.output_interfaces[source_2.medium].temperature_max == 55.0

    # second timestep, source 1 is too cold to supply demand, source 2 can't fully
    # supply demand
    EnergySystems.reset(demand)
    EnergySystems.reset(source_1)
    EnergySystems.reset(source_2)
    EnergySystems.reset(bus)

    source_1.constant_temperature = 50.0

    EnergySystems.control(demand, components, simulation_parameters)
    EnergySystems.control(source_1, components, simulation_parameters)
    EnergySystems.control(source_2, components, simulation_parameters)
    EnergySystems.control(bus, components, simulation_parameters)

    EnergySystems.process(demand, simulation_parameters)
    EnergySystems.process(bus, simulation_parameters)

    EnergySystems.process(source_1, simulation_parameters)

    @test source_1.output_interfaces[source_1.medium].balance == 0.0
    @test source_1.output_interfaces[source_1.medium].temperature_max == 50.0

    EnergySystems.process(source_2, simulation_parameters)

    @test source_2.output_interfaces[source_2.medium].balance == 500.0
    @test source_2.output_interfaces[source_2.medium].temperature_max == 55.0

    blnc = EnergySystems.balance(bus)
    @test blnc == -500.0

    EnergySystems.distribute!(bus)
end

@testset "many_to_one" begin
    test_many_to_one()
end
