using Debugger
using Test
using Resie
using Resie.EnergySystems
using Resie.Profiles

EnergySystems.set_timestep(900)

function test_many_plus_storage_to_one()
    components_config = Dict{String,Any}(
        "TST_SRC_01" => Dict{String,Any}(
            "type" => "BoundedSupply",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "output_refs" => ["TST_BUS_TH_01"],
            "constant_power" => 2000,
            "constant_temperature" => 55,
        ),
        "TST_SRC_02" => Dict{String,Any}(
            "type" => "BoundedSupply",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "output_refs" => ["TST_BUS_TH_01"],
            "constant_power" => 2000,
            "constant_temperature" => 55,
        ),
        "TST_DEM_01" => Dict{String,Any}(
            "type" => "Demand",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "output_refs" => [],
            "constant_demand" => 4000,
            "constant_temperature" => 55,
        ),
        "TST_TES_01" => Dict{String,Any}(
            "type" => "Storage",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "output_refs" => ["TST_BUS_TH_01"],
            "capacity" => 10000,
            "load" => 5000,
        ),
        "TST_BUS_TH_01" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "connections" => Dict{String,Any}(
                "input_order" => [
                    "TST_SRC_01",
                    "TST_SRC_02",
                    "TST_TES_01",
                ],
                "output_order" => [
                    "TST_DEM_01",
                    "TST_TES_01",
                ],
                "energy_flow" => [
                    [1,1],
                    [1,1],
                    [1,0],
                ],
            ),
        ),
    )

    simulation_parameters = Dict{String,Any}(
        "time_step_seconds" => 900,
        "time" => 0,
    )

    components = Resie.load_components(components_config, simulation_parameters)
    source_1 = components["TST_SRC_01"]
    source_2 = components["TST_SRC_02"]
    demand = components["TST_DEM_01"]
    bus = components["TST_BUS_TH_01"]
    storage = components["TST_TES_01"]

    # time step 1: both sources can cover the demand, no storage needed

    EnergySystems.reset(demand)
    EnergySystems.reset(bus)
    EnergySystems.reset(storage)
    EnergySystems.reset(source_1)
    EnergySystems.reset(source_2)

    EnergySystems.control(demand, components, simulation_parameters)
    EnergySystems.control(bus, components, simulation_parameters)
    EnergySystems.control(storage, components, simulation_parameters)
    EnergySystems.control(source_1, components, simulation_parameters)
    EnergySystems.control(source_2, components, simulation_parameters)

    EnergySystems.process(demand, simulation_parameters)
    @test demand.input_interfaces[demand.medium].balance == -1000.0
    @test demand.input_interfaces[demand.medium].temperature == 55.0

    EnergySystems.process(bus, simulation_parameters)

    EnergySystems.process(source_1, simulation_parameters)
    @test source_1.output_interfaces[source_1.medium].balance == 500.0
    @test source_1.output_interfaces[source_1.medium].temperature == 55.0

    EnergySystems.process(source_2, simulation_parameters)
    @test source_2.output_interfaces[source_2.medium].balance == 500.0
    @test source_2.output_interfaces[source_2.medium].temperature == 55.0

    EnergySystems.process(storage, simulation_parameters)
    @test storage.output_interfaces[storage.medium].balance == 0.0
    EnergySystems.load(storage, simulation_parameters)
    @test storage.input_interfaces[storage.medium].balance == 0.0

    blnc = EnergySystems.balance(bus)
    @test blnc == 0.0

    EnergySystems.distribute!(bus)

    # time step 2: the sources can't cover the demand, storage has to provide the rest

    EnergySystems.reset(demand)
    EnergySystems.reset(bus)
    EnergySystems.reset(storage)
    EnergySystems.reset(source_1)
    EnergySystems.reset(source_2)

    demand.constant_demand = 6000

    EnergySystems.control(demand, components, simulation_parameters)
    EnergySystems.control(bus, components, simulation_parameters)
    EnergySystems.control(storage, components, simulation_parameters)
    EnergySystems.control(source_1, components, simulation_parameters)
    EnergySystems.control(source_2, components, simulation_parameters)

    EnergySystems.process(demand, simulation_parameters)
    @test demand.input_interfaces[demand.medium].balance == -1500.0
    @test demand.input_interfaces[demand.medium].temperature == 55.0

    EnergySystems.process(bus, simulation_parameters)

    EnergySystems.process(source_1, simulation_parameters)
    @test source_1.output_interfaces[source_1.medium].balance == 500.0
    @test source_1.output_interfaces[source_1.medium].temperature == 55.0

    EnergySystems.process(source_2, simulation_parameters)
    @test source_2.output_interfaces[source_2.medium].balance == 500.0
    @test source_2.output_interfaces[source_2.medium].temperature == 55.0

    EnergySystems.process(storage, simulation_parameters)
    @test storage.output_interfaces[storage.medium].balance == 500.0
    EnergySystems.load(storage, simulation_parameters)
    @test storage.input_interfaces[storage.medium].balance == 0.0

    blnc = EnergySystems.balance(bus)
    @test blnc == 0.0
    @test storage.load == 4500

    EnergySystems.distribute!(bus)

    # time step 3: the sources can cover the demand and have some production capacity left,
    # which is used to load the storage

    EnergySystems.reset(demand)
    EnergySystems.reset(bus)
    EnergySystems.reset(storage)
    EnergySystems.reset(source_1)
    EnergySystems.reset(source_2)

    demand.constant_demand = 3000

    EnergySystems.control(demand, components, simulation_parameters)
    EnergySystems.control(bus, components, simulation_parameters)
    EnergySystems.control(storage, components, simulation_parameters)
    EnergySystems.control(source_1, components, simulation_parameters)
    EnergySystems.control(source_2, components, simulation_parameters)

    EnergySystems.process(demand, simulation_parameters)
    @test demand.input_interfaces[demand.medium].balance == -750.0
    @test demand.input_interfaces[demand.medium].temperature == 55.0

    EnergySystems.process(bus, simulation_parameters)

    EnergySystems.process(source_1, simulation_parameters)
    @test source_1.output_interfaces[source_1.medium].balance == 500.0
    @test source_1.output_interfaces[source_1.medium].temperature == 55.0

    EnergySystems.process(source_2, simulation_parameters)
    @test source_2.output_interfaces[source_2.medium].balance == 500.0
    @test source_2.output_interfaces[source_2.medium].temperature == 55.0

    EnergySystems.process(storage, simulation_parameters)
    @test storage.output_interfaces[storage.medium].balance == 0.0
    EnergySystems.load(storage, simulation_parameters)
    @test storage.input_interfaces[storage.medium].balance == -250.0

    blnc = EnergySystems.balance(bus)
    @test blnc == 0.0
    @test storage.load == 4750

    EnergySystems.distribute!(bus)

    @test demand.input_interfaces[demand.medium].balance == 0.0
    @test source_1.output_interfaces[source_1.medium].balance == 0.0
    @test source_2.output_interfaces[source_2.medium].balance == 0.0
    @test storage.output_interfaces[storage.medium].balance == 0.0
    @test storage.input_interfaces[storage.medium].balance == 0.0
end

@testset "many_plus_storage_to_one" begin
    test_many_plus_storage_to_one()
end