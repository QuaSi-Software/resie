using Debugger
using Test
using Resie
using Resie.EnergySystems
using Resie.Profiles

EnergySystems.set_timestep(900)

function test_one_plus_storage_to_many()
    components_config = Dict{String,Any}(
        "TST_SRC_01" => Dict{String,Any}(
            "type" => "BoundedSupply",
            "medium" => "m_h_w_ht1",
            "output_refs" => ["TST_BUS_TH_01"],
            "constant_power" => 4000,
            "constant_temperature" => 55,
        ),
        "TST_DEM_01" => Dict{String,Any}(
            "type" => "Demand",
            "medium" => "m_h_w_ht1",
            "output_refs" => [],
            "constant_demand" => 2000,
            "constant_temperature" => 55,
        ),
        "TST_DEM_02" => Dict{String,Any}(
            "type" => "Demand",
            "medium" => "m_h_w_ht1",
            "output_refs" => [],
            "constant_demand" => 2000,
            "constant_temperature" => 55,
        ),
        "TST_TES_01" => Dict{String,Any}(
            "type" => "Storage",
            "medium" => "m_h_w_ht1",
            "output_refs" => ["TST_BUS_TH_01"],
            "capacity" => 10000,
            "load" => 5000,
        ),
        "TST_BUS_TH_01" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_h_w_ht1",
            "connections" => Dict{String,Any}(
                "input_order" => ["TST_SRC_01",
                                  "TST_TES_01"],
                "output_order" => ["TST_DEM_01",
                                   "TST_DEM_02",
                                   "TST_TES_01"],
                "energy_flow" => [[1, 1, 1],
                                  [1, 1, 0]],
            ),
        ),
    )

    simulation_parameters = get_default_sim_params()

    components = Resie.load_components(components_config, simulation_parameters)
    source = components["TST_SRC_01"]
    demand_1 = components["TST_DEM_01"]
    demand_2 = components["TST_DEM_02"]
    bus = components["TST_BUS_TH_01"]
    storage = components["TST_TES_01"]

    # time step 1: the source can't cover the demands, the storage needs to fill the rest

    source.constant_power = 3000

    EnergySystems.reset(demand_1)
    EnergySystems.reset(demand_2)
    EnergySystems.reset(bus)
    EnergySystems.reset(storage)
    EnergySystems.reset(source)

    EnergySystems.control(demand_1, components, simulation_parameters)
    EnergySystems.control(demand_2, components, simulation_parameters)
    EnergySystems.control(bus, components, simulation_parameters)
    EnergySystems.control(storage, components, simulation_parameters)
    EnergySystems.control(source, components, simulation_parameters)

    EnergySystems.process(demand_1, simulation_parameters)
    @test demand_1.input_interfaces[demand_1.medium].balance == -500.0
    @test demand_1.input_interfaces[demand_1.medium].temperature_min == 55.0

    EnergySystems.process(demand_2, simulation_parameters)
    @test demand_2.input_interfaces[demand_2.medium].balance == -500.0
    @test demand_2.input_interfaces[demand_2.medium].temperature_min == 55.0

    EnergySystems.process(bus, simulation_parameters)

    EnergySystems.process(source, simulation_parameters)
    @test source.output_interfaces[source.medium].balance == 750.0
    @test source.output_interfaces[source.medium].temperature_max == 55.0

    EnergySystems.process(storage, simulation_parameters)
    @test storage.output_interfaces[storage.medium].balance == 250.0
    EnergySystems.load(storage, simulation_parameters)
    @test storage.input_interfaces[storage.medium].balance == 0.0

    blnc = EnergySystems.balance(bus)
    @test blnc == 0.0

    EnergySystems.distribute!(bus)
    @test storage.load == 4750

    # time step 2: the source covers the demands exactly, the storage does nothing

    source.constant_power = 4000

    EnergySystems.reset(demand_1)
    EnergySystems.reset(demand_2)
    EnergySystems.reset(bus)
    EnergySystems.reset(storage)
    EnergySystems.reset(source)

    EnergySystems.control(demand_1, components, simulation_parameters)
    EnergySystems.control(demand_2, components, simulation_parameters)
    EnergySystems.control(bus, components, simulation_parameters)
    EnergySystems.control(storage, components, simulation_parameters)
    EnergySystems.control(source, components, simulation_parameters)

    EnergySystems.process(demand_1, simulation_parameters)
    @test demand_1.input_interfaces[demand_1.medium].balance == -500.0
    @test demand_1.input_interfaces[demand_1.medium].temperature_min == 55.0

    EnergySystems.process(demand_2, simulation_parameters)
    @test demand_2.input_interfaces[demand_2.medium].balance == -500.0
    @test demand_2.input_interfaces[demand_2.medium].temperature_min == 55.0

    EnergySystems.process(bus, simulation_parameters)

    EnergySystems.process(source, simulation_parameters)
    @test source.output_interfaces[source.medium].balance == 1000.0
    @test source.output_interfaces[source.medium].temperature_max == 55.0

    EnergySystems.process(storage, simulation_parameters)
    @test storage.output_interfaces[storage.medium].balance == 0.0
    EnergySystems.load(storage, simulation_parameters)
    @test storage.input_interfaces[storage.medium].balance == 0.0

    blnc = EnergySystems.balance(bus)
    @test blnc == 0.0

    EnergySystems.distribute!(bus)
    @test storage.load == 4750

    # time step 3: the source can cover the demands and has excess production capacity to
    # load the storage

    source.constant_power = 5000

    EnergySystems.reset(demand_1)
    EnergySystems.reset(demand_2)
    EnergySystems.reset(bus)
    EnergySystems.reset(storage)
    EnergySystems.reset(source)

    EnergySystems.control(demand_1, components, simulation_parameters)
    EnergySystems.control(demand_2, components, simulation_parameters)
    EnergySystems.control(bus, components, simulation_parameters)
    EnergySystems.control(storage, components, simulation_parameters)
    EnergySystems.control(source, components, simulation_parameters)

    EnergySystems.process(demand_1, simulation_parameters)
    @test demand_1.input_interfaces[demand_1.medium].balance == -500.0
    @test demand_1.input_interfaces[demand_1.medium].temperature_min == 55.0

    EnergySystems.process(demand_2, simulation_parameters)
    @test demand_2.input_interfaces[demand_2.medium].balance == -500.0
    @test demand_2.input_interfaces[demand_2.medium].temperature_min == 55.0

    EnergySystems.process(bus, simulation_parameters)

    EnergySystems.process(source, simulation_parameters)
    @test source.output_interfaces[source.medium].balance == 1250.0
    @test source.output_interfaces[source.medium].temperature_max == 55.0

    EnergySystems.process(storage, simulation_parameters)
    @test storage.output_interfaces[storage.medium].balance == 0.0
    EnergySystems.load(storage, simulation_parameters)
    @test storage.input_interfaces[storage.medium].balance == -250.0

    blnc = EnergySystems.balance(bus)
    @test blnc == 0.0

    EnergySystems.distribute!(bus)
    @test storage.load == 5000
end

@testset "one_plus_storage_to_many" begin
    test_one_plus_storage_to_many()
end
