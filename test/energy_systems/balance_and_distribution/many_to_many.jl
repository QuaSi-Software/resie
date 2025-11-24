using Debugger
using Test
using Resie
using Resie.EnergySystems
using Resie.Profiles

include("../../test_util.jl")

function test_many_to_many()
    components_config = Dict{String,Any}(
        "TST_SRC_01" => Dict{String,Any}(
            "type" => "FlexibleSupply",
            "medium" => "m_h_w_ht1",
            "output_refs" => ["TST_BUS_TH_01"],
            "constant_power" => 1000,
            "constant_temperature" => 55,
        ),
        "TST_SRC_02" => Dict{String,Any}(
            "type" => "FlexibleSupply",
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
        "TST_TES_02" => Dict{String,Any}(
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
                                  "TST_TES_01",
                                  "TST_TES_02",
                                  "TST_SRC_02"],
                "output_order" => ["TST_DEM_01",
                                   "TST_DEM_02",
                                   "TST_TES_01",
                                   "TST_TES_02"],
                "energy_flow" => [[1, 1, 1, 1],
                                  [1, 0, 0, 0],
                                  [0, 1, 0, 0],
                                  [1, 1, 0, 0]],
            ),
        ),
    )

    simulation_parameters = get_default_sim_params()

    components = Resie.load_components(components_config, simulation_parameters)
    source_1 = components["TST_SRC_01"]
    source_2 = components["TST_SRC_02"]
    demand_1 = components["TST_DEM_01"]
    demand_2 = components["TST_DEM_02"]
    bus = components["TST_BUS_TH_01"]
    storage_1 = components["TST_TES_01"]
    storage_2 = components["TST_TES_02"]

    # time step 1: source 1 is not enough to meet demands and storages are empty. source 2
    # fills the rest, but is not allowed to load storages
    storage_1.load = 0
    storage_2.load = 0

    EnergySystems.reset(demand_1)
    EnergySystems.reset(demand_2)
    EnergySystems.reset(bus)
    EnergySystems.reset(storage_1)
    EnergySystems.reset(storage_2)
    EnergySystems.reset(source_1)
    EnergySystems.reset(source_2)

    EnergySystems.control(demand_1, components, simulation_parameters)
    EnergySystems.control(demand_2, components, simulation_parameters)
    EnergySystems.control(bus, components, simulation_parameters)
    EnergySystems.control(storage_1, components, simulation_parameters)
    EnergySystems.control(storage_2, components, simulation_parameters)
    EnergySystems.control(source_1, components, simulation_parameters)
    EnergySystems.control(source_2, components, simulation_parameters)

    EnergySystems.process(demand_1, simulation_parameters)
    @test demand_1.input_interfaces[demand_1.medium].balance == -500.0
    @test demand_1.input_interfaces[demand_1.medium].max_energy.temperature_min == [55.0]

    EnergySystems.process(demand_2, simulation_parameters)
    @test demand_2.input_interfaces[demand_2.medium].balance == -500.0
    @test demand_2.input_interfaces[demand_2.medium].max_energy.temperature_min == [55.0]

    EnergySystems.process(bus, simulation_parameters)

    EnergySystems.process(source_1, simulation_parameters)
    @test source_1.output_interfaces[source_1.medium].balance == 250.0
    @test source_1.output_interfaces[source_1.medium].max_energy.temperature_max == [55.0]

    EnergySystems.process(storage_1, simulation_parameters)
    @test storage_1.output_interfaces[storage_1.medium].balance == 0.0
    EnergySystems.process(storage_2, simulation_parameters)
    @test storage_2.output_interfaces[storage_2.medium].balance == 0.0

    EnergySystems.process(source_2, simulation_parameters)
    @test source_2.output_interfaces[source_2.medium].balance == 750.0
    @test source_2.output_interfaces[source_2.medium].max_energy.temperature_max == [55.0]

    EnergySystems.load(storage_1, simulation_parameters)
    @test storage_1.input_interfaces[storage_1.medium].balance == 0.0
    EnergySystems.load(storage_2, simulation_parameters)
    @test storage_2.input_interfaces[storage_2.medium].balance == 0.0

    blnc = EnergySystems.balance(bus)
    @test blnc == 0.0

    EnergySystems.distribute!(bus)

    # time step 2: demands are zero, source 1 can fill storage 1 and 2, source 2 does
    # nothing as it is not allowed to fill storages
    demand_1.constant_demand = 0
    demand_2.constant_demand = 0
    storage_1.load = 9900
    storage_2.load = 0

    EnergySystems.reset(demand_1)
    EnergySystems.reset(demand_2)
    EnergySystems.reset(bus)
    EnergySystems.reset(storage_1)
    EnergySystems.reset(storage_2)
    EnergySystems.reset(source_1)
    EnergySystems.reset(source_2)

    EnergySystems.control(demand_1, components, simulation_parameters)
    EnergySystems.control(demand_2, components, simulation_parameters)
    EnergySystems.control(bus, components, simulation_parameters)
    EnergySystems.control(storage_1, components, simulation_parameters)
    EnergySystems.control(storage_2, components, simulation_parameters)
    EnergySystems.control(source_1, components, simulation_parameters)
    EnergySystems.control(source_2, components, simulation_parameters)

    EnergySystems.process(demand_1, simulation_parameters)
    @test demand_1.input_interfaces[demand_1.medium].balance == 0.0
    @test demand_1.input_interfaces[demand_1.medium].max_energy.temperature_min == [55.0]

    EnergySystems.process(demand_2, simulation_parameters)
    @test demand_2.input_interfaces[demand_2.medium].balance == 0.0
    @test demand_2.input_interfaces[demand_2.medium].max_energy.temperature_min == [55.0]

    EnergySystems.process(bus, simulation_parameters)

    EnergySystems.process(source_1, simulation_parameters)
    @test source_1.output_interfaces[source_1.medium].balance == 250.0
    @test source_1.output_interfaces[source_1.medium].max_energy.temperature_max == [55.0]

    EnergySystems.process(storage_1, simulation_parameters)
    @test storage_1.output_interfaces[storage_1.medium].balance == 0.0
    EnergySystems.process(storage_2, simulation_parameters)
    @test storage_2.output_interfaces[storage_2.medium].balance == 0.0

    EnergySystems.process(source_2, simulation_parameters)
    @test source_2.output_interfaces[source_2.medium].balance == 0.0
    @test source_2.output_interfaces[source_2.medium].max_energy.temperature_max == [55.0]

    EnergySystems.load(storage_1, simulation_parameters)
    @test storage_1.input_interfaces[storage_1.medium].balance == -100.0
    EnergySystems.load(storage_2, simulation_parameters)
    @test storage_2.input_interfaces[storage_2.medium].balance == -150.0

    blnc = EnergySystems.balance(bus)
    @test blnc == 0.0

    EnergySystems.distribute!(bus)
end

@testset "many_to_many" begin
    test_many_to_many()
end
