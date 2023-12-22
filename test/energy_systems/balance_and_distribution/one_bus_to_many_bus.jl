using Debugger
using Test
using Resie
using Resie.EnergySystems
using Resie.Profiles

EnergySystems.set_timestep(900)

function test_one_bus_to_many_bus()
    components_config = Dict{String,Any}(
        "TST_SRC_01" => Dict{String,Any}(
            "type" => "BoundedSupply",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "strategy" => Dict{String,Any}(
                "name" => "extended_storage_control",
                "load_any_storage" => true,
            ),
            "output_refs" => ["TST_BUS_01"],
            "is_source" => true,
            "constant_power" => 400,
        ),
        "TST_BUS_01" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "output_refs" => [
                "TST_BUS_02",
                "TST_BUS_03",
            ],
            "connection_matrix" => Dict{String, Any}(
                "input_order" => [
                    "TST_SRC_01",
                ],
                "output_order" => [
                    "TST_BUS_02",
                    "TST_BUS_03",
                ],
            )
        ),
        "TST_BUS_02" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "output_refs" => [
                "TST_DEM_01",
                "TST_TES_01",
            ],
            "connection_matrix" => Dict{String, Any}(
                "input_order" => [
                    "TST_BUS_01",
                    "TST_TES_01",
                ],
                "output_order" => [
                    "TST_DEM_01",
                    "TST_TES_01",
                ],
            )
        ),
        "TST_BUS_03" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "output_refs" => [
                "TST_DEM_02",
                "TST_TES_02",
            ],
            "connection_matrix" => Dict{String, Any}(
                "input_order" => [
                    "TST_BUS_01",
                    "TST_TES_02",
                ],
                "output_order" => [
                    "TST_DEM_02",
                    "TST_TES_02",
                ],
            )
        ),
        "TST_TES_01" => Dict{String,Any}(
            "type" => "Storage",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "output_refs" => ["TST_BUS_02"],
            "capacity" => 1000,
            "load" => 500,
        ),
        "TST_TES_02" => Dict{String,Any}(
            "type" => "Storage",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "output_refs" => ["TST_BUS_03"],
            "capacity" => 1000,
            "load" => 500,
        ),
        "TST_DEM_01" => Dict{String,Any}(
            "type" => "Demand",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "output_refs" => [],
            "constant_demand" => 400.0,
            "constant_temperature" => 55.0,
        ),
        "TST_DEM_02" => Dict{String,Any}(
            "type" => "Demand",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "output_refs" => [],
            "constant_demand" => 400.0,
            "constant_temperature" => 55.0,
        ),
    )

    simulation_parameters = Dict{String,Any}(
        "time_step_seconds" => 900,
        "time" => 0,
        "epsilon" => 1e-9
    )

    components = Resie.load_components(components_config, simulation_parameters)
    demand_1 = components["TST_DEM_01"]
    demand_2 = components["TST_DEM_02"]
    storage_1 = components["TST_TES_01"]
    storage_2 = components["TST_TES_02"]
    source = components["TST_SRC_01"]
    bus_1 = components["TST_BUS_01"]
    bus_2 = components["TST_BUS_02"]
    bus_3 = components["TST_BUS_03"]


    # time step 1: source can only cover the first demand in the second bus, the second
    # demand has to be covered by the second storage in the third bus

    source.constant_power = 400

    EnergySystems.reset(demand_1)
    EnergySystems.reset(demand_2)
    EnergySystems.reset(bus_1)
    EnergySystems.reset(bus_2)
    EnergySystems.reset(bus_3)
    EnergySystems.reset(storage_1)
    EnergySystems.reset(storage_1)
    EnergySystems.reset(source)

    EnergySystems.control(demand_2, components, simulation_parameters)
    EnergySystems.control(demand_1, components, simulation_parameters)
    EnergySystems.control(bus_2, components, simulation_parameters)
    EnergySystems.control(bus_1, components, simulation_parameters)
    EnergySystems.control(bus_3, components, simulation_parameters)
    EnergySystems.control(storage_1, components, simulation_parameters)
    EnergySystems.control(storage_2, components, simulation_parameters)
    EnergySystems.control(source, components, simulation_parameters)

    EnergySystems.process(demand_2, simulation_parameters)
    EnergySystems.process(demand_1, simulation_parameters)

    exchanges = EnergySystems.balance_on(bus_2.input_interfaces[1], bus_2)
    @test EnergySystems.balance(exchanges) ≈ -100.0
    @test EnergySystems.storage_potential(exchanges) ≈ -500.0
    @test EnergySystems.energy_potential(exchanges) ≈ 0.0
    @test EnergySystems.temperature_first(exchanges) === 55.0

    exchanges = EnergySystems.balance_on(bus_3.input_interfaces[1], bus_3)
    @test EnergySystems.balance(exchanges) ≈ -100.0
    @test EnergySystems.storage_potential(exchanges) ≈ -500.0
    @test EnergySystems.energy_potential(exchanges) ≈ 0.0
    @test EnergySystems.temperature_first(exchanges) === 55.0

    exchanges = EnergySystems.balance_on(bus_1.output_interfaces[1], bus_1)
    @test EnergySystems.balance(exchanges) ≈ -100.0 # not -200 because we "call from" bus 2
    @test EnergySystems.storage_potential(exchanges) ≈ 0.0
    @test EnergySystems.energy_potential(exchanges) ≈ 100.0
    @test EnergySystems.temperature_first(exchanges) === nothing

    EnergySystems.process(bus_2, simulation_parameters)
    EnergySystems.process(bus_1, simulation_parameters)
    EnergySystems.process(bus_3, simulation_parameters)
    EnergySystems.process(storage_1, simulation_parameters)
    EnergySystems.process(storage_2, simulation_parameters)

    exchanges = EnergySystems.balance_on(bus_2.input_interfaces[1], bus_2)
    @test EnergySystems.balance(exchanges) ≈ 0.0
    @test EnergySystems.storage_potential(exchanges) ≈ -600.0
    @test EnergySystems.energy_potential(exchanges) ≈ 0.0
    @test EnergySystems.temperature_first(exchanges) === 55.0

    exchanges = EnergySystems.balance_on(bus_3.input_interfaces[1], bus_3)
    @test EnergySystems.balance(exchanges) ≈ 0.0
    @test EnergySystems.storage_potential(exchanges) ≈ -600.0
    @test EnergySystems.energy_potential(exchanges) ≈ 0.0
    @test EnergySystems.temperature_first(exchanges) === 55.0

    EnergySystems.load(storage_1, simulation_parameters)
    EnergySystems.load(storage_2, simulation_parameters)

    # there is a "bug" here where the storages don't load any energy as the source has not
    # processed yet and storages only load according to balance, not energy potential.
    # when the source processes, the balance is already zero due to the storages covering
    # demand, but storage potential is still there, so the sources produces energy which is
    # not taken in by anything

    exchanges = EnergySystems.balance_on(bus_2.input_interfaces[1], bus_2)
    @test EnergySystems.balance(exchanges) ≈ 0.0
    @test EnergySystems.storage_potential(exchanges) ≈ -600.0
    @test EnergySystems.energy_potential(exchanges) ≈ 0.0
    @test EnergySystems.temperature_first(exchanges) === 55.0

    exchanges = EnergySystems.balance_on(bus_3.input_interfaces[1], bus_3)
    @test EnergySystems.balance(exchanges) ≈ 0.0
    @test EnergySystems.storage_potential(exchanges) ≈ -600.0
    @test EnergySystems.energy_potential(exchanges) ≈ 0.0
    @test EnergySystems.temperature_first(exchanges) === 55.0

    EnergySystems.process(source, simulation_parameters)
    @test source.output_interfaces[source.medium].sum_abs_change ≈ 100.0

    EnergySystems.distribute!(bus_2)
    EnergySystems.distribute!(bus_3)
    EnergySystems.distribute!(bus_1)

    # because the energy is with the source, it's communicated to all three busses
    @test EnergySystems.balance(bus_1) ≈ 100.0
    @test EnergySystems.balance(bus_2) ≈ 100.0
    @test EnergySystems.balance(bus_3) ≈ 100.0
end

@testset "one_bus_to_many_bus" begin
    test_one_bus_to_many_bus()
end
