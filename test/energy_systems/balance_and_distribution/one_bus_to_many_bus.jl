using Debugger
using Test
using Resie
using Resie.EnergySystems
using Resie.Profiles

include("../../test_util.jl")

function test_one_bus_to_many_bus()
    components_config = Dict{String,Any}(
        "TST_SRC_01" => Dict{String,Any}(
            "type" => "BoundedSupply",
            "medium" => "m_h_w_ht1",
            "output_refs" => ["TST_BUS_01"],
            "is_source" => true,
            "constant_power" => 400,
        ),
        "TST_BUS_01" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_h_w_ht1",
            "connections" => Dict{String,Any}(
                "input_order" => ["TST_SRC_01"],
                "output_order" => ["TST_BUS_02",
                                   "TST_BUS_03"],
            ),
        ),
        "TST_BUS_02" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_h_w_ht1",
            "connections" => Dict{String,Any}(
                "input_order" => ["TST_BUS_01",
                                  "TST_TES_01"],
                "output_order" => ["TST_DEM_01",
                                   "TST_TES_01"],
            ),
        ),
        "TST_BUS_03" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_h_w_ht1",
            "connections" => Dict{String,Any}(
                "input_order" => ["TST_BUS_01",
                                  "TST_TES_02"],
                "output_order" => ["TST_DEM_02",
                                   "TST_TES_02"],
            ),
        ),
        "TST_TES_01" => Dict{String,Any}(
            "type" => "Storage",
            "medium" => "m_h_w_ht1",
            "output_refs" => ["TST_BUS_02"],
            "capacity" => 1000,
            "load" => 500,
        ),
        "TST_TES_02" => Dict{String,Any}(
            "type" => "Storage",
            "medium" => "m_h_w_ht1",
            "output_refs" => ["TST_BUS_03"],
            "capacity" => 1000,
            "load" => 500,
        ),
        "TST_DEM_01" => Dict{String,Any}(
            "type" => "Demand",
            "medium" => "m_h_w_ht1",
            "output_refs" => [],
            "constant_demand" => 400.0,
            "constant_temperature" => 55.0,
        ),
        "TST_DEM_02" => Dict{String,Any}(
            "type" => "Demand",
            "medium" => "m_h_w_ht1",
            "output_refs" => [],
            "constant_demand" => 400.0,
            "constant_temperature" => 55.0,
        ),
    )

    simulation_parameters = get_default_sim_params()

    components = Resie.load_components(components_config, simulation_parameters)
    demand_1 = components["TST_DEM_01"]
    demand_2 = components["TST_DEM_02"]
    storage_1 = components["TST_TES_01"]
    storage_2 = components["TST_TES_02"]
    source = components["TST_SRC_01"]
    bus_1 = components["TST_BUS_01"]
    bus_2 = components["TST_BUS_02"]
    bus_3 = components["TST_BUS_03"]
    bus_proxy = components["Proxy-TST_BUS_01|TST_BUS_03|TST_BUS_02"]

    # time step 1: source can only cover the first demand in the second bus, the second
    # demand has to be covered by the second storage in the third bus

    source.constant_power = 400

    EnergySystems.reset(demand_1)
    EnergySystems.reset(demand_2)
    EnergySystems.reset(bus_1)
    EnergySystems.reset(bus_2)
    EnergySystems.reset(bus_3)
    EnergySystems.reset(bus_proxy)
    EnergySystems.reset(storage_1)
    EnergySystems.reset(storage_1)
    EnergySystems.reset(source)

    EnergySystems.control(demand_2, components, simulation_parameters)
    EnergySystems.control(demand_1, components, simulation_parameters)
    EnergySystems.control(bus_2, components, simulation_parameters)
    EnergySystems.control(bus_1, components, simulation_parameters)
    EnergySystems.control(bus_3, components, simulation_parameters)
    EnergySystems.control(bus_proxy, components, simulation_parameters)
    EnergySystems.control(storage_1, components, simulation_parameters)
    EnergySystems.control(storage_2, components, simulation_parameters)
    EnergySystems.control(source, components, simulation_parameters)

    exchanges = EnergySystems.balance_on(bus_proxy.input_interfaces[1], bus_1)
    @test EnergySystems.balance(exchanges) ≈ 0.0
    @test EnergySystems.energy_potential(exchanges) ≈ -100.0  # limited by the available energy of the source
    @test EnergySystems.temp_min_highest(exchanges) === 55.0

    exchanges = EnergySystems.balance_on(bus_proxy.output_interfaces[3], bus_3)  # DEM_02
    @test EnergySystems.balance(exchanges) ≈ 0.0
    @test EnergySystems.energy_potential(exchanges) ≈ 100.0  # limited by the actual demand
    @test EnergySystems.temp_min_highest(exchanges) == 55.0
    exchanges = EnergySystems.balance_on(bus_proxy.output_interfaces[1], bus_2)  # DEM_01
    @test EnergySystems.balance(exchanges) ≈ 0.0
    @test EnergySystems.energy_potential(exchanges) ≈ 100.0  # limited by the actual demand
    @test EnergySystems.temp_min_highest(exchanges) == 55.0

    EnergySystems.process(demand_2, simulation_parameters)
    EnergySystems.process(demand_1, simulation_parameters)
    EnergySystems.process(bus_2, simulation_parameters)
    EnergySystems.process(bus_1, simulation_parameters)
    EnergySystems.process(bus_3, simulation_parameters)
    EnergySystems.process(bus_proxy, simulation_parameters)

    exchanges = EnergySystems.balance_on(bus_proxy.input_interfaces[3], bus_2)  # TES_01
    @test EnergySystems.balance(exchanges) ≈ 0.0
    @test EnergySystems.energy_potential(exchanges) ≈ 0.0
    @test EnergySystems.temp_min_highest(exchanges) === nothing

    exchanges = EnergySystems.balance_on(bus_proxy.input_interfaces[2], bus_3)  # TES_02
    @test EnergySystems.balance(exchanges) ≈ 0.0
    @test EnergySystems.energy_potential(exchanges) ≈ -100.0
    @test EnergySystems.temp_min_highest(exchanges) === 55.0

    EnergySystems.process(storage_1, simulation_parameters)
    EnergySystems.process(storage_2, simulation_parameters)

    exchanges = EnergySystems.balance_on(bus_proxy.output_interfaces[2], bus_2)  # TES_01
    @test EnergySystems.balance(exchanges) ≈ 0.0
    @test EnergySystems.energy_potential(exchanges) ≈ 0.0
    @test EnergySystems.temp_min_highest(exchanges) === nothing

    exchanges = EnergySystems.balance_on(bus_proxy.output_interfaces[4], bus_3)  # TES_02
    @test EnergySystems.balance(exchanges) ≈ 0.0
    @test EnergySystems.energy_potential(exchanges) ≈ 0.0
    @test EnergySystems.temp_min_highest(exchanges) === nothing

    EnergySystems.load(storage_1, simulation_parameters)
    EnergySystems.load(storage_2, simulation_parameters)

    exchanges = EnergySystems.balance_on(bus_proxy.input_interfaces[1], bus_2)  # SRC_01
    @test EnergySystems.balance(exchanges) ≈ 0.0
    @test EnergySystems.energy_potential(exchanges) ≈ -100.0
    @test EnergySystems.temp_min_highest(exchanges) === 55.0

    EnergySystems.process(source, simulation_parameters)
    @test source.output_interfaces[source.medium].sum_abs_change ≈ 100.0

    EnergySystems.distribute!(bus_proxy)
    EnergySystems.distribute!(bus_2)
    EnergySystems.distribute!(bus_3)
    EnergySystems.distribute!(bus_1)

    # energy has been communicated to all interfaces of all three busses
    @test EnergySystems.balance(bus_1) ≈ 0.0
    @test EnergySystems.balance(bus_2) ≈ 0.0
    @test EnergySystems.balance(bus_3) ≈ 0.0
    @test EnergySystems.balance(bus_proxy) ≈ 0.0

    @test demand_1.input_interfaces[demand_1.medium].sum_abs_change ≈ 2 * 100.0
    @test demand_2.input_interfaces[demand_2.medium].sum_abs_change ≈ 2 * 100.0
    @test storage_1.input_interfaces[storage_1.medium].sum_abs_change ≈ 0.0
    @test storage_2.input_interfaces[storage_2.medium].sum_abs_change ≈ 0.0
    @test storage_1.output_interfaces[storage_1.medium].sum_abs_change ≈ 0.0
    @test storage_2.output_interfaces[storage_2.medium].sum_abs_change ≈ 2 * 100.0
    @test source.output_interfaces[source.medium].sum_abs_change ≈ 2 * 100.0

    @test bus_2.input_interfaces[1].sum_abs_change ≈ 2 * 100.0   # from bus_1
    @test bus_2.input_interfaces[2].sum_abs_change ≈ 0.0         # from storage_1
    @test bus_2.output_interfaces[1].sum_abs_change ≈ 2 * 100.0  # to dem_1
    @test bus_2.output_interfaces[2].sum_abs_change ≈ 0.0        # to storage_1

    @test bus_3.input_interfaces[1].sum_abs_change ≈ 0.0         # from bus_1
    @test bus_3.input_interfaces[2].sum_abs_change ≈ 2 * 100.0   # from storage_2
    @test bus_3.output_interfaces[1].sum_abs_change ≈ 2 * 100.0  # to dem_2
    @test bus_3.output_interfaces[2].sum_abs_change ≈ 0.0        # to storage_2

    @test bus_1.output_interfaces[1].sum_abs_change ≈ 2 * 100.0  # to bus_1
    @test bus_1.output_interfaces[2].sum_abs_change ≈ 0.0        # to bus_2
    @test bus_1.input_interfaces[1].sum_abs_change ≈ 2 * 100.0   # from src_1
end

@testset "one_bus_to_many_bus" begin
    test_one_bus_to_many_bus()
end
