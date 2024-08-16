using Debugger
using Test
using Resie
using Resie.EnergySystems
using Resie.Profiles

EnergySystems.set_timestep(900)

function test_one_bus_to_one_bus()
    components_config = Dict{String,Any}(
        "TST_GRI_01" => Dict{String,Any}(
            "type" => "GridConnection",
            "medium" => "m_h_w_ht1",
            "output_refs" => ["TST_BUS_01"],
            "is_source" => true,
        ),
        "TST_BUS_01" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_h_w_ht1",
            "connections" => Dict{String,Any}(
                "input_order" => ["TST_GRI_01"],
                "output_order" => ["TST_BUS_02"],
            ),
        ),
        "TST_BUS_02" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_h_w_ht1",
            "connections" => Dict{String,Any}(
                "input_order" => ["TST_BUS_01"],
                "output_order" => ["TST_DEM_01"],
            ),
        ),
        "TST_DEM_01" => Dict{String,Any}(
            "type" => "Demand",
            "medium" => "m_h_w_ht1",
            "output_refs" => [],
            "constant_demand" => 300.0,
            "constant_temperature" => 55.0,
        ),
    )

    simulation_parameters = Dict{String,Any}(
        "time_step_seconds" => 900,
        "time" => 0,
        "epsilon" => 1e-9,
    )

    components = Resie.load_components(components_config, simulation_parameters)
    demand = components["TST_DEM_01"]
    grid = components["TST_GRI_01"]
    bus_1 = components["TST_BUS_01"]
    bus_2 = components["TST_BUS_02"]
    bus_proxy = components["Proxy-TST_BUS_01|TST_BUS_02"]

    EnergySystems.reset(demand)
    EnergySystems.reset(grid)
    EnergySystems.reset(bus_1)
    EnergySystems.reset(bus_2)
    EnergySystems.reset(bus_proxy)

    EnergySystems.control(demand, components, simulation_parameters)
    EnergySystems.control(grid, components, simulation_parameters)
    EnergySystems.control(bus_1, components, simulation_parameters)
    EnergySystems.control(bus_2, components, simulation_parameters)
    EnergySystems.control(bus_proxy, components, simulation_parameters)

    @test demand.input_interfaces[demand.medium].balance ≈ 0.0
    @test demand.input_interfaces[demand.medium].temperature_min === 55.0
    @test EnergySystems.balance(bus_1) ≈ 0.0
    @test EnergySystems.balance(bus_2) ≈ 0.0
    @test grid.output_interfaces[grid.medium].balance ≈ 0.0

    # demand not processed yet --> balance is zero, but energy_potential not
    exchanges = EnergySystems.balance_on(bus_1.input_interfaces[1], bus_1)
    @test EnergySystems.balance(exchanges) ≈ 0.0
    @test EnergySystems.energy_potential(exchanges) ≈ -75.0
    @test EnergySystems.temp_min_highest(exchanges) === 55.0
    exchanges = EnergySystems.balance_on(bus_proxy.input_interfaces[1], bus_1)
    @test EnergySystems.balance(exchanges) ≈ 0.0
    @test EnergySystems.energy_potential(exchanges) ≈ -75.0
    @test EnergySystems.temp_min_highest(exchanges) === 55.0
    exchanges = EnergySystems.balance_on(bus_2.output_interfaces[1], bus_2)
    @test EnergySystems.balance(exchanges) ≈ 0.0
    @test EnergySystems.energy_potential(exchanges) ≈ Inf
    @test EnergySystems.temp_min_highest(exchanges) === nothing
    exchanges = EnergySystems.balance_on(bus_proxy.output_interfaces[1], bus_2)
    @test EnergySystems.balance(exchanges) ≈ 0.0
    @test EnergySystems.energy_potential(exchanges) ≈ Inf
    @test EnergySystems.temp_min_highest(exchanges) === nothing

    EnergySystems.process(demand, simulation_parameters)

    @test demand.input_interfaces[demand.medium].balance ≈ -75.0
    @test demand.input_interfaces[demand.medium].temperature_min === 55.0
    @test EnergySystems.balance(bus_proxy) ≈ -75.0
    @test grid.output_interfaces[grid.medium].balance ≈ 0.0

    # demand already processed --> balance is zero in both busses as balance is only 
    # non-zero for 1-1 connections without busses. energy potential is zero for inputs 
    # to the chain as the demand is now written by set_max_energy, but is still infinite
    # for the outputs as the grid has not processed energy yet
    exchanges = EnergySystems.balance_on(bus_1.input_interfaces[1], bus_1)
    @test EnergySystems.balance(exchanges) ≈ 0.0
    @test EnergySystems.energy_potential(exchanges) === -75.0
    @test EnergySystems.temp_min_highest(exchanges) === 55.0
    exchanges = EnergySystems.balance_on(bus_proxy.input_interfaces[1], bus_2)
    @test EnergySystems.balance(exchanges) ≈ 0.0
    @test EnergySystems.energy_potential(exchanges) === -75.0
    @test EnergySystems.temp_min_highest(exchanges) === 55.0
    # The interfaces of the proxy and the orignal bus shoule be the same at this point
    # The interface between bus 1 and bus 2 is not filled before distribute()
    exchanges = EnergySystems.balance_on(bus_2.output_interfaces[1], bus_2)
    @test EnergySystems.balance(exchanges) ≈ 0.0
    @test EnergySystems.energy_potential(exchanges) ≈ Inf
    @test EnergySystems.temp_max_highest(exchanges) === nothing
    exchanges = EnergySystems.balance_on(bus_proxy.output_interfaces[1], bus_2)
    @test EnergySystems.balance(exchanges) ≈ 0.0
    @test EnergySystems.energy_potential(exchanges) ≈ Inf
    @test EnergySystems.temp_max_highest(exchanges) === nothing

    EnergySystems.process(bus_2, simulation_parameters)
    EnergySystems.process(bus_1, simulation_parameters)
    EnergySystems.process(bus_proxy, simulation_parameters)
    EnergySystems.process(grid, simulation_parameters)

    # everything processed
    # balance of balance_on() of bus is always zero.
    exchanges = EnergySystems.balance_on(bus_1.input_interfaces[1], bus_1)
    @test EnergySystems.balance(exchanges) ≈ 0.0
    @test EnergySystems.energy_potential(exchanges) ≈ -75.0  # as max_energy is set on bus_1.input_interfaces[1]
    @test EnergySystems.temp_min_highest(exchanges) === 55.0
    exchanges = EnergySystems.balance_on(bus_proxy.input_interfaces[1], bus_1)
    @test EnergySystems.balance(exchanges) ≈ 0.0
    @test EnergySystems.energy_potential(exchanges) ≈ -75.0  # as max_energy is set on bus_1.input_interfaces[1]
    @test EnergySystems.temp_min_highest(exchanges) === 55.0
    exchanges = EnergySystems.balance_on(bus_2.output_interfaces[1], bus_2)
    @test EnergySystems.balance(exchanges) ≈ 0.0
    @test EnergySystems.energy_potential(exchanges) ≈ 75.0  # as grid is produced, this is the produced energy at this point
    @test EnergySystems.temp_max_highest(exchanges) === nothing
    exchanges = EnergySystems.balance_on(bus_proxy.output_interfaces[1], bus_2)
    @test EnergySystems.balance(exchanges) ≈ 0.0
    @test EnergySystems.energy_potential(exchanges) ≈ 75.0  # as grid is produced, this is the produced energy at this point
    @test EnergySystems.temp_max_highest(exchanges) === nothing

    # balance in interfaces of components are non-zero:
    @test demand.input_interfaces[demand.medium].balance ≈ -75.0
    @test demand.input_interfaces[demand.medium].temperature_min === 55.0
    @test EnergySystems.balance(bus_1) ≈ 75.0  # no distribute!() beween busses yet
    @test EnergySystems.balance(bus_2) ≈ -75.0  # no distribute!() beween busses yet
    @test EnergySystems.balance(bus_proxy) ≈ 0.0
    @test bus_1.remainder ≈ 0.0
    @test bus_2.remainder ≈ 0.0
    @test bus_proxy.remainder ≈ 0.0
    @test grid.output_interfaces[grid.medium].balance ≈ 75.0
    @test bus_proxy.balance_table[1, 1] ≈ 75.0 # energy transferred from source 1 to sink 1 ...
    @test bus_proxy.balance_table[1, 2] ≈ 55.0 # with temperature of 55.0 °C

    EnergySystems.distribute!(bus_proxy)
    EnergySystems.distribute!(bus_2)
    EnergySystems.distribute!(bus_1)

    @test demand.input_interfaces[demand.medium].balance ≈ 0.0
    @test demand.input_interfaces[demand.medium].temperature_min === 55.0
    @test demand.input_interfaces[demand.medium].sum_abs_change ≈ 150.0
    @test EnergySystems.balance(bus_1) ≈ 0.0
    @test EnergySystems.balance(bus_2) ≈ 0.0
    @test EnergySystems.balance(bus_proxy) ≈ 0.0
    @test bus_1.remainder ≈ 0.0
    @test bus_2.remainder ≈ 0.0
    @test grid.output_interfaces[grid.medium].balance ≈ 0.0
    @test grid.output_interfaces[grid.medium].sum_abs_change ≈ 150.0
    @test bus_1.output_interfaces[1].balance ≈ 0.0
    @test bus_1.output_interfaces[1].sum_abs_change ≈ 150.0
    @test bus_2.input_interfaces[1].balance ≈ 0.0
    @test bus_2.input_interfaces[1].sum_abs_change ≈ 150.0
end

@testset "one_bus_to_one_bus" begin
    test_one_bus_to_one_bus()
end
