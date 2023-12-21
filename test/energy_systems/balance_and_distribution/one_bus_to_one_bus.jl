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
            "control_refs" => [],
            "output_refs" => ["TST_BUS_01"],
            "is_source" => true,
        ),
        "TST_BUS_01" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "output_refs" => ["TST_BUS_02"],
            "connection_matrix" => Dict{String, Any}(
                "input_order" => [
                    "TST_GRI_01",
                ],
                "output_order" => [
                    "TST_BUS_02",
                ],
            )
        ),
        "TST_BUS_02" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "output_refs" => ["TST_DEM_01"],
            "connection_matrix" => Dict{String, Any}(
                "input_order" => [
                    "TST_BUS_01",
                ],
                "output_order" => [
                    "TST_DEM_01",
                ],
            )
        ),
        "TST_DEM_01" => Dict{String,Any}(
            "type" => "Demand",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "output_refs" => [],
            "constant_demand" => 300.0,
            "constant_temperature" => 55.0,
        ),
    )

    simulation_parameters = Dict{String,Any}(
        "time_step_seconds" => 900,
        "time" => 0,
        "epsilon" => 1e-9
    )

    components = Resie.load_components(components_config, simulation_parameters)
    demand = components["TST_DEM_01"]
    grid = components["TST_GRI_01"]
    bus_1 = components["TST_BUS_01"]
    bus_2 = components["TST_BUS_02"]

    EnergySystems.reset(demand)
    EnergySystems.reset(grid)
    EnergySystems.reset(bus_1)
    EnergySystems.reset(bus_2)

    EnergySystems.control(demand, components, simulation_parameters)
    EnergySystems.control(grid, components, simulation_parameters)
    EnergySystems.control(bus_1, components, simulation_parameters)
    EnergySystems.control(bus_2, components, simulation_parameters)

    @test demand.input_interfaces[demand.medium].balance ≈ 0.0
    @test demand.input_interfaces[demand.medium].temperature === 55.0
    @test EnergySystems.balance(bus_1) ≈ 0.0
    @test EnergySystems.balance(bus_2) ≈ 0.0
    @test grid.output_interfaces[grid.medium].balance ≈ 0.0

    # demand not processed yet --> balance is zero, but energy_potential not
    exchanges = EnergySystems.balance_on(bus_1.input_interfaces[1], bus_1)
    @test EnergySystems.balance(exchanges) ≈ 0.0
    @test EnergySystems.storage_potential(exchanges) ≈ 0.0
    @test EnergySystems.energy_potential(exchanges) ≈ -75.0
    @test EnergySystems.temperature_first(exchanges) === 55.0
    exchanges = EnergySystems.balance_on(bus_2.input_interfaces[1], bus_2)
    @test EnergySystems.balance(exchanges) ≈ 0.0
    @test EnergySystems.storage_potential(exchanges) ≈ 0.0
    @test EnergySystems.energy_potential(exchanges) ≈ -75.0
    @test EnergySystems.temperature_first(exchanges) === 55.0
    exchanges = EnergySystems.balance_on(bus_1.output_interfaces[1], bus_1)
    @test EnergySystems.balance(exchanges) ≈ 0.0
    @test EnergySystems.storage_potential(exchanges) ≈ 0.0
    @test EnergySystems.energy_potential(exchanges) ≈ Inf
    @test EnergySystems.temperature_first(exchanges) === nothing
    exchanges = EnergySystems.balance_on(bus_2.output_interfaces[1], bus_2)
    @test EnergySystems.balance(exchanges) ≈ 0.0
    @test EnergySystems.storage_potential(exchanges) ≈ 0.0
    @test EnergySystems.energy_potential(exchanges) ≈ Inf
    @test EnergySystems.temperature_first(exchanges) === 55.0

    EnergySystems.process(demand, simulation_parameters)

    @test demand.input_interfaces[demand.medium].balance ≈ -75.0
    @test demand.input_interfaces[demand.medium].temperature === 55.0
    @test EnergySystems.balance(bus_1) ≈ -75.0
    @test EnergySystems.balance(bus_2) ≈ -75.0
    @test grid.output_interfaces[grid.medium].balance ≈ 0.0

    # demand already processed --> balance is not zero in both busses. energy potential
    # is zero for inputs to the chain as the demand is now written, but is still infinite
    # as the grid has not processed energy yet
    exchanges = EnergySystems.balance_on(bus_1.input_interfaces[1], bus_1)
    @test EnergySystems.balance(exchanges) ≈ -75.0
    @test EnergySystems.storage_potential(exchanges) ≈ 0.0
    @test EnergySystems.energy_potential(exchanges) ≈ 0.0
    @test EnergySystems.temperature_first(exchanges) === 55.0
    exchanges = EnergySystems.balance_on(bus_2.input_interfaces[1], bus_2)
    @test EnergySystems.balance(exchanges) ≈ -75.0
    @test EnergySystems.storage_potential(exchanges) ≈ 0.0
    @test EnergySystems.energy_potential(exchanges) ≈ 0.0
    @test EnergySystems.temperature_first(exchanges) === 55.0
    # from the perspective of the interface between bus 1 and 2, bus 1 has a balance of
    # zero as the demand does not "backflow" from bus 2 to 1
    exchanges = EnergySystems.balance_on(bus_1.output_interfaces[1], bus_1)
    @test EnergySystems.balance(exchanges) ≈ 0.0
    @test EnergySystems.storage_potential(exchanges) ≈ 0.0
    @test EnergySystems.energy_potential(exchanges) ≈ Inf
    @test EnergySystems.temperature_first(exchanges) === nothing
    exchanges = EnergySystems.balance_on(bus_2.output_interfaces[1], bus_2)
    @test EnergySystems.balance(exchanges) ≈ -75.0
    @test EnergySystems.storage_potential(exchanges) ≈ 0.0
    @test EnergySystems.energy_potential(exchanges) ≈ Inf
    @test EnergySystems.temperature_first(exchanges) === 55.0

    EnergySystems.process(bus_2, simulation_parameters)
    EnergySystems.process(bus_1, simulation_parameters)
    EnergySystems.process(grid, simulation_parameters)

    # everything processed --> energy_potential should be zero!
    exchanges = EnergySystems.balance_on(bus_1.input_interfaces[1], bus_1)
    @test EnergySystems.balance(exchanges) ≈ 0.0
    @test EnergySystems.storage_potential(exchanges) ≈ 0.0
    @test EnergySystems.energy_potential(exchanges) ≈ 0.0
    @test EnergySystems.temperature_first(exchanges) === 55.0
    exchanges = EnergySystems.balance_on(bus_2.input_interfaces[1], bus_2)
    @test EnergySystems.balance(exchanges) ≈ -75.0
    @test EnergySystems.storage_potential(exchanges) ≈ 0.0
    @test EnergySystems.energy_potential(exchanges) ≈ 0.0
    @test EnergySystems.temperature_first(exchanges) === 55.0
    exchanges = EnergySystems.balance_on(bus_1.output_interfaces[1], bus_1)
    @test EnergySystems.balance(exchanges) ≈ 75.0
    @test EnergySystems.storage_potential(exchanges) ≈ 0.0
    @test EnergySystems.energy_potential(exchanges) ≈ 0.0
    @test EnergySystems.temperature_first(exchanges) === nothing
    exchanges = EnergySystems.balance_on(bus_2.output_interfaces[1], bus_2)
    @test EnergySystems.balance(exchanges) ≈ 0.0
    @test EnergySystems.storage_potential(exchanges) ≈ 0.0
    @test EnergySystems.energy_potential(exchanges) ≈ 0.0
    @test EnergySystems.temperature_first(exchanges) === 55.0

    @test demand.input_interfaces[demand.medium].balance ≈ -75.0
    @test demand.input_interfaces[demand.medium].temperature === 55.0
    @test EnergySystems.balance(bus_1) ≈ 0.0
    @test EnergySystems.balance(bus_2) ≈ 0.0
    @test bus_1.remainder ≈ 0.0
    @test bus_2.remainder ≈ 0.0
    @test grid.output_interfaces[grid.medium].balance ≈ 75.0

    EnergySystems.distribute!(bus_2)
    EnergySystems.distribute!(bus_1)

    @test demand.input_interfaces[demand.medium].balance ≈ 0.0
    @test demand.input_interfaces[demand.medium].temperature === 55.0
    @test demand.input_interfaces[demand.medium].sum_abs_change ≈ 150.0
    @test EnergySystems.balance(bus_1) ≈ 0.0
    @test EnergySystems.balance(bus_2) ≈ 0.0
    @test bus_1.remainder ≈ 0.0
    @test bus_2.remainder ≈ 0.0
    @test grid.output_interfaces[grid.medium].balance ≈ 0.0
    @test grid.output_interfaces[grid.medium].sum_abs_change ≈ 150.0
end

@testset "one_bus_to_one_bus" begin
    test_one_bus_to_one_bus()
end
