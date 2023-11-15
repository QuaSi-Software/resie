using Debugger
using Test
using Resie
using Resie.EnergySystems
using Resie.Profiles

EnergySystems.set_timestep(900)

function test_busses_communicate_demand()
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
            "energy_profile_file_path" => "./profiles/tests/demand_heating_energy.prf",
            "temperature_profile_file_path" => "./profiles/tests/demand_heating_temperature.prf",
            "scale" => 1000
        ),
    )
    components = Resie.load_components(components_config)
    demand = components["TST_DEM_01"]
    grid = components["TST_GRI_01"]
    bus_1 = components["TST_BUS_01"]
    bus_2 = components["TST_BUS_02"]

    simulation_parameters = Dict{String,Any}(
        "time_step_seconds" => 900,
        "time" => 0,
        "epsilon" => 1e-9
    )

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
    exchange = EnergySystems.balance_on(bus_1.input_interfaces[1], bus_1)
    @test exchange.balance ≈ 0.0
    @test exchange.storage_potential ≈ 0.0
    @test exchange.energy_potential ≈ -75.0
    @test exchange.temperature === 55.0
    exchange = EnergySystems.balance_on(bus_2.input_interfaces[1], bus_2)
    @test exchange.balance ≈ 0.0
    @test exchange.storage_potential ≈ 0.0
    @test exchange.energy_potential ≈ -75.0
    @test exchange.temperature === 55.0
    exchange = EnergySystems.balance_on(bus_1.output_interfaces[1], bus_1)
    @test exchange.balance ≈ 0.0
    @test exchange.storage_potential ≈ 0.0
    @test exchange.energy_potential ≈ Inf
    @test exchange.temperature === nothing
    exchange = EnergySystems.balance_on(bus_2.output_interfaces[1], bus_2)
    @test exchange.balance ≈ 0.0
    @test exchange.storage_potential ≈ 0.0
    @test exchange.energy_potential ≈ Inf
    @test exchange.temperature === 55.0 # is not nothing as temperature is given in interface from demand

    EnergySystems.process(demand, simulation_parameters)

    @test demand.input_interfaces[demand.medium].balance ≈ -75.0
    @test demand.input_interfaces[demand.medium].temperature === 55.0
    @test EnergySystems.balance(bus_1) ≈ -75.0
    @test EnergySystems.balance(bus_2) ≈ -75.0
    @test grid.output_interfaces[grid.medium].balance ≈ 0.0

    # demand already processed --> balance is not zero, but energy_potential
    exchange = EnergySystems.balance_on(bus_1.input_interfaces[1], bus_1)
    @test exchange.balance ≈ -75.0
    @test exchange.storage_potential ≈ 0.0
    @test exchange.energy_potential ≈ 0.0
    @test exchange.temperature === 55.0
    exchange = EnergySystems.balance_on(bus_2.input_interfaces[1], bus_2)
    @test exchange.balance ≈ -75.0
    @test exchange.storage_potential ≈ 0.0
    @test exchange.energy_potential ≈ 0.0
    @test exchange.temperature === 55.0
    exchange = EnergySystems.balance_on(bus_1.output_interfaces[1], bus_1)
    @test exchange.balance ≈ -75.0
    @test exchange.storage_potential ≈ 0.0
    @test exchange.energy_potential ≈ Inf
    @test exchange.temperature === nothing
    exchange = EnergySystems.balance_on(bus_2.output_interfaces[1], bus_2)
    @test exchange.balance ≈ -75.0
    @test exchange.storage_potential ≈ 0.0
    @test exchange.energy_potential ≈ 0.0 # balance in current interface already written!
    @test exchange.temperature === 55.0 # is not nothing as temperature is given in interface from demand

    EnergySystems.process(bus_2, simulation_parameters)
    EnergySystems.process(bus_1, simulation_parameters)
    EnergySystems.process(grid, simulation_parameters)

    # everything processed --> energy_potential should be zero!
    exchange = EnergySystems.balance_on(bus_1.input_interfaces[1], bus_1)
    @test exchange.balance ≈ 0.0
    @test exchange.storage_potential ≈ 0.0
    @test exchange.energy_potential ≈ 0.0
    @test exchange.temperature === 55.0
    exchange = EnergySystems.balance_on(bus_2.input_interfaces[1], bus_2)
    @test exchange.balance ≈ 0.0
    @test exchange.storage_potential ≈ 0.0
    @test exchange.energy_potential ≈ 0.0
    @test exchange.temperature === 55.0
    exchange = EnergySystems.balance_on(bus_1.output_interfaces[1], bus_1)
    @test exchange.balance ≈ 0.0
    @test exchange.storage_potential ≈ 0.0
    @test exchange.energy_potential ≈ 0.0
    @test exchange.temperature === nothing
    exchange = EnergySystems.balance_on(bus_2.output_interfaces[1], bus_2)
    @test exchange.balance ≈ 0.0
    @test exchange.storage_potential ≈ 0.0
    @test exchange.energy_potential ≈ 0.0
    @test exchange.temperature === 55.0 # is not nothing as temperature is given in interface from demand

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

@testset "busses_communicate_demand" begin
    test_busses_communicate_demand()
end

function test_demand_over_busses_supply_is_transformer()
    components_config = Dict{String,Any}(
        "TST_GRI_01" => Dict{String,Any}(
            "type" => "GridConnection",
            "medium" => "m_c_g_natgas",
            "control_refs" => [],
            "output_refs" => ["TST_GBO_01"],
            "is_source" => true,
        ),
        "TST_GBO_01" => Dict{String,Any}(
            "type" => "FuelBoiler",
            "m_fuel_in" => "m_c_g_natgas",
            "control_refs" => ["TST_BUS_01"],
            "output_refs" => ["TST_BUS_01"],
            "strategy" => Dict{String,Any}(
                "name" => "demand_driven",
            ),
            "power" => 10000
        ),
        "TST_BUS_01" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "output_refs" => ["TST_BUS_02", "TST_BUS_03"],
            "connection_matrix" => Dict{String, Any}(
                "input_order" => [
                    "TST_GBO_01",
                ],
                "output_order" => [
                    "TST_BUS_02",
                    "TST_BUS_03"
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
        "TST_BUS_03" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "output_refs" => ["TST_DEM_02"],
            "connection_matrix" => Dict{String, Any}(
                "input_order" => [
                    "TST_BUS_01",
                ],
                "output_order" => [
                    "TST_DEM_02",
                ],
            )
        ),
        "TST_DEM_01" => Dict{String,Any}(
            "type" => "Demand",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "output_refs" => [],
            "static_load" => 1000,
            "static_temperature" => 60,
            "scale" => 1
        ),
        "TST_DEM_02" => Dict{String,Any}(
            "type" => "Demand",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "output_refs" => [],
            "static_load" => 1000,
            "static_temperature" => 60,
            "scale" => 1
        ),
    )
    components = Resie.load_components(components_config)
    demand_1 = components["TST_DEM_01"]
    demand_2 = components["TST_DEM_02"]
    grid = components["TST_GRI_01"]
    boiler = components["TST_GBO_01"]
    bus_1 = components["TST_BUS_01"]
    bus_2 = components["TST_BUS_02"]
    bus_3 = components["TST_BUS_03"]

    simulation_parameters = Dict{String,Any}(
        "time_step_seconds" => 900,
        "time" => 0,
        "epsilon" => 1e-9
    )

    # first timestep, all works as expected, all demands can be met

    EnergySystems.reset(demand_2)
    EnergySystems.reset(demand_1)
    EnergySystems.reset(bus_2)
    EnergySystems.reset(bus_3)
    EnergySystems.reset(bus_1)
    EnergySystems.reset(boiler)
    EnergySystems.reset(grid)

    EnergySystems.control(demand_2, components, simulation_parameters)
    EnergySystems.control(demand_1, components, simulation_parameters)
    EnergySystems.control(boiler, components, simulation_parameters)
    EnergySystems.control(bus_2, components, simulation_parameters)
    EnergySystems.control(bus_1, components, simulation_parameters)
    EnergySystems.control(bus_3, components, simulation_parameters)
    EnergySystems.control(grid, components, simulation_parameters)

    @test demand_1.input_interfaces[demand_1.medium].balance ≈ 0.0
    @test demand_1.input_interfaces[demand_1.medium].temperature === 60.0
    @test demand_2.input_interfaces[demand_2.medium].balance ≈ 0.0
    @test demand_2.input_interfaces[demand_2.medium].temperature === 60.0
    @test EnergySystems.balance(bus_1) ≈ 0.0
    @test EnergySystems.balance(bus_2) ≈ 0.0
    @test EnergySystems.balance(bus_3) ≈ 0.0

    # demand not processed yet --> balance is zero, but energy_potential not
    # input interfaces
    exchange = EnergySystems.balance_on(bus_1.input_interfaces[1], bus_1)
    @test exchange.balance ≈ 0.0
    @test exchange.storage_potential ≈ 0.0
    @test exchange.energy_potential ≈ -2000.0
    @test exchange.temperature === 60.0
    exchange = EnergySystems.balance_on(bus_2.input_interfaces[1], bus_2)
    @test exchange.balance ≈ 0.0
    @test exchange.storage_potential ≈ 0.0
    @test exchange.energy_potential ≈ -1000.0
    @test exchange.temperature === 60.0
    exchange = EnergySystems.balance_on(bus_3.input_interfaces[1], bus_3)
    @test exchange.balance ≈ 0.0
    @test exchange.storage_potential ≈ 0.0
    @test exchange.energy_potential ≈ -1000.0
    @test exchange.temperature === 60.0

    # output interfaces
    exchange = EnergySystems.balance_on(bus_2.input_interfaces[1], bus_1)
    @test exchange.balance ≈ 0.0
    @test exchange.storage_potential ≈ 0.0
    @test exchange.energy_potential ≈ 0.0 # gasboiler is not doing potential() step, therefore max_energy is not written
    @test exchange.temperature === nothing # is nothing as no temperature is written in interface and gasboilder has nothing as temperature
    exchange = EnergySystems.balance_on(bus_3.input_interfaces[1], bus_1)
    @test exchange.balance ≈ 0.0
    @test exchange.storage_potential ≈ 0.0
    @test exchange.energy_potential ≈ 0.0 # gasboiler is not doing potential() step, therefore max_energy is not written
    @test exchange.temperature ===  nothing # is nothing as no temperature is written in interface and gasboilder has nothing as temperature
    exchange = EnergySystems.balance_on(bus_2.output_interfaces[1], bus_2)
    @test exchange.balance ≈ 0.0
    @test exchange.storage_potential ≈ 0.0
    @test exchange.energy_potential ≈ 0.0 # gasboiler is not doing potential() step, therefore max_energy is not written
    @test exchange.temperature === 60.0
    exchange = EnergySystems.balance_on(bus_3.output_interfaces[1], bus_3)
    @test exchange.balance ≈ 0.0
    @test exchange.storage_potential ≈ 0.0
    @test exchange.energy_potential ≈ 0.0 # gasboiler is not doing potential() step, therefore max_energy is not written
    @test exchange.temperature === 60.0

    EnergySystems.process(demand_2, simulation_parameters)
    EnergySystems.process(demand_1, simulation_parameters)

    # demand already processed --> balance is not zero anymore, but energy_potential
    # input interfaces
    exchange = EnergySystems.balance_on(bus_1.input_interfaces[1], bus_1)
    @test exchange.balance ≈ -2000.0
    @test exchange.storage_potential ≈ 0.0
    @test exchange.energy_potential ≈ 0.0
    @test exchange.temperature === 60.0
    exchange = EnergySystems.balance_on(bus_2.input_interfaces[1], bus_2)
    @test exchange.balance ≈ -1000.0
    @test exchange.storage_potential ≈ 0.0
    @test exchange.energy_potential ≈ 0.0
    @test exchange.temperature === 60.0
    exchange = EnergySystems.balance_on(bus_3.input_interfaces[1], bus_3)
    @test exchange.balance ≈ -1000.0
    @test exchange.storage_potential ≈ 0.0 
    @test exchange.energy_potential ≈ 0.0
    @test exchange.temperature === 60.0

     # output interfaces
     exchange = EnergySystems.balance_on(bus_2.input_interfaces[1], bus_1)
     @test exchange.balance ≈ -2000.0
     @test exchange.storage_potential ≈ 0.0
     @test exchange.energy_potential ≈ 0.0 
     @test exchange.temperature === nothing # is nothing as no temperature is written in interface and gasboilder has nothing as temperature
     exchange = EnergySystems.balance_on(bus_3.input_interfaces[1], bus_1)
     @test exchange.balance ≈ -2000.0
     @test exchange.storage_potential ≈ 0.0
     @test exchange.energy_potential ≈ 0.0
     @test exchange.temperature === nothing # is nothing as no temperature is written in interface and gasboilder has nothing as temperature
     exchange = EnergySystems.balance_on(bus_2.output_interfaces[1], bus_2)
     @test exchange.balance ≈ -1000.0
     @test exchange.storage_potential ≈ 0.0
     @test exchange.energy_potential ≈ 0.0
     @test exchange.temperature === 60.0
     exchange = EnergySystems.balance_on(bus_3.output_interfaces[1], bus_3)
     @test exchange.balance ≈ -1000.0
     @test exchange.storage_potential ≈ 0.0
     @test exchange.energy_potential ≈ 0.0
     @test exchange.temperature === 60.0 

    @test demand_1.input_interfaces[demand_1.medium].balance ≈ -1000.0
    @test demand_1.input_interfaces[demand_1.medium].temperature === 60.0
    @test demand_2.input_interfaces[demand_2.medium].balance ≈ -1000.0
    @test demand_2.input_interfaces[demand_2.medium].temperature === 60.0
    @test EnergySystems.balance(bus_3) ≈ -1000.0
    @test EnergySystems.balance(bus_2) ≈ -1000.0
    @test EnergySystems.balance(bus_1) ≈ -2000.0

    EnergySystems.process(bus_2, simulation_parameters)
    EnergySystems.process(bus_1, simulation_parameters)
    EnergySystems.process(bus_3, simulation_parameters)
    EnergySystems.process(boiler, simulation_parameters)

    # gasboiler already processed --> balance and energy_potential is zero
    # input interfaces
    exchange = EnergySystems.balance_on(bus_1.input_interfaces[1], bus_1)
    @test exchange.balance ≈ 0.0
    @test exchange.storage_potential ≈ 0.0
    @test exchange.energy_potential ≈ 0.0
    @test exchange.temperature === 60.0
    exchange = EnergySystems.balance_on(bus_2.input_interfaces[1], bus_2)
    @test exchange.balance ≈ 0.0
    @test exchange.storage_potential ≈ 0.0
    @test exchange.energy_potential ≈ 0.0
    @test exchange.temperature === 60.0
    exchange = EnergySystems.balance_on(bus_3.input_interfaces[1], bus_3)
    @test exchange.balance ≈ 0.0
    @test exchange.storage_potential ≈ 0.0 
    @test exchange.energy_potential ≈ 0.0
    @test exchange.temperature === 60.0

     # output interfaces
     exchange = EnergySystems.balance_on(bus_2.input_interfaces[1], bus_1)
     @test exchange.balance ≈ 0.0
     @test exchange.storage_potential ≈ 0.0
     @test exchange.energy_potential ≈ 0.0 
     @test exchange.temperature === nothing # is nothing as no temperature is written in interface and gasboilder has nothing as temperature
     exchange = EnergySystems.balance_on(bus_3.input_interfaces[1], bus_1)
     @test exchange.balance ≈ 0.0
     @test exchange.storage_potential ≈ 0.0
     @test exchange.energy_potential ≈ 0.0 
     @test exchange.temperature === nothing # is nothing as no temperature is written in interface and gasboilder has nothing as temperature
     exchange = EnergySystems.balance_on(bus_2.output_interfaces[1], bus_2)
     @test exchange.balance ≈ 0.0
     @test exchange.storage_potential ≈ 0.0
     @test exchange.energy_potential ≈ 0.0 
     @test exchange.temperature === 60.0
     exchange = EnergySystems.balance_on(bus_3.output_interfaces[1], bus_3)
     @test exchange.balance ≈ 0.0
     @test exchange.storage_potential ≈ 0.0
     @test exchange.energy_potential ≈ 0.0 
     @test exchange.temperature === 60.0 

    @test EnergySystems.balance(bus_3) ≈ 0.0
    @test EnergySystems.balance(bus_2) ≈ 0.0
    @test EnergySystems.balance(bus_1) ≈ 0.0
    @test boiler.input_interfaces[boiler.m_fuel_in].balance ≈ -2000.0

    EnergySystems.process(grid, simulation_parameters)

    @test demand_1.input_interfaces[demand_1.medium].balance ≈ -1000.0
    @test demand_1.input_interfaces[demand_1.medium].temperature === 60.0
    @test demand_2.input_interfaces[demand_2.medium].balance ≈ -1000.0
    @test demand_2.input_interfaces[demand_2.medium].temperature === 60.0
    @test EnergySystems.balance(bus_3) ≈ 0.0
    @test EnergySystems.balance(bus_2) ≈ 0.0
    @test EnergySystems.balance(bus_1) ≈ 0.0
    @test grid.output_interfaces[boiler.m_fuel_in].balance ≈ 0.0

    EnergySystems.distribute!(bus_2)
    EnergySystems.distribute!(bus_3)
    EnergySystems.distribute!(bus_1)

    @test demand_1.input_interfaces[demand_1.medium].balance ≈ 0.0
    @test demand_1.input_interfaces[demand_1.medium].temperature === 60.0
    @test demand_1.input_interfaces[demand_1.medium].sum_abs_change ≈ 2000.0
    @test demand_2.input_interfaces[demand_2.medium].balance ≈ 0.0
    @test demand_2.input_interfaces[demand_2.medium].temperature === 60.0
    @test demand_2.input_interfaces[demand_2.medium].sum_abs_change ≈ 2000.0

    @test EnergySystems.balance(bus_1) ≈ 0.0
    @test EnergySystems.balance(bus_2) ≈ 0.0
    @test EnergySystems.balance(bus_3) ≈ 0.0
    @test bus_1.output_interfaces[1].sum_abs_change ≈ 2000.0
    @test bus_1.output_interfaces[2].sum_abs_change ≈ 2000.0

    # second timestep, there's not enough supply to meet demand, bus 2 has priority
    # over bus 3
    # not tests for balance_on() of busses here as for this function the missmatch of 
    # demand and supply is not important

    boiler.power = 6000

    EnergySystems.reset(demand_2)
    EnergySystems.reset(demand_1)
    EnergySystems.reset(bus_2)
    EnergySystems.reset(bus_3)
    EnergySystems.reset(bus_1)
    EnergySystems.reset(boiler)
    EnergySystems.reset(grid)

    EnergySystems.control(demand_2, components, simulation_parameters)
    EnergySystems.control(demand_1, components, simulation_parameters)
    EnergySystems.control(bus_2, components, simulation_parameters)
    EnergySystems.control(bus_1, components, simulation_parameters)
    EnergySystems.control(bus_3, components, simulation_parameters)
    EnergySystems.control(boiler, components, simulation_parameters)
    EnergySystems.control(grid, components, simulation_parameters)

    @test demand_1.input_interfaces[demand_1.medium].balance ≈ 0.0
    @test demand_1.input_interfaces[demand_1.medium].temperature === 60.0
    @test demand_2.input_interfaces[demand_2.medium].balance ≈ 0.0
    @test demand_2.input_interfaces[demand_2.medium].temperature === 60.0
    @test EnergySystems.balance(bus_1) ≈ 0.0
    @test EnergySystems.balance(bus_2) ≈ 0.0
    @test EnergySystems.balance(bus_3) ≈ 0.0

    EnergySystems.process(demand_2, simulation_parameters)
    EnergySystems.process(demand_1, simulation_parameters)

    @test demand_1.input_interfaces[demand_1.medium].balance ≈ -1000.0
    @test demand_1.input_interfaces[demand_1.medium].temperature === 60.0
    @test demand_2.input_interfaces[demand_2.medium].balance ≈ -1000.0
    @test demand_2.input_interfaces[demand_2.medium].temperature === 60.0
    @test EnergySystems.balance(bus_3) ≈ -1000.0
    @test EnergySystems.balance(bus_2) ≈ -1000.0
    @test EnergySystems.balance(bus_1) ≈ -2000.0

    EnergySystems.process(bus_2, simulation_parameters)
    EnergySystems.process(bus_1, simulation_parameters)
    EnergySystems.process(bus_3, simulation_parameters)
    EnergySystems.process(boiler, simulation_parameters)

    # busses don't consider output priority in the balance() function, so bus_2 also
    # thinks it has a negative balance even though it will later, in distribute(), be
    # prefered over bus_3 due to the output priorities
    @test EnergySystems.balance(bus_3) ≈ -500.0
    @test EnergySystems.balance(bus_2) ≈ -500.0
    @test EnergySystems.balance(bus_1) ≈ -500.0
    @test boiler.input_interfaces[boiler.m_fuel_in].balance ≈ -1500.0

    EnergySystems.process(grid, simulation_parameters)

    @test demand_1.input_interfaces[demand_1.medium].balance ≈ -1000.0
    @test demand_1.input_interfaces[demand_1.medium].temperature === 60.0
    @test demand_2.input_interfaces[demand_2.medium].balance ≈ -1000.0
    @test demand_2.input_interfaces[demand_2.medium].temperature === 60.0
    @test grid.output_interfaces[:m_c_g_natgas].balance ≈ 0.0

    EnergySystems.distribute!(bus_2)
    EnergySystems.distribute!(bus_3)
    EnergySystems.distribute!(bus_1)

    @test demand_1.input_interfaces[demand_1.medium].balance ≈ 0.0
    @test demand_1.input_interfaces[demand_1.medium].temperature === 60.0
    @test demand_1.input_interfaces[demand_1.medium].sum_abs_change ≈ 2000.0
    @test demand_2.input_interfaces[demand_2.medium].balance ≈ 0.0
    @test demand_2.input_interfaces[demand_2.medium].temperature === 60.0
    @test demand_2.input_interfaces[demand_2.medium].sum_abs_change ≈ 2000.0

    @test EnergySystems.balance(bus_1) ≈ -500.0
    @test EnergySystems.balance(bus_2) ≈ 0.0
    @test EnergySystems.balance(bus_3) ≈ 0.0
    @test bus_1.output_interfaces[1].balance ≈ 0.0
    @test bus_1.output_interfaces[1].sum_abs_change ≈ 2000.0
    @test bus_1.output_interfaces[2].balance ≈ -500.0
    @test bus_1.output_interfaces[2].sum_abs_change ≈ 1500.0
    @test bus_1.remainder ≈ 0.0
    @test bus_2.remainder ≈ 0.0
    @test bus_3.remainder ≈ 0.0
end

@testset "demand_over_busses_supply_is_transformer" begin
    test_demand_over_busses_supply_is_transformer()
end

function test_busses_communicate_storage_potential()
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
            "output_refs" => ["TST_BUS_02", "TST_BFT_01"],
            "connection_matrix" => Dict{String, Any}(
                "input_order" => [
                    "TST_BFT_01",
                    "TST_GRI_01"
                ],
                "output_order" => [
                    "TST_BUS_02",
                    "TST_BFT_01"
                ],
            )
        ),
        "TST_BFT_01" => Dict{String,Any}(
            "type" => "BufferTank",
            "control_refs" => [],
            "output_refs" => [
                "TST_BUS_01"
            ],
            "capacity" => 40000,
            "load" => 20000,
            "high_temperature" => 55.0
        ),
        "TST_BUS_02" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "output_refs" => ["TST_DEM_01", "TST_BFT_02"],
            "connection_matrix" => Dict{String, Any}(
                "input_order" => [
                    "TST_BFT_02",
                    "TST_BUS_01"
                ],
                "output_order" => [
                    "TST_DEM_01",
                    "TST_BFT_02"
                ],
            )
        ),
        "TST_BFT_02" => Dict{String,Any}(
            "type" => "BufferTank",
            "control_refs" => [],
            "output_refs" => [
                "TST_BUS_02"
            ],
            "capacity" => 20000,
            "load" => 10000,
            "high_temperature" => 55.0
        ),
        "TST_DEM_01" => Dict{String,Any}(
            "type" => "Demand",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "output_refs" => [],
            "energy_profile_file_path" => "./profiles/tests/demand_heating_energy.prf",
            "temperature_profile_file_path" => "./profiles/tests/demand_heating_temperature.prf",
            "scale" => 1000
        ),
    )
    components = Resie.load_components(components_config)
    demand = components["TST_DEM_01"]
    grid = components["TST_GRI_01"]
    bus_1 = components["TST_BUS_01"]
    bus_2 = components["TST_BUS_02"]
    tank_1 = components["TST_BFT_01"]
    tank_2 = components["TST_BFT_02"]

    simulation_parameters = Dict{String,Any}(
        "time_step_seconds" => 900,
        "time" => 0,
        "epsilon" => 1e-9
    )

    EnergySystems.reset(demand)
    EnergySystems.reset(bus_2)
    EnergySystems.reset(bus_1)
    EnergySystems.reset(tank_2)
    EnergySystems.reset(tank_1)
    EnergySystems.reset(grid)

    EnergySystems.control(demand, components, simulation_parameters)
    EnergySystems.control(bus_2, components, simulation_parameters)
    EnergySystems.control(tank_2, components, simulation_parameters)
    EnergySystems.control(bus_1, components, simulation_parameters)
    EnergySystems.control(tank_1, components, simulation_parameters)
    EnergySystems.control(grid, components, simulation_parameters)

    # demand not processed yet --> balance is zero, but energy_potential not
    # input interfaces
    exchange = EnergySystems.balance_on(grid.output_interfaces[grid.medium], bus_1)
    @test exchange.balance ≈ 0.0
    @test exchange.storage_potential ≈ -30000
    @test exchange.energy_potential ≈ -75.0
    @test exchange.temperature === 55.0
    interface_between_bus1_bus2 = bus_2.input_interfaces[1].source.uac == "TST_BUS_01" ? bus_2.input_interfaces[1] : bus_2.input_interfaces[2]
    exchange = EnergySystems.balance_on(interface_between_bus1_bus2, bus_2)
    @test exchange.balance ≈ 0.0
    @test exchange.storage_potential ≈ -10000
    @test exchange.energy_potential ≈ -75.0
    @test exchange.temperature === 55.0

    # output interfaces
    exchange = EnergySystems.balance_on(interface_between_bus1_bus2, bus_1)
    @test exchange.balance ≈ 0.0
    @test exchange.storage_potential ≈ 20000
    @test exchange.energy_potential ≈ Inf
    @test exchange.temperature === 55.0
    exchange = EnergySystems.balance_on(demand.input_interfaces[demand.medium], bus_2)
    @test exchange.balance ≈ 0.0
    @test exchange.storage_potential ≈ 30000
    @test exchange.energy_potential ≈ Inf
    @test exchange.temperature === 55.0

    EnergySystems.process(demand, simulation_parameters)

    exchange = EnergySystems.balance_on(
        tank_2.output_interfaces[tank_2.medium], bus_2
    )
    @test exchange.balance ≈ -75.0
    @test exchange.storage_potential ≈ 0.0
    @test exchange.energy_potential ≈ 0.0
    @test exchange.temperature === 55.0

    exchange = EnergySystems.balance_on(
        tank_1.output_interfaces[tank_1.medium], bus_1
    )
    @test exchange.balance ≈ -75.0
    @test exchange.storage_potential ≈ -10000.0
    @test exchange.energy_potential ≈ 0.0
    @test exchange.temperature === 55.0

    EnergySystems.process(bus_2, simulation_parameters)
    EnergySystems.process(tank_2, simulation_parameters)
    EnergySystems.process(bus_1, simulation_parameters)
    EnergySystems.process(tank_1, simulation_parameters)
    EnergySystems.load(tank_2, simulation_parameters)
    EnergySystems.load(tank_1, simulation_parameters)
    EnergySystems.process(grid, simulation_parameters)
    EnergySystems.distribute!(bus_2)
    EnergySystems.distribute!(bus_1)

    # a peculiar thing happens after distribute! has been called on a bus: when a system
    # interface that previously held a demand at a not-nothing temperature has been matched
    # by a corresponding supply, calling balance_on on the bus now returns a temperature of
    # nothing, even if the system interface still has the same temperature. this happens
    # because a balance of 0 is not considered a demand and is not considered for "the
    # highest demand temperature on the bus". this behaviour is not wrong, but unintuitive

    exchange = EnergySystems.balance_on(
        tank_2.output_interfaces[tank_2.medium], bus_2
    )
    @test exchange.balance ≈ 0.0
    @test exchange.storage_potential ≈ 0.0
    @test exchange.energy_potential ≈ 0.0
    @test exchange.temperature === 55.0

    exchange = EnergySystems.balance_on(
        tank_1.output_interfaces[tank_1.medium], bus_1
    )
    @test exchange.balance ≈ 0.0
    @test exchange.storage_potential ≈ -10075.0
    @test exchange.energy_potential ≈ 0.0
    @test exchange.temperature === 55.0
end

@testset "busses_communicate_storage_potential" begin
    test_busses_communicate_storage_potential()
end