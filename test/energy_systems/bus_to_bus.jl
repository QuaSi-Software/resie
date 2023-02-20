watt_to_wh = function (watts :: Float64)
    watts * 900 / 3600.0
end

@testset "busses_communicate_demand" begin
    systems_config = Dict{String, Any}(
        "TST_GRI_01" => Dict{String, Any}(
            "type" => "GridConnection",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "production_refs" => ["TST_BUS_01"],
            "is_source" => true,
        ),
        "TST_BUS_01" => Dict{String, Any}(
            "type" => "Bus",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "production_refs" => ["TST_BUS_02"],
            "input_priorities" => ["TST_GRI_01"]
        ),
        "TST_BUS_02" => Dict{String, Any}(
            "type" => "Bus",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "production_refs" => ["TST_DEM_01"],
            "input_priorities" => ["TST_BUS_01"]
        ),
        "TST_DEM_01" => Dict{String, Any}(
            "type" => "Demand",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "production_refs" => [],
            "energy_profile_file_path" => "./profiles/tests/demand_heating_energy.prf",
            "temperature_profile_file_path" => "./profiles/tests/demand_heating_temperature.prf",
            "scale" => 1000
        ),
    )
    systems = Resie.load_systems(systems_config)
    demand = systems["TST_DEM_01"]
    grid = systems["TST_GRI_01"]
    bus_1 = systems["TST_BUS_01"]
    bus_2 = systems["TST_BUS_02"]

    simulation_parameters = Dict{String, Any}(
        "time_step_seconds" => 900,
        "time" => 0,
    )

    EnergySystems.reset(demand)
    EnergySystems.reset(grid)
    EnergySystems.reset(bus_1)
    EnergySystems.reset(bus_2)

    EnergySystems.control(demand, systems, simulation_parameters)
    EnergySystems.control(grid, systems, simulation_parameters)
    EnergySystems.control(bus_1, systems, simulation_parameters)
    EnergySystems.control(bus_2, systems, simulation_parameters)

    @test demand.input_interfaces[EnergySystems.m_h_w_ht1].balance == 0.0
    @test demand.input_interfaces[EnergySystems.m_h_w_ht1].temperature == 55.0
    @test EnergySystems.balance(bus_1) == 0.0
    @test EnergySystems.balance(bus_2) == 0.0
    @test grid.output_interfaces[EnergySystems.m_h_w_ht1].balance == 0.0

    EnergySystems.produce(demand, simulation_parameters, watt_to_wh)

    @test demand.input_interfaces[EnergySystems.m_h_w_ht1].balance == -75.0
    @test demand.input_interfaces[EnergySystems.m_h_w_ht1].temperature == 55.0
    @test EnergySystems.balance(bus_1) == -75.0
    @test EnergySystems.balance(bus_2) == -75.0
    @test grid.output_interfaces[EnergySystems.m_h_w_ht1].balance == 0.0

    EnergySystems.produce(bus_2, simulation_parameters, watt_to_wh)
    EnergySystems.produce(bus_1, simulation_parameters, watt_to_wh)
    EnergySystems.produce(grid, simulation_parameters, watt_to_wh)

    @test demand.input_interfaces[EnergySystems.m_h_w_ht1].balance == -75.0
    @test demand.input_interfaces[EnergySystems.m_h_w_ht1].temperature == 55.0
    @test EnergySystems.balance(bus_1) == 0.0
    @test EnergySystems.balance(bus_2) == 0.0
    @test bus_1.remainder == 0.0
    @test bus_2.remainder == 0.0
    @test grid.output_interfaces[EnergySystems.m_h_w_ht1].balance == 75.0

    EnergySystems.distribute!(bus_2)
    EnergySystems.distribute!(bus_1)

    @test demand.input_interfaces[EnergySystems.m_h_w_ht1].balance == 0.0
    @test demand.input_interfaces[EnergySystems.m_h_w_ht1].temperature == 55.0
    @test demand.input_interfaces[EnergySystems.m_h_w_ht1].sum_abs_change == 150.0
    @test EnergySystems.balance(bus_1) == 0.0
    @test EnergySystems.balance(bus_2) == 0.0
    @test bus_1.remainder == 75.0
    @test bus_2.remainder == -75.0
    @test grid.output_interfaces[EnergySystems.m_h_w_ht1].balance == 0.0
    @test grid.output_interfaces[EnergySystems.m_h_w_ht1].sum_abs_change == 150.0
end

@testset "demand_over_busses_supply_is_transformer" begin
    systems_config = Dict{String, Any}(
        "TST_GRI_01" => Dict{String, Any}(
            "type" => "GridConnection",
            "medium" => "m_c_g_natgas",
            "control_refs" => [],
            "production_refs" => ["TST_GBO_01"],
            "is_source" => true,
        ),
        "TST_GBO_01" => Dict{String, Any}(
            "type" => "GasBoiler",
            "control_refs" => ["TST_BUS_01"],
            "production_refs" => ["TST_BUS_01"],
            "strategy" => {
                "name" => "demand_driven",
            },
            "power" => 10000
        ),
        "TST_BUS_01" => Dict{String, Any}(
            "type" => "Bus",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "production_refs" => ["TST_BUS_02"],
            "input_priorities" => ["TST_GRI_01"]
        ),
        "TST_BUS_02" => Dict{String, Any}(
            "type" => "Bus",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "production_refs" => ["TST_DEM_01"],
            "input_priorities" => ["TST_BUS_01"]
        ),
        "TST_DEM_01" => Dict{String, Any}(
            "type" => "Demand",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "production_refs" => [],
            "energy_profile_file_path" => "./profiles/tests/demand_heating_energy.prf",
            "temperature_profile_file_path" => "./profiles/tests/demand_heating_temperature.prf",
            "scale" => 1000
        ),
    )
    systems = Resie.load_systems(systems_config)
    demand = systems["TST_DEM_01"]
    grid = systems["TST_GRI_01"]
    bus_1 = systems["TST_BUS_01"]
    bus_2 = systems["TST_BUS_02"]

    simulation_parameters = Dict{String, Any}(
        "time_step_seconds" => 900,
        "time" => 0,
    )

    EnergySystems.reset(demand)
    EnergySystems.reset(grid)
    EnergySystems.reset(bus_1)
    EnergySystems.reset(bus_2)

    EnergySystems.control(demand, systems, simulation_parameters)
    EnergySystems.control(grid, systems, simulation_parameters)
    EnergySystems.control(bus_1, systems, simulation_parameters)
    EnergySystems.control(bus_2, systems, simulation_parameters)

    @test demand.input_interfaces[EnergySystems.m_h_w_ht1].balance == 0.0
    @test demand.input_interfaces[EnergySystems.m_h_w_ht1].temperature == 55.0
    @test EnergySystems.balance(bus_1) == 0.0
    @test EnergySystems.balance(bus_2) == 0.0
    @test grid.output_interfaces[EnergySystems.m_h_w_ht1].balance == 0.0

    EnergySystems.produce(demand, simulation_parameters, watt_to_wh)

    @test demand.input_interfaces[EnergySystems.m_h_w_ht1].balance == -75.0
    @test demand.input_interfaces[EnergySystems.m_h_w_ht1].temperature == 55.0
    @test EnergySystems.balance(bus_1) == -75.0
    @test EnergySystems.balance(bus_2) == -75.0
    @test grid.output_interfaces[EnergySystems.m_h_w_ht1].balance == 0.0

    EnergySystems.produce(bus_2, simulation_parameters, watt_to_wh)
    EnergySystems.produce(bus_1, simulation_parameters, watt_to_wh)
    EnergySystems.produce(grid, simulation_parameters, watt_to_wh)

    @test demand.input_interfaces[EnergySystems.m_h_w_ht1].balance == -75.0
    @test demand.input_interfaces[EnergySystems.m_h_w_ht1].temperature == 55.0
    @test EnergySystems.balance(bus_1) == 0.0
    @test EnergySystems.balance(bus_2) == 0.0
    @test bus_1.remainder == 0.0
    @test bus_2.remainder == 0.0
    @test grid.output_interfaces[EnergySystems.m_h_w_ht1].balance == 75.0

    EnergySystems.distribute!(bus_2)
    EnergySystems.distribute!(bus_1)

    @test demand.input_interfaces[EnergySystems.m_h_w_ht1].balance == 0.0
    @test demand.input_interfaces[EnergySystems.m_h_w_ht1].temperature == 55.0
    @test demand.input_interfaces[EnergySystems.m_h_w_ht1].sum_abs_change == 150.0
    @test EnergySystems.balance(bus_1) == 0.0
    @test EnergySystems.balance(bus_2) == 0.0
    @test bus_1.remainder == 75.0
    @test bus_2.remainder == -75.0
    @test grid.output_interfaces[EnergySystems.m_h_w_ht1].balance == 0.0
    @test grid.output_interfaces[EnergySystems.m_h_w_ht1].sum_abs_change == 150.0
end

@testset "busses_communicate_storage_potential" begin
    systems_config = Dict{String, Any}(
        "TST_GRI_01" => Dict{String, Any}(
            "type" => "GridConnection",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "production_refs" => ["TST_BUS_01"],
            "is_source" => true,
        ),
        "TST_BUS_01" => Dict{String, Any}(
            "type" => "Bus",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "production_refs" => ["TST_BUS_02", "TST_BFT_01"],
            "input_priorities" => ["TST_BFT_01", "TST_GRI_01"]
        ),
        "TST_BFT_01" => Dict{String, Any}(
            "type" => "BufferTank",
            "control_refs" => [],
            "production_refs" => [
                "TST_BUS_01"
            ],
            "capacity" => 40000,
            "load" => 20000
        ),
        "TST_BUS_02" => Dict{String, Any}(
            "type" => "Bus",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "production_refs" => ["TST_DEM_01", "TST_BFT_02"],
            "input_priorities" => ["TST_BFT_02", "TST_BUS_01"]
        ),
        "TST_BFT_02" => Dict{String, Any}(
            "type" => "BufferTank",
            "control_refs" => [],
            "production_refs" => [
                "TST_BUS_02"
            ],
            "capacity" => 20000,
            "load" => 10000
        ),
        "TST_DEM_01" => Dict{String, Any}(
            "type" => "Demand",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "production_refs" => [],
            "energy_profile_file_path" => "./profiles/tests/demand_heating_energy.prf",
            "temperature_profile_file_path" => "./profiles/tests/demand_heating_temperature.prf",
            "scale" => 1000
        ),
    )
    systems = Resie.load_systems(systems_config)
    demand = systems["TST_DEM_01"]
    grid = systems["TST_GRI_01"]
    bus_1 = systems["TST_BUS_01"]
    bus_2 = systems["TST_BUS_02"]
    tank_1 = systems["TST_BFT_01"]
    tank_2 = systems["TST_BFT_02"]

    simulation_parameters = Dict{String, Any}(
        "time_step_seconds" => 900,
        "time" => 0,
    )

    EnergySystems.reset(demand)
    EnergySystems.reset(bus_2)
    EnergySystems.reset(bus_1)
    EnergySystems.reset(tank_2)
    EnergySystems.reset(tank_1)
    EnergySystems.reset(grid)

    EnergySystems.control(demand, systems, simulation_parameters)
    EnergySystems.control(bus_2, systems, simulation_parameters)
    EnergySystems.control(tank_2, systems, simulation_parameters)
    EnergySystems.control(bus_1, systems, simulation_parameters)
    EnergySystems.control(tank_1, systems, simulation_parameters)
    EnergySystems.control(grid, systems, simulation_parameters)

    EnergySystems.produce(demand, simulation_parameters, watt_to_wh)

    balance, potential, temperature = EnergySystems.balance_on(
        tank_2.output_interfaces[EnergySystems.m_h_w_ht1], bus_2
    )
    @test balance == -75.0
    @test potential == -10000.0
    @test temperature == 55.0

    balance, potential, temperature = EnergySystems.balance_on(
        tank_1.output_interfaces[EnergySystems.m_h_w_ht1], bus_1
    )
    @test balance == -75.0
    @test potential == -30000.0
    @test temperature == 55.0

    EnergySystems.produce(bus_2, simulation_parameters, watt_to_wh)
    EnergySystems.produce(tank_2, simulation_parameters, watt_to_wh)
    EnergySystems.produce(bus_1, simulation_parameters, watt_to_wh)
    EnergySystems.produce(tank_1, simulation_parameters, watt_to_wh)
    EnergySystems.load(tank_2, simulation_parameters, watt_to_wh)
    EnergySystems.load(tank_1, simulation_parameters, watt_to_wh)
    EnergySystems.produce(grid, simulation_parameters, watt_to_wh)
    EnergySystems.distribute!(bus_2)
    EnergySystems.distribute!(bus_1)

    # a peculiar thing happens after distribute! has been called on a bus: when a system
    # interface that previously held a demand at a not-nothing temperature has been matched
    # by a corresponding supply, calling balance_on on the bus now returns a temperature of
    # nothing, even if the system interface still has the same temperature. this happens
    # because a balance of 0 is not considered a demand and is not considered for "the
    # highest demand temperature on the bus". this behaviour is not wrong, but unintuitive

    balance, potential, temperature = EnergySystems.balance_on(
        tank_2.output_interfaces[EnergySystems.m_h_w_ht1], bus_2
    )
    @test balance == 0.0
    @test potential == -10075.0
    @test temperature === nothing

    balance, potential, temperature = EnergySystems.balance_on(
        tank_1.output_interfaces[EnergySystems.m_h_w_ht1], bus_1
    )
    @test balance == 0.0
    @test potential == -30075.0
    @test temperature === nothing
end