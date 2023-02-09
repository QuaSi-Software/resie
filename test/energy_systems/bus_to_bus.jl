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