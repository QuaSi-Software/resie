using Debugger
using Test
using Resie.EnergySystems
using Resie.Profiles

watt_to_wh = function (watts :: Float64)
    watts * 900 / 3600.0
end

function test_primary_producer_can_load_storage()
    systems_config = Dict{String, Any}(
        "TST_GRI_01" => Dict{String, Any}(
            "type" => "GridConnection",
            "medium" => "m_c_g_natgas",
            "control_refs" => [],
            "production_refs" => ["TST_GBO_01"],
            "is_source" => true,
        ),
        "TST_GRI_02" => Dict{String, Any}(
            "type" => "GridConnection",
            "medium" => "m_c_g_natgas",
            "control_refs" => [],
            "production_refs" => ["TST_GBO_02"],
            "is_source" => true,
        ),
        "TST_GBO_01" => Dict{String, Any}(
            "type" => "GasBoiler",
            "control_refs" => ["TST_BFT_01"],
            "production_refs" => [
                "TST_BUS_01"
            ],
            "strategy" => {
                "name" => "storage_driven",
                "high_threshold" => 0.5,
                "low_threshold" => 0.1
            },
            "power" => 10000
        ),
        "TST_GBO_02" => Dict{String, Any}(
            "type" => "GasBoiler",
            "control_refs" => ["TST_BFT_01"],
            "production_refs" => [
                "TST_BUS_01"
            ],
            "strategy" => Dict{String, Any}(
                "name" => "demand_driven"
            ),
            "power" => 40000
        ),
        "TST_BUS_01" => Dict{String, Any}(
            "type" => "Bus",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "production_refs" => ["TST_DEM_01", "TST_BFT_01"],
            "input_priorities" => ["TST_GBO_01", "TST_BFT_01", "TST_GBO_02"]
        ),
        "TST_BFT_01" => Dict{String, Any}(
            "type" => "BufferTank",
            "control_refs" => [],
            "production_refs" => [
                "TST_BUS_01"
            ],
            "capacity" => 40000,
            "load" => 0
        ),
        "TST_DEM_01" => Dict{String, Any}(
            "type" => "Demand",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "production_refs" => [],
            "energy_profile_file_path" => "./profiles/tests/demand_heating_energy.prf",
            "temperature_profile_file_path" => "./profiles/tests/demand_heating_temperature.prf",
            "scale" => 1000,
            "static_load" => 1000,
            "static_temperature" => 60
        ),
    )

    systems = Resie.load_systems(systems_config)
    demand = systems["TST_DEM_01"]
    grid_1 = systems["TST_GRI_01"]
    grid_2 = systems["TST_GRI_02"]
    bus = systems["TST_BUS_01"]
    tank = systems["TST_BFT_01"]
    boiler_1 = systems["TST_GBO_01"]
    boiler_2 = systems["TST_GBO_02"]

    simulation_parameters = Dict{String, Any}(
        "time_step_seconds" => 900,
        "time" => 0,
    )

    EnergySystems.reset(grid_1)
    EnergySystems.reset(grid_2)
    EnergySystems.reset(boiler_1)
    EnergySystems.reset(boiler_2)
    EnergySystems.reset(bus)
    EnergySystems.reset(tank)
    EnergySystems.reset(demand)

    EnergySystems.control(demand, systems, simulation_parameters)
    EnergySystems.control(bus, systems, simulation_parameters)
    EnergySystems.control(boiler_1, systems, simulation_parameters)
    EnergySystems.control(tank, systems, simulation_parameters)
    EnergySystems.control(boiler_2, systems, simulation_parameters)
    EnergySystems.control(grid_1, systems, simulation_parameters)
    EnergySystems.control(grid_2, systems, simulation_parameters)

    demand.input_interfaces[EnergySystems.m_h_w_ht1].bala

    @test boiler_1.controller.state_machine.state == 2
    @test boiler_2.controller.state_machine.state == 1

    EnergySystems.produce(demand, systems, simulation_parameters)
    EnergySystems.produce(bus, systems, simulation_parameters)
    EnergySystems.produce(boiler_1, systems, simulation_parameters)
    EnergySystems.produce(tank, systems, simulation_parameters)
    EnergySystems.produce(boiler_2, systems, simulation_parameters)
    EnergySystems.load(tank, simulation_parameters, watt_to_wh)
    EnergySystems.produce(grid_1, systems, simulation_parameters)
    EnergySystems.produce(grid_2, systems, simulation_parameters)
    EnergySystems.distribute!(bus)

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

@testset "primary_producer_can_load_storage" begin
    test_primary_producer_can_load_storage()
end