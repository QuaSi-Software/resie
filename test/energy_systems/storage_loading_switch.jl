using Debugger
using Test
using Resie
using Resie.EnergySystems
using Resie.Profiles

watt_to_wh = function (watts::Float64)
    watts * 900 / 3600.0
end

function test_primary_producer_can_load_storage()
    systems_config = Dict{String,Any}(
        "TST_GRI_01" => Dict{String,Any}(
            "type" => "GridConnection",
            "medium" => "m_c_g_natgas",
            "control_refs" => [],
            "production_refs" => ["TST_GBO_01"],
            "is_source" => true,
        ),
        "TST_GRI_02" => Dict{String,Any}(
            "type" => "GridConnection",
            "medium" => "m_c_g_natgas",
            "control_refs" => [],
            "production_refs" => ["TST_GBO_02"],
            "is_source" => true,
        ),
        "TST_GBO_01" => Dict{String,Any}(
            "type" => "GasBoiler",
            "control_refs" => ["TST_BFT_01"],
            "production_refs" => [
                "TST_BUS_01"
            ],
            "strategy" => Dict{String,Any}(
                "name" => "storage_driven",
                "high_threshold" => 0.5,
                "low_threshold" => 0.1
            ),
            "power" => 10000
        ),
        "TST_GBO_02" => Dict{String,Any}(
            "type" => "GasBoiler",
            "control_refs" => ["TST_BFT_01"],
            "production_refs" => [
                "TST_BUS_01"
            ],
            "strategy" => Dict{String,Any}(
                "name" => "demand_driven"
            ),
            "power" => 40000
        ),
        "TST_BUS_01" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "production_refs" => ["TST_DEM_01", "TST_BFT_01"],
            "connection_matrix" => Dict{String, Any}(
                "input_order" => [
                    "TST_GBO_01",
                    "TST_BFT_01",
                    "TST_GBO_02"
                ],
                "output_order" => [
                    "TST_DEM_01",
                    "TST_BFT_01"
                ],
                "storage_loading" => [
                    [1, 1],
                    [1, 0],
                    [1, 0]
                ]
            )
        ),
        "TST_BFT_01" => Dict{String,Any}(
            "type" => "BufferTank",
            "control_refs" => [],
            "production_refs" => [
                "TST_BUS_01"
            ],
            "capacity" => 40000,
            "load" => 0
        ),
        "TST_DEM_01" => Dict{String,Any}(
            "type" => "Demand",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "production_refs" => [],
            "energy_profile_file_path" => "./profiles/tests/demand_heating_energy.prf",
            "temperature_profile_file_path" => "./profiles/tests/demand_heating_temperature.prf",
            "scale" => 1,
            "static_load" => 5000,
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

    simulation_parameters = Dict{String,Any}(
        "time_step_seconds" => 900,
        "time" => 0,
    )

    # first timestep, demand is higher than primary producer can provide, so both should
    # operate. however as the secondary should not load the storage, it only covers the
    # remaining demand

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

    @test boiler_1.controller.state_machine.state == 2
    @test boiler_2.controller.state_machine.state == 1

    EnergySystems.produce(demand, simulation_parameters, watt_to_wh)
    EnergySystems.produce(bus, simulation_parameters, watt_to_wh)

    EnergySystems.produce(boiler_1, simulation_parameters, watt_to_wh)
    @test boiler_1.output_interfaces[boiler_1.m_heat_out].balance == 2500.0

    EnergySystems.produce(tank, simulation_parameters, watt_to_wh)
    @test tank.output_interfaces[tank.medium].sum_abs_change == 0.0

    EnergySystems.produce(boiler_2, simulation_parameters, watt_to_wh)
    @test boiler_2.output_interfaces[boiler_2.m_heat_out].balance == 2500.0

    EnergySystems.load(tank, simulation_parameters, watt_to_wh)
    @test tank.input_interfaces[tank.medium].sum_abs_change == 0.0

    EnergySystems.produce(grid_2, simulation_parameters, watt_to_wh)
    EnergySystems.produce(grid_1, simulation_parameters, watt_to_wh)
    EnergySystems.distribute!(bus)

    InterfaceInfo = EnergySystems.balance_on(
        boiler_1.output_interfaces[boiler_1.m_heat_out], bus
    )
    @test InterfaceInfo.balance == 0.0
    @test InterfaceInfo.storage_potential == -40000.0

    InterfaceInfo = EnergySystems.balance_on(
        boiler_2.output_interfaces[boiler_2.m_heat_out], bus
    )
    @test InterfaceInfo.balance == 0.0
    @test InterfaceInfo.storage_potential == 0.0

    # in the second timestep, demand is lower than primary producer can provide, so the
    # excess can go into storage. because the secondary should not load the storage, it
    # should be inactive during this time step

    EnergySystems.reset(grid_1)
    EnergySystems.reset(grid_2)
    EnergySystems.reset(boiler_1)
    EnergySystems.reset(boiler_2)
    EnergySystems.reset(bus)
    EnergySystems.reset(tank)
    EnergySystems.reset(demand)

    demand.static_load = 1500

    EnergySystems.control(demand, systems, simulation_parameters)
    EnergySystems.control(bus, systems, simulation_parameters)
    EnergySystems.control(boiler_1, systems, simulation_parameters)
    EnergySystems.control(tank, systems, simulation_parameters)
    EnergySystems.control(boiler_2, systems, simulation_parameters)
    EnergySystems.control(grid_1, systems, simulation_parameters)
    EnergySystems.control(grid_2, systems, simulation_parameters)

    @test boiler_1.controller.state_machine.state == 2
    @test boiler_2.controller.state_machine.state == 1

    EnergySystems.produce(demand, simulation_parameters, watt_to_wh)
    EnergySystems.produce(bus, simulation_parameters, watt_to_wh)

    EnergySystems.produce(boiler_1, simulation_parameters, watt_to_wh)
    @test boiler_1.output_interfaces[boiler_1.m_heat_out].balance == 2500.0

    EnergySystems.produce(tank, simulation_parameters, watt_to_wh)
    @test tank.output_interfaces[tank.medium].sum_abs_change == 0.0

    EnergySystems.produce(boiler_2, simulation_parameters, watt_to_wh)
    @test boiler_2.output_interfaces[boiler_2.m_heat_out].sum_abs_change == 0.0

    EnergySystems.load(tank, simulation_parameters, watt_to_wh)
    @test tank.load == 1000.0

    EnergySystems.produce(grid_2, simulation_parameters, watt_to_wh)
    EnergySystems.produce(grid_1, simulation_parameters, watt_to_wh)
    EnergySystems.distribute!(bus)

    @test tank.input_interfaces[tank.medium].sum_abs_change == 2000.0
    @test boiler_1.output_interfaces[boiler_2.m_heat_out].sum_abs_change == 5000.0

    InterfaceInfo = EnergySystems.balance_on(
        boiler_1.output_interfaces[boiler_1.m_heat_out], bus
    )
    @test InterfaceInfo.balance == 0.0
    @test InterfaceInfo.storage_potential == -39000.0

    InterfaceInfo = EnergySystems.balance_on(
        boiler_2.output_interfaces[boiler_2.m_heat_out], bus
    )
    @test InterfaceInfo.balance == 0.0
    @test InterfaceInfo.storage_potential == 0.0
end

@testset "primary_producer_can_load_storage" begin
    test_primary_producer_can_load_storage()
end