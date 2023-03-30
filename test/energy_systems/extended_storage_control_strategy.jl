using Debugger
using Test
using Resie
using Resie.EnergySystems
using Resie.Profiles

watt_to_wh = function (watts::Float64)
    watts * 900 / 3600.0
end

function test_extended_storage_control_strategy_allow_loading_by_storage()
    systems_config = Dict{String,Any}(
        "TST_GRI_01" => Dict{String,Any}(
            "type" => "GridConnection",
            "medium" => "m_c_g_natgas",
            "control_refs" => [],
            "production_refs" => ["TST_GBO_01"],
            "is_source" => true,
        ),
        "TST_GBO_01" => Dict{String,Any}(
            "type" => "GasBoiler",
            "control_refs" => [],
            "production_refs" => [
                "TST_BUS_01"
            ],
            "strategy" => Dict{String,Any}(
                "name" => "demand_driven",
                "load_storages" => false
            ),
            "power" => 10000
        ),
        "TST_BUS_01" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "production_refs" => ["TST_BUS_02", "TST_BFT_01"],
            "connection_matrix" => Dict{String, Any}(
                "input_order" => [
                    "TST_GBO_01",
                    "TST_BFT_01"
                ],
                "output_order" => [
                    "TST_BUS_02",
                    "TST_BFT_01"
                ]
            )
        ),
        "TST_BUS_02" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "production_refs" => ["TST_DEM_01", "TST_BFT_02"],
            "connection_matrix" => Dict{String, Any}(
                "input_order" => [
                    "TST_BUS_01",
                    "TST_BFT_02"
                ],
                "output_order" => [
                    "TST_DEM_01",
                    "TST_BFT_02"
                ]
            )
        ),
        "TST_BFT_01" => Dict{String,Any}(
            "type" => "BufferTank",
            "control_refs" => [],
            "production_refs" => [
                "TST_BUS_01"
            ],
            "strategy" => Dict{String,Any}(
                "name" => "extended_storage_control",
                "load_any_storage" => true
            ),
            "capacity" => 40000,
            "load" => 20000
        ),
        "TST_BFT_02" => Dict{String,Any}(
            "type" => "BufferTank",
            "control_refs" => [],
            "production_refs" => [
                "TST_BUS_02"
            ],
            "capacity" => 40000,
            "load" => 20000
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
    grid = systems["TST_GRI_01"]
    bus_1 = systems["TST_BUS_01"]
    bus_2 = systems["TST_BUS_02"]
    storage_1 = systems["TST_BFT_01"]
    storage_2 = systems["TST_BFT_02"]
    boiler = systems["TST_GBO_01"]

    simulation_parameters = Dict{String,Any}(
        "time_step_seconds" => 900,
        "time" => 0,
        "epsilon" => 1e-9
    )

    # Gasboiler IS NOT allowed to load storages, but storage_1 IS allowed to load storage_2.
    # As result the gasboilder should be provide max energy, the rest of the demand should be 
    # supplied by storage_2 which will be refilled by storage_1

    EnergySystems.reset(grid)
    EnergySystems.reset(bus_1)
    EnergySystems.reset(bus_2)
    EnergySystems.reset(storage_1)
    EnergySystems.reset(storage_2)
    EnergySystems.reset(boiler)
    EnergySystems.reset(demand)

    EnergySystems.control(demand, systems, simulation_parameters)
    EnergySystems.control(bus_2, systems, simulation_parameters)
    EnergySystems.control(bus_1, systems, simulation_parameters)
    EnergySystems.control(boiler, systems, simulation_parameters)
    EnergySystems.control(storage_2, systems, simulation_parameters)
    EnergySystems.control(storage_1, systems, simulation_parameters)
    EnergySystems.control(grid, systems, simulation_parameters)

    EnergySystems.produce(demand, simulation_parameters, watt_to_wh)
    EnergySystems.produce(bus_2, simulation_parameters, watt_to_wh)
    EnergySystems.produce(bus_1, simulation_parameters, watt_to_wh)

    EnergySystems.produce(boiler, simulation_parameters, watt_to_wh)
    @test boiler.output_interfaces[boiler.m_heat_out].balance == 2500.0

    EnergySystems.produce(storage_2, simulation_parameters, watt_to_wh)
    @test storage_2.output_interfaces[storage_2.medium].balance == 2500.0
    @test storage_2.load == 20000.0 - 2500.0

    EnergySystems.produce(storage_1, simulation_parameters, watt_to_wh)
    @test storage_1.output_interfaces[storage_1.medium].balance == 20000.0
    @test storage_1.load == 0

    EnergySystems.load(storage_2, simulation_parameters, watt_to_wh)
    @test storage_2.input_interfaces[storage_2.medium].balance == -20000.0
    @test storage_2.load == 40000.0 - 2500.0

    EnergySystems.load(storage_1, simulation_parameters, watt_to_wh)
    @test storage_1.input_interfaces[storage_1.medium].balance == 0.0
    @test storage_1.load == 0

    EnergySystems.produce(grid, simulation_parameters, watt_to_wh)
    EnergySystems.distribute!(bus_2)
    EnergySystems.distribute!(bus_1)
    @test storage_2.input_interfaces[storage_2.medium].balance == 0.0
  
end

function test_extended_storage_control_strategy_deny_loading_by_storage()
    systems_config = Dict{String,Any}(
        "TST_GRI_01" => Dict{String,Any}(
            "type" => "GridConnection",
            "medium" => "m_c_g_natgas",
            "control_refs" => [],
            "production_refs" => ["TST_GBO_01"],
            "is_source" => true,
        ),
        "TST_GBO_01" => Dict{String,Any}(
            "type" => "GasBoiler",
            "control_refs" => [],
            "production_refs" => [
                "TST_BUS_01"
            ],
            "strategy" => Dict{String,Any}(
                "name" => "demand_driven",
                "load_storages" => false
            ),
            "power" => 10000
        ),
        "TST_BUS_01" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "production_refs" => ["TST_BUS_02", "TST_BFT_01"],
            "connection_matrix" => Dict{String, Any}(
                "input_order" => [
                    "TST_GBO_01",
                    "TST_BFT_01"
                ],
                "output_order" => [
                    "TST_BUS_02",
                    "TST_BFT_01"
                ]
            )
        ),
        "TST_BUS_02" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "production_refs" => ["TST_DEM_01", "TST_BFT_02"],
            "connection_matrix" => Dict{String, Any}(
                "input_order" => [
                    "TST_BUS_01",
                    "TST_BFT_02"
                ],
                "output_order" => [
                    "TST_DEM_01",
                    "TST_BFT_02"
                ]
            )
        ),
        "TST_BFT_01" => Dict{String,Any}(
            "type" => "BufferTank",
            "control_refs" => [],
            "production_refs" => [
                "TST_BUS_01"
            ],
            "strategy" => Dict{String,Any}(
                "name" => "extended_storage_control"
            ),
            "capacity" => 40000,
            "load" => 20000
        ),
        "TST_BFT_02" => Dict{String,Any}(
            "type" => "BufferTank",
            "control_refs" => [],
            "production_refs" => [
                "TST_BUS_02"
            ],
            "capacity" => 40000,
            "load" => 20000
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
    grid = systems["TST_GRI_01"]
    bus_1 = systems["TST_BUS_01"]
    bus_2 = systems["TST_BUS_02"]
    storage_1 = systems["TST_BFT_01"]
    storage_2 = systems["TST_BFT_02"]
    boiler = systems["TST_GBO_01"]

    simulation_parameters = Dict{String,Any}(
        "time_step_seconds" => 900,
        "time" => 0,
        "epsilon" => 1e-9
    )

    # Gasboiler IS NOT allowed to load storages, storage_1 IS NOT allowed to load storage_2.
    # As result the gasboilder should be provide max energy, the rest of the demand should be 
    # supplied by storage_2 and storage_1 should do nothing as all demand can be supplied by
    # storage_2

    EnergySystems.reset(grid)
    EnergySystems.reset(bus_1)
    EnergySystems.reset(bus_2)
    EnergySystems.reset(storage_1)
    EnergySystems.reset(storage_2)
    EnergySystems.reset(boiler)
    EnergySystems.reset(demand)

    EnergySystems.control(demand, systems, simulation_parameters)
    EnergySystems.control(bus_2, systems, simulation_parameters)
    EnergySystems.control(bus_1, systems, simulation_parameters)
    EnergySystems.control(boiler, systems, simulation_parameters)
    EnergySystems.control(storage_2, systems, simulation_parameters)
    EnergySystems.control(storage_1, systems, simulation_parameters)
    EnergySystems.control(grid, systems, simulation_parameters)

    EnergySystems.produce(demand, simulation_parameters, watt_to_wh)
    EnergySystems.produce(bus_2, simulation_parameters, watt_to_wh)
    EnergySystems.produce(bus_1, simulation_parameters, watt_to_wh)

    EnergySystems.produce(boiler, simulation_parameters, watt_to_wh)
    @test boiler.output_interfaces[boiler.m_heat_out].balance == 2500.0

    EnergySystems.produce(storage_2, simulation_parameters, watt_to_wh)
    @test storage_2.output_interfaces[storage_2.medium].balance == 2500.0
    @test storage_2.load == 20000.0 - 2500.0

    EnergySystems.produce(storage_1, simulation_parameters, watt_to_wh)
    @test storage_1.output_interfaces[storage_1.medium].balance == 0.0
    @test storage_1.load == 20000.0

    EnergySystems.load(storage_2, simulation_parameters, watt_to_wh)
    @test storage_2.input_interfaces[storage_2.medium].balance == 0.0
    @test storage_2.load == 20000.0 - 2500.0

    EnergySystems.load(storage_1, simulation_parameters, watt_to_wh)
    @test storage_1.input_interfaces[storage_1.medium].balance == 0.0
    @test storage_1.load == 20000.0

    EnergySystems.produce(grid, simulation_parameters, watt_to_wh)
    EnergySystems.distribute!(bus_2)
    EnergySystems.distribute!(bus_1)
  
end

function test_extended_storage_control_strategy_allow_loading_by_storage_and_gasboiler()
    systems_config = Dict{String,Any}(
        "TST_GRI_01" => Dict{String,Any}(
            "type" => "GridConnection",
            "medium" => "m_c_g_natgas",
            "control_refs" => [],
            "production_refs" => ["TST_GBO_01"],
            "is_source" => true,
        ),
        "TST_GBO_01" => Dict{String,Any}(
            "type" => "GasBoiler",
            "control_refs" => [],
            "production_refs" => [
                "TST_BUS_01"
            ],
            "strategy" => Dict{String,Any}(
                "name" => "demand_driven"
            ),
            "power" => 10000
        ),
        "TST_BUS_01" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "production_refs" => ["TST_BUS_02", "TST_BFT_01"],
            "connection_matrix" => Dict{String, Any}(
                "input_order" => [
                    "TST_GBO_01",
                    "TST_BFT_01"
                ],
                "output_order" => [
                    "TST_BUS_02",
                    "TST_BFT_01"
                ]
            )
        ),
        "TST_BUS_02" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "production_refs" => ["TST_DEM_01", "TST_BFT_02"],
            "connection_matrix" => Dict{String, Any}(
                "input_order" => [
                    "TST_BUS_01",
                    "TST_BFT_02"
                ],
                "output_order" => [
                    "TST_DEM_01",
                    "TST_BFT_02"
                ]
            )
        ),
        "TST_BFT_01" => Dict{String,Any}(
            "type" => "BufferTank",
            "control_refs" => [],
            "production_refs" => [
                "TST_BUS_01"
            ],
            "strategy" => Dict{String,Any}(
                "name" => "extended_storage_control",
                "load_any_storage" => true
            ),
            "capacity" => 40000,
            "load" => 20000
        ),
        "TST_BFT_02" => Dict{String,Any}(
            "type" => "BufferTank",
            "control_refs" => [],
            "production_refs" => [
                "TST_BUS_02"
            ],
            "strategy" => Dict{String,Any}(
                "name" => "extended_storage_control",
                "load_any_storage" => true
            ),
            "capacity" => 40000,
            "load" => 20000
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
    grid = systems["TST_GRI_01"]
    bus_1 = systems["TST_BUS_01"]
    bus_2 = systems["TST_BUS_02"]
    storage_1 = systems["TST_BFT_01"]
    storage_2 = systems["TST_BFT_02"]
    boiler = systems["TST_GBO_01"]

    simulation_parameters = Dict{String,Any}(
        "time_step_seconds" => 900,
        "time" => 0,
        "epsilon" => 1e-9
    )

    # timestep 1:
    # Gasboiler IS  allowed to load storages, and storage_1 IS allowed to load storage_2.
    # The demand is higher than the gasboiler max_energy.
    # As result the gasboilder should be provide max energy, the rest of the demand should be 
    # supplied by storage_2 which will be refilled by storage_1

    # Also "load_any_storage" is activated at storage_2 which should not have any effect on the 
    # results as storage_2 has no storage to load downstream.

    EnergySystems.reset(grid)
    EnergySystems.reset(bus_1)
    EnergySystems.reset(bus_2)
    EnergySystems.reset(storage_1)
    EnergySystems.reset(storage_2)
    EnergySystems.reset(boiler)
    EnergySystems.reset(demand)

    EnergySystems.control(demand, systems, simulation_parameters)
    EnergySystems.control(bus_2, systems, simulation_parameters)
    EnergySystems.control(bus_1, systems, simulation_parameters)
    EnergySystems.control(boiler, systems, simulation_parameters)
    EnergySystems.control(storage_2, systems, simulation_parameters)
    EnergySystems.control(storage_1, systems, simulation_parameters)
    EnergySystems.control(grid, systems, simulation_parameters)

    EnergySystems.produce(demand, simulation_parameters, watt_to_wh)
    EnergySystems.produce(bus_2, simulation_parameters, watt_to_wh)
    EnergySystems.produce(bus_1, simulation_parameters, watt_to_wh)

    EnergySystems.produce(boiler, simulation_parameters, watt_to_wh)
    @test boiler.output_interfaces[boiler.m_heat_out].balance == 2500.0

    EnergySystems.produce(storage_2, simulation_parameters, watt_to_wh)
    @test storage_2.output_interfaces[storage_2.medium].balance == 2500.0
    @test storage_2.load == 20000.0 - 2500.0

    EnergySystems.produce(storage_1, simulation_parameters, watt_to_wh)
    @test storage_1.output_interfaces[storage_1.medium].balance == 20000.0
    @test storage_1.load == 0

    EnergySystems.load(storage_2, simulation_parameters, watt_to_wh)
    @test storage_2.input_interfaces[storage_2.medium].balance == -20000.0
    @test storage_2.load == 40000.0 - 2500.0

    EnergySystems.load(storage_1, simulation_parameters, watt_to_wh)
    @test storage_1.input_interfaces[storage_1.medium].balance == 0.0
    @test storage_1.load == 0

    EnergySystems.produce(grid, simulation_parameters, watt_to_wh)
    EnergySystems.distribute!(bus_2)
    EnergySystems.distribute!(bus_1)
    @test storage_2.input_interfaces[storage_2.medium].balance == 0.0
  
    # timestep 2:
    # Gasboiler IS  allowed to load storages, and storage_1 IS allowed to load storage_2.
    # Now the demand is lower than the gasboiler max_energy.
    # As result the gasboilder should provide max energy, supply demand and load storage_2.
    # As storage_1 is futher upstream and therfore calculated later, is will not be loaded.

    EnergySystems.reset(grid)
    EnergySystems.reset(bus_1)
    EnergySystems.reset(bus_2)
    EnergySystems.reset(storage_1)
    EnergySystems.reset(storage_2)
    EnergySystems.reset(boiler)
    EnergySystems.reset(demand)

    demand.static_load = 2000

    EnergySystems.control(demand, systems, simulation_parameters)
    EnergySystems.control(bus_2, systems, simulation_parameters)
    EnergySystems.control(bus_1, systems, simulation_parameters)
    EnergySystems.control(boiler, systems, simulation_parameters)
    EnergySystems.control(storage_2, systems, simulation_parameters)
    EnergySystems.control(storage_1, systems, simulation_parameters)
    EnergySystems.control(grid, systems, simulation_parameters)

    EnergySystems.produce(demand, simulation_parameters, watt_to_wh)
    EnergySystems.produce(bus_2, simulation_parameters, watt_to_wh)
    EnergySystems.produce(bus_1, simulation_parameters, watt_to_wh)

    EnergySystems.produce(boiler, simulation_parameters, watt_to_wh)
    @test boiler.output_interfaces[boiler.m_heat_out].balance == 2500.0

    EnergySystems.produce(storage_2, simulation_parameters, watt_to_wh)
    @test storage_2.output_interfaces[storage_2.medium].balance == 0.0
    @test storage_2.load == 40000.0 - 2500.0

    EnergySystems.produce(storage_1, simulation_parameters, watt_to_wh)
    @test storage_1.output_interfaces[storage_1.medium].balance == 0.0
    @test storage_1.load == 0

    EnergySystems.load(storage_2, simulation_parameters, watt_to_wh)
    @test storage_2.input_interfaces[storage_2.medium].balance == -500.0
    @test storage_2.load == 40000.0 - 2500.0 + 500.0

    EnergySystems.load(storage_1, simulation_parameters, watt_to_wh)
    @test storage_1.input_interfaces[storage_1.medium].balance == 0.0
    @test storage_1.load == 0

    EnergySystems.produce(grid, simulation_parameters, watt_to_wh)
    EnergySystems.distribute!(bus_2)
    EnergySystems.distribute!(bus_1)
    @test storage_2.input_interfaces[storage_2.medium].balance == 0.0

end

@testset "extended_storage_control_strategy" begin
    test_extended_storage_control_strategy_allow_loading_by_storage()
    test_extended_storage_control_strategy_deny_loading_by_storage()
    test_extended_storage_control_strategy_allow_loading_by_storage_and_gasboiler()

end