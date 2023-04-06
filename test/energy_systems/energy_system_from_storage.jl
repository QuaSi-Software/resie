using Debugger
using Test
using Resie
using Resie.EnergySystems
using Resie.Profiles

watt_to_wh = function (watts::Float64)
    watts * 900 / 3600.0
end

function test_run_energy_system_from_storage()
    systems_config = Dict{String,Any}(
        "TST_DEM_01" => Dict{String,Any}(
            "type" => "Demand",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "production_refs" => [],
            "energy_profile_file_path" => "./profiles/tests/demand_heating_energy.prf",
            "temperature_profile_file_path" => "./profiles/tests/demand_heating_temperature.prf",
            "scale" => 1500
        ),
        "TST_BUS_01" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_h_w_lt1",
            "control_refs" => [],
            "production_refs" => ["TST_HP_01", "TST_BFT_01"],
            "connection_matrix" => Dict{String, Any}(
                "input_order" => [
                    "TST_BFT_01"
                ],
                "output_order" => [
                    "TST_HP_01",
                    "TST_BFT_01"
                ],
            )
        ),
        "TST_BFT_01" => Dict{String,Any}(
            "type" => "BufferTank",
            "medium" => "m_h_w_lt1",
            "control_refs" => [],
            "production_refs" => [
                "TST_BUS_01"
            ],
            "capacity" => 40000,
            "load" => 30000,
            "high_temperature" => 35.0
        ),
        "TST_GRI_01" => Dict{String,Any}(
            "type" => "GridConnection",
            "medium" => "m_e_ac_230v",
            "control_refs" => [],
            "production_refs" => ["TST_HP_01"],
            "is_source" => true,
        ),
        "TST_HP_01" => Dict{String,Any}(
            "type" => "HeatPump",
            "control_refs" => ["TST_DEM_01"],
            "production_refs" => ["TST_DEM_01"],
            "strategy" => Dict{String,Any}(
                "name" => "demand_driven",
                "unload_storages" => true
            ),
            "power" => 12000,
            "fixed_cop" => 3.0,
            "min_power_fraction" => 0.0
        ),
    )
    systems = Resie.load_systems(systems_config)
    heat_pump = systems["TST_HP_01"]
    hheat_demand = systems["TST_DEM_01"]
    power_grid = systems["TST_GRI_01"]
    lheat_storage = systems["TST_BFT_01"]
    lheat_bus = systems["TST_BUS_01"]

    simulation_parameters = Dict{String,Any}(
        "time_step_seconds" => 900,
        "time" => 0,
        "epsilon" => 1e-9
    )

    @test heat_pump.controller.state_machine.state == 1

    # first time step: storage is full to power heatpump
    
    for unit in values(systems)
        EnergySystems.reset(unit)
    end

    EnergySystems.control(hheat_demand, systems, simulation_parameters)
    EnergySystems.control(heat_pump, systems, simulation_parameters)
    EnergySystems.control(power_grid, systems, simulation_parameters)
    EnergySystems.control(lheat_storage, systems, simulation_parameters)
    EnergySystems.control(lheat_bus, systems, simulation_parameters)

    hheat_demand.load = 800
    hheat_demand.temperature = 45.0
    hheat_demand.input_interfaces[hheat_demand.medium].temperature = 45.0

    EnergySystems.produce(hheat_demand, simulation_parameters, watt_to_wh)
    @test hheat_demand.input_interfaces[hheat_demand.medium].balance ≈ -800
    @test hheat_demand.input_interfaces[hheat_demand.medium].temperature === 45.0

    # demand not produced yet --> balance is zero, but energy_potential not
    # input interfaces
    exchange = EnergySystems.balance_on(heat_pump.input_interfaces[lheat_bus.medium], lheat_bus)
    @test exchange.balance ≈ 0.0
    @test exchange.storage_potential ≈ 30000
    @test exchange.energy_potential ≈ 0.0
    @test exchange.temperature === 35.0

    EnergySystems.produce(heat_pump, simulation_parameters, watt_to_wh)
    @test heat_pump.output_interfaces[heat_pump.m_heat_out].balance ≈ 0
    @test heat_pump.output_interfaces[heat_pump.m_heat_out].sum_abs_change ≈ 1600
    @test heat_pump.output_interfaces[heat_pump.m_heat_out].temperature === 45.0
    @test heat_pump.input_interfaces[heat_pump.m_el_in].balance ≈ -800/3
    @test heat_pump.input_interfaces[heat_pump.m_el_in].temperature === nothing
    @test heat_pump.input_interfaces[heat_pump.m_heat_in].balance ≈ -800*2/3
    @test heat_pump.input_interfaces[heat_pump.m_heat_in].temperature === nothing

    EnergySystems.produce(lheat_bus, simulation_parameters, watt_to_wh)
    EnergySystems.produce(lheat_storage, simulation_parameters, watt_to_wh)
    EnergySystems.load(lheat_storage, simulation_parameters, watt_to_wh)

    @test lheat_storage.output_interfaces[lheat_storage.medium].balance ≈ 800*2/3
    @test lheat_storage.output_interfaces[lheat_storage.medium].sum_abs_change ≈ 800*2/3
    @test lheat_storage.output_interfaces[lheat_storage.medium].temperature === 35.0

    EnergySystems.produce(power_grid, simulation_parameters, watt_to_wh)
    @test power_grid.output_interfaces[power_grid.medium].balance ≈ 0
    @test power_grid.output_interfaces[power_grid.medium].sum_abs_change ≈ 2*800/3
    @test power_grid.output_interfaces[power_grid.medium].temperature === nothing

    # second step: storage is nearly empty, operation of heatpump is limited

    for unit in values(systems)
        EnergySystems.reset(unit)
    end

    EnergySystems.control(hheat_demand, systems, simulation_parameters)
    EnergySystems.control(heat_pump, systems, simulation_parameters)
    EnergySystems.control(power_grid, systems, simulation_parameters)
    EnergySystems.control(lheat_storage, systems, simulation_parameters)
    EnergySystems.control(lheat_bus, systems, simulation_parameters)

    lheat_storage.load = 100.0

    hheat_demand.load = 800
    hheat_demand.temperature = 45.0
    hheat_demand.input_interfaces[hheat_demand.medium].temperature = 45.0

    EnergySystems.produce(hheat_demand, simulation_parameters, watt_to_wh)
    @test hheat_demand.input_interfaces[hheat_demand.medium].balance ≈ -800
    @test hheat_demand.input_interfaces[hheat_demand.medium].temperature === 45.0

    # demand not produced yet --> balance is zero, but energy_potential not
    # input interfaces
    exchange = EnergySystems.balance_on(heat_pump.input_interfaces[lheat_bus.medium], lheat_bus)
    @test exchange.balance ≈ 0.0
    @test exchange.storage_potential ≈ 100
    @test exchange.energy_potential ≈ 0.0
    @test exchange.temperature === 35.0

    EnergySystems.produce(heat_pump, simulation_parameters, watt_to_wh)
    @test heat_pump.output_interfaces[heat_pump.m_heat_out].balance ≈ -800 + 100*3/2
    @test heat_pump.output_interfaces[heat_pump.m_heat_out].sum_abs_change ≈ 800+100*3/2
    @test heat_pump.output_interfaces[heat_pump.m_heat_out].temperature === 45.0
    @test heat_pump.input_interfaces[heat_pump.m_el_in].balance ≈ -(100*3/2)/3
    @test heat_pump.input_interfaces[heat_pump.m_el_in].temperature === nothing
    @test heat_pump.input_interfaces[heat_pump.m_heat_in].balance ≈ -100
    @test heat_pump.input_interfaces[heat_pump.m_heat_in].temperature === nothing

    EnergySystems.produce(lheat_bus, simulation_parameters, watt_to_wh)
    EnergySystems.produce(lheat_storage, simulation_parameters, watt_to_wh)
    EnergySystems.load(lheat_storage, simulation_parameters, watt_to_wh)

    @test lheat_storage.output_interfaces[lheat_storage.medium].balance ≈ 100
    @test lheat_storage.output_interfaces[lheat_storage.medium].sum_abs_change ≈ 100
    @test lheat_storage.output_interfaces[lheat_storage.medium].temperature === 35.0

    EnergySystems.produce(power_grid, simulation_parameters, watt_to_wh)
    @test power_grid.output_interfaces[power_grid.medium].balance ≈ 0
    @test power_grid.output_interfaces[power_grid.medium].sum_abs_change ≈ 2*(100*3/2)/3
    @test power_grid.output_interfaces[power_grid.medium].temperature === nothing

end

function test_run_energy_system_from_storage_denied()
    systems_config = Dict{String,Any}(
        "TST_DEM_01" => Dict{String,Any}(
            "type" => "Demand",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "production_refs" => [],
            "energy_profile_file_path" => "./profiles/tests/demand_heating_energy.prf",
            "temperature_profile_file_path" => "./profiles/tests/demand_heating_temperature.prf",
            "scale" => 1500
        ),
        "TST_BUS_01" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_h_w_lt1",
            "control_refs" => [],
            "production_refs" => ["TST_HP_01", "TST_BFT_01"],
            "connection_matrix" => Dict{String, Any}(
                "input_order" => [
                    "TST_BFT_01"
                ],
                "output_order" => [
                    "TST_HP_01",
                    "TST_BFT_01"
                ],
            )
        ),
        "TST_BFT_01" => Dict{String,Any}(
            "type" => "BufferTank",
            "medium" => "m_h_w_lt1",
            "control_refs" => [],
            "production_refs" => [
                "TST_BUS_01"
            ],
            "capacity" => 40000,
            "load" => 30000,
            "high_temperature" => 35
        ),
        "TST_GRI_01" => Dict{String,Any}(
            "type" => "GridConnection",
            "medium" => "m_e_ac_230v",
            "control_refs" => [],
            "production_refs" => ["TST_HP_01"],
            "is_source" => true,
        ),
        "TST_HP_01" => Dict{String,Any}(
            "type" => "HeatPump",
            "control_refs" => ["TST_DEM_01"],
            "production_refs" => ["TST_DEM_01"],
            "strategy" => Dict{String,Any}(
                "name" => "demand_driven",
                "unload_storages" => false
            ),
            "power" => 12000,
            "fixed_cop" => 3.0,
            "min_power_fraction" => 0.0
        ),
    )
    systems = Resie.load_systems(systems_config)
    heat_pump = systems["TST_HP_01"]
    hheat_demand = systems["TST_DEM_01"]
    power_grid = systems["TST_GRI_01"]
    lheat_storage = systems["TST_BFT_01"]
    lheat_bus = systems["TST_BUS_01"]

    simulation_parameters = Dict{String,Any}(
        "time_step_seconds" => 900,
        "time" => 0,
        "epsilon" => 1e-9
    )

    @test heat_pump.controller.state_machine.state == 1

    # first time step: storage is full to power heatpump, but heatpump unloading storages is test_run_energy_system_from_storage_denied
    # not energy should be transferred at all.
    
    for unit in values(systems)
        EnergySystems.reset(unit)
    end

    EnergySystems.control(hheat_demand, systems, simulation_parameters)
    EnergySystems.control(heat_pump, systems, simulation_parameters)
    EnergySystems.control(power_grid, systems, simulation_parameters)
    EnergySystems.control(lheat_storage, systems, simulation_parameters)
    EnergySystems.control(lheat_bus, systems, simulation_parameters)

    hheat_demand.load = 800
    hheat_demand.temperature = 45.0
    hheat_demand.input_interfaces[hheat_demand.medium].temperature = 45.0

    EnergySystems.produce(hheat_demand, simulation_parameters, watt_to_wh)
    @test hheat_demand.input_interfaces[hheat_demand.medium].balance ≈ -800
    @test hheat_demand.input_interfaces[hheat_demand.medium].temperature === 45.0

    # demand not produced yet --> balance is zero, but energy_potential not
    # input interfaces
    exchange = EnergySystems.balance_on(heat_pump.input_interfaces[lheat_bus.medium], lheat_bus)
    @test exchange.balance ≈ 0.0
    @test exchange.storage_potential ≈ 30000
    @test exchange.energy_potential ≈ 0.0
    @test exchange.temperature === 35.0

    EnergySystems.produce(heat_pump, simulation_parameters, watt_to_wh)
    @test heat_pump.output_interfaces[heat_pump.m_heat_out].balance ≈ -800
    @test heat_pump.output_interfaces[heat_pump.m_heat_out].sum_abs_change ≈ 800
    @test heat_pump.output_interfaces[heat_pump.m_heat_out].temperature === 45.0
    @test heat_pump.input_interfaces[heat_pump.m_el_in].balance ≈ 0.0
    @test heat_pump.input_interfaces[heat_pump.m_el_in].temperature === nothing
    @test heat_pump.input_interfaces[heat_pump.m_heat_in].balance ≈ 0.0
    @test heat_pump.input_interfaces[heat_pump.m_heat_in].temperature === nothing

    EnergySystems.produce(lheat_bus, simulation_parameters, watt_to_wh)
    EnergySystems.produce(lheat_storage, simulation_parameters, watt_to_wh)
    EnergySystems.load(lheat_storage, simulation_parameters, watt_to_wh)

    @test lheat_storage.output_interfaces[lheat_storage.medium].balance ≈ 0.0
    @test lheat_storage.output_interfaces[lheat_storage.medium].sum_abs_change ≈ 0.0
    @test lheat_storage.output_interfaces[lheat_storage.medium].temperature === 35.0

    EnergySystems.produce(power_grid, simulation_parameters, watt_to_wh)
    @test power_grid.output_interfaces[power_grid.medium].balance ≈ 0
    @test power_grid.output_interfaces[power_grid.medium].sum_abs_change ≈ 0.0
    @test power_grid.output_interfaces[power_grid.medium].temperature === nothing

end

@testset "run_energy_system_from_storage" begin
    test_run_energy_system_from_storage()
    test_run_energy_system_from_storage_denied()
end