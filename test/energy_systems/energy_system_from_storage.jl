using Debugger
using Test
using Resie
using Resie.EnergySystems
using Resie.Profiles

EnergySystems.set_timestep(900)

function test_run_energy_system_from_storage()
    components_config = Dict{String,Any}(
        "TST_DEM_01" => Dict{String,Any}(
            "type" => "Demand",
            "medium" => "m_h_w_ht1",
            "output_refs" => [],
            "energy_profile_file_path" => "./profiles/tests/demand_heating_energy.prf",
            "temperature_profile_file_path" => "./profiles/tests/demand_heating_temperature.prf",
            "scale" => 1500,
        ),
        "TST_BUS_01" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_h_w_lt1",
            "connections" => Dict{String,Any}(
                "input_order" => ["TST_BFT_01"],
                "output_order" => ["TST_HP_01",
                                   "TST_BFT_01"],
            ),
        ),
        "TST_BFT_01" => Dict{String,Any}(
            "type" => "BufferTank",
            "medium" => "m_h_w_lt1",
            "output_refs" => ["TST_BUS_01"],
            "capacity" => 40000,
            "load" => 30000,
            "high_temperature" => 35.0,
        ),
        "TST_GRI_01" => Dict{String,Any}(
            "type" => "GridConnection",
            "medium" => "m_e_ac_230v",
            "output_refs" => ["TST_HP_01"],
            "is_source" => true,
        ),
        "TST_HP_01" => Dict{String,Any}(
            "type" => "HeatPump",
            "output_refs" => ["TST_DEM_01"],
            "control_parameters" => Dict{String,Any}(
                "unload_storages m_e_ac_230v" => true,
            ),
            "m_el_in" => "m_e_ac_230v",
            "power_th" => 12000,
            "constant_cop" => 3.0,
            "min_power_fraction" => 0.0,
        ),
    )

    simulation_parameters = get_default_sim_params()

    components = Resie.load_components(components_config, simulation_parameters)
    heat_pump = components["TST_HP_01"]
    hheat_demand = components["TST_DEM_01"]
    power_grid = components["TST_GRI_01"]
    lheat_storage = components["TST_BFT_01"]
    lheat_bus = components["TST_BUS_01"]

    # first time step: storage is full to power heatpump

    for unit in values(components)
        EnergySystems.reset(unit)
    end

    EnergySystems.control(heat_pump, components, simulation_parameters)
    EnergySystems.control(power_grid, components, simulation_parameters)
    EnergySystems.control(lheat_storage, components, simulation_parameters)
    EnergySystems.control(lheat_bus, components, simulation_parameters)

    hheat_demand.constant_demand = 800 * 4
    hheat_demand.constant_temperature = 45.0
    EnergySystems.control(hheat_demand, components, simulation_parameters)

    EnergySystems.process(hheat_demand, simulation_parameters)
    @test hheat_demand.input_interfaces[hheat_demand.medium].balance ≈ -800
    @test hheat_demand.input_interfaces[hheat_demand.medium].temperature_min === 45.0

    # demand not processed yet --> balance is zero, but energy_potential not
    # input interfaces
    exchanges = EnergySystems.balance_on(heat_pump.input_interfaces[lheat_bus.medium], lheat_bus)
    @test EnergySystems.balance(exchanges) ≈ 0.0
    @test EnergySystems.energy_potential(exchanges) ≈ 30000
    @test EnergySystems.temp_max_highest(exchanges) === 35.0

    EnergySystems.process(heat_pump, simulation_parameters)
    @test heat_pump.output_interfaces[heat_pump.m_heat_out].balance ≈ 0
    @test heat_pump.output_interfaces[heat_pump.m_heat_out].sum_abs_change ≈ 1600
    @test heat_pump.output_interfaces[heat_pump.m_heat_out].temperature_min ≈ 45.0
    @test heat_pump.output_interfaces[heat_pump.m_heat_out].temperature_max === nothing
    @test heat_pump.input_interfaces[heat_pump.m_el_in].balance ≈ -800 / 3
    @test heat_pump.input_interfaces[heat_pump.m_el_in].temperature_min === nothing
    @test heat_pump.input_interfaces[heat_pump.m_heat_in].balance ≈ -800 * 2 / 3
    @test heat_pump.input_interfaces[heat_pump.m_heat_in].temperature_min === nothing

    EnergySystems.process(lheat_bus, simulation_parameters)
    EnergySystems.process(lheat_storage, simulation_parameters)
    EnergySystems.load(lheat_storage, simulation_parameters)

    @test lheat_storage.output_interfaces[lheat_storage.medium].balance ≈ 800 * 2 / 3
    @test lheat_storage.output_interfaces[lheat_storage.medium].sum_abs_change ≈ 800 * 2 / 3
    @test lheat_storage.output_interfaces[lheat_storage.medium].temperature_max ≈ 35.0

    EnergySystems.process(power_grid, simulation_parameters)
    @test power_grid.output_interfaces[power_grid.medium].balance ≈ 0
    @test power_grid.output_interfaces[power_grid.medium].sum_abs_change ≈ 2 * 800 / 3
    @test power_grid.output_interfaces[power_grid.medium].temperature_max === nothing

    # second step: storage is nearly empty, operation of heatpump is limited

    for unit in values(components)
        EnergySystems.reset(unit)
    end

    EnergySystems.control(heat_pump, components, simulation_parameters)
    EnergySystems.control(power_grid, components, simulation_parameters)
    EnergySystems.control(lheat_bus, components, simulation_parameters)

    lheat_storage.load = 100.0
    EnergySystems.control(lheat_storage, components, simulation_parameters)

    hheat_demand.constant_demand = 800 * 4
    hheat_demand.constant_temperature = 45.0
    EnergySystems.control(hheat_demand, components, simulation_parameters)

    EnergySystems.process(hheat_demand, simulation_parameters)
    @test hheat_demand.input_interfaces[hheat_demand.medium].balance ≈ -800
    @test hheat_demand.input_interfaces[hheat_demand.medium].temperature_min ≈ 45.0

    # demand not processed yet --> balance is zero, but energy_potential not
    # input interfaces
    exchanges = EnergySystems.balance_on(heat_pump.input_interfaces[lheat_bus.medium], lheat_bus)
    @test EnergySystems.balance(exchanges) ≈ 0.0
    @test EnergySystems.energy_potential(exchanges) ≈ 100
    @test EnergySystems.temp_max_highest(exchanges) === 35.0

    EnergySystems.process(heat_pump, simulation_parameters)
    @test heat_pump.output_interfaces[heat_pump.m_heat_out].balance ≈ -800 + 100 * 3 / 2
    @test heat_pump.output_interfaces[heat_pump.m_heat_out].sum_abs_change ≈ 800 + 100 * 3 / 2
    @test heat_pump.output_interfaces[heat_pump.m_heat_out].temperature_min ≈ 45.0
    @test heat_pump.output_interfaces[heat_pump.m_heat_out].temperature_max === nothing
    @test heat_pump.input_interfaces[heat_pump.m_el_in].balance ≈ -(100 * 3 / 2) / 3
    @test heat_pump.input_interfaces[heat_pump.m_el_in].temperature_min === nothing
    @test heat_pump.input_interfaces[heat_pump.m_heat_in].balance ≈ -100
    @test heat_pump.input_interfaces[heat_pump.m_heat_in].temperature_min === nothing

    EnergySystems.process(lheat_bus, simulation_parameters)
    EnergySystems.process(lheat_storage, simulation_parameters)
    EnergySystems.load(lheat_storage, simulation_parameters)

    @test lheat_storage.output_interfaces[lheat_storage.medium].balance ≈ 100
    @test lheat_storage.output_interfaces[lheat_storage.medium].sum_abs_change ≈ 100
    @test lheat_storage.output_interfaces[lheat_storage.medium].temperature_max === 35.0

    EnergySystems.process(power_grid, simulation_parameters)
    @test power_grid.output_interfaces[power_grid.medium].balance ≈ 0
    @test power_grid.output_interfaces[power_grid.medium].sum_abs_change ≈ 2 * (100 * 3 / 2) / 3
    @test power_grid.output_interfaces[power_grid.medium].temperature_max === nothing
end

function test_run_energy_system_from_storage_denied()
    components_config = Dict{String,Any}(
        "TST_DEM_01" => Dict{String,Any}(
            "type" => "Demand",
            "medium" => "m_h_w_ht1",
            "output_refs" => [],
            "energy_profile_file_path" => "./profiles/tests/demand_heating_energy.prf",
            "temperature_profile_file_path" => "./profiles/tests/demand_heating_temperature.prf",
            "scale" => 1500,
        ),
        "TST_BUS_01" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_h_w_lt1",
            "connections" => Dict{String,Any}(
                "input_order" => ["TST_BFT_01"],
                "output_order" => ["TST_HP_01",
                                   "TST_BFT_01"],
            ),
        ),
        "TST_BFT_01" => Dict{String,Any}(
            "type" => "BufferTank",
            "medium" => "m_h_w_lt1",
            "output_refs" => ["TST_BUS_01"],
            "capacity" => 40000,
            "load" => 30000,
            "high_temperature" => 35,
        ),
        "TST_GRI_01" => Dict{String,Any}(
            "type" => "GridConnection",
            "medium" => "m_e_ac_230v",
            "output_refs" => ["TST_HP_01"],
            "is_source" => true,
        ),
        "TST_HP_01" => Dict{String,Any}(
            "type" => "HeatPump",
            "output_refs" => ["TST_DEM_01"],
            "control_parameters" => Dict{String,Any}(
                "unload_storages m_h_w_lt1" => false,
            ),
            "m_el_in" => "m_e_ac_230v",
            "power_th" => 12000,
            "constant_cop" => 3.0,
            "min_power_fraction" => 0.0,
        ),
    )

    simulation_parameters = get_default_sim_params()

    components = Resie.load_components(components_config, simulation_parameters)
    heat_pump = components["TST_HP_01"]
    hheat_demand = components["TST_DEM_01"]
    power_grid = components["TST_GRI_01"]
    lheat_storage = components["TST_BFT_01"]
    lheat_bus = components["TST_BUS_01"]

    # first time step: storage is full to power heatpump, but heatpump unloading storages is test_run_energy_system_from_storage_denied
    # not energy should be transferred at all.

    for unit in values(components)
        EnergySystems.reset(unit)
    end

    EnergySystems.control(heat_pump, components, simulation_parameters)
    EnergySystems.control(power_grid, components, simulation_parameters)
    EnergySystems.control(lheat_storage, components, simulation_parameters)
    EnergySystems.control(lheat_bus, components, simulation_parameters)

    hheat_demand.constant_demand = 800 * 4
    hheat_demand.constant_temperature = 45.0
    EnergySystems.control(hheat_demand, components, simulation_parameters)

    EnergySystems.process(hheat_demand, simulation_parameters)
    @test hheat_demand.input_interfaces[hheat_demand.medium].balance ≈ -800
    @test hheat_demand.input_interfaces[hheat_demand.medium].temperature_min === 45.0

    # demand not processed yet --> balance is zero, but energy_potential not
    # input interfaces
    exchanges = EnergySystems.balance_on(heat_pump.input_interfaces[lheat_bus.medium], lheat_bus)
    @test EnergySystems.balance(exchanges) ≈ 0.0
    @test EnergySystems.energy_potential(exchanges) ≈ 0
    @test EnergySystems.temp_max_highest(exchanges) === nothing

    EnergySystems.process(heat_pump, simulation_parameters)
    @test heat_pump.output_interfaces[heat_pump.m_heat_out].balance ≈ -800
    @test heat_pump.output_interfaces[heat_pump.m_heat_out].sum_abs_change ≈ 800
    @test heat_pump.output_interfaces[heat_pump.m_heat_out].temperature_min ≈ 45.0
    @test heat_pump.output_interfaces[heat_pump.m_heat_out].temperature_max === nothing
    @test heat_pump.input_interfaces[heat_pump.m_el_in].balance ≈ 0.0
    @test heat_pump.input_interfaces[heat_pump.m_el_in].temperature_min === nothing
    @test heat_pump.input_interfaces[heat_pump.m_heat_in].balance ≈ 0.0
    @test heat_pump.input_interfaces[heat_pump.m_heat_in].temperature_min === nothing

    EnergySystems.process(lheat_bus, simulation_parameters)
    EnergySystems.process(lheat_storage, simulation_parameters)
    EnergySystems.load(lheat_storage, simulation_parameters)

    @test lheat_storage.output_interfaces[lheat_storage.medium].balance ≈ 0.0
    @test lheat_storage.output_interfaces[lheat_storage.medium].sum_abs_change ≈ 0.0
    @test lheat_storage.output_interfaces[lheat_storage.medium].temperature_max === 35.0

    EnergySystems.process(power_grid, simulation_parameters)
    @test power_grid.output_interfaces[power_grid.medium].balance ≈ 0
    @test power_grid.output_interfaces[power_grid.medium].sum_abs_change ≈ 0.0
    @test power_grid.output_interfaces[power_grid.medium].temperature_max === nothing
end

@testset "run_energy_system_from_storage" begin
    test_run_energy_system_from_storage()
    test_run_energy_system_from_storage_denied()
end
