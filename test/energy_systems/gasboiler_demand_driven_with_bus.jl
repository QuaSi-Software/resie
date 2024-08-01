using Debugger
using Test
using Resie
using Resie.EnergySystems
using Resie.Profiles

EnergySystems.set_timestep(900)

function test_gasboiler_demand_driven_with_bus()
   components_config = Dict{String,Any}(
        "TST_DEM_01" => Dict{String,Any}(
            "type" => "Demand",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "output_refs" => [],
            "energy_profile_file_path" => "./profiles/tests/demand_heating_energy.prf",
            "temperature_profile_file_path" => "./profiles/tests/demand_heating_temperature.prf",
            "scale" => 1500
        ),
        "TST_BUS_01" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "connections" => Dict{String,Any}(
                "input_order" => ["TST_GB_01"],
                "output_order" => ["TST_DEM_01"]
            )
        ),
        "TST_GRI_01" => Dict{String,Any}(
            "type" => "GridConnection",
            "medium" => "m_c_g_natgas",
            "control_refs" => [],
            "output_refs" => ["TST_GB_01"],
            "is_source" => true,
        ),
        "TST_GB_01" => Dict{String,Any}(
            "type" => "FuelBoiler",
            "m_fuel_in" => "m_c_g_natgas",
            "control_refs" => ["TST_BUS_01"],
            "output_refs" => ["TST_BUS_01"],
            "strategy" => Dict{String,Any}(
                "name" => "demand_driven",
            ),
            "power_th" => 12000,
            "efficiency_fuel_in" => "const:1.0",
        ),
    )

    simulation_parameters = Dict{String,Any}(
        "time_step_seconds" => 900,
        "time" => 0,
        "epsilon" => 1e-9,
        "is_first_timestep" => true
    )

    components = Resie.load_components(components_config, simulation_parameters)
    gasboiler = components["TST_GB_01"]
    grid = components["TST_GRI_01"]
    demand = components["TST_DEM_01"]
    bus = components["TST_BUS_01"]

    # test if correct demand is processed by gasboiler to make sure that the information of
    # demand is transported through Bus.

    # first time step: demand is exactly the max power of GasBoiler

    for unit in values(components)
        EnergySystems.reset(unit)
    end

    EnergySystems.control(demand, components, simulation_parameters)
    @test demand.demand == EnergySystems.get_max_energy(demand.input_interfaces[demand.medium].max_energy)
    @test demand.temperature == demand.input_interfaces[demand.medium].temperature_min
    EnergySystems.reset(demand)
    demand.constant_demand = 12000
    demand.constant_temperature = 85
    EnergySystems.control(demand, components, simulation_parameters)

    EnergySystems.control(bus, components, simulation_parameters)
    EnergySystems.control(gasboiler, components, simulation_parameters)
    EnergySystems.control(grid, components, simulation_parameters)
    @test EnergySystems.get_max_energy(grid.output_interfaces[grid.medium].max_energy) == Inf

    # no processing so far, balance is zero, energy_potential is non-zero
    exchanges = EnergySystems.balance_on(gasboiler.output_interfaces[bus.medium], bus)
    @test EnergySystems.balance(exchanges) ≈ 0.0
    @test EnergySystems.energy_potential(exchanges) ≈ -12000/4

    # no processing so far, balance is zero, energy_potential is zero
    exchanges = EnergySystems.balance_on(demand.input_interfaces[bus.medium], bus)
    @test EnergySystems.balance(exchanges) ≈ 0.0
    @test EnergySystems.energy_potential(exchanges) ≈ Inf

    EnergySystems.process(demand, simulation_parameters)
    @test demand.input_interfaces[demand.medium].balance ≈ -12000/4
    @test demand.input_interfaces[demand.medium].temperature_min ≈ 85

    # demand was processed --> energy_potential shows the possible energy to deliver, balance is zero
    exchanges = EnergySystems.balance_on(gasboiler.output_interfaces[bus.medium], bus)
    @test EnergySystems.balance(exchanges) ≈ 0.0  # balance of bus is always zero
    @test EnergySystems.energy_potential(exchanges) ≈ -12000.0/4

    EnergySystems.process(bus, simulation_parameters)
    EnergySystems.process(gasboiler, simulation_parameters)
    @test gasboiler.output_interfaces[gasboiler.m_heat_out].balance ≈ 12000/4
    @test gasboiler.input_interfaces[gasboiler.m_fuel_in].balance ≈ -12000/4

    EnergySystems.process(grid, simulation_parameters)
    @test grid.output_interfaces[grid.medium].balance ≈ 0
    @test grid.output_interfaces[grid.medium].sum_abs_change ≈ 2*12000/4

    EnergySystems.distribute!(bus)
    @test gasboiler.output_interfaces[gasboiler.m_heat_out].balance ≈ 0
    @test gasboiler.output_interfaces[gasboiler.m_heat_out].sum_abs_change ≈ 2*12000/4
    @test demand.input_interfaces[demand.medium].balance ≈ 0
    @test demand.input_interfaces[demand.medium].sum_abs_change ≈ 2*12000/4
    @test bus.remainder ≈ 0

    # second step: demand is above max power of GasBoiler

    for unit in values(components)
        EnergySystems.reset(unit)
    end

    EnergySystems.control(demand, components, simulation_parameters)
    @test demand.demand == EnergySystems.get_max_energy(demand.input_interfaces[demand.medium].max_energy)
    @test demand.temperature == demand.input_interfaces[demand.medium].temperature_min
    EnergySystems.reset(demand)
    demand.constant_demand = 15000
    demand.constant_temperature = 85
    EnergySystems.control(demand, components, simulation_parameters)

    EnergySystems.control(bus, components, simulation_parameters)
    EnergySystems.control(gasboiler, components, simulation_parameters)
    EnergySystems.control(grid, components, simulation_parameters)
    @test EnergySystems.get_max_energy(grid.output_interfaces[grid.medium].max_energy) == Inf

    # no processing so far, balance is zero, energy_potential is non-zero
    exchanges = EnergySystems.balance_on(gasboiler.output_interfaces[bus.medium], bus)
    @test EnergySystems.balance(exchanges) ≈ 0.0
    @test EnergySystems.energy_potential(exchanges) ≈ -15000/4

    # no processing so far, balance is zero, energy_potential is zero
    exchanges = EnergySystems.balance_on(demand.input_interfaces[bus.medium], bus)
    @test EnergySystems.balance(exchanges) ≈ 0.0
    @test EnergySystems.energy_potential(exchanges) ≈ Inf

    EnergySystems.process(demand, simulation_parameters)
    @test demand.input_interfaces[demand.medium].balance ≈ -15000/4
    @test demand.input_interfaces[demand.medium].temperature_min ≈ 85

    # demand was processed --> energy_potential shows the possible energy to deliver, balance is zero
    exchanges = EnergySystems.balance_on(gasboiler.output_interfaces[bus.medium], bus)
    @test EnergySystems.balance(exchanges) ≈ 0.0  # balance of bus is always zero
    @test EnergySystems.energy_potential(exchanges) ≈ -15000.0/4

    EnergySystems.process(bus, simulation_parameters)
    EnergySystems.process(gasboiler, simulation_parameters)
    @test gasboiler.output_interfaces[gasboiler.m_heat_out].balance ≈ 12000/4
    @test gasboiler.input_interfaces[gasboiler.m_fuel_in].balance ≈ -12000/4

    EnergySystems.process(grid, simulation_parameters)
    @test grid.output_interfaces[grid.medium].balance ≈ 0
    @test grid.output_interfaces[grid.medium].sum_abs_change ≈ 2*12000/4

    EnergySystems.distribute!(bus)
    @test gasboiler.output_interfaces[gasboiler.m_heat_out].balance ≈ 0
    @test gasboiler.output_interfaces[gasboiler.m_heat_out].sum_abs_change ≈ 2*12000/4
    @test demand.input_interfaces[demand.medium].balance ≈ 0 #-15000/4 + 12000/4   # is Bus intended to supply all demand???
    @test demand.input_interfaces[demand.medium].sum_abs_change ≈ 15000/4 + 15000/4
    @test bus.remainder ≈ -15000/4 + 12000/4

    # third step: demand is below max power of GasBoiler

    for unit in values(components)
        EnergySystems.reset(unit)
    end

    EnergySystems.control(demand, components, simulation_parameters)
    @test demand.demand == EnergySystems.get_max_energy(demand.input_interfaces[demand.medium].max_energy)
    @test demand.temperature == demand.input_interfaces[demand.medium].temperature_min
    EnergySystems.reset(demand)
    demand.constant_demand = 10000
    demand.constant_temperature = 85
    EnergySystems.control(demand, components, simulation_parameters)

    EnergySystems.control(bus, components, simulation_parameters)
    EnergySystems.control(gasboiler, components, simulation_parameters)
    EnergySystems.control(grid, components, simulation_parameters)
    @test EnergySystems.get_max_energy(grid.output_interfaces[grid.medium].max_energy) == Inf

    # no processing so far, balance is zero, energy_potential is non-zero
    exchanges = EnergySystems.balance_on(gasboiler.output_interfaces[bus.medium], bus)
    @test EnergySystems.balance(exchanges) ≈ 0.0
    @test EnergySystems.energy_potential(exchanges) ≈ -10000/4

    # no processing so far, balance is zero, energy_potential is zero
    exchanges = EnergySystems.balance_on(demand.input_interfaces[bus.medium], bus)
    @test EnergySystems.balance(exchanges) ≈ 0.0
    @test EnergySystems.energy_potential(exchanges) ≈ Inf
    EnergySystems.process(demand, simulation_parameters)
    @test demand.input_interfaces[demand.medium].balance ≈ -10000/4
    @test demand.input_interfaces[demand.medium].temperature_min ≈ 85

    # demand was processed --> energy_potential shows the possible energy to deliver, balance is zero
    exchanges = EnergySystems.balance_on(gasboiler.output_interfaces[bus.medium], bus)
    @test EnergySystems.balance(exchanges) ≈ 0.0  # balance of bus is always zero
    @test EnergySystems.energy_potential(exchanges) ≈ -10000.0/4

    EnergySystems.process(bus, simulation_parameters)
    EnergySystems.process(gasboiler, simulation_parameters)
    @test gasboiler.output_interfaces[gasboiler.m_heat_out].balance ≈ 10000/4
    @test gasboiler.input_interfaces[gasboiler.m_fuel_in].balance ≈ -10000/4

    EnergySystems.process(grid, simulation_parameters)
    @test grid.output_interfaces[grid.medium].balance ≈ 0
    @test grid.output_interfaces[grid.medium].sum_abs_change ≈ 2*10000/4

    EnergySystems.distribute!(bus)
    @test gasboiler.output_interfaces[gasboiler.m_heat_out].balance ≈ 0
    @test gasboiler.output_interfaces[gasboiler.m_heat_out].sum_abs_change ≈ 2*10000/4
    @test demand.input_interfaces[demand.medium].balance ≈ 0
    @test demand.input_interfaces[demand.medium].sum_abs_change ≈ 2*10000/4
    @test bus.remainder ≈ 0

end

function test_gasboiler_demand_driven_without_bus()
    components_config = Dict{String,Any}(
        "TST_DEM_01" => Dict{String,Any}(
            "type" => "Demand",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "output_refs" => [],
            "energy_profile_file_path" => "./profiles/tests/demand_heating_energy.prf",
            "temperature_profile_file_path" => "./profiles/tests/demand_heating_temperature.prf",
            "scale" => 1500
        ),
        "TST_GRI_01" => Dict{String,Any}(
            "type" => "GridConnection",
            "medium" => "m_c_g_natgas",
            "control_refs" => [],
            "output_refs" => ["TST_GB_01"],
            "is_source" => true,
        ),
        "TST_GB_01" => Dict{String,Any}(
            "type" => "FuelBoiler",
            "m_fuel_in" => "m_c_g_natgas",
            "control_refs" => ["TST_DEM_01"],
            "output_refs" => ["TST_DEM_01"],
            "strategy" => Dict{String,Any}(
                "name" => "demand_driven",
            ),
            "power_th" => 12000,
            "efficiency_fuel_in" => "const:1.0",
        ),
    )

    simulation_parameters = Dict{String,Any}(
        "time_step_seconds" => 900,
        "time" => 0,
        "epsilon" => 1e-9,
        "is_first_timestep" => true
    )

    components = Resie.load_components(components_config, simulation_parameters)
    gasboiler = components["TST_GB_01"]
    grid = components["TST_GRI_01"]
    demand = components["TST_DEM_01"]

    # test if correct demand is processed by gasboiler to make sure that the information of
    # demand is transported through Interface.

    # first time step: demand is exactly the max power of GasBoiler

    for unit in values(components)
        EnergySystems.reset(unit)
    end

    EnergySystems.control(demand, components, simulation_parameters)
    @test demand.demand == EnergySystems.get_max_energy(demand.input_interfaces[demand.medium].max_energy)
    @test demand.temperature == demand.input_interfaces[demand.medium].temperature_min
    EnergySystems.reset(demand)
    demand.constant_demand = 12000
    demand.constant_temperature = 85
    EnergySystems.control(demand, components, simulation_parameters)

    EnergySystems.control(gasboiler, components, simulation_parameters)
    EnergySystems.control(grid, components, simulation_parameters)
    @test EnergySystems.get_max_energy(grid.output_interfaces[grid.medium].max_energy) == Inf

    # no processing so far, balance is zero, energy_potential is non-zero
    exchanges = EnergySystems.balance_on(gasboiler.output_interfaces[gasboiler.m_heat_out], gasboiler)
    @test EnergySystems.balance(exchanges) ≈ 0.0
    @test EnergySystems.energy_potential(exchanges) ≈ 12000/4

    # no processing so far, balance is zero, energy_potential is non-zero
    exchanges = EnergySystems.balance_on(gasboiler.input_interfaces[gasboiler.m_fuel_in], gasboiler)
    @test EnergySystems.balance(exchanges) ≈ 0.0
    @test EnergySystems.energy_potential(exchanges) == -Inf

    EnergySystems.process(demand, simulation_parameters)
    @test demand.input_interfaces[demand.medium].balance ≈ -12000/4
    @test demand.input_interfaces[demand.medium].temperature_min == 85

    # demand was processed --> energy_potential should be zero, but not the balance
    exchanges = EnergySystems.balance_on(gasboiler.output_interfaces[gasboiler.m_heat_out], gasboiler)
    @test EnergySystems.balance(exchanges) ≈ -12000/4
    @test EnergySystems.energy_potential(exchanges) ≈ 0.0

    EnergySystems.process(gasboiler, simulation_parameters)
    @test gasboiler.output_interfaces[gasboiler.m_heat_out].balance ≈ 0
    @test gasboiler.output_interfaces[gasboiler.m_heat_out].sum_abs_change ≈ 2*12000/4
    @test gasboiler.input_interfaces[gasboiler.m_fuel_in].balance ≈ -12000/4
    @test demand.input_interfaces[demand.medium].balance ≈ 0
    @test demand.input_interfaces[demand.medium].sum_abs_change ≈ 2*12000/4

    # demand and gasboiler were processed --> energy_potential should be zero, but not the balance
    exchanges = EnergySystems.balance_on(gasboiler.input_interfaces[gasboiler.m_fuel_in], gasboiler)
    @test EnergySystems.balance(exchanges) ≈ -12000/4
    @test EnergySystems.energy_potential(exchanges) ≈ 0.0

    EnergySystems.process(grid, simulation_parameters)
    @test grid.output_interfaces[grid.medium].balance ≈ 0
    @test grid.output_interfaces[grid.medium].sum_abs_change ≈ 2*12000/4

    # second step: demand is above max power of GasBoiler

    for unit in values(components)
        EnergySystems.reset(unit)
    end

    EnergySystems.control(demand, components, simulation_parameters)
    @test demand.demand == EnergySystems.get_max_energy(demand.input_interfaces[demand.medium].max_energy)
    @test demand.temperature == demand.input_interfaces[demand.medium].temperature_min
    EnergySystems.reset(demand)
    demand.constant_demand = 15000
    demand.constant_temperature = 85
    EnergySystems.control(demand, components, simulation_parameters)

    EnergySystems.control(gasboiler, components, simulation_parameters)
    EnergySystems.control(grid, components, simulation_parameters)
    @test EnergySystems.get_max_energy(grid.output_interfaces[grid.medium].max_energy) == Inf

    # no processing so far, balance is zero, energy_potential is non-zero
    exchanges = EnergySystems.balance_on(gasboiler.output_interfaces[gasboiler.m_heat_out], gasboiler)
    @test EnergySystems.balance(exchanges) ≈ 0.0
    @test EnergySystems.energy_potential(exchanges) ≈ 15000/4

    # no processing so far, balance is zero, energy_potential is non-zero
    exchanges = EnergySystems.balance_on(gasboiler.input_interfaces[gasboiler.m_fuel_in], gasboiler)
    @test EnergySystems.balance(exchanges) ≈ 0.0
    @test EnergySystems.energy_potential(exchanges) == -Inf

    EnergySystems.process(demand, simulation_parameters)
    @test demand.input_interfaces[demand.medium].balance ≈ -15000/4
    @test demand.input_interfaces[demand.medium].temperature_min == 85

    # demand was processed --> energy_potential should be zero, but not the balance
    exchanges = EnergySystems.balance_on(gasboiler.output_interfaces[gasboiler.m_heat_out], gasboiler)
    @test EnergySystems.balance(exchanges) ≈ -15000/4
    @test EnergySystems.energy_potential(exchanges) ≈ 0.0

    EnergySystems.process(gasboiler, simulation_parameters)
    @test gasboiler.output_interfaces[gasboiler.m_heat_out].balance ≈ -15000/4 + 12000/4
    @test gasboiler.output_interfaces[gasboiler.m_heat_out].sum_abs_change ≈ 12000/4 + 15000/4
    @test gasboiler.input_interfaces[gasboiler.m_fuel_in].balance ≈ -12000/4
    @test demand.input_interfaces[demand.medium].balance ≈ -15000/4 + 12000/4
    @test demand.input_interfaces[demand.medium].sum_abs_change ≈ 15000/4 + 12000/4

    # demand and gasboiler were processed --> energy_potential should be zero, but not the balance
    exchanges = EnergySystems.balance_on(gasboiler.input_interfaces[gasboiler.m_fuel_in], gasboiler)
    @test EnergySystems.balance(exchanges) ≈ -12000/4
    @test EnergySystems.energy_potential(exchanges) ≈ 0.0

    EnergySystems.process(grid, simulation_parameters)
    @test grid.output_interfaces[grid.medium].balance ≈ 0
    @test grid.output_interfaces[grid.medium].sum_abs_change ≈ 2*12000/4

    # third step: demand is below max power of GasBoiler

    for unit in values(components)
        EnergySystems.reset(unit)
    end

    EnergySystems.control(demand, components, simulation_parameters)
    @test demand.demand == EnergySystems.get_max_energy(demand.input_interfaces[demand.medium].max_energy)
    @test demand.temperature == demand.input_interfaces[demand.medium].temperature_min
    EnergySystems.reset(demand)
    demand.constant_demand = 10000
    demand.constant_temperature = 85
    EnergySystems.control(demand, components, simulation_parameters)

    EnergySystems.control(gasboiler, components, simulation_parameters)
    EnergySystems.control(grid, components, simulation_parameters)
    @test EnergySystems.get_max_energy(grid.output_interfaces[grid.medium].max_energy) == Inf

    # no processing so far, balance is zero, energy_potential is non-zero
    exchanges = EnergySystems.balance_on(gasboiler.output_interfaces[gasboiler.m_heat_out], gasboiler)
    @test EnergySystems.balance(exchanges) ≈ 0.0
    @test EnergySystems.energy_potential(exchanges) ≈ 10000/4

    # no processing so far, balance is zero, energy_potential is non-zero
    exchanges = EnergySystems.balance_on(gasboiler.input_interfaces[gasboiler.m_fuel_in], gasboiler)
    @test EnergySystems.balance(exchanges) ≈ 0.0
    @test EnergySystems.energy_potential(exchanges) == -Inf

    EnergySystems.process(demand, simulation_parameters)
    @test demand.input_interfaces[demand.medium].balance ≈ -10000/4
    @test demand.input_interfaces[demand.medium].temperature_min == 85

    # demand was processed --> energy_potential should be zero, but not the balance
    exchanges = EnergySystems.balance_on(gasboiler.output_interfaces[gasboiler.m_heat_out], gasboiler)
    @test EnergySystems.balance(exchanges) ≈ -10000/4
    @test EnergySystems.energy_potential(exchanges) ≈ 0.0

    EnergySystems.process(gasboiler, simulation_parameters)
    @test gasboiler.output_interfaces[gasboiler.m_heat_out].balance ≈ 0
    @test gasboiler.output_interfaces[gasboiler.m_heat_out].sum_abs_change ≈ 2*10000/4
    @test gasboiler.input_interfaces[gasboiler.m_fuel_in].balance ≈ -10000/4
    @test demand.input_interfaces[demand.medium].balance ≈ 0
    @test demand.input_interfaces[demand.medium].sum_abs_change ≈ 2*10000/4

    # demand and gasboiler were processed --> energy_potential should be zero, but not the balance
    exchanges = EnergySystems.balance_on(gasboiler.input_interfaces[gasboiler.m_fuel_in], gasboiler)
    @test EnergySystems.balance(exchanges) ≈ -10000/4
    @test EnergySystems.energy_potential(exchanges) ≈ 0.0

    EnergySystems.process(grid, simulation_parameters)
    @test grid.output_interfaces[grid.medium].balance ≈ 0
    @test grid.output_interfaces[grid.medium].sum_abs_change ≈ 2*10000/4

end


@testset "test_gasboiler_demand_driven_with_bus" begin
    test_gasboiler_demand_driven_with_bus()
    test_gasboiler_demand_driven_without_bus()
end