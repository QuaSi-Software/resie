using Debugger
using Test
using Resie
using Resie.EnergySystems
using Resie.Profiles

EnergySystems.set_timestep(900)

function test_gas_boiler_demand_driven_plr()
    components_config = Dict{String,Any}(
        "TST_DEM_01" => Dict{String,Any}(
            "type" => "Demand",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "output_refs" => [],
            "static_demand" => 1000,
            "scale" => 1,
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
            "power" => 4000,
            "is_plr_dependant" => true,
            "max_thermal_efficiency" => 1.0,
        ),
    )
    components = Resie.load_components(components_config)
    gasboiler = components["TST_GB_01"]
    grid = components["TST_GRI_01"]
    demand = components["TST_DEM_01"]

    simulation_parameters = Dict{String,Any}(
        "time_step_seconds" => 900,
        "time" => 0,
        "epsilon" => 1e-9
    )

    # test if the thermal efficiency for a gas boiler in demand-driven mode is calculated
    # correctly depending on the part-load-ratio (PLR)

    # first time step: demand is exactly the max power of GasBoiler

    for unit in values(components)
        EnergySystems.reset(unit)
    end

    expected_efficiency = 1.0 

    EnergySystems.control(demand, components, simulation_parameters)
    EnergySystems.control(gasboiler, components, simulation_parameters)
    EnergySystems.control(grid, components, simulation_parameters)

    EnergySystems.process(demand, simulation_parameters)
    @test demand.input_interfaces[demand.medium].balance ≈ -1000

    EnergySystems.process(gasboiler, simulation_parameters)
    @test gasboiler.output_interfaces[gasboiler.m_heat_out].balance ≈ 0
    @test gasboiler.input_interfaces[gasboiler.m_fuel_in].balance ≈ -1000 / expected_efficiency
    @test demand.input_interfaces[demand.medium].balance ≈ 0

    EnergySystems.process(grid, simulation_parameters)
    @test grid.output_interfaces[grid.medium].balance ≈ 0

    # second time step: demand is half the max power of GasBoiler

    for unit in values(components)
        EnergySystems.reset(unit)
    end

    demand.static_demand = 500
    expected_efficiency = 0.7440249999999999

    EnergySystems.control(demand, components, simulation_parameters)
    EnergySystems.control(gasboiler, components, simulation_parameters)
    EnergySystems.control(grid, components, simulation_parameters)

    EnergySystems.process(demand, simulation_parameters)
    @test demand.input_interfaces[demand.medium].balance ≈ -500

    EnergySystems.process(gasboiler, simulation_parameters)
    @test gasboiler.output_interfaces[gasboiler.m_heat_out].balance ≈ 0
    @test gasboiler.input_interfaces[gasboiler.m_fuel_in].balance ≈ -500 / expected_efficiency
    @test demand.input_interfaces[demand.medium].balance ≈ 0

    EnergySystems.process(grid, simulation_parameters)
    @test grid.output_interfaces[grid.medium].balance ≈ 0
end

function test_gas_boiler_supply_driven_plr()
    components_config = Dict{String,Any}(
        "TST_DEM_01" => Dict{String,Any}(
            "type" => "BoundedSink",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "output_refs" => [],
            "max_power_profile_file_path" => "./profiles/tests/source_heat_max_power.prf",
            "scale" => 1.0,
        ),
        "TST_GRI_01" => Dict{String,Any}(
            "type" => "FixedSupply",
            "medium" => "m_c_g_natgas",
            "control_refs" => [],
            "output_refs" => ["TST_GB_01"],
            "is_source" => true,
            "energy_profile_file_path" => "./profiles/tests/demand_heating_energy.prf",
            "supply" => 0,
            "scale" => 1.0,
        ),
        "TST_GB_01" => Dict{String,Any}(
            "type" => "FuelBoiler",
            "m_fuel_in" => "m_c_g_natgas",
            "control_refs" => ["TST_DEM_01"],
            "output_refs" => ["TST_DEM_01"],
            "strategy" => Dict{String,Any}(
                "name" => "supply_driven",
            ),
            "power" => 4000,
            "is_plr_dependant" => true,
            "max_thermal_efficiency" => 1.0,
        ),
    )
    components = Resie.load_components(components_config)
    gasboiler = components["TST_GB_01"]
    grid = components["TST_GRI_01"]
    demand = components["TST_DEM_01"]

    simulation_parameters = Dict{String,Any}(
        "time_step_seconds" => 900,
        "time" => 0,
        "epsilon" => 1e-9
    )

    # test if the thermal efficiency for a gas boiler in supply-driven mode is calculated
    # correctly depending on the part-load-ratio (PLR)

    # first time step: demand is exactly the max power of GasBoiler

    for unit in values(components)
        EnergySystems.reset(unit)
    end

    expected_efficiency = 1.0

    EnergySystems.control(grid, components, simulation_parameters)
    grid.supply = 1000 / expected_efficiency
    grid.output_interfaces[grid.medium].max_energy = grid.supply
    EnergySystems.control(gasboiler, components, simulation_parameters)
    EnergySystems.control(demand, components, simulation_parameters)
    demand.input_interfaces[demand.medium].max_energy = 1000 * 10
    demand.max_energy = 1000 * 10

    EnergySystems.process(grid, simulation_parameters)
    @test grid.output_interfaces[grid.medium].balance ≈ 1000 / expected_efficiency

    EnergySystems.process(gasboiler, simulation_parameters)
    @test gasboiler.output_interfaces[gasboiler.m_heat_out].balance ≈ 1000
    @test gasboiler.input_interfaces[gasboiler.m_fuel_in].balance ≈ 0
    @test demand.input_interfaces[demand.medium].balance ≈ 1000

    EnergySystems.process(demand, simulation_parameters)
    @test demand.input_interfaces[demand.medium].balance ≈ 0

    # second time step: demand is half the max power of GasBoiler

    for unit in values(components)
        EnergySystems.reset(unit)
    end

    expected_efficiency = 0.744025

    EnergySystems.control(grid, components, simulation_parameters)
    grid.supply = 500 / expected_efficiency
    grid.output_interfaces[grid.medium].max_energy = grid.supply
    EnergySystems.control(gasboiler, components, simulation_parameters)
    EnergySystems.control(demand, components, simulation_parameters)
    demand.input_interfaces[demand.medium].max_energy = 500 * 10
    demand.max_energy = 500 * 10

    EnergySystems.process(grid, simulation_parameters)
    @test grid.output_interfaces[grid.medium].balance ≈ 500 / expected_efficiency
    @test gasboiler.output_interfaces[gasboiler.m_heat_out].balance ≈ 0

    EnergySystems.process(gasboiler, simulation_parameters)
    @test gasboiler.output_interfaces[gasboiler.m_heat_out].balance ≈ 500
    @test gasboiler.input_interfaces[gasboiler.m_fuel_in].balance ≈ 0
    @test demand.input_interfaces[demand.medium].balance ≈ 500

    EnergySystems.process(demand, simulation_parameters)
    @test gasboiler.output_interfaces[gasboiler.m_heat_out].balance ≈ 0
    @test demand.input_interfaces[demand.medium].balance ≈ 0
end

@testset "test_gas_boiler_demand_driven_plr" begin
    test_gas_boiler_demand_driven_plr()
    test_gas_boiler_supply_driven_plr()
end