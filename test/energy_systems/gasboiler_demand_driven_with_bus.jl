using Debugger
using Test
using Resie
using Resie.EnergySystems
using Resie.Profiles

watt_to_wh = function (watts::Float64)
    watts * 900 / 3600.0
end

function test_gasboiler_demand_driven_with_bus()
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
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "production_refs" => ["TST_DEM_01"],
            "connection_matrix" => Dict{String,Any}(
                "input_order" => ["TST_GB_01"],
                "output_order" => ["TST_DEM_01"]
            )
        ),
        "TST_GRI_01" => Dict{String,Any}(
            "type" => "GridConnection",
            "medium" => "m_c_g_natgas",
            "control_refs" => [],
            "production_refs" => ["TST_GB_01"],
            "is_source" => true,
        ),
        "TST_GB_01" => Dict{String,Any}(
            "type" => "GasBoiler",
            "control_refs" => ["TST_BUS_01"],
            "production_refs" => ["TST_BUS_01"],
            "strategy" => Dict{String,Any}(
                "name" => "demand_driven",
            ),
            "power" => 12000
        ),
    )
    systems = Resie.load_systems(systems_config)
    gasboiler = systems["TST_GB_01"]
    grid = systems["TST_GRI_01"]
    demand = systems["TST_DEM_01"]
    bus = systems["TST_BUS_01"]

    simulation_parameters = Dict{String,Any}(
        "time_step_seconds" => 900,
        "time" => 0,
        "epsilon" => 1e-9
    )

    # test if correct demand is produced by gasboiler to make sure that the information of
    # demand is transported through Bus.

    # first time step: demand is exactly the max power of GasBoiler  

    for unit in values(systems)
        EnergySystems.reset(unit)
    end

    EnergySystems.control(demand, systems, simulation_parameters)
    @test demand.load == demand.input_interfaces[demand.medium].max_energy
    @test demand.temperature == demand.input_interfaces[demand.medium].temperature
    demand.load = 12000/4
    demand.temperature = 85
    demand.input_interfaces[demand.medium].max_energy = 12000/4
    demand.input_interfaces[demand.medium].temperature = 85

    EnergySystems.control(bus, systems, simulation_parameters)
    EnergySystems.control(gasboiler, systems, simulation_parameters)
    EnergySystems.control(grid, systems, simulation_parameters)
    @test grid.output_interfaces[grid.medium].max_energy == Inf

    # no production so far, balance is zero, energy_potential is non-zero
    exchange = EnergySystems.balance_on(gasboiler.output_interfaces[bus.medium], bus)
    @test exchange.balance ≈ 0.0
    @test exchange.storage_potential ≈ 0.0
    @test exchange.energy_potential ≈ -12000/4
    
    # no production so far, balance is zero, energy_potential is zero
    exchange = EnergySystems.balance_on(demand.input_interfaces[bus.medium], bus)
    @test exchange.balance ≈ 0.0
    @test exchange.storage_potential ≈ 0.0
    @test exchange.energy_potential ≈ 0.0

    EnergySystems.produce(demand, simulation_parameters, watt_to_wh)
    @test demand.input_interfaces[demand.medium].balance ≈ -12000/4
    @test demand.input_interfaces[demand.medium].temperature ≈ 85

    # demand was produced --> energy_potential should be zero, but not the balance
    exchange = EnergySystems.balance_on(gasboiler.output_interfaces[bus.medium], bus)
    @test exchange.balance ≈ -12000/4  
    @test exchange.storage_potential ≈ 0.0
    @test exchange.energy_potential ≈ 0.0

    EnergySystems.produce(bus, simulation_parameters, watt_to_wh)
    EnergySystems.produce(gasboiler, simulation_parameters, watt_to_wh)
    @test gasboiler.output_interfaces[gasboiler.m_heat_out].balance ≈ 12000/4
    @test gasboiler.input_interfaces[gasboiler.m_gas_in].balance ≈ -12000/4

    EnergySystems.produce(grid, simulation_parameters, watt_to_wh)
    @test grid.output_interfaces[grid.medium].balance ≈ 0
    @test grid.output_interfaces[grid.medium].sum_abs_change ≈ 2*12000/4

    EnergySystems.distribute!(bus)
    @test gasboiler.output_interfaces[gasboiler.m_heat_out].balance ≈ 0
    @test gasboiler.output_interfaces[gasboiler.m_heat_out].sum_abs_change ≈ 2*12000/4
    @test demand.input_interfaces[demand.medium].balance ≈ 0
    @test demand.input_interfaces[demand.medium].sum_abs_change ≈ 2*12000/4
    @test bus.remainder ≈ 0

    # second step: demand is above max power of GasBoiler 

    for unit in values(systems)
        EnergySystems.reset(unit)
    end

    EnergySystems.control(demand, systems, simulation_parameters)
    @test demand.load == demand.input_interfaces[demand.medium].max_energy
    @test demand.temperature == demand.input_interfaces[demand.medium].temperature
    demand.load = 15000/4
    demand.temperature = 85
    demand.input_interfaces[demand.medium].max_energy = 15000/4
    demand.input_interfaces[demand.medium].temperature = 85

    EnergySystems.control(bus, systems, simulation_parameters)
    EnergySystems.control(gasboiler, systems, simulation_parameters)
    EnergySystems.control(grid, systems, simulation_parameters)
    @test grid.output_interfaces[grid.medium].max_energy == Inf

    # no production so far, balance is zero, energy_potential is non-zero
    exchange = EnergySystems.balance_on(gasboiler.output_interfaces[bus.medium], bus)
    @test exchange.balance ≈ 0.0
    @test exchange.storage_potential ≈ 0.0
    @test exchange.energy_potential ≈ -15000/4
    
    # no production so far, balance is zero, energy_potential is zero
    exchange = EnergySystems.balance_on(demand.input_interfaces[bus.medium], bus)
    @test exchange.balance ≈ 0.0
    @test exchange.storage_potential ≈ 0.0
    @test exchange.energy_potential ≈ 0.0

    EnergySystems.produce(demand, simulation_parameters, watt_to_wh)
    @test demand.input_interfaces[demand.medium].balance ≈ -15000/4
    @test demand.input_interfaces[demand.medium].temperature ≈ 85

    # demand was produced --> energy_potential should be zero, but not the balance
    exchange = EnergySystems.balance_on(gasboiler.output_interfaces[bus.medium], bus)
    @test exchange.balance ≈ -15000/4  
    @test exchange.storage_potential ≈ 0.0
    @test exchange.energy_potential ≈ 0.0

    EnergySystems.produce(bus, simulation_parameters, watt_to_wh)
    EnergySystems.produce(gasboiler, simulation_parameters, watt_to_wh)
    @test gasboiler.output_interfaces[gasboiler.m_heat_out].balance ≈ 12000/4
    @test gasboiler.input_interfaces[gasboiler.m_gas_in].balance ≈ -12000/4

    EnergySystems.produce(grid, simulation_parameters, watt_to_wh)
    @test grid.output_interfaces[grid.medium].balance ≈ 0
    @test grid.output_interfaces[grid.medium].sum_abs_change ≈ 2*12000/4

    EnergySystems.distribute!(bus)
    @test gasboiler.output_interfaces[gasboiler.m_heat_out].balance ≈ 0
    @test gasboiler.output_interfaces[gasboiler.m_heat_out].sum_abs_change ≈ 2*12000/4
    @test demand.input_interfaces[demand.medium].balance ≈ 0 #-15000/4 + 12000/4   # is Bus intended to supply all demand???
    @test demand.input_interfaces[demand.medium].sum_abs_change ≈ 15000/4 + 15000/4
    @test bus.remainder ≈ -15000/4 + 12000/4

    # third step: demand is below max power of GasBoiler 

    for unit in values(systems)
        EnergySystems.reset(unit)
    end

    EnergySystems.control(demand, systems, simulation_parameters)
    @test demand.load == demand.input_interfaces[demand.medium].max_energy
    @test demand.temperature == demand.input_interfaces[demand.medium].temperature
    demand.load = 10000/4
    demand.temperature = 85
    demand.input_interfaces[demand.medium].max_energy = 10000/4
    demand.input_interfaces[demand.medium].temperature = 85

    EnergySystems.control(bus, systems, simulation_parameters)
    EnergySystems.control(gasboiler, systems, simulation_parameters)
    EnergySystems.control(grid, systems, simulation_parameters)
    @test grid.output_interfaces[grid.medium].max_energy == Inf

    # no production so far, balance is zero, energy_potential is non-zero
    exchange = EnergySystems.balance_on(gasboiler.output_interfaces[bus.medium], bus)
    @test exchange.balance ≈ 0.0
    @test exchange.storage_potential ≈ 0.0
    @test exchange.energy_potential ≈ -10000/4
    
    # no production so far, balance is zero, energy_potential is zero
    exchange = EnergySystems.balance_on(demand.input_interfaces[bus.medium], bus)
    @test exchange.balance ≈ 0.0
    @test exchange.storage_potential ≈ 0.0
    @test exchange.energy_potential ≈ 0.0

    EnergySystems.produce(demand, simulation_parameters, watt_to_wh)
    @test demand.input_interfaces[demand.medium].balance ≈ -10000/4
    @test demand.input_interfaces[demand.medium].temperature ≈ 85

    # demand was produced --> energy_potential should be zero, but not the balance
    exchange = EnergySystems.balance_on(gasboiler.output_interfaces[bus.medium], bus)
    @test exchange.balance ≈ -10000/4  
    @test exchange.storage_potential ≈ 0.0
    @test exchange.energy_potential ≈ 0.0

    EnergySystems.produce(bus, simulation_parameters, watt_to_wh)
    EnergySystems.produce(gasboiler, simulation_parameters, watt_to_wh)
    @test gasboiler.output_interfaces[gasboiler.m_heat_out].balance ≈ 10000/4
    @test gasboiler.input_interfaces[gasboiler.m_gas_in].balance ≈ -10000/4

    EnergySystems.produce(grid, simulation_parameters, watt_to_wh)
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
        "TST_GRI_01" => Dict{String,Any}(
            "type" => "GridConnection",
            "medium" => "m_c_g_natgas",
            "control_refs" => [],
            "production_refs" => ["TST_GB_01"],
            "is_source" => true,
        ),
        "TST_GB_01" => Dict{String,Any}(
            "type" => "GasBoiler",
            "control_refs" => ["TST_DEM_01"],
            "production_refs" => ["TST_DEM_01"],
            "strategy" => Dict{String,Any}(
                "name" => "demand_driven",
            ),
            "power" => 12000
        ),
    )
    systems = Resie.load_systems(systems_config)
    gasboiler = systems["TST_GB_01"]
    grid = systems["TST_GRI_01"]
    demand = systems["TST_DEM_01"]

    simulation_parameters = Dict{String,Any}(
        "time_step_seconds" => 900,
        "time" => 0,
        "epsilon" => 1e-9
    )

    # test if correct demand is produced by gasboiler to make sure that the information of
    # demand is transported through Interface.

    # first time step: demand is exactly the max power of GasBoiler  

    for unit in values(systems)
        EnergySystems.reset(unit)
    end

    EnergySystems.control(demand, systems, simulation_parameters)
    @test demand.load == demand.input_interfaces[demand.medium].max_energy
    @test demand.temperature == demand.input_interfaces[demand.medium].temperature
    demand.load = 12000/4
    demand.temperature = 85
    demand.input_interfaces[demand.medium].max_energy = 12000/4
    demand.input_interfaces[demand.medium].temperature = 85

    EnergySystems.control(gasboiler, systems, simulation_parameters)
    EnergySystems.control(grid, systems, simulation_parameters)
    @test grid.output_interfaces[grid.medium].max_energy == Inf

    # no production so far, balance is zero, energy_potential is non-zero
    exchange = EnergySystems.balance_on(gasboiler.output_interfaces[gasboiler.m_heat_out], gasboiler)
    @test exchange.balance ≈ 0.0
    @test exchange.storage_potential ≈ 0.0
    @test exchange.energy_potential ≈ 12000/4
    
    # no production so far, balance is zero, energy_potential is non-zero
    exchange = EnergySystems.balance_on(gasboiler.input_interfaces[gasboiler.m_gas_in], gasboiler)
    @test exchange.balance ≈ 0.0
    @test exchange.storage_potential ≈ 0.0
    @test exchange.energy_potential == -Inf

    EnergySystems.produce(demand, simulation_parameters, watt_to_wh)
    @test demand.input_interfaces[demand.medium].balance ≈ -12000/4
    @test demand.input_interfaces[demand.medium].temperature == 85

    # demand was produced --> energy_potential should be zero, but not the balance
    exchange = EnergySystems.balance_on(gasboiler.output_interfaces[gasboiler.m_heat_out], gasboiler)
    @test exchange.balance ≈ -12000/4  
    @test exchange.storage_potential ≈ 0.0
    @test exchange.energy_potential ≈ 0.0

    EnergySystems.produce(gasboiler, simulation_parameters, watt_to_wh)
    @test gasboiler.output_interfaces[gasboiler.m_heat_out].balance ≈ 0
    @test gasboiler.output_interfaces[gasboiler.m_heat_out].sum_abs_change ≈ 2*12000/4
    @test gasboiler.input_interfaces[gasboiler.m_gas_in].balance ≈ -12000/4
    @test demand.input_interfaces[demand.medium].balance ≈ 0
    @test demand.input_interfaces[demand.medium].sum_abs_change ≈ 2*12000/4

    # demand and gasboiler were produced --> energy_potential should be zero, but not the balance
    exchange = EnergySystems.balance_on(gasboiler.input_interfaces[gasboiler.m_gas_in], gasboiler)
    @test exchange.balance ≈ -12000/4
    @test exchange.storage_potential ≈ 0.0
    @test exchange.energy_potential ≈ 0.0

    EnergySystems.produce(grid, simulation_parameters, watt_to_wh)
    @test grid.output_interfaces[grid.medium].balance ≈ 0
    @test grid.output_interfaces[grid.medium].sum_abs_change ≈ 2*12000/4

    # second step: demand is above max power of GasBoiler 

    for unit in values(systems)
        EnergySystems.reset(unit)
    end

    EnergySystems.control(demand, systems, simulation_parameters)
    @test demand.load == demand.input_interfaces[demand.medium].max_energy
    @test demand.temperature == demand.input_interfaces[demand.medium].temperature
    demand.load = 15000/4
    demand.temperature = 85
    demand.input_interfaces[demand.medium].max_energy = 15000/4
    demand.input_interfaces[demand.medium].temperature = 85

    EnergySystems.control(gasboiler, systems, simulation_parameters)
    EnergySystems.control(grid, systems, simulation_parameters)
    @test grid.output_interfaces[grid.medium].max_energy == Inf

    # no production so far, balance is zero, energy_potential is non-zero
    exchange = EnergySystems.balance_on(gasboiler.output_interfaces[gasboiler.m_heat_out], gasboiler)
    @test exchange.balance ≈ 0.0
    @test exchange.storage_potential ≈ 0.0
    @test exchange.energy_potential ≈ 15000/4
    
    # no production so far, balance is zero, energy_potential is non-zero
    exchange = EnergySystems.balance_on(gasboiler.input_interfaces[gasboiler.m_gas_in], gasboiler)
    @test exchange.balance ≈ 0.0
    @test exchange.storage_potential ≈ 0.0
    @test exchange.energy_potential == -Inf

    EnergySystems.produce(demand, simulation_parameters, watt_to_wh)
    @test demand.input_interfaces[demand.medium].balance ≈ -15000/4
    @test demand.input_interfaces[demand.medium].temperature == 85

    # demand was produced --> energy_potential should be zero, but not the balance
    exchange = EnergySystems.balance_on(gasboiler.output_interfaces[gasboiler.m_heat_out], gasboiler)
    @test exchange.balance ≈ -15000/4  
    @test exchange.storage_potential ≈ 0.0
    @test exchange.energy_potential ≈ 0.0

    EnergySystems.produce(gasboiler, simulation_parameters, watt_to_wh)
    @test gasboiler.output_interfaces[gasboiler.m_heat_out].balance ≈ -15000/4 + 12000/4
    @test gasboiler.output_interfaces[gasboiler.m_heat_out].sum_abs_change ≈ 12000/4 + 15000/4
    @test gasboiler.input_interfaces[gasboiler.m_gas_in].balance ≈ -12000/4
    @test demand.input_interfaces[demand.medium].balance ≈ -15000/4 + 12000/4   
    @test demand.input_interfaces[demand.medium].sum_abs_change ≈ 15000/4 + 12000/4

    # demand and gasboiler were produced --> energy_potential should be zero, but not the balance
    exchange = EnergySystems.balance_on(gasboiler.input_interfaces[gasboiler.m_gas_in], gasboiler)
    @test exchange.balance ≈ -12000/4
    @test exchange.storage_potential ≈ 0.0
    @test exchange.energy_potential ≈ 0.0

    EnergySystems.produce(grid, simulation_parameters, watt_to_wh)
    @test grid.output_interfaces[grid.medium].balance ≈ 0
    @test grid.output_interfaces[grid.medium].sum_abs_change ≈ 2*12000/4

    # third step: demand is below max power of GasBoiler 

    for unit in values(systems)
        EnergySystems.reset(unit)
    end

    EnergySystems.control(demand, systems, simulation_parameters)
    @test demand.load == demand.input_interfaces[demand.medium].max_energy
    @test demand.temperature == demand.input_interfaces[demand.medium].temperature
    demand.load = 10000/4
    demand.temperature = 85
    demand.input_interfaces[demand.medium].max_energy = 10000/4
    demand.input_interfaces[demand.medium].temperature = 85

    EnergySystems.control(gasboiler, systems, simulation_parameters)
    EnergySystems.control(grid, systems, simulation_parameters)
    @test grid.output_interfaces[grid.medium].max_energy == Inf

    # no production so far, balance is zero, energy_potential is non-zero
    exchange = EnergySystems.balance_on(gasboiler.output_interfaces[gasboiler.m_heat_out], gasboiler)
    @test exchange.balance ≈ 0.0
    @test exchange.storage_potential ≈ 0.0
    @test exchange.energy_potential ≈ 10000/4
    
    # no production so far, balance is zero, energy_potential is non-zero
    exchange = EnergySystems.balance_on(gasboiler.input_interfaces[gasboiler.m_gas_in], gasboiler)
    @test exchange.balance ≈ 0.0
    @test exchange.storage_potential ≈ 0.0
    @test exchange.energy_potential == -Inf

    EnergySystems.produce(demand, simulation_parameters, watt_to_wh)
    @test demand.input_interfaces[demand.medium].balance ≈ -10000/4
    @test demand.input_interfaces[demand.medium].temperature == 85

    # demand was produced --> energy_potential should be zero, but not the balance
    exchange = EnergySystems.balance_on(gasboiler.output_interfaces[gasboiler.m_heat_out], gasboiler)
    @test exchange.balance ≈ -10000/4  
    @test exchange.storage_potential ≈ 0.0
    @test exchange.energy_potential ≈ 0.0

    EnergySystems.produce(gasboiler, simulation_parameters, watt_to_wh)
    @test gasboiler.output_interfaces[gasboiler.m_heat_out].balance ≈ 0
    @test gasboiler.output_interfaces[gasboiler.m_heat_out].sum_abs_change ≈ 2*10000/4 
    @test gasboiler.input_interfaces[gasboiler.m_gas_in].balance ≈ -10000/4
    @test demand.input_interfaces[demand.medium].balance ≈ 0
    @test demand.input_interfaces[demand.medium].sum_abs_change ≈ 2*10000/4 

    # demand and gasboiler were produced --> energy_potential should be zero, but not the balance
    exchange = EnergySystems.balance_on(gasboiler.input_interfaces[gasboiler.m_gas_in], gasboiler)
    @test exchange.balance ≈ -10000/4
    @test exchange.storage_potential ≈ 0.0
    @test exchange.energy_potential ≈ 0.0

    EnergySystems.produce(grid, simulation_parameters, watt_to_wh)
    @test grid.output_interfaces[grid.medium].balance ≈ 0
    @test grid.output_interfaces[grid.medium].sum_abs_change ≈ 2*10000/4

end


@testset "test_gasboiler_demand_driven_with_bus" begin
    test_gasboiler_demand_driven_with_bus()
    test_gasboiler_demand_driven_without_bus()
end