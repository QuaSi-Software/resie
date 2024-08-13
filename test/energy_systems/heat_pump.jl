using Debugger
using Test
using Resie
using Resie.EnergySystems
using Resie.Profiles

EnergySystems.set_timestep(900)

function get_config_heat_pump_1S1D()
    return Dict{String,Any}(
        "TST_DEM_01" => Dict{String,Any}(
            "type" => "Demand",
            "medium" => "m_h_w_ht1",
            "output_refs" => [],
            "energy_profile_file_path" => "./profiles/tests/demand_heating_energy.prf",
            "temperature_profile_file_path" => "./profiles/tests/demand_heating_temperature.prf",
            "scale" => 1500,
        ),
        "TST_SRC_01" => Dict{String,Any}(
            "type" => "BoundedSupply",
            "medium" => "m_h_w_lt1",
            "output_refs" => ["TST_HP_01"],
            "max_power_profile_file_path" => "./profiles/tests/demand_heating_energy.prf",
            "temperature_profile_file_path" => "./profiles/tests/demand_heating_temperature.prf",
            "scale" => 6000,
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
            "power_th" => 12000,
        ),
    )
end

function test_heat_pump_one_source_dynamic_cop()
    components_config = get_config_heat_pump_1S1D()

    simulation_parameters = Dict{String,Any}(
        "time_step_seconds" => 900,
        "time" => 0,
        "epsilon" => 1e-9,
        "is_first_timestep" => true,
    )

    components = Resie.load_components(components_config, simulation_parameters)
    heat_pump = components["TST_HP_01"]
    source = components["TST_SRC_01"]
    demand = components["TST_DEM_01"]
    grid = components["TST_GRI_01"]

    # first time step: demand is below max power of source (adjusted for additional input
    # of electricity), small delta T leads to high COP = 12.725999999999999

    for unit in values(components)
        EnergySystems.reset(unit)
    end

    demand.constant_demand = 900.0 * 4
    demand.constant_temperature = 45.0
    EnergySystems.control(demand, components, simulation_parameters)

    @test heat_pump.input_interfaces[heat_pump.m_heat_in].temperature_min === nothing

    source.constant_power = 5000.0
    source.constant_temperature = 35
    EnergySystems.control(source, components, simulation_parameters)

    source.output_interfaces[source.medium].temperature_max = 35
    source.output_interfaces[source.medium].max_energy.max_energy[1] = 5000 / 4

    EnergySystems.control(heat_pump, components, simulation_parameters)
    EnergySystems.control(grid, components, simulation_parameters)

    EnergySystems.process(demand, simulation_parameters)
    @test demand.input_interfaces[demand.medium].balance ≈ -900
    @test demand.input_interfaces[demand.medium].temperature_min == 45

    EnergySystems.process(heat_pump, simulation_parameters)
    @test heat_pump.output_interfaces[heat_pump.m_heat_out].balance ≈ 0
    @test heat_pump.output_interfaces[heat_pump.m_heat_out].sum_abs_change ≈ 1800
    @test heat_pump.output_interfaces[heat_pump.m_heat_out].temperature_min == 45
    @test heat_pump.input_interfaces[heat_pump.m_el_in].balance ≈ -70.7213578500708
    @test heat_pump.input_interfaces[heat_pump.m_el_in].temperature_min === nothing
    @test heat_pump.input_interfaces[heat_pump.m_heat_in].balance ≈ -829.2786421499292
    @test heat_pump.input_interfaces[heat_pump.m_heat_in].temperature_max == 35

    EnergySystems.process(source, simulation_parameters)
    @test source.output_interfaces[source.medium].balance ≈ 0
    @test source.output_interfaces[source.medium].sum_abs_change ≈ 1658.5572842998583
    @test source.output_interfaces[source.medium].temperature_max == 35

    EnergySystems.process(grid, simulation_parameters)
    @test grid.output_interfaces[grid.medium].balance ≈ 0
    @test grid.output_interfaces[grid.medium].sum_abs_change ≈ 141.4427157001416
    @test grid.output_interfaces[grid.medium].temperature_max === nothing

    # second step: demand is above max power of source, big delta T leads to low COP = 3.4814999999999996

    for unit in values(components)
        EnergySystems.reset(unit)
    end

    demand.constant_demand = 2100 * 4
    demand.constant_temperature = 75.0
    EnergySystems.control(demand, components, simulation_parameters)

    @test heat_pump.input_interfaces[heat_pump.m_heat_in].temperature_min === nothing

    source.constant_power = 500 * 4
    source.constant_temperature = 35.0
    EnergySystems.control(source, components, simulation_parameters)

    source.output_interfaces[source.medium].temperature_max = 35
    source.output_interfaces[source.medium].max_energy.max_energy[1] = 500.0

    EnergySystems.control(heat_pump, components, simulation_parameters)
    EnergySystems.control(grid, components, simulation_parameters)

    EnergySystems.process(demand, simulation_parameters)
    @test demand.input_interfaces[demand.medium].balance ≈ -2100
    @test demand.input_interfaces[demand.medium].temperature_min == 75

    EnergySystems.process(heat_pump, simulation_parameters)
    @test heat_pump.output_interfaces[heat_pump.m_heat_out].balance ≈
          -2100 + 500 * (3.4814999999999996 / (3.4814999999999996 - 1))
    @test heat_pump.output_interfaces[heat_pump.m_heat_out].sum_abs_change ≈
          2100 + 500 * (3.4814999999999996 / (3.4814999999999996 - 1))
    @test heat_pump.output_interfaces[heat_pump.m_heat_out].temperature_min == 75
    @test heat_pump.input_interfaces[heat_pump.m_el_in].balance ≈
          -(500 * (3.4814999999999996 / (3.4814999999999996 - 1)) - 500)
    @test heat_pump.input_interfaces[heat_pump.m_el_in].temperature_min === nothing
    @test heat_pump.input_interfaces[heat_pump.m_heat_in].balance ≈ -500
    @test heat_pump.input_interfaces[heat_pump.m_heat_in].temperature_max == 35

    EnergySystems.process(source, simulation_parameters)
    @test source.output_interfaces[source.medium].balance ≈ 0
    @test source.output_interfaces[source.medium].sum_abs_change ≈ 1000
    @test source.output_interfaces[source.medium].temperature_max == 35

    EnergySystems.process(grid, simulation_parameters)
    @test grid.output_interfaces[grid.medium].balance ≈ 0
    @test grid.output_interfaces[grid.medium].sum_abs_change ≈
          2 * (500 * (3.4814999999999996 / (3.4814999999999996 - 1)) - 500)
    @test grid.output_interfaces[grid.medium].temperature_max === nothing
end

@testset "heat_pump_one_source_dynamic_cop" begin
    test_heat_pump_one_source_dynamic_cop()
end

function get_config_heat_pump_2S2D()
    return Dict{String,Any}(
        "TST_DEM_01" => Dict{String,Any}(
            "type" => "Demand",
            "medium" => "m_h_w_ht1",
            "output_refs" => [],
            "constant_demand" => 8000,
            "constant_temperature" => 70,
        ),
        "TST_DEM_02" => Dict{String,Any}(
            "type" => "Demand",
            "medium" => "m_h_w_ht1",
            "output_refs" => [],
            "constant_demand" => 4000,
            "constant_temperature" => 50,
        ),
        "TST_SRC_01" => Dict{String,Any}(
            "type" => "BoundedSupply",
            "medium" => "m_h_w_lt1",
            "output_refs" => ["TST_BUS_01"],
            "constant_power" => 4000,
            "constant_temperature" => 40,
        ),
        "TST_SRC_02" => Dict{String,Any}(
            "type" => "BoundedSupply",
            "medium" => "m_h_w_lt1",
            "output_refs" => ["TST_BUS_01"],
            "constant_power" => 8000,
            "constant_temperature" => 20,
        ),
        "TST_BUS_01" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_h_w_lt1",
            "connections" => Dict{String,Any}(
                "input_order" => ["TST_SRC_01",
                                  "TST_SRC_02"],
                "output_order" => ["TST_HP_01"],
            ),
        ),
        "TST_BUS_02" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_h_w_ht1",
            "connections" => Dict{String,Any}(
                "input_order" => ["TST_HP_01"],
                "output_order" => ["TST_DEM_01",
                                   "TST_DEM_02"],
            ),
        ),
        "TST_GRI_01" => Dict{String,Any}(
            "type" => "GridConnection",
            "medium" => "m_e_ac_230v",
            "output_refs" => ["TST_HP_01"],
            "is_source" => true,
        ),
        "TST_HP_01" => Dict{String,Any}(
            "type" => "HeatPump",
            "output_refs" => ["TST_BUS_02"],
            "power_th" => 12000,
        ),
    )
end

function test_heat_pump_2S2D_constant_cop()
    components_config = get_config_heat_pump_2S2D()

    simulation_parameters = Dict{String,Any}(
        "time_step_seconds" => 900,
        "time" => 0,
        "epsilon" => 1e-9,
        "is_first_timestep" => true,
    )

    components = Resie.load_components(components_config, simulation_parameters)
    heat_pump = components["TST_HP_01"]
    source_1 = components["TST_SRC_01"]
    source_2 = components["TST_SRC_02"]
    demand_1 = components["TST_DEM_01"]
    demand_2 = components["TST_DEM_01"]
    grid = components["TST_GRI_01"]
    bus_1 = components["TST_BUS_01"]
    bus_2 = components["TST_BUS_02"]

    heat_pump.constant_cop = 3.0
    heat_pump.bypass_cop = 10.0

    # first time step: no bypass, regular calculation

    for unit in values(components)
        EnergySystems.reset(unit)
    end
    for unit in values(components)
        EnergySystems.control(unit, components, simulation_parameters)
    end

    EnergySystems.process(demand_1, simulation_parameters)
    EnergySystems.process(demand_2, simulation_parameters)
    EnergySystems.process(heat_pump, simulation_parameters)
    EnergySystems.process(source_1, simulation_parameters)
    EnergySystems.process(source_2, simulation_parameters)
    EnergySystems.process(grid, simulation_parameters)

    EnergySystems.distribute!(bus_1)
    EnergySystems.distribute!(bus_2)

    @test heat_pump.input_interfaces[heat_pump.m_heat_in].sum_abs_change ≈ 2000.0 * 2
    @test heat_pump.input_interfaces[heat_pump.m_el_in].sum_abs_change ≈ 1000.0 * 2
    @test heat_pump.output_interfaces[heat_pump.m_heat_out].sum_abs_change ≈ 3000.0 * 2

    # second time step: bypass would be applicable, but is ignored with constant cop
    source_1.constant_temperature = 71.0

    for unit in values(components)
        EnergySystems.reset(unit)
    end
    for unit in values(components)
        EnergySystems.control(unit, components, simulation_parameters)
    end

    EnergySystems.process(demand_1, simulation_parameters)
    EnergySystems.process(demand_2, simulation_parameters)
    EnergySystems.process(heat_pump, simulation_parameters)
    EnergySystems.process(source_1, simulation_parameters)
    EnergySystems.process(source_2, simulation_parameters)
    EnergySystems.process(grid, simulation_parameters)

    EnergySystems.distribute!(bus_1)
    EnergySystems.distribute!(bus_2)

    @test heat_pump.input_interfaces[heat_pump.m_heat_in].sum_abs_change ≈ 2000.0 * 2
    @test heat_pump.input_interfaces[heat_pump.m_el_in].sum_abs_change ≈ 1000.0 * 2
    @test heat_pump.output_interfaces[heat_pump.m_heat_out].sum_abs_change ≈ 3000.0 * 2
end

@testset "heat_pump_2S2D_constant_cop" begin
    test_heat_pump_2S2D_constant_cop()
end

function test_heat_pump_2S2D_dynamic_cop()
    components_config = get_config_heat_pump_2S2D()

    simulation_parameters = Dict{String,Any}(
        "time_step_seconds" => 900,
        "time" => 0,
        "epsilon" => 1e-9,
        "is_first_timestep" => true,
    )

    components = Resie.load_components(components_config, simulation_parameters)
    heat_pump = components["TST_HP_01"]
    source_1 = components["TST_SRC_01"]
    source_2 = components["TST_SRC_02"]
    demand_1 = components["TST_DEM_01"]
    demand_2 = components["TST_DEM_01"]
    grid = components["TST_GRI_01"]
    bus_1 = components["TST_BUS_01"]
    bus_2 = components["TST_BUS_02"]

    heat_pump.bypass_cop = 10.0

    # first time step: no bypass, regular calculation
    for unit in values(components)
        EnergySystems.reset(unit)
    end
    for unit in values(components)
        EnergySystems.control(unit, components, simulation_parameters)
    end

    EnergySystems.process(demand_1, simulation_parameters)
    EnergySystems.process(demand_2, simulation_parameters)
    EnergySystems.process(heat_pump, simulation_parameters)
    EnergySystems.process(source_1, simulation_parameters)
    EnergySystems.process(source_2, simulation_parameters)
    EnergySystems.process(grid, simulation_parameters)

    EnergySystems.distribute!(bus_1)
    EnergySystems.distribute!(bus_2)

    # cop_1_1 = 4.5753333333333
    # cop_1_2 = 12.926
    # cop_2_1 = 2.7452
    # cop_2_2 = 4.3086666666666
    @test heat_pump.input_interfaces[heat_pump.m_heat_in].sum_abs_change ≈ 4187.2924962864345
    @test heat_pump.input_interfaces[heat_pump.m_el_in].sum_abs_change ≈ 1812.707530713566
    @test heat_pump.output_interfaces[heat_pump.m_heat_out].sum_abs_change ≈ 3000.0 * 2
    @test heat_pump.output_interfaces[heat_pump.m_heat_out].sum_abs_change /
          heat_pump.input_interfaces[heat_pump.m_el_in].sum_abs_change ≈ 3.309965837477

    # second time step: bypass should reduce required electricity
    source_1.constant_temperature = 71.0

    for unit in values(components)
        EnergySystems.reset(unit)
    end
    for unit in values(components)
        EnergySystems.control(unit, components, simulation_parameters)
    end

    EnergySystems.process(demand_1, simulation_parameters)
    EnergySystems.process(demand_2, simulation_parameters)
    EnergySystems.process(heat_pump, simulation_parameters)
    EnergySystems.process(source_1, simulation_parameters)
    EnergySystems.process(source_2, simulation_parameters)
    EnergySystems.process(grid, simulation_parameters)

    EnergySystems.distribute!(bus_1)
    EnergySystems.distribute!(bus_2)

    @test heat_pump.input_interfaces[heat_pump.m_heat_in].sum_abs_change ≈ 4401.638415335049
    @test heat_pump.input_interfaces[heat_pump.m_el_in].sum_abs_change ≈ 1598.3615846649507
    @test heat_pump.output_interfaces[heat_pump.m_heat_out].sum_abs_change ≈ 3000.0 * 2
    @test heat_pump.output_interfaces[heat_pump.m_heat_out].sum_abs_change /
          heat_pump.input_interfaces[heat_pump.m_el_in].sum_abs_change ≈ 3.753843972206
end

@testset "heat_pump_2S2D_dynamic_cop" begin
    test_heat_pump_2S2D_dynamic_cop()
end

function test_heat_pump_2S2D_reorder_inputs()
    components_config = get_config_heat_pump_2S2D()

    components_config["TST_HP_01"]["control_modules"] = [Dict{String,Any}(
                                                             "name" => "temperature_sorting",
                                                         )]

    simulation_parameters = Dict{String,Any}(
        "time_step_seconds" => 900,
        "time" => 0,
        "epsilon" => 1e-9,
        "is_first_timestep" => true,
    )

    components = Resie.load_components(components_config, simulation_parameters)
    heat_pump = components["TST_HP_01"]
    source_1 = components["TST_SRC_01"]
    source_2 = components["TST_SRC_02"]
    demand_1 = components["TST_DEM_01"]
    demand_2 = components["TST_DEM_01"]
    grid = components["TST_GRI_01"]
    bus_1 = components["TST_BUS_01"]
    bus_2 = components["TST_BUS_02"]

    heat_pump.constant_cop = 3.0
    source_1.constant_power = 40000

    # first time step: highest priority source also has highest temperature and cover all demand
    for unit in values(components)
        EnergySystems.reset(unit)
    end
    for unit in values(components)
        EnergySystems.control(unit, components, simulation_parameters)
    end

    EnergySystems.process(demand_1, simulation_parameters)
    EnergySystems.process(demand_2, simulation_parameters)
    EnergySystems.process(heat_pump, simulation_parameters)
    EnergySystems.process(source_1, simulation_parameters)
    EnergySystems.process(source_2, simulation_parameters)
    EnergySystems.process(grid, simulation_parameters)

    EnergySystems.distribute!(bus_1)
    EnergySystems.distribute!(bus_2)

    @test heat_pump.input_interfaces[heat_pump.m_heat_in].sum_abs_change ≈ 2000.0 * 2
    @test heat_pump.input_interfaces[heat_pump.m_el_in].sum_abs_change ≈ 1000.0 * 2
    @test heat_pump.output_interfaces[heat_pump.m_heat_out].sum_abs_change ≈ 3000.0 * 2
    @test source_1.output_interfaces[source_1.medium].sum_abs_change ≈ 2000.0 * 2
    @test source_2.output_interfaces[source_2.medium].sum_abs_change ≈ 0.0

    # second time step: source 2, having a lower priority, has the highest temperature and
    # can cover all demand
    source_2.constant_power = 40000
    source_2.constant_temperature = 41.0

    for unit in values(components)
        EnergySystems.reset(unit)
    end
    for unit in values(components)
        EnergySystems.control(unit, components, simulation_parameters)
    end

    EnergySystems.process(demand_1, simulation_parameters)
    EnergySystems.process(demand_2, simulation_parameters)
    EnergySystems.process(heat_pump, simulation_parameters)
    EnergySystems.process(source_1, simulation_parameters)
    EnergySystems.process(source_2, simulation_parameters)
    EnergySystems.process(grid, simulation_parameters)

    EnergySystems.distribute!(bus_1)
    EnergySystems.distribute!(bus_2)

    @test heat_pump.input_interfaces[heat_pump.m_heat_in].sum_abs_change ≈ 2000.0 * 2
    @test heat_pump.input_interfaces[heat_pump.m_el_in].sum_abs_change ≈ 1000.0 * 2
    @test heat_pump.output_interfaces[heat_pump.m_heat_out].sum_abs_change ≈ 3000.0 * 2
    @test source_1.output_interfaces[source_1.medium].sum_abs_change ≈ 0.0
    @test source_2.output_interfaces[source_2.medium].sum_abs_change ≈ 2000.0 * 2
end

@testset "heat_pump_2S2D_reorder_inputs" begin
    test_heat_pump_2S2D_reorder_inputs()
end
