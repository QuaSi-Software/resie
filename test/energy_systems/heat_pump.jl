using Debugger
using Test
using Resie
using Resie.EnergySystems
using Resie.Profiles

include("../test_util.jl")

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
    simulation_parameters = get_default_sim_params()

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

function test_heat_pump_1S1D_icing_losses()
    components_config = get_config_heat_pump_1S1D()
    simulation_parameters = get_default_sim_params()

    components = Resie.load_components(components_config, simulation_parameters)
    heat_pump = components["TST_HP_01"]
    source = components["TST_SRC_01"]
    demand = components["TST_DEM_01"]
    grid = components["TST_GRI_01"]

    heat_pump.consider_icing = true

    # first time step: source is warm enough to not incur icing losses

    for unit in values(components)
        EnergySystems.reset(unit)
    end

    demand.constant_demand = 900.0 * 4
    demand.constant_temperature = 45.0

    source.constant_power = 1000.0 * 4
    source.constant_temperature = 20.0

    EnergySystems.control(demand, components, simulation_parameters)
    EnergySystems.control(source, components, simulation_parameters)
    EnergySystems.control(heat_pump, components, simulation_parameters)
    EnergySystems.control(grid, components, simulation_parameters)
    EnergySystems.process(demand, simulation_parameters)
    EnergySystems.process(heat_pump, simulation_parameters)
    EnergySystems.process(source, simulation_parameters)
    EnergySystems.process(grid, simulation_parameters)

    heat_out = heat_pump.output_interfaces[heat_pump.m_heat_out].sum_abs_change * 0.5
    el_in = heat_pump.input_interfaces[heat_pump.m_el_in].sum_abs_change * 0.5
    heat_in = heat_pump.input_interfaces[heat_pump.m_heat_in].sum_abs_change * 0.5

    @test heat_out ≈ 900
    @test el_in + heat_in ≈ heat_out
    # COP without icing losses would be 5.0904
    @test heat_out / el_in ≈ 5.090384423755175

    # second step: source is in a range where icing losses are significant

    for unit in values(components)
        EnergySystems.reset(unit)
    end

    demand.constant_demand = 900.0 * 4
    demand.constant_temperature = 45.0

    source.constant_power = 1000.0 * 4
    source.constant_temperature = 5.0

    EnergySystems.control(demand, components, simulation_parameters)
    EnergySystems.control(source, components, simulation_parameters)
    EnergySystems.control(heat_pump, components, simulation_parameters)
    EnergySystems.control(grid, components, simulation_parameters)
    EnergySystems.process(demand, simulation_parameters)
    EnergySystems.process(heat_pump, simulation_parameters)
    EnergySystems.process(source, simulation_parameters)
    EnergySystems.process(grid, simulation_parameters)

    heat_out = heat_pump.output_interfaces[heat_pump.m_heat_out].sum_abs_change * 0.5
    el_in = heat_pump.input_interfaces[heat_pump.m_el_in].sum_abs_change * 0.5
    heat_in = heat_pump.input_interfaces[heat_pump.m_heat_in].sum_abs_change * 0.5

    @test heat_out ≈ 900
    @test el_in + heat_in ≈ heat_out
    # COP without icing losses would be 3.1815
    @test heat_out / el_in ≈ 2.7993295246351666
end

@testset "test_heat_pump_1S1D_icing_losses" begin
    test_heat_pump_1S1D_icing_losses()
end

function get_config_heat_pump_1S1D_infinities(; inf_as_source::Bool=true)
    base_dict = Dict{String,Any}(
        "TST_GRI_01" => Dict{String,Any}(
            "type" => "GridConnection",
            "medium" => "m_e_ac_230v",
            "output_refs" => ["TST_HP_01"],
            "is_source" => true,
        ),
        "TST_HP_01" => Dict{String,Any}(
            "type" => "HeatPump",
            "output_refs" => ["TST_DEM_01"],
            "power_th" => 2500 * 4,
        ),
    )

    if inf_as_source
        base_dict["TST_HP_01"]["input_temperature"] = 30

        extended = Dict{String,Any}(
            "TST_DEM_01" => Dict{String,Any}(
                "type" => "Demand",
                "medium" => "m_h_w_ht1",
                "output_refs" => [],
                "constant_demand" => 2000 * 4,
                "constant_temperature" => 60,
            ),
            "TST_SRC_01" => Dict{String,Any}(
                "type" => "GridConnection",
                "medium" => "m_h_w_lt1",
                "output_refs" => ["TST_HP_01"],
                "is_source" => true,
            ),
        )
    else
        base_dict["TST_HP_01"]["output_temperature"] = 60

        extended = Dict{String,Any}(
            "TST_DEM_01" => Dict{String,Any}(
                "type" => "GridConnection",
                "medium" => "m_h_w_lt1",
                "output_refs" => ["TST_HP_01"],
                "is_source" => false,
            ),
            "TST_SRC_01" => Dict{String,Any}(
                "type" => "FixedSupply",
                "medium" => "m_h_w_lt1",
                "output_refs" => ["TST_HP_01"],
                "constant_supply" => 2000 * 4,
                "constant_temperature" => 30,
            ),
        )
    end

    return Base.merge(base_dict, extended)
end

function test_heat_pump_1S1D_infinite_input()
    components_config = get_config_heat_pump_1S1D_infinities(; inf_as_source=true)
    simulation_parameters = get_default_sim_params()

    components = Resie.load_components(components_config, simulation_parameters)
    heat_pump = components["TST_HP_01"]
    source = components["TST_SRC_01"]
    demand = components["TST_DEM_01"]
    grid = components["TST_GRI_01"]

    for unit in values(components)
        EnergySystems.reset(unit)
    end

    EnergySystems.control(demand, components, simulation_parameters)
    EnergySystems.control(source, components, simulation_parameters)
    EnergySystems.control(heat_pump, components, simulation_parameters)
    EnergySystems.control(grid, components, simulation_parameters)

    EnergySystems.process(demand, simulation_parameters)
    EnergySystems.process(heat_pump, simulation_parameters)
    EnergySystems.process(source, simulation_parameters)
    EnergySystems.process(grid, simulation_parameters)

    el_in = heat_pump.input_interfaces[heat_pump.m_el_in].sum_abs_change * 0.5
    heat_in = heat_pump.input_interfaces[heat_pump.m_heat_in].sum_abs_change * 0.5
    heat_out = heat_pump.output_interfaces[heat_pump.m_heat_out].sum_abs_change * 0.5
    cop = heat_out / el_in
    @test el_in ≈ 450.24763619991
    @test heat_in ≈ 1549.7523638
    @test heat_out ≈ 2000
    @test cop ≈ 4.442
end

@testset "test_heat_pump_1S1D_infinite_output" begin
    test_heat_pump_1S1D_infinite_input()
end

function test_heat_pump_1S1D_infinite_output()
    components_config = get_config_heat_pump_1S1D_infinities(; inf_as_source=false)
    simulation_parameters = get_default_sim_params()

    components = Resie.load_components(components_config, simulation_parameters)
    heat_pump = components["TST_HP_01"]
    source = components["TST_SRC_01"]
    demand = components["TST_DEM_01"]
    grid = components["TST_GRI_01"]

    for unit in values(components)
        EnergySystems.reset(unit)
    end

    EnergySystems.control(demand, components, simulation_parameters)
    EnergySystems.control(source, components, simulation_parameters)
    EnergySystems.control(heat_pump, components, simulation_parameters)
    EnergySystems.control(grid, components, simulation_parameters)

    EnergySystems.process(demand, simulation_parameters)
    EnergySystems.process(heat_pump, simulation_parameters)
    EnergySystems.process(source, simulation_parameters)
    EnergySystems.process(grid, simulation_parameters)

    el_in = heat_pump.input_interfaces[heat_pump.m_el_in].sum_abs_change * 0.5
    heat_in = heat_pump.input_interfaces[heat_pump.m_heat_in].sum_abs_change * 0.5
    heat_out = heat_pump.output_interfaces[heat_pump.m_heat_out].sum_abs_change * 0.5
    cop = heat_out / el_in
    @test el_in ≈ 581.0575246949
    @test heat_in ≈ 2000
    @test heat_out ≈ 2581.0575246949
    @test cop ≈ 4.442
end

@testset "test_heat_pump_1S1D_infinite_input" begin
    test_heat_pump_1S1D_infinite_input()
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
    simulation_parameters = get_default_sim_params()

    components = Resie.load_components(components_config, simulation_parameters)
    heat_pump = components["TST_HP_01"]
    source_1 = components["TST_SRC_01"]
    source_2 = components["TST_SRC_02"]
    demand_1 = components["TST_DEM_01"]
    demand_2 = components["TST_DEM_02"]
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
    simulation_parameters = get_default_sim_params()

    components = Resie.load_components(components_config, simulation_parameters)
    heat_pump = components["TST_HP_01"]
    source_1 = components["TST_SRC_01"]
    source_2 = components["TST_SRC_02"]
    demand_1 = components["TST_DEM_01"]
    demand_2 = components["TST_DEM_02"]
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
    @test heat_pump.input_interfaces[heat_pump.m_heat_in].sum_abs_change ≈ 4451.65614361108
    @test heat_pump.input_interfaces[heat_pump.m_el_in].sum_abs_change ≈ 1548.3438563889197
    @test heat_pump.output_interfaces[heat_pump.m_heat_out].sum_abs_change ≈ 3000.0 * 2
    @test heat_pump.output_interfaces[heat_pump.m_heat_out].sum_abs_change /
          heat_pump.input_interfaces[heat_pump.m_el_in].sum_abs_change ≈ 3.8751082165904203

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

    @test heat_pump.input_interfaces[heat_pump.m_heat_in].sum_abs_change ≈ 4666.002062659695
    @test heat_pump.input_interfaces[heat_pump.m_el_in].sum_abs_change ≈ 1333.9979373403046
    @test heat_pump.output_interfaces[heat_pump.m_heat_out].sum_abs_change ≈ 3000.0 * 2
    @test heat_pump.output_interfaces[heat_pump.m_heat_out].sum_abs_change /
          heat_pump.input_interfaces[heat_pump.m_el_in].sum_abs_change ≈ 4.497758078968747
end

@testset "heat_pump_2S2D_dynamic_cop" begin
    test_heat_pump_2S2D_dynamic_cop()
end

function test_heat_pump_2S2D_reorder_inputs()
    components_config = get_config_heat_pump_2S2D()

    components_config["TST_HP_01"]["control_modules"] = [Dict{String,Any}(
                                                             "name" => "temperature_sorting",
                                                         )]
    simulation_parameters = get_default_sim_params()

    components = Resie.load_components(components_config, simulation_parameters)
    heat_pump = components["TST_HP_01"]
    source_1 = components["TST_SRC_01"]
    source_2 = components["TST_SRC_02"]
    demand_1 = components["TST_DEM_01"]
    demand_2 = components["TST_DEM_02"]
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

function test_heat_pump_2S2D_min_power()
    components_config = get_config_heat_pump_2S2D()

    components_config["TST_HP_01"]["power_th"] = 28000
    components_config["TST_HP_01"]["max_power_function"] = "const:1.0"
    components_config["TST_HP_01"]["min_power_function"] = "const:0.35"

    simulation_parameters = get_default_sim_params()

    components = Resie.load_components(components_config, simulation_parameters)
    heat_pump = components["TST_HP_01"]
    source_1 = components["TST_SRC_01"]
    source_2 = components["TST_SRC_02"]
    demand_1 = components["TST_DEM_01"]
    demand_2 = components["TST_DEM_02"]
    grid = components["TST_GRI_01"]
    bus_1 = components["TST_BUS_01"]
    bus_2 = components["TST_BUS_02"]

    # first time step: slices can be "slowed down" enough to compensate for high power

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

    energy_full_power = EnergySystems.watt_to_wh(heat_pump.design_power_th)
    produced_heat = heat_pump.output_interfaces[heat_pump.m_heat_out].sum_abs_change * 0.5
    @test produced_heat ≈ 3000.0
    @test produced_heat / energy_full_power > 0.35

    # second time step: min power of each slice is so high that they can't be "slowed down"
    # enough to compensate
    heat_pump.min_power_function = (x, y) -> 0.8

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

    energy_full_power = EnergySystems.watt_to_wh(heat_pump.design_power_th)
    produced_heat = heat_pump.output_interfaces[heat_pump.m_heat_out].sum_abs_change * 0.5
    @test produced_heat ≈ 0.0
    @test produced_heat / energy_full_power ≈ 0.0
end

@testset "heat_pump_2S2D_min_power" begin
    test_heat_pump_2S2D_min_power()
end

function test_heat_pump_2S2D_optimising_slices()
    components_config = get_config_heat_pump_2S2D()

    components_config["TST_HP_01"]["power_th"] = 28000
    components_config["TST_HP_01"]["max_power_function"] = "const:1.0"
    components_config["TST_HP_01"]["min_power_function"] = "const:0.2"
    components_config["TST_HP_01"]["min_power_fraction"] = 0.2
    components_config["TST_HP_01"]["cop_function"] = "carnot:0.4:poly:-1.93407,1.53407,0.9"
    components_config["TST_HP_01"]["optimise_slice_dispatch"] = true
    components_config["TST_HP_01"]["optimal_plr"] = 0.4

    simulation_parameters = get_default_sim_params()
    eps = simulation_parameters["epsilon"]

    components = Resie.load_components(components_config, simulation_parameters)
    heat_pump = components["TST_HP_01"]
    source_1 = components["TST_SRC_01"]
    source_2 = components["TST_SRC_02"]
    demand_1 = components["TST_DEM_01"]
    demand_2 = components["TST_DEM_02"]
    grid = components["TST_GRI_01"]
    bus_1 = components["TST_BUS_01"]
    bus_2 = components["TST_BUS_02"]

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

    # without optimisation, at full power each slice, energies would be:
    # 1403.5371253472942
    # 1596.4628746527055
    # 3000.0

    @test heat_pump.input_interfaces[heat_pump.m_el_in].sum_abs_change * 0.5 < 1403.5371253472942 - eps
    @test heat_pump.input_interfaces[heat_pump.m_heat_in].sum_abs_change * 0.5 > 1596.4628746527055 + eps
    @test heat_pump.output_interfaces[heat_pump.m_heat_out].sum_abs_change * 0.5 ≈ 3000.0
end

@testset "heat_pump_2S2D_optimising_slices" begin
    test_heat_pump_2S2D_optimising_slices()
end