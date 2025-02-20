using Debugger
using Test
using Resie
using Resie.EnergySystems
using Resie.Profiles

include("../test_util.jl")

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
            "power_losses_factor" => 1.0,
            "heat_losses_factor" => 1.0,
            "min_power_function" => "const:0.0",
        ),
    )
end

function test_heat_pump_one_source_dynamic_cop()
    components_config = get_config_heat_pump_1S1D()
    simulation_parameters = get_default_sim_params()

    components = Resie.load_components(components_config, simulation_parameters)
    setup_mock_run!(components, simulation_parameters)

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
    setup_mock_run!(components, simulation_parameters)

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
            "power_losses_factor" => 1.0,
            "heat_losses_factor" => 1.0,
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
    setup_mock_run!(components, simulation_parameters)

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

function test_heat_pump_1S1D_losses()
    components_config = get_config_heat_pump_1S1D()
    simulation_parameters = get_default_sim_params()
    components_config["TST_HP_01"]["cop_function"] = "const:2.0"
    delete!(components_config["TST_SRC_01"], "max_power_profile_file_path")
    components_config["TST_SRC_01"]["constant_power"] = 400

    components = Resie.load_components(components_config, simulation_parameters)
    setup_mock_run!(components, simulation_parameters)

    heat_pump = components["TST_HP_01"]
    source = components["TST_SRC_01"]
    demand = components["TST_DEM_01"]
    grid = components["TST_GRI_01"]

    heat_pump.power_losses_factor = 0.97
    heat_pump.heat_losses_factor = 0.95

    for unit in values(components)
        EnergySystems.reset(unit)
    end

    # first timestep, demand is well below max power of source
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

    @test el_in ≈ 56.25 / 0.97
    @test heat_in ≈ 56.25 / 0.95
    @test heat_out ≈ 112.5
    @test heat_pump.cop ≈ 2.0
    @test heat_pump.effective_cop ≈ 112.5 / (56.25 / 0.97)
    @test heat_pump.losses_heat ≈ (1.0 / 0.95 - 1.0) * 56.25
    @test heat_pump.losses_power ≈ (1.0 / 0.97 - 1.0) * 56.25
    @test heat_pump.losses ≈ (1.0 / 0.97 - 1.0) * 56.25 + (1.0 / 0.95 - 1.0) * 56.25

    # second timestep, max power of source is just above what heat input would be without
    # considering losses, but just below the required heat input including losses
    EnergySystems.reset(heat_pump)
    @test heat_pump.losses == 0.0
    @test heat_pump.losses_heat == 0.0
    @test heat_pump.losses_power == 0.0
    EnergySystems.reset(demand)
    EnergySystems.reset(source)
    EnergySystems.reset(grid)

    source.constant_power = 57.75 * 4
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
    heat_out_balance = heat_pump.output_interfaces[heat_pump.m_heat_out].balance
    heat_out = heat_pump.output_interfaces[heat_pump.m_heat_out].sum_abs_change - 112.5

    # source can only produce 57.75 Wh, which is 54.8625 after losses
    @test el_in ≈ 54.8625 / 0.97
    @test heat_in ≈ 54.8625 / 0.95
    @test heat_out ≈ 54.8625 * 2.0
    @test heat_out_balance ≈ -(112.5 - 54.8625 * 2.0)
    @test heat_pump.cop ≈ 2.0
    @test heat_pump.effective_cop ≈ (54.8625 * 2.0) / (54.8625 / 0.97)
    @test heat_pump.losses_heat ≈ (1.0 / 0.95 - 1.0) * 54.8625
    @test heat_pump.losses_power ≈ (1.0 / 0.97 - 1.0) * 54.8625
    @test heat_pump.losses ≈ (1.0 / 0.95 - 1.0) * 54.8625 + (1.0 / 0.97 - 1.0) * 54.8625
end

@testset "test_heat_pump_1S1D_losses" begin
    test_heat_pump_1S1D_losses()
end

function test_heat_pump_1S1D_infinite_output()
    components_config = get_config_heat_pump_1S1D_infinities(; inf_as_source=false)
    simulation_parameters = get_default_sim_params()

    components = Resie.load_components(components_config, simulation_parameters)
    setup_mock_run!(components, simulation_parameters)

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
            "power_losses_factor" => 1.0,
            "heat_losses_factor" => 1.0,
        ),
    )
end

function test_heat_pump_2S2D_constant_cop()
    components_config = get_config_heat_pump_2S2D()
    simulation_parameters = get_default_sim_params()

    components = Resie.load_components(components_config, simulation_parameters)
    setup_mock_run!(components, simulation_parameters)

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

function test_heat_pump_2S2D_losses()
    components_config = get_config_heat_pump_2S2D()
    simulation_parameters = get_default_sim_params()
    components_config["TST_HP_01"]["power_losses_factor"] = 0.97
    components_config["TST_HP_01"]["heat_losses_factor"] = 0.95
    components_config["TST_HP_01"]["power_th"] = 16000

    components = Resie.load_components(components_config, simulation_parameters)
    setup_mock_run!(components, simulation_parameters)

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

    el_in = heat_pump.input_interfaces[heat_pump.m_el_in].sum_abs_change * 0.5
    heat_in = heat_pump.input_interfaces[heat_pump.m_heat_in].sum_abs_change * 0.5
    heat_out = heat_pump.output_interfaces[heat_pump.m_heat_out].sum_abs_change * 0.5

    # without losses:
    # el_in = 774,1719281944598
    # heat_in = 2225,82807180554
    # heat_out = 3000
    # cop = effective_cop = 3,87510821659042
    #
    # with losses:
    # el_in 807,7268744618796 (+4,334301599604284 %)
    # heat_in 2333,163086075765 (+4,8222509020275 %)
    # heat_out = 3000
    # cop = 3,8289966608020816 (-1,189942401890145 %)
    # effective_cop = 3,714126760978019 (-4,15424412983344442 %)
    # losses_power (actual) = 24,231806233856332
    # losses_heat (actual) = 116,6581543037887
    # losses_power (as diff) = 33,55494626741972
    # losses_heat (as diff) = 107,3350142702252
    #
    # what we suspect happens here is that, because the three slices are different in
    # efficiency, including losses shifts the actually processed energy slightly towards the
    # other slices. namely, it requires a bit less power and utilises more lower temperature
    # heat. the COP however is calculated based on the unchanged heat output, which is why
    # it is lower when losses are included. the inner COP with losses is also smaller than
    # the effective COP without losses, but larger than the effective COP with losses.
    @test el_in ≈ 807.7268744618796
    @test heat_in ≈ 2333.163086075765
    @test heat_out ≈ 3000
    @test heat_pump.cop > 3.714126760978019 && heat_pump.cop < 3.87510821659042
    @test heat_pump.effective_cop ≈ 3.714126760978019
    @test heat_pump.losses_heat > 2333.163086075765 - 2225.82807180554
    @test heat_pump.losses_power < 807.7268744618796 - 774.1719281944598
    @test heat_pump.time_active ≈ 0.75
    @test heat_pump.avg_plr ≈ 1.0
end

@testset "heat_pump_2S2D_losses" begin
    test_heat_pump_2S2D_losses()
end

function test_heat_pump_2S2D_dynamic_cop()
    components_config = get_config_heat_pump_2S2D()
    simulation_parameters = get_default_sim_params()

    components = Resie.load_components(components_config, simulation_parameters)
    setup_mock_run!(components, simulation_parameters)

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
    setup_mock_run!(components, simulation_parameters)

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
    setup_mock_run!(components, simulation_parameters)

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

    energy_full_power = simulation_parameters["watt_to_wh"](heat_pump.design_power_th)
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

    energy_full_power = simulation_parameters["watt_to_wh"](heat_pump.design_power_th)
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
    components_config["TST_HP_01"]["plf_function"] = "poly:-1.93407,1.53407,0.9"
    components_config["TST_HP_01"]["cop_function"] = "carnot:0.4"
    components_config["TST_HP_01"]["optimise_slice_dispatch"] = true
    components_config["TST_HP_01"]["optimal_plr"] = 0.4

    simulation_parameters = get_default_sim_params()
    eps = simulation_parameters["epsilon"]

    components = Resie.load_components(components_config, simulation_parameters)
    setup_mock_run!(components, simulation_parameters)

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
    # fully optimal would be 1.0 for time_active and 0.4 for plr
    @test heat_pump.time_active ≈ 0.9222481721394582
    @test heat_pump.avg_plr ≈ 0.4653057040302939
end

@testset "heat_pump_2S2D_optimising_slices" begin
    test_heat_pump_2S2D_optimising_slices()
end

function test_heat_pump_2S2D_auto_calculate_optimal_plr()
    components_config = get_config_heat_pump_2S2D()

    components_config["TST_HP_01"]["power_th"] = 28000
    components_config["TST_HP_01"]["max_power_function"] = "const:1.0"
    components_config["TST_HP_01"]["min_power_function"] = "const:0.2"
    components_config["TST_HP_01"]["min_power_fraction"] = 0.2
    components_config["TST_HP_01"]["plf_function"] = "poly:-1.93407,1.53407,0.9"
    components_config["TST_HP_01"]["cop_function"] = "carnot:0.4"
    components_config["TST_HP_01"]["optimise_slice_dispatch"] = true

    simulation_parameters = get_default_sim_params()
    components = Resie.load_components(components_config, simulation_parameters)
    heat_pump = components["TST_HP_01"]
    @test heat_pump.optimal_plr > 0.39659112648456 - 0.001
    @test heat_pump.optimal_plr < 0.39659112648456 + 0.001
end

@testset "heat_pump_2S2D_auto_calculate_optimal_plr" begin
    test_heat_pump_2S2D_auto_calculate_optimal_plr()
end
