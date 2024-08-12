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

function test_heat_pump_demand_one_source_dynamic_cop()
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

@testset "heat_pump_demand_one_source_dynamic_cop" begin
    test_heat_pump_demand_one_source_dynamic_cop()
end
