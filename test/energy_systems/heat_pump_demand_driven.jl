watt_to_wh = function (watts :: Float64)
    watts * 900 / 3600.0
end

@testset "heat_pump_demand_driven_correct_order" begin
    systems_config = Dict{String, Any}(
        "TST_DEM_01" => Dict{String, Any}(
            "type" => "Demand",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "production_refs" => [],
            "energy_profile_file_path" => "./profiles/tests/demand_heating_energy.prf",
            "temperature_profile_file_path" => "./profiles/tests/demand_heating_temperature.prf",
            "scale" => 1500
        ),
        "TST_SRC_01" => Dict{String, Any}(
            "type" => "DispatchableSupply",
            "medium" => "m_h_w_lt1",
            "control_refs" => [],
            "production_refs" => ["TST_HP_01"],
            "max_power_profile_file_path" => "./profiles/tests/demand_heating_energy.prf",
            "temperature_profile_file_path" => "./profiles/tests/demand_heating_temperature.prf",
            "scale" => 6000
        ),
        "TST_GRI_01" => Dict{String, Any}(
            "type" => "GridConnection",
            "medium" => "m_e_ac_230v",
            "control_refs" => [],
            "production_refs" => ["TST_HP_01"],
            "is_source" => true,
        ),
        "TST_HP_01" => Dict{String, Any}(
            "type" => "HeatPump",
            "control_refs" => ["TST_DEM_01"],
            "production_refs" => ["TST_DEM_01"],
            "strategy" => Dict{String, Any}(
                "name" => "demand_driven",
            ),
            "power" => 12000,
            "cop" => 3.0
        ),
    )
    systems = Bran.load_systems(systems_config)
    heat_pump = systems["TST_HP_01"]
    source = systems["TST_SRC_01"]
    demand = systems["TST_DEM_01"]
    grid = systems["TST_GRI_01"]

    simulation_parameters = Dict{String, Any}(
        "time_step_seconds" => 900,
        "time" => 0,
    )

    @test heat_pump.controller.state_machine.state == 1

    # first time step: demand is below max power of source (adjusted for additional input
    # of electricity), small delta T leads to high COP = 4.594631

    for unit in values(systems)
        EnergySystems.reset(unit)
    end

    EnergySystems.control(demand, systems, simulation_parameters)

    demand.load = 900
    demand.temperature = 45

    EnergySystems.control(source, systems, simulation_parameters)

    source.max_power = 2500
    source.temperature = 35

    EnergySystems.control(heat_pump, systems, simulation_parameters)
    EnergySystems.control(grid, systems, simulation_parameters)

    EnergySystems.produce(demand, simulation_parameters, watt_to_wh)
    @test demand.input_interfaces[EnergySystems.m_h_w_ht1].balance ≈ -900
    @test demand.input_interfaces[EnergySystems.m_h_w_ht1].temperature == 45

    EnergySystems.produce(heat_pump, simulation_parameters, watt_to_wh)
    @test heat_pump.output_interfaces[EnergySystems.m_h_w_ht1].balance ≈ 0
    @test heat_pump.output_interfaces[EnergySystems.m_h_w_ht1].sum_abs_change ≈ 1800
    @test heat_pump.output_interfaces[EnergySystems.m_h_w_ht1].temperature == 45
    @test heat_pump.input_interfaces[EnergySystems.m_e_ac_230v].balance ≈ -300
    @test heat_pump.input_interfaces[EnergySystems.m_e_ac_230v].temperature === nothing
    @test heat_pump.input_interfaces[EnergySystems.m_h_w_lt1].balance ≈ -600
    @test heat_pump.input_interfaces[EnergySystems.m_h_w_lt1].temperature === nothing

    EnergySystems.produce(source, simulation_parameters, watt_to_wh)
    @test source.output_interfaces[EnergySystems.m_h_w_lt1].balance ≈ 0
    @test source.output_interfaces[EnergySystems.m_h_w_lt1].sum_abs_change ≈ 1200
    @test source.output_interfaces[EnergySystems.m_h_w_lt1].temperature == 35

    EnergySystems.produce(grid, simulation_parameters, watt_to_wh)
    @test grid.output_interfaces[EnergySystems.m_e_ac_230v].balance ≈ 0
    @test grid.output_interfaces[EnergySystems.m_e_ac_230v].sum_abs_change ≈ 600
    @test grid.output_interfaces[EnergySystems.m_e_ac_230v].temperature === nothing

    # second step: demand is above max power of source (adjusted for additional input
    # of electricity), big delta T leads to low COP = 1.326097

    for unit in values(systems)
        EnergySystems.reset(unit)
    end

    EnergySystems.control(demand, systems, simulation_parameters)

    demand.load = 2100
    demand.temperature = 75

    EnergySystems.control(source, systems, simulation_parameters)

    source.max_power = 2000
    source.temperature = 35

    EnergySystems.control(heat_pump, systems, simulation_parameters)
    EnergySystems.control(grid, systems, simulation_parameters)

    EnergySystems.produce(demand, simulation_parameters, watt_to_wh)
    @test demand.input_interfaces[EnergySystems.m_h_w_ht1].balance ≈ -2100
    @test demand.input_interfaces[EnergySystems.m_h_w_ht1].temperature == 75

    EnergySystems.produce(heat_pump, simulation_parameters, watt_to_wh)
    @test heat_pump.output_interfaces[EnergySystems.m_h_w_ht1].balance ≈ 0
    @test heat_pump.output_interfaces[EnergySystems.m_h_w_ht1].sum_abs_change ≈ 4200
    @test heat_pump.output_interfaces[EnergySystems.m_h_w_ht1].temperature == 75
    @test heat_pump.input_interfaces[EnergySystems.m_e_ac_230v].balance ≈ -700
    @test heat_pump.input_interfaces[EnergySystems.m_e_ac_230v].temperature === nothing
    @test heat_pump.input_interfaces[EnergySystems.m_h_w_lt1].balance ≈ -1400
    @test heat_pump.input_interfaces[EnergySystems.m_h_w_lt1].temperature === nothing

    EnergySystems.produce(source, simulation_parameters, watt_to_wh)
    @test source.output_interfaces[EnergySystems.m_h_w_lt1].balance ≈ -900
    @test source.output_interfaces[EnergySystems.m_h_w_lt1].sum_abs_change ≈ 1900
    @test source.output_interfaces[EnergySystems.m_h_w_lt1].temperature == 35

    EnergySystems.produce(grid, simulation_parameters, watt_to_wh)
    @test grid.output_interfaces[EnergySystems.m_e_ac_230v].balance ≈ 0
    @test grid.output_interfaces[EnergySystems.m_e_ac_230v].sum_abs_change ≈ 1400
    @test grid.output_interfaces[EnergySystems.m_e_ac_230v].temperature === nothing
end