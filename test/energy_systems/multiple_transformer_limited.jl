using Debugger
using Test
using Resie
using Resie.EnergySystems
using Resie.Profiles

include("../test_util.jl")

function test_multiple_transformer_with_limitations()
    components_config = Dict{String,Any}(
        "TST_DEM_heat_01" => Dict{String,Any}(
            "type" => "Demand",
            "medium" => "m_h_w_ht1",
            "output_refs" => [],
            "energy_profile_file_path" => "./profiles/tests/demand_heating_energy.prf",
            "temperature_profile_file_path" => "./profiles/tests/demand_heating_temperature.prf",
            "scale" => 1500,
        ),
        "TST_DEM_H2_01" => Dict{String,Any}(
            "type" => "Demand",
            "medium" => "m_c_g_h2",
            "output_refs" => [],
            "energy_profile_file_path" => "./profiles/tests/demand_h2.prf",
            "scale" => 3000,
        ),
        "TST_GRI_O2_01" => Dict{String,Any}(
            "type" => "GridConnection",
            "medium" => "m_c_g_o2",
            "output_refs" => [],
            "is_source" => false,
        ),
        "TST_GRI_el_01" => Dict{String,Any}(
            "type" => "GridConnection",
            "medium" => "m_e_ac_230v",
            "output_refs" => ["TST_ELY_01"],
            "is_source" => true,
        ),
        "TST_GRI_el_02" => Dict{String,Any}(
            "type" => "GridConnection",
            "medium" => "m_e_ac_230v",
            "output_refs" => ["TST_HP_01"],
            "is_source" => true,
        ),
        "TST_ELY_01" => Dict{String,Any}(
            "type" => "Electrolyser",
            "output_refs" => ["TST_HP_01", "TST_DEM_H2_01", "TST_GRI_O2_01"],
            "power_el" => 4000,
            "output_temperature_ht" => 55.0,
            "m_heat_ht_out" => "m_h_w_lt1",
            "m_heat_lt_out" => "m_h_w_lt2",
            "nr_switchable_units" => 1,
            "dispatch_strategy" => "all_equal",
            "heat_lt_is_usable" => false,
            "min_power_fraction" => 0.0,
            "min_power_fraction_total" => 0.0,
            "linear_interface" => "el_in",
            "efficiency_el_in" => "const:1.0",
            "efficiency_h2_out" => "const:0.6",
            "efficiency_h2_out_lossless" => "const:0.6",
            "efficiency_o2_out" => "const:0.6",
            "efficiency_heat_ht_out" => "const:0.4",
        ),
        "TST_HP_01" => Dict{String,Any}(
            "type" => "HeatPump",
            "output_refs" => ["TST_DEM_heat_01"],
            "power_th" => 2240,
            "cop_function" => "const:3.5",
            "min_power_function" => "const:0.0",
        ),
    )

    simulation_parameters = get_default_sim_params()

    components = Resie.load_components(components_config, simulation_parameters)
    heat_pump = components["TST_HP_01"]
    electrolyser = components["TST_ELY_01"]
    grid_el1 = components["TST_GRI_el_01"]
    grid_el2 = components["TST_GRI_el_02"]
    demand_h2 = components["TST_DEM_H2_01"]
    demand_heat = components["TST_DEM_heat_01"]
    grid_o2 = components["TST_GRI_O2_01"]

    for unit in values(components)
        EnergySystems.reset(unit)
    end

    # first time step: demand is met perfectly
    demand_h2.constant_demand = 2400
    demand_heat.constant_demand = 2240
    EnergySystems.control(demand_heat, components, simulation_parameters)
    EnergySystems.control(demand_h2, components, simulation_parameters)
    EnergySystems.control(heat_pump, components, simulation_parameters)
    EnergySystems.control(electrolyser, components, simulation_parameters)
    EnergySystems.control(grid_el1, components, simulation_parameters)
    EnergySystems.control(grid_el2, components, simulation_parameters)
    EnergySystems.control(grid_o2, components, simulation_parameters)

    EnergySystems.process(demand_heat, simulation_parameters)
    EnergySystems.process(demand_h2, simulation_parameters)
    @test demand_heat.input_interfaces[demand_heat.medium].balance ≈ -2240 / 4
    @test demand_h2.input_interfaces[demand_h2.medium].balance ≈ -2400 / 4

    EnergySystems.potential(heat_pump, simulation_parameters)
    @test EnergySystems.get_max_energy(heat_pump.input_interfaces[heat_pump.m_heat_in].max_energy) ≈
          2240 / 4 * ((3.5 - 1) / 3.5)
    @test EnergySystems.get_max_energy(heat_pump.input_interfaces[heat_pump.m_el_in].max_energy) ≈ 2240 / 4 / 3.5
    @test EnergySystems.get_max_energy(heat_pump.output_interfaces[heat_pump.m_heat_out].max_energy) ≈ 2240 / 4

    exchanges = EnergySystems.balance_on(heat_pump.input_interfaces[heat_pump.m_heat_in], heat_pump)
    @test EnergySystems.balance(exchanges) ≈ 0.0
    @test EnergySystems.energy_potential(exchanges) ≈ -2400 / 4 / 0.6 * 0.4

    EnergySystems.potential(electrolyser, simulation_parameters)
    @test EnergySystems.get_max_energy(electrolyser.output_interfaces[electrolyser.m_heat_ht_out].max_energy) ≈
          2400 / 4 / 0.6 * 0.4
    @test EnergySystems.get_max_energy(electrolyser.output_interfaces[electrolyser.m_h2_out].max_energy) ≈ 2400 / 4
    @test EnergySystems.get_max_energy(electrolyser.output_interfaces[electrolyser.m_o2_out].max_energy) ≈ 2400 / 4
    @test EnergySystems.get_max_energy(electrolyser.input_interfaces[electrolyser.m_el_in].max_energy) ≈ 2400 / 4 / 0.6

    exchanges = EnergySystems.balance_on(heat_pump.input_interfaces[heat_pump.m_heat_in], heat_pump)
    @test EnergySystems.balance(exchanges) ≈ 0.0
    @test EnergySystems.energy_potential(exchanges) ≈ -2400 / 4 / 0.6 * 0.4

    exchanges = EnergySystems.balance_on(electrolyser.output_interfaces[electrolyser.m_heat_ht_out], electrolyser)
    @test EnergySystems.balance(exchanges) ≈ 0.0
    @test EnergySystems.energy_potential(exchanges) ≈ 2400 / 4 / 0.6 * 0.4

    EnergySystems.process(electrolyser, simulation_parameters)
    @test electrolyser.output_interfaces[electrolyser.m_heat_ht_out].balance ≈ 1600 / 4
    @test electrolyser.output_interfaces[electrolyser.m_heat_ht_out].sum_abs_change ≈ 1600 / 4
    @test electrolyser.output_interfaces[electrolyser.m_h2_out].balance ≈ 0
    @test electrolyser.output_interfaces[electrolyser.m_h2_out].sum_abs_change ≈ 2 * 2400 / 4
    @test electrolyser.output_interfaces[electrolyser.m_o2_out].balance ≈ 2400 / 4
    @test electrolyser.input_interfaces[electrolyser.m_el_in].balance ≈ -4000 / 4

    EnergySystems.process(heat_pump, simulation_parameters)
    @test heat_pump.output_interfaces[heat_pump.m_heat_out].balance ≈ 0
    @test heat_pump.output_interfaces[heat_pump.m_heat_out].sum_abs_change ≈ 2 * 2240 / 4
    @test heat_pump.input_interfaces[heat_pump.m_el_in].balance ≈ -640 / 4
    @test heat_pump.input_interfaces[heat_pump.m_heat_in].balance ≈ 0

    EnergySystems.process(grid_el1, simulation_parameters)
    EnergySystems.process(grid_el2, simulation_parameters)
    EnergySystems.process(grid_o2, simulation_parameters)
    @test grid_o2.input_interfaces[grid_o2.medium].balance ≈ 0
    @test grid_o2.input_interfaces[grid_o2.medium].sum_abs_change ≈ 2 * 2400 / 4
    @test grid_el1.output_interfaces[grid_el1.medium].balance ≈ 0
    @test grid_el1.output_interfaces[grid_el1.medium].sum_abs_change ≈ 2 * 4000 / 4
    @test grid_el2.output_interfaces[grid_el2.medium].balance ≈ 0
    @test grid_el2.output_interfaces[grid_el2.medium].sum_abs_change ≈ 2 * 640 / 4

    # second timestep: demand of h2 is limiting the state of the electrolyser. Therefore, the demand of waste heat in
    # the interface between electrolyser and heat pump can not be met. As a result, the state of the heat pump
    # should also be less that nominal state.
    # 12 test are failing at the moment. Reasons:
    # - implemented stratey of electrolyse is not considering demand of H2
    # - calculation of state of one energy system is only done once. If the required demand are not satisfied, this has
    #   no effect on the actual state of this energy system!
    #  --> electrolyser is running on 100% instead of 50% due to control stragety
    #  --> if that is fixed, the heat pump is running on 100% instead of 50% due to lack of information that electrolyser
    #      can only run 50% due to a limitation in hydrogen demand

    for unit in values(components)
        EnergySystems.reset(unit)
    end

    demand_h2.constant_demand = 0.5 * 2400  # reducing h2 demand by half
    demand_heat.constant_demand = 2240  # same heat demand as bevore
    EnergySystems.control(demand_heat, components, simulation_parameters)
    EnergySystems.control(demand_h2, components, simulation_parameters)
    EnergySystems.control(heat_pump, components, simulation_parameters)
    EnergySystems.control(electrolyser, components, simulation_parameters)
    EnergySystems.control(grid_el1, components, simulation_parameters)
    EnergySystems.control(grid_el2, components, simulation_parameters)
    EnergySystems.control(grid_o2, components, simulation_parameters)

    EnergySystems.process(demand_heat, simulation_parameters)
    EnergySystems.process(demand_h2, simulation_parameters)
    @test demand_heat.input_interfaces[demand_heat.medium].balance ≈ -2240 / 4
    @test demand_h2.input_interfaces[demand_h2.medium].balance ≈ -0.5 * 2400 / 4

    EnergySystems.potential(heat_pump, simulation_parameters)
    @test EnergySystems.get_max_energy(heat_pump.input_interfaces[heat_pump.m_heat_in].max_energy) ≈
          2240 / 4 * ((3.5 - 1) / 3.5)
    @test EnergySystems.get_max_energy(heat_pump.input_interfaces[heat_pump.m_el_in].max_energy) ≈ 2240 / 4 / 3.5
    @test EnergySystems.get_max_energy(heat_pump.output_interfaces[heat_pump.m_heat_out].max_energy) ≈ 2240 / 4

    exchanges = EnergySystems.balance_on(heat_pump.input_interfaces[heat_pump.m_heat_in], heat_pump)
    @test EnergySystems.balance(exchanges) ≈ 0.0
    @test EnergySystems.energy_potential(exchanges) ≈ -2400 / 4 / 0.6 * 0.4

    EnergySystems.potential(electrolyser, simulation_parameters)
    @test EnergySystems.get_max_energy(electrolyser.output_interfaces[electrolyser.m_heat_ht_out].max_energy) ≈
          0.5 * 2400 / 4 / 0.6 * 0.4
    @test EnergySystems.get_max_energy(electrolyser.output_interfaces[electrolyser.m_h2_out].max_energy) ≈
          0.5 * 2400 / 4
    @test EnergySystems.get_max_energy(electrolyser.output_interfaces[electrolyser.m_o2_out].max_energy) ≈
          0.5 * 2400 / 4
    @test EnergySystems.get_max_energy(electrolyser.input_interfaces[electrolyser.m_el_in].max_energy) ≈
          0.5 * 2400 / 4 / 0.6

    exchanges = EnergySystems.balance_on(heat_pump.input_interfaces[heat_pump.m_heat_in], heat_pump)
    @test EnergySystems.balance(exchanges) ≈ 0.0
    @test EnergySystems.energy_potential(exchanges) ≈ -0.5 * 2400 / 4 / 0.6 * 0.4

    exchanges = EnergySystems.balance_on(electrolyser.output_interfaces[electrolyser.m_heat_ht_out], electrolyser)
    @test EnergySystems.balance(exchanges) ≈ 0.0
    @test EnergySystems.energy_potential(exchanges) ≈ 0.5 * 2400 / 4 / 0.6 * 0.4

    EnergySystems.process(electrolyser, simulation_parameters)
    EnergySystems.process(heat_pump, simulation_parameters)

    @test heat_pump.output_interfaces[heat_pump.m_heat_out].balance ≈ -0.5 * 2240 / 4
    @test heat_pump.output_interfaces[heat_pump.m_heat_out].sum_abs_change ≈ 560 + 560 / 2
    @test heat_pump.input_interfaces[heat_pump.m_el_in].balance ≈ -0.5 * 640 / 4
    @test heat_pump.input_interfaces[heat_pump.m_heat_in].balance ≈ 0  # process of electrolyser already done
    @test heat_pump.input_interfaces[heat_pump.m_heat_in].sum_abs_change ≈ 0.5 * 2 * 1600 / 4 # 200 are transferred between ely and hp

    @test electrolyser.output_interfaces[electrolyser.m_heat_ht_out].balance ≈ 0 #actually the same as two test above
    @test electrolyser.output_interfaces[electrolyser.m_heat_ht_out].sum_abs_change ≈ 0.5 * 2 * 1600 / 4  #actually the same as two test above

    @test electrolyser.output_interfaces[electrolyser.m_h2_out].balance ≈ 0
    @test electrolyser.output_interfaces[electrolyser.m_h2_out].sum_abs_change ≈ 2 * 0.5 * 2400 / 4
    @test electrolyser.output_interfaces[electrolyser.m_o2_out].balance ≈ 0.5 * 2400 / 4
    @test electrolyser.input_interfaces[electrolyser.m_el_in].balance ≈ -0.5 * 4000 / 4

    EnergySystems.process(grid_el1, simulation_parameters)
    EnergySystems.process(grid_el2, simulation_parameters)
    EnergySystems.process(grid_o2, simulation_parameters)
    @test grid_o2.input_interfaces[grid_o2.medium].balance ≈ 0
    @test grid_o2.input_interfaces[grid_o2.medium].sum_abs_change ≈ 2 * 0.5 * 2400 / 4
    @test grid_el1.output_interfaces[grid_el1.medium].balance ≈ 0
    @test grid_el1.output_interfaces[grid_el1.medium].sum_abs_change ≈ 2 * 0.5 * 4000 / 4
    @test grid_el2.output_interfaces[grid_el2.medium].balance ≈ 0
    @test grid_el2.output_interfaces[grid_el2.medium].sum_abs_change ≈ 2 * 0.5 * 640 / 4
end

@testset "multiple_transformer_with_limitations" begin
    test_multiple_transformer_with_limitations()
end
