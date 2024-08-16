using Debugger
using Test
using Resie
using Resie.EnergySystems

function test_distance_from_sink()
    components_config = Dict{String,Any}(
        "TST_BUS_01" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_h_w_ht1",
            "connections" => Dict{String,Any}(
                "input_order" => [],
                "output_order" => ["TST_BUS_02", "TST_BUS_03"],
            ),
        ),
        "TST_BUS_02" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_h_w_ht1",
            "connections" => Dict{String,Any}(
                "input_order" => ["TST_BUS_01"],
                "output_order" => ["TST_BUS_04", "TST_BUS_05"],
            ),
        ),
        "TST_BUS_03" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_h_w_ht1",
            "connections" => Dict{String,Any}(
                "input_order" => ["TST_BUS_01"],
                "output_order" => [],
            ),
        ),
        "TST_BUS_04" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_h_w_ht1",
            "connections" => Dict{String,Any}(
                "input_order" => ["TST_BUS_02"],
                "output_order" => [],
            ),
        ),
        "TST_BUS_05" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_h_w_ht1",
            "connections" => Dict{String,Any}(
                "input_order" => ["TST_BUS_02"],
                "output_order" => ["TST_BUS_06"],
            ),
        ),
        "TST_BUS_06" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_h_w_ht1",
            "connections" => Dict{String,Any}(
                "input_order" => ["TST_BUS_05"],
                "output_order" => [],
            ),
        ),
    )

    simulation_parameters = Dict{String,Any}(
        "time_step_seconds" => 900,
        "time" => 0,
        "epsilon" => 1e-9,
        "is_first_timestep" => true,
    )

    components = Resie.load_components(components_config, simulation_parameters)

    @test Resie.distance_to_sink(components["TST_BUS_01"], EnergySystems.sf_bus, [], "") == 3
    @test Resie.distance_to_sink(components["TST_BUS_02"], EnergySystems.sf_bus, [], "") == 2
    @test Resie.distance_to_sink(components["TST_BUS_03"], EnergySystems.sf_bus, [], "") == 0
    @test Resie.distance_to_sink(components["TST_BUS_04"], EnergySystems.sf_bus, [], "") == 0
    @test Resie.distance_to_sink(components["TST_BUS_05"], EnergySystems.sf_bus, [], "") == 1
    @test Resie.distance_to_sink(components["TST_BUS_06"], EnergySystems.sf_bus, [], "") == 0
end

@testset "distance_from_sink" begin
    test_distance_from_sink()
end

function test_iteration_order()
    components_config = Dict{String,Any}(
        "TST_BUS_01" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_h_w_ht1",
            "connections" => Dict{String,Any}(
                "input_order" => [],
                "output_order" => ["TST_BUS_02", "TST_BUS_03"],
            ),
        ),
        "TST_BUS_02" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_h_w_ht1",
            "connections" => Dict{String,Any}(
                "input_order" => ["TST_BUS_01"],
                "output_order" => ["TST_BUS_04", "TST_BUS_05"],
            ),
        ),
        "TST_BUS_03" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_h_w_ht1",
            "connections" => Dict{String,Any}(
                "input_order" => ["TST_BUS_01"],
                "output_order" => [],
            ),
        ),
        "TST_BUS_04" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_h_w_ht1",
            "connections" => Dict{String,Any}(
                "input_order" => ["TST_BUS_02"],
                "output_order" => [],
            ),
        ),
        "TST_BUS_05" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_h_w_ht1",
            "connections" => Dict{String,Any}(
                "input_order" => ["TST_BUS_02"],
                "output_order" => ["TST_BUS_06"],
            ),
        ),
        "TST_BUS_06" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_h_w_ht1",
            "connections" => Dict{String,Any}(
                "input_order" => ["TST_BUS_06"],
                "output_order" => [],
            ),
        ),
    )

    simulation_parameters = Dict{String,Any}(
        "time_step_seconds" => 900,
        "time" => 0,
        "epsilon" => 1e-9,
        "is_first_timestep" => true,
    )

    components = Resie.load_components(components_config, simulation_parameters)

    chain = [components["TST_BUS_01"],
             components["TST_BUS_02"],
             components["TST_BUS_03"],
             components["TST_BUS_04"],
             components["TST_BUS_05"],
             components["TST_BUS_06"]]
    expected = [components["TST_BUS_03"],
                components["TST_BUS_04"],
                components["TST_BUS_06"],
                components["TST_BUS_05"],
                components["TST_BUS_02"],
                components["TST_BUS_01"]]
    calculated = Resie.iterate_chain(chain, EnergySystems.sf_bus; reverse=false)
    for i in 1:6
        @test "#$i: " * calculated[i].uac == "#$i: " * expected[i].uac
    end
end

@testset "iteration_order" begin
    test_iteration_order()
end

function test_find_chains()
    components_config = Dict{String,Any}(
        "TST_BUS_01" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_h_w_lt1",
            "connections" => Dict{String,Any}(
                "input_order" => [],
                "output_order" => ["TST_BUS_02"],
            ),
        ),
        "TST_BUS_02" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_h_w_lt1",
            "connections" => Dict{String,Any}(
                "input_order" => ["TST_BUS_01"],
                "output_order" => ["TST_HTP_01", "TST_DEM_01"],
            ),
        ),
        "TST_BUS_03" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_h_w_ht1",
            "connections" => Dict{String,Any}(
                "input_order" => ["TST_HTP_01"],
                "output_order" => ["TST_BUS_04", "TST_BUS_05"],
            ),
        ),
        "TST_BUS_04" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_h_w_ht1",
            "connections" => Dict{String,Any}(
                "input_order" => ["TST_BUS_03"],
                "output_order" => [],
            ),
        ),
        "TST_BUS_05" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_h_w_ht1",
            "connections" => Dict{String,Any}(
                "input_order" => ["TST_BUS_03"],
                "output_order" => [],
            ),
        ),
        "TST_GRI_01" => Dict{String,Any}(
            "type" => "GridConnection",
            "medium" => "m_e_ac_230v",
            "output_refs" => ["TST_HTP_01"],
            "is_source" => true,
        ),
        "TST_HTP_01" => Dict{String,Any}(
            "type" => "HeatPump",
            "output_refs" => ["TST_BUS_03"],
            "power_th" => 12000,
            "constant_cop" => 3.0,
        ),
        "TST_DEM_01" => Dict{String,Any}(
            "type" => "Demand",
            "medium" => "m_h_w_lt1",
            "energy_profile_file_path" => "./profiles/tests/demand_heating_energy.prf",
            "temperature_profile_file_path" => "./profiles/tests/demand_heating_temperature.prf",
            "scale" => 1000,
        ),
    )

    simulation_parameters = Dict{String,Any}(
        "time_step_seconds" => 900,
        "time" => 0,
        "epsilon" => 1e-9,
        "is_first_timestep" => true,
    )

    components = Resie.load_components(components_config, simulation_parameters)

    expected_1 = Set([components["TST_BUS_01"],
                      components["TST_BUS_02"]])
    expected_2 = Set([components["TST_BUS_03"],
                      components["TST_BUS_04"],
                      components["TST_BUS_05"]])

    chains = Resie.find_chains([u for u in values(components) if !startswith(u.uac, "Proxy")],
                               EnergySystems.sf_bus)

    @test length(chains) == 2
    @test chains[1] == expected_1
    @test chains[2] == expected_2
end

function test_find_indirect_chains()
    components_config = Dict{String,Any}(
        "TST_DEM_01" => Dict{String,Any}(
            "type" => "Demand",
            "medium" => "m_h_w_ht1",
            "output_refs" => [],
            "energy_profile_file_path" => "./profiles/examples/general/dem_heat_nrg_var_hi-amp.prf",
            "constant_temperature" => 70,
            "scale" => 500,
        ),
        "TST_DEM_02" => Dict{String,Any}(
            "type" => "Demand",
            "medium" => "m_h_w_ht1",
            "output_refs" => [],
            "energy_profile_file_path" => "./profiles/examples/general/dem_heat_nrg_var_hi-amp.prf",
            "constant_temperature" => 90,
            "scale" => 500,
        ),
        "TST_DEM_03" => Dict{String,Any}(
            "type" => "Demand",
            "medium" => "m_h_w_ht1",
            "output_refs" => [],
            "energy_profile_file_path" => "./profiles/examples/general/dem_heat_nrg_var_hi-amp.prf",
            "constant_temperature" => 95,
            "scale" => 500,
        ),
        "TST_GRI_H2" => Dict{String,Any}(
            "type" => "GridConnection",
            "medium" => "m_c_g_h2",
            "input_refs" => ["TST_01_ELY_01"],
            "is_source" => false,
        ),
        "TST_GRI_O2" => Dict{String,Any}(
            "type" => "GridConnection",
            "medium" => "m_c_g_o2",
            "input_refs" => ["TST_01_ELY_01"],
            "is_source" => false,
        ),
        "TST_SRC_01" => Dict{String,Any}(
            "type" => "BoundedSupply",
            "medium" => "m_h_w_lt1",
            "output_refs" => ["TST_BUS_00"],
            "max_power_profile_file_path" => "./profiles/examples/general/src_heat_maxpow_var_lo-amp.prf",
            "temperature_profile_file_path" => "./profiles/examples/general/src_heat_temp_var_avg25.prf",
            "scale" => 1000,
        ),
        "TST_GRI_00" => Dict{String,Any}(
            "type" => "GridConnection",
            "medium" => "m_e_ac_230v",
            "output_refs" => ["TST_01_ELY_01"],
            "is_source" => true,
        ),
        "TST_GRI_01" => Dict{String,Any}(
            "type" => "GridConnection",
            "medium" => "m_e_ac_230v",
            "output_refs" => ["TST_HP_01"],
            "is_source" => true,
        ),
        "TST_SRC_1b" => Dict{String,Any}(
            "type" => "BoundedSupply",
            "medium" => "m_e_ac_230v",
            "output_refs" => ["TST_HP_01b"],
            "is_source" => true,
            "constant_power" => 400,
        ),
        "TST_GRI_02" => Dict{String,Any}(
            "type" => "GridConnection",
            "medium" => "m_e_ac_230v",
            "output_refs" => ["TST_HP_02"],
            "is_source" => true,
        ),
        "TST_GRI_03" => Dict{String,Any}(
            "type" => "GridConnection",
            "medium" => "m_e_ac_230v",
            "output_refs" => ["TST_HP_03"],
            "is_source" => true,
        ),
        "TST_GRI_04" => Dict{String,Any}(
            "type" => "GridConnection",
            "medium" => "m_e_ac_230v",
            "output_refs" => ["TST_HP_04"],
            "is_source" => true,
        ),
        "TST_HP_01" => Dict{String,Any}(
            "type" => "HeatPump",
            "output_refs" => ["TST_HP_01b"],
            "power_th" => 1200,
            "output_temperature" => 60,
            "min_power_fraction" => 0.0,
        ),
        "TST_HP_01b" => Dict{String,Any}(
            "type" => "HeatPump",
            "m_heat_in" => "m_h_w_ht1",
            "output_refs" => ["TST_BUS_01"],
            "power_th" => 9000,
            "min_power_fraction" => 0.0,
        ),
        "TST_HP_02" => Dict{String,Any}(
            "type" => "HeatPump",
            "output_refs" => ["TST_BUS_01"],
            "power_th" => 9000,
            "min_power_fraction" => 0.0,
        ),
        "TST_HP_03" => Dict{String,Any}(
            "type" => "HeatPump",
            "m_heat_in" => "m_h_w_ht1",
            "output_refs" => ["TST_DEM_02"],
            "power_th" => 9000,
            "min_power_fraction" => 0.0,
        ),
        "TST_HP_04" => Dict{String,Any}(
            "type" => "HeatPump",
            "m_heat_in" => "m_h_w_ht1",
            "output_refs" => ["TST_DEM_03"],
            "power_th" => 9000,
            "min_power_fraction" => 0.0,
        ),
        "TST_01_ELY_01" => Dict{String,Any}(
            "type" => "Electrolyser",
            "output_refs" => ["TST_BUS_00",
                              "TST_GRI_H2",
                              "TST_GRI_O2"],
            "m_heat_ht_out" => "m_h_w_lt1",
            "output_temperature" => 45,
            "power_el" => 40000,
            "heat_lt_is_usable" => false,
            "nr_switchable_units" => 4,
            "dispatch_strategy" => "equal_with_mpf",
            "min_power_fraction_total" => 0.0,
            "min_power_fraction" => 0.0,
        ),
        "TST_BUS_00" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_h_w_lt1",
            "connections" => Dict{String,Any}(
                "input_order" => ["TST_01_ELY_01", "TST_SRC_01"],
                "output_order" => ["TST_HP_02", "TST_HP_01"],
                "energy_flow" => [[1, 1],
                                  [1, 0]],
            ),
        ),
        "TST_BUS_01" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_h_w_ht1",
            "connections" => Dict{String,Any}(
                "input_order" => ["TST_HP_01b", "TST_HP_02"],
                "output_order" => ["TST_HP_03", "TST_HP_04", "TST_DEM_01"],
                "energy_flow" => [[1, 1, 1],
                                  [1, 1, 1]],
            ),
        ),
    )

    simulation_parameters = Dict{String,Any}(
        "time_step_seconds" => 900,
        "time" => 0,
        "epsilon" => 1e-9,
        "is_first_timestep" => true,
    )

    components = Resie.load_components(components_config, simulation_parameters)

    expected_bus1 = Set([components["TST_BUS_00"]])
    expected_bus2 = Set([components["TST_BUS_01"]])
    expected_transformer = Set([components["TST_HP_01"],
                                components["TST_HP_01b"],
                                components["TST_HP_02"],
                                components["TST_HP_03"],
                                components["TST_HP_04"],
                                components["TST_01_ELY_01"]])

    chains_bus = Resie.find_chains([u for u in values(components) if !startswith(u.uac, "Proxy")],
                                   EnergySystems.sf_bus)

    chains_transformer = Resie.find_chains([u for u in values(components) if !startswith(u.uac, "Proxy")],
                                           EnergySystems.sf_transformer;
                                           direct_connection_only=false)

    @test length(chains_bus) == 2
    @test chains_bus[1] == expected_bus1
    @test chains_bus[2] == expected_bus2
    @test length(chains_transformer) == 1
    @test chains_transformer[1] == expected_transformer

    # test distance to sink, considers energy flow matrix of bus
    @test Resie.distance_to_sink(components["TST_HP_01"], EnergySystems.sf_transformer, [], "") == 3
    @test Resie.distance_to_sink(components["TST_HP_01b"], EnergySystems.sf_transformer, [], "") == 2
    @test Resie.distance_to_sink(components["TST_HP_02"], EnergySystems.sf_transformer, [], "") == 2
    @test Resie.distance_to_sink(components["TST_HP_03"], EnergySystems.sf_transformer, [], "") == 0
    @test Resie.distance_to_sink(components["TST_HP_04"], EnergySystems.sf_transformer, [], "") == 0
    @test Resie.distance_to_sink(components["TST_01_ELY_01"], EnergySystems.sf_transformer, [], "") == 5

    # test function to determine if transformers are at interfaces, considers energy flow matrix of bus
    @test Resie.check_interface_for_transformer(components["TST_HP_02"].output_interfaces[components["TST_HP_02"].m_heat_out],
                                                "output") == true
    @test Resie.check_interface_for_transformer(components["TST_HP_01b"].output_interfaces[components["TST_HP_01b"].m_heat_out],
                                                "output") == true
    @test Resie.check_interface_for_transformer(components["TST_HP_01b"].input_interfaces[components["TST_HP_01b"].m_heat_in],
                                                "input") == true
    @test Resie.check_interface_for_transformer(components["TST_01_ELY_01"].output_interfaces[components["TST_01_ELY_01"].m_heat_ht_out],
                                                "output") == true
    @test Resie.check_interface_for_transformer(components["TST_HP_03"].output_interfaces[components["TST_HP_03"].m_heat_out],
                                                "output") == false
    @test Resie.check_interface_for_transformer(components["TST_HP_03"].input_interfaces[components["TST_HP_03"].m_heat_in],
                                                "input") == true
    @test Resie.check_interface_for_transformer(components["TST_HP_04"].input_interfaces[components["TST_HP_04"].m_heat_in],
                                                "input") == true
end

function test_find_indirect_chains_denied()
    # This is the same energy system as in test_find_indirect_chains(), but the connection matrix
    # of the two busses denies the connection to and from HP1/HP1b to other transformers. 
    # Therefore, here two independent transformer chains are detected instead of one.
    components_config = Dict{String,Any}(
        "TST_DEM_01" => Dict{String,Any}(
            "type" => "Demand",
            "medium" => "m_h_w_ht1",
            "output_refs" => [],
            "energy_profile_file_path" => "./profiles/examples/general/dem_heat_nrg_var_hi-amp.prf",
            "constant_temperature" => 70,
            "scale" => 500,
        ),
        "TST_DEM_02" => Dict{String,Any}(
            "type" => "Demand",
            "medium" => "m_h_w_ht1",
            "output_refs" => [],
            "energy_profile_file_path" => "./profiles/examples/general/dem_heat_nrg_var_hi-amp.prf",
            "constant_temperature" => 90,
            "scale" => 500,
        ),
        "TST_DEM_03" => Dict{String,Any}(
            "type" => "Demand",
            "medium" => "m_h_w_ht1",
            "output_refs" => [],
            "energy_profile_file_path" => "./profiles/examples/general/dem_heat_nrg_var_hi-amp.prf",
            "constant_temperature" => 95,
            "scale" => 500,
        ),
        "TST_GRI_H2" => Dict{String,Any}(
            "type" => "GridConnection",
            "medium" => "m_c_g_h2",
            "input_refs" => ["TST_01_ELY_01"],
            "is_source" => false,
        ),
        "TST_GRI_O2" => Dict{String,Any}(
            "type" => "GridConnection",
            "medium" => "m_c_g_o2",
            "input_refs" => ["TST_01_ELY_01"],
            "is_source" => false,
        ),
        "TST_SRC_01" => Dict{String,Any}(
            "type" => "BoundedSupply",
            "medium" => "m_h_w_lt1",
            "output_refs" => ["TST_BUS_00"],
            "max_power_profile_file_path" => "./profiles/examples/general/src_heat_maxpow_var_lo-amp.prf",
            "temperature_profile_file_path" => "./profiles/examples/general/src_heat_temp_var_avg25.prf",
            "scale" => 1000,
        ),
        "TST_GRI_00" => Dict{String,Any}(
            "type" => "GridConnection",
            "medium" => "m_e_ac_230v",
            "output_refs" => ["TST_01_ELY_01"],
            "is_source" => true,
        ),
        "TST_GRI_01" => Dict{String,Any}(
            "type" => "GridConnection",
            "medium" => "m_e_ac_230v",
            "output_refs" => ["TST_HP_01"],
            "is_source" => true,
        ),
        "TST_SRC_1b" => Dict{String,Any}(
            "type" => "BoundedSupply",
            "medium" => "m_e_ac_230v",
            "output_refs" => ["TST_HP_01b"],
            "is_source" => true,
            "constant_power" => 400,
        ),
        "TST_GRI_02" => Dict{String,Any}(
            "type" => "GridConnection",
            "medium" => "m_e_ac_230v",
            "output_refs" => ["TST_HP_02"],
            "is_source" => true,
        ),
        "TST_GRI_03" => Dict{String,Any}(
            "type" => "GridConnection",
            "medium" => "m_e_ac_230v",
            "output_refs" => ["TST_HP_03"],
            "is_source" => true,
        ),
        "TST_GRI_04" => Dict{String,Any}(
            "type" => "GridConnection",
            "medium" => "m_e_ac_230v",
            "output_refs" => ["TST_HP_04"],
            "is_source" => true,
        ),
        "TST_HP_01" => Dict{String,Any}(
            "type" => "HeatPump",
            "output_refs" => ["TST_HP_01b"],
            "power_th" => 1200,
            "output_temperature" => 60,
            "min_power_fraction" => 0.0,
        ),
        "TST_HP_01b" => Dict{String,Any}(
            "type" => "HeatPump",
            "m_heat_in" => "m_h_w_ht1",
            "output_refs" => ["TST_BUS_01"],
            "power_th" => 9000,
            "min_power_fraction" => 0.0,
        ),
        "TST_HP_02" => Dict{String,Any}(
            "type" => "HeatPump",
            "output_refs" => ["TST_BUS_01"],
            "power_th" => 9000,
            "min_power_fraction" => 0.0,
        ),
        "TST_HP_03" => Dict{String,Any}(
            "type" => "HeatPump",
            "m_heat_in" => "m_h_w_ht1",
            "output_refs" => ["TST_DEM_02"],
            "power_th" => 9000,
            "min_power_fraction" => 0.0,
        ),
        "TST_HP_04" => Dict{String,Any}(
            "type" => "HeatPump",
            "m_heat_in" => "m_h_w_ht1",
            "output_refs" => ["TST_DEM_03"],
            "power_th" => 9000,
            "min_power_fraction" => 0.0,
        ),
        "TST_01_ELY_01" => Dict{String,Any}(
            "type" => "Electrolyser",
            "output_refs" => ["TST_BUS_00",
                              "TST_GRI_H2",
                              "TST_GRI_O2"],
            "m_heat_ht_out" => "m_h_w_lt1",
            "output_temperature" => 45,
            "power_el" => 40000,
            "heat_lt_is_usable" => false,
            "nr_switchable_units" => 4,
            "dispatch_strategy" => "equal_with_mpf",
            "min_power_fraction_total" => 0.0,
            "min_power_fraction" => 0.0,
        ),
        "TST_BUS_00" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_h_w_lt1",
            "connections" => Dict{String,Any}(
                "input_order" => ["TST_01_ELY_01", "TST_SRC_01"],
                "output_order" => ["TST_HP_02", "TST_HP_01"],
                "energy_flow" => [[1, 0],
                                  [1, 1]],
            ),
        ),
        "TST_BUS_01" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_h_w_ht1",
            "connections" => Dict{String,Any}(
                "input_order" => ["TST_HP_01b", "TST_HP_02"],
                "output_order" => ["TST_HP_03", "TST_HP_04", "TST_DEM_01"],
                "energy_flow" => [[0, 0, 1],
                                  [1, 1, 1]],
            ),
        ),
    )

    simulation_parameters = Dict{String,Any}(
        "time_step_seconds" => 900,
        "time" => 0,
        "epsilon" => 1e-9,
        "is_first_timestep" => true,
    )

    components = Resie.load_components(components_config, simulation_parameters)

    expected_bus1 = Set([components["TST_BUS_00"]])
    expected_bus2 = Set([components["TST_BUS_01"]])
    expected_transformer1 = Set([components["TST_HP_02"],
                                 components["TST_HP_03"],
                                 components["TST_HP_04"],
                                 components["TST_01_ELY_01"]])
    expected_transformer2 = Set([components["TST_HP_01"],
                                 components["TST_HP_01b"]])

    chains_bus = Resie.find_chains([u for u in values(components) if !startswith(u.uac, "Proxy")],
                                   EnergySystems.sf_bus)

    chains_transformer = Resie.find_chains([u for u in values(components) if !startswith(u.uac, "Proxy")],
                                           EnergySystems.sf_transformer;
                                           direct_connection_only=false)

    @test length(chains_bus) == 2
    @test chains_bus[1] == expected_bus1
    @test chains_bus[2] == expected_bus2
    @test length(chains_transformer) == 2
    @test chains_transformer[1] == expected_transformer2
    @test chains_transformer[2] == expected_transformer1

    # test distance to sink, considers energy flow matrix of bus
    @test Resie.distance_to_sink(components["TST_HP_01"], EnergySystems.sf_transformer, [], "") == 1
    @test Resie.distance_to_sink(components["TST_HP_01b"], EnergySystems.sf_transformer, [], "") == 0
    @test Resie.distance_to_sink(components["TST_HP_02"], EnergySystems.sf_transformer, [], "") == 2
    @test Resie.distance_to_sink(components["TST_HP_03"], EnergySystems.sf_transformer, [], "") == 0
    @test Resie.distance_to_sink(components["TST_HP_04"], EnergySystems.sf_transformer, [], "") == 0
    @test Resie.distance_to_sink(components["TST_01_ELY_01"], EnergySystems.sf_transformer, [], "") == 4

    # test function to determine if transformers are at interfaces, considers energy flow matrix of bus
    @test Resie.check_interface_for_transformer(components["TST_HP_02"].output_interfaces[components["TST_HP_02"].m_heat_out],
                                                "output") == true
    @test Resie.check_interface_for_transformer(components["TST_HP_01b"].output_interfaces[components["TST_HP_01b"].m_heat_out],
                                                "output") == false
    @test Resie.check_interface_for_transformer(components["TST_HP_01b"].input_interfaces[components["TST_HP_01b"].m_heat_in],
                                                "input") == true
    @test Resie.check_interface_for_transformer(components["TST_01_ELY_01"].output_interfaces[components["TST_01_ELY_01"].m_heat_ht_out],
                                                "output") == true
    @test Resie.check_interface_for_transformer(components["TST_HP_03"].output_interfaces[components["TST_HP_03"].m_heat_out],
                                                "output") == false
    @test Resie.check_interface_for_transformer(components["TST_HP_03"].input_interfaces[components["TST_HP_03"].m_heat_in],
                                                "input") == true
    @test Resie.check_interface_for_transformer(components["TST_HP_04"].input_interfaces[components["TST_HP_04"].m_heat_in],
                                                "input") == true
end

function test_find_indirect_chains_denied2()
    # This is the same energy system as in test_find_indirect_chains(), but the connection matrix
    # of the two busses denies the connection to HP4 from other transformers. 
    # Therefore, here two independent transformer chains are detected instead of one. One chain
    # is only HP4.
    components_config = Dict{String,Any}(
        "TST_DEM_01" => Dict{String,Any}(
            "type" => "Demand",
            "medium" => "m_h_w_ht1",
            "output_refs" => [],
            "energy_profile_file_path" => "./profiles/examples/general/dem_heat_nrg_var_hi-amp.prf",
            "constant_temperature" => 70,
            "scale" => 500,
        ),
        "TST_DEM_02" => Dict{String,Any}(
            "type" => "Demand",
            "medium" => "m_h_w_ht1",
            "output_refs" => [],
            "energy_profile_file_path" => "./profiles/examples/general/dem_heat_nrg_var_hi-amp.prf",
            "constant_temperature" => 90,
            "scale" => 500,
        ),
        "TST_DEM_03" => Dict{String,Any}(
            "type" => "Demand",
            "medium" => "m_h_w_ht1",
            "output_refs" => [],
            "energy_profile_file_path" => "./profiles/examples/general/dem_heat_nrg_var_hi-amp.prf",
            "constant_temperature" => 95,
            "scale" => 500,
        ),
        "TST_GRI_H2" => Dict{String,Any}(
            "type" => "GridConnection",
            "medium" => "m_c_g_h2",
            "input_refs" => ["TST_01_ELY_01"],
            "is_source" => false,
        ),
        "TST_GRI_O2" => Dict{String,Any}(
            "type" => "GridConnection",
            "medium" => "m_c_g_o2",
            "input_refs" => ["TST_01_ELY_01"],
            "is_source" => false,
        ),
        "TST_SRC_01" => Dict{String,Any}(
            "type" => "BoundedSupply",
            "medium" => "m_h_w_lt1",
            "output_refs" => ["TST_BUS_00"],
            "max_power_profile_file_path" => "./profiles/examples/general/src_heat_maxpow_var_lo-amp.prf",
            "temperature_profile_file_path" => "./profiles/examples/general/src_heat_temp_var_avg25.prf",
            "scale" => 1000,
        ),
        "TST_GRI_00" => Dict{String,Any}(
            "type" => "GridConnection",
            "medium" => "m_e_ac_230v",
            "output_refs" => ["TST_01_ELY_01"],
            "is_source" => true,
        ),
        "TST_GRI_01" => Dict{String,Any}(
            "type" => "GridConnection",
            "medium" => "m_e_ac_230v",
            "output_refs" => ["TST_HP_01"],
            "is_source" => true,
        ),
        "TST_SRC_1b" => Dict{String,Any}(
            "type" => "BoundedSupply",
            "medium" => "m_e_ac_230v",
            "output_refs" => ["TST_HP_01b"],
            "is_source" => true,
            "constant_power" => 400,
        ),
        "TST_GRI_02" => Dict{String,Any}(
            "type" => "GridConnection",
            "medium" => "m_e_ac_230v",
            "output_refs" => ["TST_HP_02"],
            "is_source" => true,
        ),
        "TST_GRI_03" => Dict{String,Any}(
            "type" => "GridConnection",
            "medium" => "m_e_ac_230v",
            "output_refs" => ["TST_HP_03"],
            "is_source" => true,
        ),
        "TST_GRI_04" => Dict{String,Any}(
            "type" => "GridConnection",
            "medium" => "m_e_ac_230v",
            "output_refs" => ["TST_HP_04"],
            "is_source" => true,
        ),
        "TST_HP_01" => Dict{String,Any}(
            "type" => "HeatPump",
            "output_refs" => ["TST_HP_01b"],
            "power_th" => 1200,
            "output_temperature" => 60,
            "min_power_fraction" => 0.0,
        ),
        "TST_HP_01b" => Dict{String,Any}(
            "type" => "HeatPump",
            "m_heat_in" => "m_h_w_ht1",
            "output_refs" => ["TST_BUS_01"],
            "power_th" => 9000,
            "min_power_fraction" => 0.0,
        ),
        "TST_HP_02" => Dict{String,Any}(
            "type" => "HeatPump",
            "output_refs" => ["TST_BUS_01"],
            "power_th" => 9000,
            "min_power_fraction" => 0.0,
        ),
        "TST_HP_03" => Dict{String,Any}(
            "type" => "HeatPump",
            "m_heat_in" => "m_h_w_ht1",
            "output_refs" => ["TST_DEM_02"],
            "power_th" => 9000,
            "min_power_fraction" => 0.0,
        ),
        "TST_HP_04" => Dict{String,Any}(
            "type" => "HeatPump",
            "m_heat_in" => "m_h_w_ht1",
            "output_refs" => ["TST_DEM_03"],
            "power_th" => 9000,
            "min_power_fraction" => 0.0,
        ),
        "TST_01_ELY_01" => Dict{String,Any}(
            "type" => "Electrolyser",
            "output_refs" => ["TST_BUS_00",
                              "TST_GRI_H2",
                              "TST_GRI_O2"],
            "m_heat_ht_out" => "m_h_w_lt1",
            "output_temperature" => 45,
            "power_el" => 40000,
            "heat_lt_is_usable" => false,
            "nr_switchable_units" => 4,
            "dispatch_strategy" => "equal_with_mpf",
            "min_power_fraction_total" => 0.0,
            "min_power_fraction" => 0.0,
        ),
        "TST_BUS_00" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_h_w_lt1",
            "connections" => Dict{String,Any}(
                "input_order" => ["TST_01_ELY_01", "TST_SRC_01"],
                "output_order" => ["TST_HP_02", "TST_HP_01"],
                "energy_flow" => [[1, 1],
                                  [1, 1]],
            ),
        ),
        "TST_BUS_01" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_h_w_ht1",
            "connections" => Dict{String,Any}(
                "input_order" => ["TST_HP_01b", "TST_HP_02"],
                "output_order" => ["TST_HP_03", "TST_HP_04", "TST_DEM_01"],
                "energy_flow" => [[1, 0, 1],
                                  [1, 0, 1]],
            ),
        ),
    )

    simulation_parameters = Dict{String,Any}(
        "time_step_seconds" => 900,
        "time" => 0,
        "epsilon" => 1e-9,
        "is_first_timestep" => true,
    )

    components = Resie.load_components(components_config, simulation_parameters)

    expected_bus1 = Set([components["TST_BUS_00"]])
    expected_bus2 = Set([components["TST_BUS_01"]])
    expected_transformer1 = Set([components["TST_HP_01"],
                                 components["TST_HP_01b"],
                                 components["TST_HP_02"],
                                 components["TST_HP_03"],
                                 components["TST_01_ELY_01"]])
    expected_transformer2 = Set([components["TST_HP_04"]])

    chains_bus = Resie.find_chains([u for u in values(components) if !startswith(u.uac, "Proxy")],
                                   EnergySystems.sf_bus)

    chains_transformer = Resie.find_chains([u for u in values(components) if !startswith(u.uac, "Proxy")],
                                           EnergySystems.sf_transformer;
                                           direct_connection_only=false)

    @test length(chains_bus) == 2
    @test chains_bus[1] == expected_bus1
    @test chains_bus[2] == expected_bus2
    @test length(chains_transformer) == 2
    @test chains_transformer[1] == expected_transformer1
    @test chains_transformer[2] == expected_transformer2

    # test distance to sink, considers energy flow matrix of bus
    @test Resie.distance_to_sink(components["TST_HP_01"], EnergySystems.sf_transformer, [], "") == 3
    @test Resie.distance_to_sink(components["TST_HP_01b"], EnergySystems.sf_transformer, [], "") == 2
    @test Resie.distance_to_sink(components["TST_HP_02"], EnergySystems.sf_transformer, [], "") == 2
    @test Resie.distance_to_sink(components["TST_HP_03"], EnergySystems.sf_transformer, [], "") == 0
    @test Resie.distance_to_sink(components["TST_HP_04"], EnergySystems.sf_transformer, [], "") == 0
    @test Resie.distance_to_sink(components["TST_01_ELY_01"], EnergySystems.sf_transformer, [], "") == 5

    # test function to determine if transformers are at interfaces, considers energy flow matrix of bus
    @test Resie.check_interface_for_transformer(components["TST_HP_02"].output_interfaces[components["TST_HP_02"].m_heat_out],
                                                "output") == true
    @test Resie.check_interface_for_transformer(components["TST_HP_01b"].output_interfaces[components["TST_HP_01b"].m_heat_out],
                                                "output") == true
    @test Resie.check_interface_for_transformer(components["TST_HP_01b"].input_interfaces[components["TST_HP_01b"].m_heat_in],
                                                "input") == true
    @test Resie.check_interface_for_transformer(components["TST_01_ELY_01"].output_interfaces[components["TST_01_ELY_01"].m_heat_ht_out],
                                                "output") == true
    @test Resie.check_interface_for_transformer(components["TST_HP_03"].output_interfaces[components["TST_HP_03"].m_heat_out],
                                                "output") == false
    @test Resie.check_interface_for_transformer(components["TST_HP_03"].input_interfaces[components["TST_HP_03"].m_heat_in],
                                                "input") == true
    @test Resie.check_interface_for_transformer(components["TST_HP_04"].input_interfaces[components["TST_HP_04"].m_heat_in],
                                                "input") == false
end

@testset "find_chains" begin
    test_find_chains()
    test_find_indirect_chains()
    test_find_indirect_chains_denied()
    test_find_indirect_chains_denied2()
end

function test_find_indexes()
    # normal case
    steps = [[2, ("TST_BUS_02", EnergySystems.s_control)],
             [3, ("TST_BUS_03", EnergySystems.s_control)],
             [1, ("TST_BUS_01", EnergySystems.s_control)]]
    own_idx, target_idx = Resie.find_indexes(steps,
                                             ("TST_BUS_02", EnergySystems.s_control),
                                             ("TST_BUS_01", EnergySystems.s_control))
    @test own_idx == 1
    @test target_idx == 3

    # can't find either step due to wrong UAC in one case and wrong step in the other
    steps = [[2, ("TST_BUS_02", EnergySystems.s_control)],
             [3, ("TST_BUS_03", EnergySystems.s_control)],
             [1, ("TST_BUS_01", EnergySystems.s_control)]]
    own_idx, target_idx = Resie.find_indexes(steps,
                                             ("TST_BUS_02", EnergySystems.s_process),
                                             ("TST_BSS_01", EnergySystems.s_control))
    @test own_idx === nothing
    @test target_idx === nothing
end

@testset "find_indexes" begin
    test_find_indexes()
end

function test_place_one_higher()
    # normal test case
    steps = [[2, ("TST_BUS_02", EnergySystems.s_control)],
             [3, ("TST_BUS_03", EnergySystems.s_control)],
             [1, ("TST_BUS_01", EnergySystems.s_control)]]
    Resie.place_one_higher!(steps,
                            ("TST_BUS_02", EnergySystems.s_control),
                            ("TST_BUS_01", EnergySystems.s_control))
    @test steps[1] == [2, ("TST_BUS_02", EnergySystems.s_control)]
    @test steps[2] == [4, ("TST_BUS_03", EnergySystems.s_control)]
    @test steps[3] == [3, ("TST_BUS_01", EnergySystems.s_control)]

    # priority already is higher and force was not used => no change
    steps = [[2, ("TST_BUS_02", EnergySystems.s_control)],
             [3, ("TST_BUS_03", EnergySystems.s_control)],
             [4, ("TST_BUS_01", EnergySystems.s_control)]]
    Resie.place_one_higher!(steps,
                            ("TST_BUS_02", EnergySystems.s_control),
                            ("TST_BUS_01", EnergySystems.s_control))
    @test steps[1] == [2, ("TST_BUS_02", EnergySystems.s_control)]
    @test steps[2] == [3, ("TST_BUS_03", EnergySystems.s_control)]
    @test steps[3] == [4, ("TST_BUS_01", EnergySystems.s_control)]

    # priority already is higher and force was used => works as normal
    steps = [[2, ("TST_BUS_02", EnergySystems.s_control)],
             [3, ("TST_BUS_03", EnergySystems.s_control)],
             [4, ("TST_BUS_01", EnergySystems.s_control)]]
    Resie.place_one_higher!(steps,
                            ("TST_BUS_02", EnergySystems.s_control),
                            ("TST_BUS_01", EnergySystems.s_control);
                            force=true)
    @test steps[1] == [2, ("TST_BUS_02", EnergySystems.s_control)]
    @test steps[2] == [4, ("TST_BUS_03", EnergySystems.s_control)]
    @test steps[3] == [3, ("TST_BUS_01", EnergySystems.s_control)]
end

@testset "place_one_higher" begin
    test_place_one_higher()
end

function test_place_one_lower()
    # normal test case
    steps = [[2, ("TST_BUS_02", EnergySystems.s_control)],
             [3, ("TST_BUS_03", EnergySystems.s_control)],
             [1, ("TST_BUS_01", EnergySystems.s_control)]]
    Resie.place_one_lower!(steps,
                           ("TST_BUS_02", EnergySystems.s_control),
                           ("TST_BUS_03", EnergySystems.s_control))
    @test steps[1] == [2, ("TST_BUS_02", EnergySystems.s_control)]
    @test steps[2] == [1, ("TST_BUS_03", EnergySystems.s_control)]
    @test steps[3] == [0, ("TST_BUS_01", EnergySystems.s_control)]

    # priority already is lower and force was not used => no change
    steps = [[2, ("TST_BUS_02", EnergySystems.s_control)],
             [3, ("TST_BUS_03", EnergySystems.s_control)],
             [1, ("TST_BUS_01", EnergySystems.s_control)]]
    Resie.place_one_lower!(steps,
                           ("TST_BUS_03", EnergySystems.s_control),
                           ("TST_BUS_01", EnergySystems.s_control))
    @test steps[1] == [2, ("TST_BUS_02", EnergySystems.s_control)]
    @test steps[2] == [3, ("TST_BUS_03", EnergySystems.s_control)]
    @test steps[3] == [1, ("TST_BUS_01", EnergySystems.s_control)]

    # priority already is lower and force was used => works as normal
    steps = [[2, ("TST_BUS_02", EnergySystems.s_control)],
             [3, ("TST_BUS_03", EnergySystems.s_control)],
             [1, ("TST_BUS_01", EnergySystems.s_control)]]
    Resie.place_one_lower!(steps,
                           ("TST_BUS_03", EnergySystems.s_control),
                           ("TST_BUS_01", EnergySystems.s_control);
                           force=true)
    @test steps[1] == [1, ("TST_BUS_02", EnergySystems.s_control)]
    @test steps[2] == [3, ("TST_BUS_03", EnergySystems.s_control)]
    @test steps[3] == [2, ("TST_BUS_01", EnergySystems.s_control)]
end

@testset "place_one_lower" begin
    test_place_one_lower()
end
