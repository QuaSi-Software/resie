using Debugger
using Test
using Resie
using Resie.EnergySystems

include("../test_util.jl")

function test_data_input_priorities()
    components_config = Dict{String,Any}(
        "TST_BUS_01" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_h_w_ht1",
            "connections" => Dict{String,Any}(
                "input_order" => [],
                "output_order" => ["TST_BUS_03"],
            ),
        ),
        "TST_BUS_02" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_h_w_ht1",
            "connections" => Dict{String,Any}(
                "input_order" => [],
                "output_order" => ["TST_BUS_03"],
            ),
        ),
        "TST_BUS_03" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_h_w_ht1",
            "connections" => Dict{String,Any}(
                "input_order" => ["TST_BUS_02",
                                  "TST_BUS_01"],
                "output_order" => [],
            ),
        ),
    )

    simulation_parameters = get_default_sim_params()

    components = Resie.load_components(components_config, simulation_parameters)
    by_function = Resie.categorize_by_function(components)
    return components, by_function
end

function test_input_priorities_reordered_inputs()
    components, by_function = test_data_input_priorities()
    steps = [[100, ("TST_BUS_03", EnergySystems.s_control)],
             [99, ("TST_BUS_01", EnergySystems.s_control)],
             [98, ("TST_BUS_02", EnergySystems.s_control)],
             [97, ("TST_BUS_03", EnergySystems.s_process)],
             [96, ("TST_BUS_01", EnergySystems.s_process)],
             [95, ("TST_BUS_02", EnergySystems.s_process)]]
    expected = [[100, ("TST_BUS_03", EnergySystems.s_control)],
                [99, ("TST_BUS_01", EnergySystems.s_control)],
                [98, ("TST_BUS_02", EnergySystems.s_control)],
                [97, ("TST_BUS_03", EnergySystems.s_process)],
                [96, ("TST_BUS_01", EnergySystems.s_process)],
                [95, ("TST_BUS_02", EnergySystems.s_process)]]
    Resie.reorder_for_input_priorities(steps, components, by_function)
    @test pwc_steps_astr(expected, steps) == ""
end

@testset "input_priorities_reordered_inputs" begin
    test_input_priorities_reordered_inputs()
end

function test_input_priorities_no_change()
    components, by_function = test_data_input_priorities()
    steps = [[100, ("TST_BUS_03", EnergySystems.s_control)],
             [99, ("TST_BUS_02", EnergySystems.s_control)],
             [98, ("TST_BUS_01", EnergySystems.s_control)],
             [97, ("TST_BUS_03", EnergySystems.s_process)],
             [96, ("TST_BUS_02", EnergySystems.s_process)],
             [95, ("TST_BUS_01", EnergySystems.s_process)]]
    expected = [[100, ("TST_BUS_03", EnergySystems.s_control)],
                [99, ("TST_BUS_02", EnergySystems.s_control)],
                [98, ("TST_BUS_01", EnergySystems.s_control)],
                [97, ("TST_BUS_03", EnergySystems.s_process)],
                [96, ("TST_BUS_02", EnergySystems.s_process)],
                [95, ("TST_BUS_01", EnergySystems.s_process)]]
    Resie.reorder_for_input_priorities(steps, components, by_function)
    @test pwc_steps_astr(expected, steps) == ""
end

@testset "input_priorities_no_change" begin
    test_input_priorities_no_change()
end

function test_data_busses_distribute()
    #                ------------
    #                |   Bus 2  |   ------------
    #                ------------   |   Bus 5  |
    #              /              / ------------
    # ------------   ------------
    # |   Bus 1  |---|   Bus 3  |
    # ------------   ------------
    #              \              \ ------------
    #                ------------   |   Bus 6  |
    #                |   Bus 4  |   ------------
    #                ------------
    components_config = Dict{String,Any}(
        "TST_BUS_01" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_h_w_ht1",
            "connections" => Dict{String,Any}(
                "input_order" => [],
                "output_order" => ["TST_BUS_04",
                                   "TST_BUS_03",
                                   "TST_BUS_02"],
            ),
        ),
        "TST_BUS_02" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_h_w_ht1",
            "connections" => Dict{String,Any}(
                "input_order" => ["TST_BUS_01"],
                "output_order" => [],
            ),
        ),
        "TST_BUS_03" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_h_w_ht1",
            "connections" => Dict{String,Any}(
                "input_order" => ["TST_BUS_01"],
                "output_order" => ["TST_BUS_06",
                                   "TST_BUS_05"],
            ),
        ),
        "TST_BUS_04" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_h_w_ht1",
            "connections" => Dict{String,Any}(
                "input_order" => ["TST_BUS_01"],
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
        "TST_BUS_06" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_h_w_ht1",
            "connections" => Dict{String,Any}(
                "input_order" => ["TST_BUS_03"],
                "output_order" => [],
            ),
        ),
    )

    simulation_parameters = get_default_sim_params()

    components = Resie.load_components(components_config, simulation_parameters)
    by_function = Resie.categorize_by_function(components)
    return components, by_function
end

function test_busses_distribution_no_change()
    components, by_function = test_data_busses_distribute()
    steps = [[100, ("TST_BUS_04", EnergySystems.s_distribute)],
             [99, ("TST_BUS_06", EnergySystems.s_distribute)],
             [98, ("TST_BUS_05", EnergySystems.s_distribute)],
             [97, ("TST_BUS_03", EnergySystems.s_distribute)],
             [96, ("TST_BUS_02", EnergySystems.s_distribute)],
             [95, ("TST_BUS_01", EnergySystems.s_distribute)]]
    expected = [[100, ("TST_BUS_04", EnergySystems.s_distribute)],
                [99, ("TST_BUS_06", EnergySystems.s_distribute)],
                [98, ("TST_BUS_05", EnergySystems.s_distribute)],
                [97, ("TST_BUS_03", EnergySystems.s_distribute)],
                [96, ("TST_BUS_02", EnergySystems.s_distribute)],
                [95, ("TST_BUS_01", EnergySystems.s_distribute)]]
    Resie.reorder_distribution_of_busses(steps, components, by_function)
    @test pwc_steps_astr(expected, steps) == ""
end

@testset "busses_distribution_no_change" begin
    test_busses_distribution_no_change()
end

function test_busses_distribution_reorder_steps()
    components, by_function = test_data_busses_distribute()
    steps = [[100, ("TST_BUS_01", EnergySystems.s_distribute)],
             [99, ("TST_BUS_02", EnergySystems.s_distribute)],
             [98, ("TST_BUS_03", EnergySystems.s_distribute)],
             [97,
              ("Proxy-TST_BUS_01|TST_BUS_04|Proxy-TST_BUS_03|TST_BUS_05|TST_BUS_06|TST_BUS_02",
               EnergySystems.s_distribute)],
             [96, ("TST_BUS_04", EnergySystems.s_distribute)],
             [95, ("TST_BUS_05", EnergySystems.s_distribute)],
             [94, ("TST_BUS_06", EnergySystems.s_distribute)]]
    expected = [[94, ("TST_BUS_01", EnergySystems.s_distribute)],
                [95, ("TST_BUS_02", EnergySystems.s_distribute)],
                [96, ("TST_BUS_03", EnergySystems.s_distribute)],
                [97,
                 ("Proxy-TST_BUS_01|TST_BUS_04|Proxy-TST_BUS_03|TST_BUS_05|TST_BUS_06|TST_BUS_02",
                  EnergySystems.s_distribute)],
                [93, ("TST_BUS_04", EnergySystems.s_distribute)],
                [92, ("TST_BUS_05", EnergySystems.s_distribute)],
                [91, ("TST_BUS_06", EnergySystems.s_distribute)]]
    Resie.reorder_distribution_of_busses(steps, components, by_function)
    @test pwc_steps_astr(expected, steps) == ""
end

@testset "busses_distribution_reorder_steps" begin
    test_busses_distribution_reorder_steps()
end

function test_data_storage_loading()
    components_config = Dict{String,Any}(
        "TST_BUS_01" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_h_w_ht1",
            "connections" => Dict{String,Any}(
                "input_order" => ["TST_BFT_01"],
                "output_order" => ["TST_BUS_02",
                                   "TST_BFT_01",
                                   "TST_BUS_03"],
            ),
        ),
        "TST_BUS_02" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_h_w_ht1",
            "connections" => Dict{String,Any}(
                "input_order" => ["TST_BUS_01",
                                  "TST_BFT_02"],
                "output_order" => ["TST_BFT_02"],
            ),
        ),
        "TST_BUS_03" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_h_w_ht1",
            "connections" => Dict{String,Any}(
                "input_order" => ["TST_BUS_01",
                                  "TST_BFT_03"],
                "output_order" => ["TST_BFT_03"],
            ),
        ),
        "TST_BFT_01" => Dict{String,Any}(
            "type" => "BufferTank",
            "output_refs" => ["TST_BUS_01"],
            "model_type" => "ideally_stratified",
            "capacity" => 40000,
            "initial_load" => 0,
        ),
        "TST_BFT_02" => Dict{String,Any}(
            "type" => "BufferTank",
            "output_refs" => ["TST_BUS_02"],
            "model_type" => "ideally_stratified",
            "capacity" => 40000,
            "initial_load" => 0,
        ),
        "TST_BFT_03" => Dict{String,Any}(
            "type" => "BufferTank",
            "output_refs" => ["TST_BUS_03"],
            "model_type" => "ideally_stratified",
            "capacity" => 40000,
            "initial_load" => 0,
        ),
    )

    simulation_parameters = get_default_sim_params()

    components = Resie.load_components(components_config, simulation_parameters)
    by_function = Resie.categorize_by_function(components)
    return components, by_function
end

function test_storage_loading_no_change()
    components, by_function = test_data_storage_loading()
    steps = [[100, ("TST_BFT_01", EnergySystems.s_process)],
             [99, ("TST_BFT_03", EnergySystems.s_process)],
             [98, ("TST_BFT_02", EnergySystems.s_process)],
             [97, ("TST_BFT_02", EnergySystems.s_load)],
             [96, ("TST_BFT_01", EnergySystems.s_load)],
             [95, ("TST_BFT_03", EnergySystems.s_load)]]
    expected = [[100, ("TST_BFT_01", EnergySystems.s_process)],
                [99, ("TST_BFT_03", EnergySystems.s_process)],
                [98, ("TST_BFT_02", EnergySystems.s_process)],
                [97, ("TST_BFT_02", EnergySystems.s_load)],
                [96, ("TST_BFT_01", EnergySystems.s_load)],
                [95, ("TST_BFT_03", EnergySystems.s_load)]]
    Resie.reorder_storage_loading(steps, components, by_function)
    @test pwc_steps_astr(expected, steps) == ""
end

@testset "storage_loading_no_change" begin
    test_storage_loading_no_change()
end

function test_storage_loading_reorder_steps_1()
    components, by_function = test_data_storage_loading()
    steps = [[100, ("TST_BFT_01", EnergySystems.s_process)],
             [99, ("TST_BFT_02", EnergySystems.s_process)],
             [98, ("TST_BFT_03", EnergySystems.s_process)],
             [97, ("TST_BFT_01", EnergySystems.s_load)],
             [96, ("TST_BFT_02", EnergySystems.s_load)],
             [95, ("TST_BFT_03", EnergySystems.s_load)]]
    expected = [[100, ("TST_BFT_01", EnergySystems.s_process)],
                [99, ("TST_BFT_02", EnergySystems.s_process)],
                [98, ("TST_BFT_03", EnergySystems.s_process)],
                [95, ("TST_BFT_01", EnergySystems.s_load)],
                [96, ("TST_BFT_02", EnergySystems.s_load)],
                [94, ("TST_BFT_03", EnergySystems.s_load)]]
    Resie.reorder_storage_loading(steps, components, by_function)
    @test pwc_steps_astr(expected, steps) == ""
end

function test_storage_loading_reorder_steps_2()
    components, by_function = test_data_storage_loading()
    steps = [[100, ("TST_BFT_03", EnergySystems.s_process)],
             [99, ("TST_BFT_02", EnergySystems.s_process)],
             [98, ("TST_BFT_01", EnergySystems.s_process)],
             [97, ("TST_BFT_03", EnergySystems.s_load)],
             [96, ("TST_BFT_02", EnergySystems.s_load)],
             [95, ("TST_BFT_01", EnergySystems.s_load)]]
    expected = [[100, ("TST_BFT_03", EnergySystems.s_process)],
                [99, ("TST_BFT_02", EnergySystems.s_process)],
                [98, ("TST_BFT_01", EnergySystems.s_process)],
                [94, ("TST_BFT_03", EnergySystems.s_load)],
                [96, ("TST_BFT_02", EnergySystems.s_load)],
                [95, ("TST_BFT_01", EnergySystems.s_load)]]
    Resie.reorder_storage_loading(steps, components, by_function)
    @test pwc_steps_astr(expected, steps) == ""
end

function test_storage_loading_reorder_steps_3()
    components, by_function = test_data_storage_loading()
    steps = [[100, ("TST_BFT_03", EnergySystems.s_process)],
             [99, ("TST_BFT_01", EnergySystems.s_process)],
             [98, ("TST_BFT_02", EnergySystems.s_process)],
             [97, ("TST_BFT_03", EnergySystems.s_load)],
             [96, ("TST_BFT_01", EnergySystems.s_load)],
             [95, ("TST_BFT_02", EnergySystems.s_load)]]
    expected = [[100, ("TST_BFT_03", EnergySystems.s_process)],
                [99, ("TST_BFT_01", EnergySystems.s_process)],
                [98, ("TST_BFT_02", EnergySystems.s_process)],
                [93, ("TST_BFT_03", EnergySystems.s_load)],
                [94, ("TST_BFT_01", EnergySystems.s_load)],
                [95, ("TST_BFT_02", EnergySystems.s_load)]]
    Resie.reorder_storage_loading(steps, components, by_function)
    @test pwc_steps_astr(expected, steps) == ""
end

@testset "storage_loading_reorder_steps" begin
    test_storage_loading_reorder_steps_1()
    test_storage_loading_reorder_steps_2()
    test_storage_loading_reorder_steps_3()
end

function test_data_storage_loading_with_matrix()
    components_config = Dict{String,Any}(
        "TST_GRI_01" => Dict{String,Any}(
            "type" => "GridConnection",
            "medium" => "m_h_w_ht1",
            "output_refs" => ["TST_BUS_01"],
            "is_source" => true,
        ),
        "TST_BUS_01" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_h_w_ht1",
            "connections" => Dict{String,Any}(
                "input_order" => ["TST_BFT_01",
                                  "TST_GRI_01"],
                "output_order" => ["TST_BUS_02",
                                   "TST_BFT_01",
                                   "TST_BUS_03"],
                "energy_flow" => [[1, 0, 1],
                                  [1, 0, 1]],
            ),
        ),
        "TST_BUS_02" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_h_w_ht1",
            "connections" => Dict{String,Any}(
                "input_order" => ["TST_BFT_02",
                                  "TST_BUS_01"],
                "output_order" => ["TST_BFT_02"],
                "energy_flow" => [[0],
                                  [0]],
            ),
        ),
        "TST_BUS_03" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_h_w_ht1",
            "connections" => Dict{String,Any}(
                "input_order" => ["TST_BFT_03",
                                  "TST_BUS_01"],
                "output_order" => ["TST_BFT_03"],
                "energy_flow" => [[0],
                                  [1]],
            ),
        ),
        "TST_BFT_01" => Dict{String,Any}(
            "type" => "BufferTank",
            "output_refs" => ["TST_BUS_01"],
            "model_type" => "ideally_stratified",
            "capacity" => 40000,
            "initial_load" => 0,
        ),
        "TST_BFT_02" => Dict{String,Any}(
            "type" => "BufferTank",
            "output_refs" => ["TST_BUS_02"],
            "model_type" => "ideally_stratified",
            "capacity" => 40000,
            "initial_load" => 0,
        ),
        "TST_BFT_03" => Dict{String,Any}(
            "type" => "BufferTank",
            "output_refs" => ["TST_BUS_03"],
            "model_type" => "ideally_stratified",
            "capacity" => 40000,
            "initial_load" => 0,
        ),
    )

    simulation_parameters = get_default_sim_params()

    components = Resie.load_components(components_config, simulation_parameters)
    by_function = Resie.categorize_by_function(components)
    return components, by_function
end

function test_storage_loading_reorder_steps_with_matrix_1()
    components, by_function = test_data_storage_loading_with_matrix()
    steps = [[100, ("TST_BFT_01", EnergySystems.s_process)],
             [99, ("TST_BFT_02", EnergySystems.s_process)],
             [98, ("TST_BFT_03", EnergySystems.s_process)],
             [97, ("TST_BFT_01", EnergySystems.s_load)],
             [96, ("TST_BFT_02", EnergySystems.s_load)],
             [95, ("TST_BFT_03", EnergySystems.s_load)]]
    expected = [[100, ("TST_BFT_01", EnergySystems.s_process)],
                [99, ("TST_BFT_02", EnergySystems.s_process)],
                [98, ("TST_BFT_03", EnergySystems.s_process)],
                [95, ("TST_BFT_01", EnergySystems.s_load)],
                [96, ("TST_BFT_02", EnergySystems.s_load)],
                [94, ("TST_BFT_03", EnergySystems.s_load)]]
    Resie.reorder_storage_loading(steps, components, by_function)
    @test pwc_steps_astr(expected, steps) == ""
end

function test_storage_loading_reorder_steps_with_matrix_2()
    components, by_function = test_data_storage_loading_with_matrix()
    steps = [[100, ("TST_BFT_03", EnergySystems.s_process)],
             [99, ("TST_BFT_02", EnergySystems.s_process)],
             [98, ("TST_BFT_01", EnergySystems.s_process)],
             [97, ("TST_BFT_03", EnergySystems.s_load)],
             [96, ("TST_BFT_02", EnergySystems.s_load)],
             [95, ("TST_BFT_01", EnergySystems.s_load)]]
    expected = [[100, ("TST_BFT_03", EnergySystems.s_process)],
                [99, ("TST_BFT_02", EnergySystems.s_process)],
                [98, ("TST_BFT_01", EnergySystems.s_process)],
                [94, ("TST_BFT_03", EnergySystems.s_load)],
                [96, ("TST_BFT_02", EnergySystems.s_load)],
                [95, ("TST_BFT_01", EnergySystems.s_load)]]
    Resie.reorder_storage_loading(steps, components, by_function)
    @test pwc_steps_astr(expected, steps) == ""
end

function test_storage_loading_reorder_steps_with_matrix_3()
    components, by_function = test_data_storage_loading_with_matrix()
    steps = [[100, ("TST_BFT_03", EnergySystems.s_process)],
             [99, ("TST_BFT_01", EnergySystems.s_process)],
             [98, ("TST_BFT_02", EnergySystems.s_process)],
             [97, ("TST_BFT_03", EnergySystems.s_load)],
             [96, ("TST_BFT_01", EnergySystems.s_load)],
             [95, ("TST_BFT_02", EnergySystems.s_load)]]
    expected = [[100, ("TST_BFT_03", EnergySystems.s_process)],
                [99, ("TST_BFT_01", EnergySystems.s_process)],
                [98, ("TST_BFT_02", EnergySystems.s_process)],
                [93, ("TST_BFT_03", EnergySystems.s_load)],
                [94, ("TST_BFT_01", EnergySystems.s_load)],
                [95, ("TST_BFT_02", EnergySystems.s_load)]]
    Resie.reorder_storage_loading(steps, components, by_function)
    @test pwc_steps_astr(expected, steps) == ""
end

@testset "storage_loading_reorder_steps_with_matrix" begin
    test_storage_loading_reorder_steps_with_matrix_1()
    test_storage_loading_reorder_steps_with_matrix_2()
    test_storage_loading_reorder_steps_with_matrix_3()
end
