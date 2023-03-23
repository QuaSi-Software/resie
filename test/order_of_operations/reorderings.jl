using Debugger
using Test
using Resie
using Resie.EnergySystems

include("../test_util.jl")

function test_data_input_priorities()
    systems_config = Dict{String,Any}(
        "TST_BUS_01" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "production_refs" => ["TST_BUS_03"],
        ),
        "TST_BUS_02" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "production_refs" => ["TST_BUS_03"],
        ),
        "TST_BUS_03" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "production_refs" => [],
            "connection_matrix" => Dict{String, Any}(
                "input_order" => [
                    "TST_BUS_02",
                    "TST_BUS_01",
                ],
                "output_order" => [],
            )
        ),
    )
    systems = Resie.load_systems(systems_config)
    by_function = Resie.categorize_by_function(systems)
    return systems, by_function
end

function test_input_priorities_reordered_inputs()
    systems, by_function = test_data_input_priorities()
    steps = [
        [100, ("TST_BUS_03", EnergySystems.s_control)],
        [99, ("TST_BUS_01", EnergySystems.s_control)],
        [98, ("TST_BUS_02", EnergySystems.s_control)],
        [97, ("TST_BUS_03", EnergySystems.s_produce)],
        [96, ("TST_BUS_01", EnergySystems.s_produce)],
        [95, ("TST_BUS_02", EnergySystems.s_produce)],
    ]
    expected = [
        [100, ("TST_BUS_03", EnergySystems.s_control)],
        [97, ("TST_BUS_01", EnergySystems.s_control)],
        [98, ("TST_BUS_02", EnergySystems.s_control)],
        [96, ("TST_BUS_03", EnergySystems.s_produce)],
        [93, ("TST_BUS_01", EnergySystems.s_produce)],
        [94, ("TST_BUS_02", EnergySystems.s_produce)],
    ]
    Resie.reorder_for_input_priorities(steps, systems, by_function)
    @test pwc_steps_astr(expected, steps) == ""
end

@testset "input_priorities_reordered_inputs" begin
    test_input_priorities_reordered_inputs()
end

function test_input_priorities_no_change()
    systems, by_function = test_data_input_priorities()
    steps = [
        [100, ("TST_BUS_03", EnergySystems.s_control)],
        [99, ("TST_BUS_02", EnergySystems.s_control)],
        [98, ("TST_BUS_01", EnergySystems.s_control)],
        [97, ("TST_BUS_03", EnergySystems.s_produce)],
        [96, ("TST_BUS_02", EnergySystems.s_produce)],
        [95, ("TST_BUS_01", EnergySystems.s_produce)],
    ]
    expected = [
        [100, ("TST_BUS_03", EnergySystems.s_control)],
        [99, ("TST_BUS_02", EnergySystems.s_control)],
        [98, ("TST_BUS_01", EnergySystems.s_control)],
        [97, ("TST_BUS_03", EnergySystems.s_produce)],
        [96, ("TST_BUS_02", EnergySystems.s_produce)],
        [95, ("TST_BUS_01", EnergySystems.s_produce)],
    ]
    Resie.reorder_for_input_priorities(steps, systems, by_function)
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
    systems_config = Dict{String,Any}(
        "TST_BUS_01" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "production_refs" => [
                "TST_BUS_04",
                "TST_BUS_03",
                "TST_BUS_02",
            ],
            "connection_matrix" => Dict{String, Any}(
                "input_order" => [],
                "output_order" => [
                    "TST_BUS_04",
                    "TST_BUS_03",
                    "TST_BUS_02",
                ],
            )
        ),
        "TST_BUS_02" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "production_refs" => [],
        ),
        "TST_BUS_03" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "production_refs" => [
                "TST_BUS_06",
                "TST_BUS_05",
            ],
            "connection_matrix" => Dict{String, Any}(
                "input_order" => [
                    "TST_BUS_01",
                ],
                "output_order" => [
                    "TST_BUS_06",
                    "TST_BUS_05",
                ],
            )
        ),
        "TST_BUS_04" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "production_refs" => [],
        ),
        "TST_BUS_05" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "production_refs" => [],
        ),
        "TST_BUS_06" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "production_refs" => [],
        ),
    )
    systems = Resie.load_systems(systems_config)
    by_function = Resie.categorize_by_function(systems)
    return systems, by_function
end

function test_busses_distribution_no_change()
    systems, by_function = test_data_busses_distribute()
    steps = [
        [100, ("TST_BUS_04", EnergySystems.s_distribute)],
        [99, ("TST_BUS_06", EnergySystems.s_distribute)],
        [98, ("TST_BUS_05", EnergySystems.s_distribute)],
        [97, ("TST_BUS_03", EnergySystems.s_distribute)],
        [96, ("TST_BUS_02", EnergySystems.s_distribute)],
        [95, ("TST_BUS_01", EnergySystems.s_distribute)],
    ]
    expected = [
        [100, ("TST_BUS_04", EnergySystems.s_distribute)],
        [99, ("TST_BUS_06", EnergySystems.s_distribute)],
        [98, ("TST_BUS_05", EnergySystems.s_distribute)],
        [97, ("TST_BUS_03", EnergySystems.s_distribute)],
        [96, ("TST_BUS_02", EnergySystems.s_distribute)],
        [95, ("TST_BUS_01", EnergySystems.s_distribute)],
    ]
    Resie.reorder_distribution_of_busses(steps, systems, by_function)
    @test pwc_steps_astr(expected, steps) == ""
end

@testset "busses_distribution_no_change" begin
    test_busses_distribution_no_change()
end

function test_busses_distribution_reorder_steps()
    systems, by_function = test_data_busses_distribute()
    steps = [
        [100, ("TST_BUS_01", EnergySystems.s_distribute)],
        [99, ("TST_BUS_02", EnergySystems.s_distribute)],
        [98, ("TST_BUS_03", EnergySystems.s_distribute)],
        [97, ("TST_BUS_04", EnergySystems.s_distribute)],
        [96, ("TST_BUS_05", EnergySystems.s_distribute)],
        [95, ("TST_BUS_06", EnergySystems.s_distribute)],
    ]
    expected = [
        [100, ("TST_BUS_01", EnergySystems.s_distribute)],
        [101, ("TST_BUS_02", EnergySystems.s_distribute)],
        [102, ("TST_BUS_03", EnergySystems.s_distribute)],
        [105, ("TST_BUS_04", EnergySystems.s_distribute)],
        [103, ("TST_BUS_05", EnergySystems.s_distribute)],
        [104, ("TST_BUS_06", EnergySystems.s_distribute)],
    ]
    Resie.reorder_distribution_of_busses(steps, systems, by_function)
    @test pwc_steps_astr(expected, steps) == ""
end

@testset "busses_distribution_reorder_steps" begin
    test_busses_distribution_reorder_steps()
end