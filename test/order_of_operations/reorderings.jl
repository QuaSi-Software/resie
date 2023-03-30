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

function test_data_storage_loading()
    systems_config = Dict{String,Any}(
        "TST_BUS_01" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "production_refs" => [
                "TST_BUS_02",
                "TST_BFT_01",
                "TST_BUS_03",
            ],
            "connection_matrix" => Dict{String, Any}(
                "input_order" => [],
                "output_order" => [
                    "TST_BUS_02",
                    "TST_BFT_01",
                    "TST_BUS_03",
                ],
            )
        ),
        "TST_BUS_02" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "production_refs" => ["TST_BFT_02"],
        ),
        "TST_BUS_03" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "production_refs" => ["TST_BFT_03"],
        ),
        "TST_BFT_01" => Dict{String,Any}(
            "type" => "BufferTank",
            "control_refs" => [],
            "production_refs" => [
                "TST_BUS_01"
            ],
            "capacity" => 40000,
            "load" => 0
        ),
        "TST_BFT_02" => Dict{String,Any}(
            "type" => "BufferTank",
            "control_refs" => [],
            "production_refs" => [
                "TST_BUS_02"
            ],
            "capacity" => 40000,
            "load" => 0
        ),
        "TST_BFT_03" => Dict{String,Any}(
            "type" => "BufferTank",
            "control_refs" => [],
            "production_refs" => [
                "TST_BUS_03"
            ],
            "capacity" => 40000,
            "load" => 0
        ),
    )
    systems = Resie.load_systems(systems_config)
    by_function = Resie.categorize_by_function(systems)
    return systems, by_function
end

function test_find_storages_ordered()
    systems, _ = test_data_storage_loading()

    expected = [
        systems["TST_BFT_02"],
        systems["TST_BFT_01"],
        systems["TST_BFT_03"],
    ]
    result, limits = Resie.find_storages_ordered(systems["TST_BUS_01"], systems, nothing)
    @test pwc_units_astr(expected, result) == ""

    expected = [
        systems["TST_BFT_03"],
        systems["TST_BFT_01"],
        systems["TST_BFT_02"],
    ]
    result, limits = Resie.find_storages_ordered(systems["TST_BUS_01"], systems, nothing, reverse=true)
    @test pwc_units_astr(expected, result) == ""
end

@testset "find_storages_ordered" begin
    test_find_storages_ordered()
end

function test_storage_loading_no_change()
    systems, by_function = test_data_storage_loading()
    steps = [
        [100, ("TST_BFT_02", EnergySystems.s_produce)],
        [99, ("TST_BFT_01", EnergySystems.s_produce)],
        [98, ("TST_BFT_03", EnergySystems.s_produce)],
        [97, ("TST_BFT_02", EnergySystems.s_load)],
        [96, ("TST_BFT_01", EnergySystems.s_load)],
        [95, ("TST_BFT_03", EnergySystems.s_load)],
    ]
    expected = [
        [100, ("TST_BFT_02", EnergySystems.s_produce)],
        [99, ("TST_BFT_01", EnergySystems.s_produce)],
        [98, ("TST_BFT_03", EnergySystems.s_produce)],
        [97, ("TST_BFT_02", EnergySystems.s_load)],
        [96, ("TST_BFT_01", EnergySystems.s_load)],
        [95, ("TST_BFT_03", EnergySystems.s_load)],
    ]
    Resie.reorder_storage_loading(steps, systems, by_function)
    @test pwc_steps_astr(expected, steps) == ""
end

@testset "storage_loading_no_change" begin
    test_storage_loading_no_change()
end

function test_storage_loading_reorder_steps()
    systems, by_function = test_data_storage_loading()
    steps = [
        [100, ("TST_BFT_01", EnergySystems.s_produce)],
        [99, ("TST_BFT_02", EnergySystems.s_produce)],
        [98, ("TST_BFT_03", EnergySystems.s_produce)],
        [97, ("TST_BFT_01", EnergySystems.s_load)],
        [96, ("TST_BFT_02", EnergySystems.s_load)],
        [95, ("TST_BFT_03", EnergySystems.s_load)],
    ]
    expected = [
        [98, ("TST_BFT_01", EnergySystems.s_produce)],
        [99, ("TST_BFT_02", EnergySystems.s_produce)],
        [97, ("TST_BFT_03", EnergySystems.s_produce)],
        [94, ("TST_BFT_01", EnergySystems.s_load)],
        [95, ("TST_BFT_02", EnergySystems.s_load)],
        [93, ("TST_BFT_03", EnergySystems.s_load)],
    ]
    Resie.reorder_storage_loading(steps, systems, by_function)
    @test pwc_steps_astr(expected, steps) == ""
end

@testset "storage_loading_reorder_steps" begin
    test_storage_loading_reorder_steps()
end

function test_data_control_dependencies()
    # this is a strange example as it would not normally be used. this is only used here
    # to test specific characteristics of reorder_for_control_dependencies()
    systems_config = Dict{String,Any}(
        "TST_BUS_01" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "production_refs" => [
                "TST_HP_01",
                "TST_HP_02",
                "TST_BFT_01",
            ],
        ),
        "TST_HP_02" => Dict{String,Any}(
            "type" => "HeatPump",
            "control_refs" => ["TST_HP_01"],
            "production_refs" => [],
            "strategy" => Dict{String,Any}(
                "name" => "demand_driven",
            ),
            "power" => 12000,
            "fixed_cop" => 3.0
        ),
        "TST_HP_01" => Dict{String,Any}(
            "type" => "HeatPump",
            "control_refs" => ["TST_BFT_01"],
            "production_refs" => [],
            "strategy" => Dict{String,Any}(
                "name" => "storage_driven",
                "low_threshold" => 0.2,
                "high_threshold" => 0.8,
            ),
            "power" => 12000,
            "fixed_cop" => 3.0
        ),
        "TST_BFT_01" => Dict{String,Any}(
            "type" => "BufferTank",
            "control_refs" => [],
            "production_refs" => [
                "TST_BUS_01"
            ],
            "capacity" => 40000,
            "load" => 0
        ),
    )
    systems = Resie.load_systems(systems_config)
    by_function = Resie.categorize_by_function(systems)
    return systems, by_function
end

function test_control_dependencies_no_change()
    systems, by_function = test_data_control_dependencies()
    steps = [
        [102, ("TST_BUS_01", EnergySystems.s_control)],
        [100, ("TST_HP_02", EnergySystems.s_control)],
        [101, ("TST_HP_01", EnergySystems.s_control)],
        [98, ("TST_BFT_01", EnergySystems.s_control)],
        [97, ("TST_BUS_01", EnergySystems.s_produce)],
        [95, ("TST_HP_02", EnergySystems.s_produce)],
        [96, ("TST_HP_01", EnergySystems.s_produce)],
        [93, ("TST_BFT_01", EnergySystems.s_produce)],
    ]
    expected = [
        [102, ("TST_BUS_01", EnergySystems.s_control)],
        [100, ("TST_HP_02", EnergySystems.s_control)],
        [101, ("TST_HP_01", EnergySystems.s_control)],
        [98, ("TST_BFT_01", EnergySystems.s_control)],
        [97, ("TST_BUS_01", EnergySystems.s_produce)],
        [95, ("TST_HP_02", EnergySystems.s_produce)],
        [96, ("TST_HP_01", EnergySystems.s_produce)],
        [93, ("TST_BFT_01", EnergySystems.s_produce)],
    ]
    Resie.reorder_for_control_dependencies(steps, systems, by_function)
    @test pwc_steps_astr(expected, steps) == ""
end

@testset "control_dependencies_no_change" begin
    test_control_dependencies_no_change()
end

function test_control_dependencies_reorder_steps()
    systems, by_function = test_data_control_dependencies()
    steps = [
        [102, ("TST_BUS_01", EnergySystems.s_control)],
        [100, ("TST_HP_02", EnergySystems.s_control)],
        [101, ("TST_HP_01", EnergySystems.s_control)],
        [98, ("TST_BFT_01", EnergySystems.s_control)],
        [97, ("TST_BUS_01", EnergySystems.s_produce)],
        [95, ("TST_HP_02", EnergySystems.s_produce)],
        [96, ("TST_HP_01", EnergySystems.s_produce)],
        [93, ("TST_BFT_01", EnergySystems.s_produce)],
    ]
    expected = [
        [102, ("TST_BUS_01", EnergySystems.s_control)],
        [100, ("TST_HP_02", EnergySystems.s_control)],
        [101, ("TST_HP_01", EnergySystems.s_control)],
        [98, ("TST_BFT_01", EnergySystems.s_control)],
        [97, ("TST_BUS_01", EnergySystems.s_produce)],
        [95, ("TST_HP_02", EnergySystems.s_produce)],
        [96, ("TST_HP_01", EnergySystems.s_produce)],
        [93, ("TST_BFT_01", EnergySystems.s_produce)],
    ]
    Resie.reorder_for_control_dependencies(steps, systems, by_function)
    @test pwc_steps_astr(expected, steps) == ""
end

@testset "control_dependencies_reorder_steps" begin
    test_control_dependencies_reorder_steps()
end
