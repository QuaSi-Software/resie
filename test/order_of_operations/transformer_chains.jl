using Debugger
using Test
using Resie
using Resie.EnergySystems

function test_distance_from_sink()
    components_config = Dict{String,Any}(
        "TST_BUS_01" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "output_refs" => ["TST_BUS_02", "TST_BUS_03"],
        ),
        "TST_BUS_02" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "output_refs" => ["TST_BUS_04", "TST_BUS_05"],
        ),
        "TST_BUS_03" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "output_refs" => [],
        ),
        "TST_BUS_04" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "output_refs" => [],
        ),
        "TST_BUS_05" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "output_refs" => ["TST_BUS_06"],
        ),
        "TST_BUS_06" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "output_refs" => [],
        ),
    )
    components = Resie.load_components(components_config)

    @test Resie.distance_to_sink(components["TST_BUS_01"], EnergySystems.sf_bus) == 3
    @test Resie.distance_to_sink(components["TST_BUS_02"], EnergySystems.sf_bus) == 2
    @test Resie.distance_to_sink(components["TST_BUS_03"], EnergySystems.sf_bus) == 0
    @test Resie.distance_to_sink(components["TST_BUS_04"], EnergySystems.sf_bus) == 0
    @test Resie.distance_to_sink(components["TST_BUS_05"], EnergySystems.sf_bus) == 1
    @test Resie.distance_to_sink(components["TST_BUS_06"], EnergySystems.sf_bus) == 0
end

@testset "distance_from_sink" begin
    test_distance_from_sink()
end

function test_iteration_order()
    components_config = Dict{String,Any}(
        "TST_BUS_01" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "output_refs" => ["TST_BUS_02", "TST_BUS_03"],
        ),
        "TST_BUS_02" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "output_refs" => ["TST_BUS_04", "TST_BUS_05"],
        ),
        "TST_BUS_03" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "output_refs" => [],
        ),
        "TST_BUS_04" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "output_refs" => [],
        ),
        "TST_BUS_05" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "output_refs" => ["TST_BUS_06"],
        ),
        "TST_BUS_06" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "output_refs" => [],
        ),
    )
    components = Resie.load_components(components_config)

    chain = [
        components["TST_BUS_01"],
        components["TST_BUS_02"],
        components["TST_BUS_03"],
        components["TST_BUS_04"],
        components["TST_BUS_05"],
        components["TST_BUS_06"],
    ]
    expected = [
        components["TST_BUS_03"],
        components["TST_BUS_04"],
        components["TST_BUS_06"],
        components["TST_BUS_05"],
        components["TST_BUS_02"],
        components["TST_BUS_01"],
    ]
    calculated = Resie.iterate_chain(chain, EnergySystems.sf_bus, reverse=false)
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
            "control_refs" => [],
            "output_refs" => ["TST_BUS_02"],
        ),
        "TST_BUS_02" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_h_w_lt1",
            "control_refs" => [],
            "output_refs" => ["TST_HTP_01", "TST_DEM_01"],
        ),
        "TST_BUS_03" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "output_refs" => ["TST_BUS_04", "TST_BUS_05"],
        ),
        "TST_BUS_04" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "output_refs" => [],
        ),
        "TST_BUS_05" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "output_refs" => [],
        ),
        "TST_HTP_01" => Dict{String,Any}(
            "type" => "HeatPump",
            "control_refs" => [],
            "output_refs" => ["TST_BUS_03"],
            "power" => 12000,
            "fixed_cop" => 3.0
        ),
        "TST_DEM_01" => Dict{String,Any}(
            "type" => "Demand",
            "medium" => "m_h_w_lt1",
            "control_refs" => [],
            "output_refs" => [],
            "energy_profile_file_path" => "./profiles/tests/demand_heating_energy.prf",
            "temperature_profile_file_path" => "./profiles/tests/demand_heating_temperature.prf",
            "scale" => 1000
        ),
    )
    components = Resie.load_components(components_config)

    expected_1 = Set([
        components["TST_BUS_01"],
        components["TST_BUS_02"],
    ])
    expected_2 = Set([
        components["TST_BUS_03"],
        components["TST_BUS_04"],
        components["TST_BUS_05"],
    ])

    chains = Resie.find_chains([u for u in values(components)], EnergySystems.sf_bus)

    @test length(chains) == 2
    @test chains[1] == expected_1
    @test chains[2] == expected_2
end

@testset "find_chains" begin
    test_find_chains()
end

function test_find_indexes()
    # normal case
    steps = [
        [2, ("TST_BUS_02", EnergySystems.s_control)],
        [3, ("TST_BUS_03", EnergySystems.s_control)],
        [1, ("TST_BUS_01", EnergySystems.s_control)],
    ]
    own_idx, target_idx = Resie.find_indexes(
        steps,
        ("TST_BUS_02", EnergySystems.s_control),
        ("TST_BUS_01", EnergySystems.s_control)
    )
    @test own_idx == 1
    @test target_idx == 3

    # can't find either step due to wrong UAC in one case and wrong step in the other
    steps = [
        [2, ("TST_BUS_02", EnergySystems.s_control)],
        [3, ("TST_BUS_03", EnergySystems.s_control)],
        [1, ("TST_BUS_01", EnergySystems.s_control)],
    ]
    own_idx, target_idx = Resie.find_indexes(
        steps,
        ("TST_BUS_02", EnergySystems.s_process),
        ("TST_BSS_01", EnergySystems.s_control)
    )
    @test own_idx === nothing
    @test target_idx === nothing
end

@testset "find_indexes" begin
    test_find_indexes()
end

function test_place_one_higher()
    # normal test case
    steps = [
        [2, ("TST_BUS_02", EnergySystems.s_control)],
        [3, ("TST_BUS_03", EnergySystems.s_control)],
        [1, ("TST_BUS_01", EnergySystems.s_control)],
    ]
    Resie.place_one_higher!(
        steps,
        ("TST_BUS_02", EnergySystems.s_control),
        ("TST_BUS_01", EnergySystems.s_control)
    )
    @test steps[1] == [2, ("TST_BUS_02", EnergySystems.s_control)]
    @test steps[2] == [4, ("TST_BUS_03", EnergySystems.s_control)]
    @test steps[3] == [3, ("TST_BUS_01", EnergySystems.s_control)]

    # priority already is higher and force was not used => no change
    steps = [
        [2, ("TST_BUS_02", EnergySystems.s_control)],
        [3, ("TST_BUS_03", EnergySystems.s_control)],
        [4, ("TST_BUS_01", EnergySystems.s_control)],
    ]
    Resie.place_one_higher!(
        steps,
        ("TST_BUS_02", EnergySystems.s_control),
        ("TST_BUS_01", EnergySystems.s_control)
    )
    @test steps[1] == [2, ("TST_BUS_02", EnergySystems.s_control)]
    @test steps[2] == [3, ("TST_BUS_03", EnergySystems.s_control)]
    @test steps[3] == [4, ("TST_BUS_01", EnergySystems.s_control)]

    # priority already is higher and force was used => works as normal
    steps = [
        [2, ("TST_BUS_02", EnergySystems.s_control)],
        [3, ("TST_BUS_03", EnergySystems.s_control)],
        [4, ("TST_BUS_01", EnergySystems.s_control)],
    ]
    Resie.place_one_higher!(
        steps,
        ("TST_BUS_02", EnergySystems.s_control),
        ("TST_BUS_01", EnergySystems.s_control),
        force=true
    )
    @test steps[1] == [2, ("TST_BUS_02", EnergySystems.s_control)]
    @test steps[2] == [4, ("TST_BUS_03", EnergySystems.s_control)]
    @test steps[3] == [3, ("TST_BUS_01", EnergySystems.s_control)]
end

@testset "place_one_higher" begin
    test_place_one_higher()
end

function test_place_one_lower()
    # normal test case
    steps = [
        [2, ("TST_BUS_02", EnergySystems.s_control)],
        [3, ("TST_BUS_03", EnergySystems.s_control)],
        [1, ("TST_BUS_01", EnergySystems.s_control)],
    ]
    Resie.place_one_lower!(
        steps,
        ("TST_BUS_02", EnergySystems.s_control),
        ("TST_BUS_03", EnergySystems.s_control)
    )
    @test steps[1] == [2, ("TST_BUS_02", EnergySystems.s_control)]
    @test steps[2] == [1, ("TST_BUS_03", EnergySystems.s_control)]
    @test steps[3] == [0, ("TST_BUS_01", EnergySystems.s_control)]

    # priority already is lower and force was not used => no change
    steps = [
        [2, ("TST_BUS_02", EnergySystems.s_control)],
        [3, ("TST_BUS_03", EnergySystems.s_control)],
        [1, ("TST_BUS_01", EnergySystems.s_control)],
    ]
    Resie.place_one_lower!(
        steps,
        ("TST_BUS_03", EnergySystems.s_control),
        ("TST_BUS_01", EnergySystems.s_control)
    )
    @test steps[1] == [2, ("TST_BUS_02", EnergySystems.s_control)]
    @test steps[2] == [3, ("TST_BUS_03", EnergySystems.s_control)]
    @test steps[3] == [1, ("TST_BUS_01", EnergySystems.s_control)]

    # priority already is lower and force was used => works as normal
    steps = [
        [2, ("TST_BUS_02", EnergySystems.s_control)],
        [3, ("TST_BUS_03", EnergySystems.s_control)],
        [1, ("TST_BUS_01", EnergySystems.s_control)],
    ]
    Resie.place_one_lower!(
        steps,
        ("TST_BUS_03", EnergySystems.s_control),
        ("TST_BUS_01", EnergySystems.s_control),
        force=true
    )
    @test steps[1] == [1, ("TST_BUS_02", EnergySystems.s_control)]
    @test steps[2] == [3, ("TST_BUS_03", EnergySystems.s_control)]
    @test steps[3] == [2, ("TST_BUS_01", EnergySystems.s_control)]
end

@testset "place_one_lower" begin
    test_place_one_lower()
end