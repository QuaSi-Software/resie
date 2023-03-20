using Debugger
using Test
using Resie
using Resie.EnergySystems

function test_distance_from_sink()
    systems_config = Dict{String,Any}(
        "TST_BUS_01" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "production_refs" => ["TST_BUS_02", "TST_BUS_03"],
        ),
        "TST_BUS_02" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "production_refs" => ["TST_BUS_04", "TST_BUS_05"],
        ),
        "TST_BUS_03" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "production_refs" => [],
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
            "production_refs" => ["TST_BUS_06"],
        ),
        "TST_BUS_06" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "production_refs" => [],
        ),
    )
    systems = Resie.load_systems(systems_config)

    @test Resie.distance_to_sink(systems["TST_BUS_01"], EnergySystems.sf_bus) == 3
    @test Resie.distance_to_sink(systems["TST_BUS_02"], EnergySystems.sf_bus) == 2
    @test Resie.distance_to_sink(systems["TST_BUS_03"], EnergySystems.sf_bus) == 0
    @test Resie.distance_to_sink(systems["TST_BUS_04"], EnergySystems.sf_bus) == 0
    @test Resie.distance_to_sink(systems["TST_BUS_05"], EnergySystems.sf_bus) == 1
    @test Resie.distance_to_sink(systems["TST_BUS_06"], EnergySystems.sf_bus) == 0
end

@testset "distance_from_sink" begin
    test_distance_from_sink()
end

function test_iteration_order()
    systems_config = Dict{String,Any}(
        "TST_BUS_01" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "production_refs" => ["TST_BUS_02", "TST_BUS_03"],
        ),
        "TST_BUS_02" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "production_refs" => ["TST_BUS_04", "TST_BUS_05"],
        ),
        "TST_BUS_03" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "production_refs" => [],
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
            "production_refs" => ["TST_BUS_06"],
        ),
        "TST_BUS_06" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "production_refs" => [],
        ),
    )
    systems = Resie.load_systems(systems_config)

    chain = [
        systems["TST_BUS_01"],
        systems["TST_BUS_02"],
        systems["TST_BUS_03"],
        systems["TST_BUS_04"],
        systems["TST_BUS_05"],
        systems["TST_BUS_06"],
    ]
    expected = [
        systems["TST_BUS_03"],
        systems["TST_BUS_04"],
        systems["TST_BUS_06"],
        systems["TST_BUS_05"],
        systems["TST_BUS_02"],
        systems["TST_BUS_01"],
    ]
    calculated = Resie.iterate_chain(chain, EnergySystems.sf_bus, false)
    for i in 1:6
        @test "#$i: " * calculated[i].uac == "#$i: " * expected[i].uac
    end
end

@testset "iteration_order" begin
    test_iteration_order()
end

function test_find_chains()
    systems_config = Dict{String,Any}(
        "TST_BUS_01" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_h_w_lt1",
            "control_refs" => [],
            "production_refs" => ["TST_BUS_02"],
        ),
        "TST_BUS_02" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_h_w_lt1",
            "control_refs" => [],
            "production_refs" => ["TST_HTP_01", "TST_DEM_01"],
        ),
        "TST_BUS_03" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "production_refs" => ["TST_BUS_04", "TST_BUS_05"],
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
        "TST_HTP_01" => Dict{String,Any}(
            "type" => "HeatPump",
            "control_refs" => [],
            "production_refs" => ["TST_BUS_03"],
            "power" => 12000,
            "fixed_cop" => 3.0
        ),
        "TST_DEM_01" => Dict{String,Any}(
            "type" => "Demand",
            "medium" => "m_h_w_lt1",
            "control_refs" => [],
            "production_refs" => [],
            "energy_profile_file_path" => "./profiles/tests/demand_heating_energy.prf",
            "temperature_profile_file_path" => "./profiles/tests/demand_heating_temperature.prf",
            "scale" => 1000
        ),
    )
    systems = Resie.load_systems(systems_config)

    expected_1 = Set([
        systems["TST_BUS_01"],
        systems["TST_BUS_02"],
    ])
    expected_2 = Set([
        systems["TST_BUS_03"],
        systems["TST_BUS_04"],
        systems["TST_BUS_05"],
    ])

    chains = Resie.find_chains([u for u in values(systems)], EnergySystems.sf_bus)

    @test length(chains) == 2
    @test chains[1] == expected_1
    @test chains[2] == expected_2
end

@testset "find_chains" begin
    test_find_chains()
end