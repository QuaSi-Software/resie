using Debugger
using Test
using Resie
using Resie.EnergySystems
using Resie.Profiles

include("../../test_util.jl")

function test_one_to_one_via_bus()
    components_config = Dict{String,Any}(
        "TST_GRI_01" => Dict{String,Any}(
            "type" => "GridConnection",
            "medium" => "m_h_w_ht1",
            "output_refs" => ["TST_BUS_TH_01"],
            "is_source" => true,
        ),
        "TST_DEM_01" => Dict{String,Any}(
            "type" => "Demand",
            "medium" => "m_h_w_ht1",
            "output_refs" => [],
            "constant_demand" => 4000,
            "constant_temperature" => 55,
        ),
        "TST_BUS_TH_01" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_h_w_ht1",
            "connections" => Dict{String,Any}(
                "input_order" => ["TST_GRI_01"],
                "output_order" => ["TST_DEM_01"],
                "energy_flow" => [[1]],
            ),
        ),
    )

    simulation_parameters = get_default_sim_params()

    components = Resie.load_components(components_config, simulation_parameters)
    grid = components["TST_GRI_01"]
    demand = components["TST_DEM_01"]
    bus = components["TST_BUS_TH_01"]

    EnergySystems.reset(demand)
    EnergySystems.reset(grid)
    EnergySystems.reset(bus)

    @test demand.demand == 0.0
    @test demand.temperature === nothing

    EnergySystems.control(demand, components, simulation_parameters)
    EnergySystems.control(grid, components, simulation_parameters)
    EnergySystems.control(bus, components, simulation_parameters)

    @test demand.demand == 1000.0
    @test demand.temperature == 55.0

    EnergySystems.process(demand, simulation_parameters)

    @test demand.input_interfaces[demand.medium].balance == -1000.0
    @test demand.input_interfaces[demand.medium].temperature_min == 55.0

    EnergySystems.process(bus, simulation_parameters)
    EnergySystems.process(grid, simulation_parameters)

    @test grid.output_interfaces[grid.medium].balance == 1000.0
    @test grid.output_interfaces[grid.medium].temperature_max === nothing

    blnc = EnergySystems.balance(bus)
    @test blnc == 0.0

    EnergySystems.distribute!(bus)

    @test demand.input_interfaces[demand.medium].balance == 0.0
    @test demand.input_interfaces[demand.medium].temperature_min === 55.0

    @test grid.output_interfaces[grid.medium].balance == 0.0
    @test grid.output_interfaces[grid.medium].temperature_max === nothing
end

@testset "one_to_one_via_bus" begin
    test_one_to_one_via_bus()
end
