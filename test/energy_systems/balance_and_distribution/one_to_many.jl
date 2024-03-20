using Debugger
using Test
using Resie
using Resie.EnergySystems
using Resie.Profiles

EnergySystems.set_timestep(900)

function test_many_to_one()
    components_config = Dict{String,Any}(
        "TST_SRC_01" => Dict{String,Any}(
            "type" => "BoundedSupply",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "output_refs" => ["TST_BUS_TH_01"],
            "constant_power" => 4000,
            "constant_temperature" => 55,
        ),
        "TST_DEM_01" => Dict{String,Any}(
            "type" => "Demand",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "output_refs" => [],
            "constant_demand" => 2000,
            "constant_temperature" => 54,
        ),
        "TST_DEM_02" => Dict{String,Any}(
            "type" => "Demand",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "output_refs" => [],
            "constant_demand" => 2000,
            "constant_temperature" => 53,
        ),
        "TST_BUS_TH_01" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "connections" => Dict{String,Any}(
                "input_order" => [
                    "TST_SRC_01",
                ],
                "output_order" => [
                    "TST_DEM_01",
                    "TST_DEM_02",
                ],
                "energy_flow" => [
                    [1, 1],
                ],
            ),
        ),
    )

    simulation_parameters = Dict{String,Any}(
        "time_step_seconds" => 900,
        "time" => 0,
        "epsilon" => 1e-9
    )

    components = Resie.load_components(components_config, simulation_parameters)
    source = components["TST_SRC_01"]
    demand_1 = components["TST_DEM_01"]
    demand_2 = components["TST_DEM_02"]
    bus = components["TST_BUS_TH_01"]


    EnergySystems.reset(demand_1)
    EnergySystems.reset(demand_2)
    EnergySystems.reset(source)
    EnergySystems.reset(bus)

    @test demand_1.demand == 0.0
    @test demand_1.temperature === nothing

    @test demand_2.demand == 0.0
    @test demand_2.temperature === nothing

    EnergySystems.control(demand_1, components, simulation_parameters)
    EnergySystems.control(demand_2, components, simulation_parameters)
    EnergySystems.control(source, components, simulation_parameters)
    EnergySystems.control(bus, components, simulation_parameters)

    @test demand_1.demand == 500.0
    @test demand_1.temperature == 54.0

    @test demand_2.demand == 500.0
    @test demand_2.temperature == 53.0

    exchanges = EnergySystems.balance_on(bus.input_interfaces[1], bus)
    @test EnergySystems.balance(exchanges) ≈ 0.0            # is always zero in exchange of bus
    @test EnergySystems.energy_potential(exchanges) ≈ -1000.0 
    @test EnergySystems.temp_min_highest(exchanges) === 54.0

    exchanges = EnergySystems.balance_on(bus.output_interfaces[1], bus)
    @test EnergySystems.balance(exchanges) ≈ 0.0            # is always zero in exchange of bus
    @test EnergySystems.energy_potential(exchanges) ≈ 500.0 # already inner_distributed and fitting the demand
    @test EnergySystems.temp_min_highest(exchanges) === nothing
    @test EnergySystems.temp_max_highest(exchanges) === 55.0

    EnergySystems.process(demand_1, simulation_parameters)

    @test demand_1.input_interfaces[demand_1.medium].balance == -500.0
    @test demand_1.input_interfaces[demand_1.medium].temperature_min == 54.0

    exchanges = EnergySystems.balance_on(bus.output_interfaces[2], bus)
    @test EnergySystems.balance(exchanges) ≈ 0.0            # is always zero in exchange of bus
    @test EnergySystems.energy_potential(exchanges) ≈ 500.0 
    @test EnergySystems.temp_min_highest(exchanges) === nothing
    @test EnergySystems.temp_max_highest(exchanges) === 55.0

    EnergySystems.process(demand_2, simulation_parameters)

    @test demand_2.input_interfaces[demand_2.medium].balance == -500.0
    @test demand_2.input_interfaces[demand_2.medium].temperature_min == 53.0

    exchanges = EnergySystems.balance_on(bus.input_interfaces[1], bus)
    @test EnergySystems.balance(exchanges) ≈ 0.0            # is always zero in exchange of bus
    @test EnergySystems.energy_potential(exchanges) ≈ -1000.0 
    @test EnergySystems.temp_min_highest(exchanges) === 54.0

    EnergySystems.process(bus, simulation_parameters)
    EnergySystems.process(source, simulation_parameters)

    @test source.output_interfaces[source.medium].balance == 1000.0
    @test source.output_interfaces[source.medium].temperature_max == 55.0

    blnc = EnergySystems.balance(bus)
    @test blnc == 0.0

    EnergySystems.distribute!(bus)

    @test demand_1.input_interfaces[demand_1.medium].balance == 0.0
    @test demand_1.input_interfaces[demand_1.medium].temperature_min === 54.0

    @test demand_2.input_interfaces[demand_2.medium].balance == 0.0
    @test demand_2.input_interfaces[demand_2.medium].temperature_min === 53.0

    @test source.output_interfaces[source.medium].balance == 0.0
    @test source.output_interfaces[source.medium].temperature_max == 55.0

end

@testset "many_to_one" begin
    test_many_to_one()
end