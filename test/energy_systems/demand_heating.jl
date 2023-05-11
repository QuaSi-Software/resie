using Debugger
using Test
using Resie
using Resie.EnergySystems
using Resie.Profiles

EnergySystems.set_timestep(900)

function test_demand_heating_temperature_values()
    components_config = Dict{String,Any}(
        "TST_GRI_01" => Dict{String,Any}(
            "type" => "GridConnection",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "output_refs" => ["TST_DEM_01"],
            "is_source" => true,
        ),
        "TST_DEM_01" => Dict{String,Any}(
            "type" => "Demand",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "output_refs" => [],
            "energy_profile_file_path" => "./profiles/tests/demand_heating_energy.prf",
            "temperature_profile_file_path" => "./profiles/tests/demand_heating_temperature.prf",
            "scale" => 1000
        ),
    )
    components = Resie.load_components(components_config)
    demand = components["TST_DEM_01"]

    simulation_parameters = Dict{String,Any}(
        "time_step_seconds" => 900,
        "time" => 0,
    )

    EnergySystems.reset(demand)

    @test demand.load == 0.0
    @test demand.temperature === nothing

    EnergySystems.control(demand, components, simulation_parameters)

    @test demand.load == 75.0
    @test demand.temperature == 55.0

    EnergySystems.process(demand, simulation_parameters)

    @test demand.input_interfaces[demand.medium].balance == -75.0
    @test demand.input_interfaces[demand.medium].temperature == 55.0
end

@testset "demand_heating_temperature_values" begin
    test_demand_heating_temperature_values()
end