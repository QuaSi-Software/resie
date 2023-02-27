using Debugger
using Test
using Resie
using Resie.EnergySystems
using Resie.Profiles

watt_to_wh = function (watts :: Float64)
    watts * 900 / 3600.0
end

function test_demand_heating_temperature_values()
    systems_config = Dict{String, Any}(
        "TST_GRI_01" => Dict{String, Any}(
            "type" => "GridConnection",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "production_refs" => ["TST_DEM_01"],
            "is_source" => true,
        ),
        "TST_DEM_01" => Dict{String, Any}(
            "type" => "Demand",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "production_refs" => [],
            "energy_profile_file_path" => "./profiles/tests/demand_heating_energy.prf",
            "temperature_profile_file_path" => "./profiles/tests/demand_heating_temperature.prf",
            "scale" => 1000
        ),
    )
    _ = Resie.load_medien( Array{Any}(undef,0) )
    systems = Resie.load_systems(systems_config)
    demand = systems["TST_DEM_01"]

    simulation_parameters = Dict{String, Any}(
        "time_step_seconds" => 900,
        "time" => 0,
    )

    EnergySystems.reset(demand)

    @test demand.load == 0.0
    @test demand.temperature === nothing

    EnergySystems.control(demand, systems, simulation_parameters)

    @test demand.load == 75.0
    @test demand.temperature == 55.0

    EnergySystems.produce(demand, simulation_parameters, watt_to_wh)

    @test demand.input_interfaces[EnergySystems.m_h_w_ht1].balance == -75.0
    @test demand.input_interfaces[EnergySystems.m_h_w_ht1].temperature == 55.0
end

@testset "demand_heating_temperature_values" begin
    test_demand_heating_temperature_values()
end