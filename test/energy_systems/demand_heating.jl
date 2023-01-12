watt_to_wh = function (watts :: Float64)
    watts * Float64(parameters["time_step_seconds"]) / 3600.0
end

@testset "demand_heating_temperature_values" begin
    systems_config = Dict{String, Any}(
        "TST_DEM_01" => Dict{String, Any}(
            "type" => "Demand",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "production_refs" => [],
            "energy_profile_file_path" => "../profiles/tests/demand_heating_energy.prf",
            "temperature_profile_file_path" => "../profiles/tests/demand_heating_temperature.prf",
            "scale" => 1000
        ),
    )
    systems = Bran.load_systems(systems_config)
    demand = systems["TST_DEM_01"]

    simulation_parameters = Dict{String, Any}(
        "time_step_seconds" => 900,
        "time" => 0,
    )

    @test demand.last_load == 0.0
    @test demand.last_temperature == 0.0

    Bran.EnergySystems.control(demand, systems, simulation_parameters)

    @test demand.last_load == 75.0
    @test demand.last_temperature == 55.0

    Bran.EnergySystems.produce(demand, simulation_parameters, watt_to_wh)

    @test demand.input_interfaces[m_h_w_ht1].balance == -75.0
    @test demand.input_interfaces[m_h_w_ht1].temperature == 55.0
end