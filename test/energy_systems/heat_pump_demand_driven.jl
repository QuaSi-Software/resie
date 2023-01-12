@testset "heat_pump_demand_driven" begin
    systems_config = Dict{String, Any}(
        "TST_DEM_01" => Dict{String, Any}(
            "type" => "Demand",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "production_refs" => [],
            "energy_profile_file_path" => "./profiles/tests/demand_heating_energy.prf",
            "temperature_profile_file_path" => "./profiles/tests/demand_heating_temperature.prf",
            "scale" => 18000
        ),
        "TST_GRI_01" => Dict{String, Any}(
            "type" => "GridConnection",
            "medium" => "m_h_w_lt1",
            "control_refs" => [],
            "production_refs" => ["TST_HP_01"],
            "is_source" => true,
        ),
        "TST_GRI_02" => Dict{String, Any}(
            "type" => "GridConnection",
            "medium" => "m_e_ac_230v",
            "control_refs" => [],
            "production_refs" => ["TST_HP_01"],
            "is_source" => true,
        ),
        "TST_HP_01" => Dict{String, Any}(
            "type" => "HeatPump",
            "control_refs" => ["TST_DEM_01"],
            "production_refs" => ["TST_DEM_01"],
            "strategy" => Dict{String, Any}(
                "name" => "demand_driven",
            ),
            "power" => 20000,
            "cop" => 3.0
        ),
    )
    systems = Bran.load_systems(systems_config)
    heat_pump = systems["TST_HP_01"]

    @test heat_pump.controller.state_machine.state == 1
end