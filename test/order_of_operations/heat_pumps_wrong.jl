using Debugger
using Test
using Resie
using Resie.EnergySystems

include("../test_util.jl")

function test_ooo_for_heat_pumps_wrong()
    components_config = Dict{String,Any}(
        "TST_DEM_01" => Dict{String,Any}(
            "type" => "Demand",
            "medium" => "m_h_w_ht1",
            "output_refs" => [],
            "energy_profile_file_path" => "./profiles/tests/demand_heating_energy.prf",
            "temperature_profile_file_path" => "./profiles/tests/demand_heating_temperature.prf",
            "scale" => 1500,
        ),
        "TST_SRC_01" => Dict{String,Any}(
            "type" => "BoundedSupply",
            "medium" => "m_h_w_lt1",
            "output_refs" => ["TST_HP_01"],
            "max_power_profile_file_path" => "./profiles/tests/demand_heating_energy.prf",
            "temperature_profile_file_path" => "./profiles/tests/demand_heating_temperature.prf",
            "scale" => 6000,
        ),
        "TST_GRI_01" => Dict{String,Any}(
            "type" => "GridConnection",
            "medium" => "m_e_ac_230v",
            "output_refs" => ["TST_HP_01"],
            "is_source" => true,
        ),
        "TST_HP_01" => Dict{String,Any}(
            "type" => "HeatPump",
            "output_refs" => ["TST_DEM_01"],
            "power_th" => 12000,
            "cop_function" => "const:3.0",
        ),
    )

    expected = [("TST_DEM_01", EnergySystems.s_reset),
                ("TST_HP_01", EnergySystems.s_reset),
                ("TST_SRC_01", EnergySystems.s_reset),
                ("TST_GRI_01", EnergySystems.s_reset),
                ("TST_DEM_01", EnergySystems.s_control),
                ("TST_HP_01", EnergySystems.s_control),
                ("TST_SRC_01", EnergySystems.s_control),
                ("TST_GRI_01", EnergySystems.s_control),
                ("TST_DEM_01", EnergySystems.s_process),
                ("TST_HP_01", EnergySystems.s_process),
                ("TST_SRC_01", EnergySystems.s_process),
                ("TST_GRI_01", EnergySystems.s_process)]

    simulation_parameters = Dict{String,Any}(
        "time_step_seconds" => 900,
        "time" => 0,
        "epsilon" => 1e-9,
    )

    components = Resie.load_components(components_config, simulation_parameters)
    ooo = Resie.calculate_order_of_operations(components)
    @test pwc_steps_astr(expected, ooo) == ""
end

@testset "ooo_for_heat_pumps_wrong" begin
    test_ooo_for_heat_pumps_wrong()
end
