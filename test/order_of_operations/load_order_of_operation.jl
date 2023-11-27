using Debugger
using Test
using Resie
using Resie.EnergySystems

include("../test_util.jl")

function load_order_of_operation()
    components_config = Dict{String,Any}(
        "TST_DEM_01" => Dict{String,Any}(
            "type"=> "Demand",
            "medium"=> "m_h_w_ht1",
            "control_refs" => [],
            "output_refs"=> [],
            "energy_profile_file_path" => "./profiles/tests/demand_heating_energy.prf",
            "temperature_profile_file_path" => "./profiles/tests/demand_heating_temperature.prf",
            "scale" => 1500
        ),
        "TST_SRC_01" => Dict{String,Any}(
            "type" => "BoundedSupply",
            "medium" => "m_h_w_lt1",
            "control_refs" => [],
            "output_refs" => [
                "TST_HP_01"
            ],
            "max_power_profile_file_path" => "./profiles/tests/source_heat_max_power.prf",
            "temperature_profile_file_path" => "./profiles/tests/source_heat_temperature.prf",
            "scale" => 6000
        ),
        "TST_GRI_01" => Dict{String,Any}(
            "type" => "GridConnection",
            "medium" => "m_e_ac_230v",
            "control_refs" => [],
            "output_refs" => [
                "TST_HP_01"
            ],
            "is_source" => true
        ),
        "TST_HP_01" => Dict{String,Any}(
            "type" => "HeatPump",
            "control_refs" => [
                "TST_DEM_01"
            ],
            "output_refs" => [
                "TST_DEM_01"
            ],
            "strategy" => Dict{String,Any}(
                "name" => "demand_driven"
            ),
            "power_th" => 9000,
            "fixed_cop" => 3.0,
            "min_power_fraction" => 0.0
        )  
    )
    order_of_operation = [
        "TST_DEM_01 s_reset",
        "TST_SRC_01 s_reset",
        "TST_GRI_01 s_reset",
        "TST_DEM_01 s_control",
        "TST_HP_01 s_control",
        "TST_SRC_01 s_control",
        "TST_GRI_01 s_control",
        "TST_DEM_01 s_process",
        "TST_HP_01 s_process",
        "TST_SRC_01 s_process",
        "TST_GRI_01 s_process",
        "TST_HP_01 s_reset",
    ]

    expected_order_from_input_file = [
        ("TST_DEM_01", EnergySystems.s_reset),
        ("TST_SRC_01", EnergySystems.s_reset),
        ("TST_GRI_01", EnergySystems.s_reset),
        ("TST_DEM_01", EnergySystems.s_control),
        ("TST_HP_01", EnergySystems.s_control),
        ("TST_SRC_01", EnergySystems.s_control),
        ("TST_GRI_01", EnergySystems.s_control),
        ("TST_DEM_01", EnergySystems.s_process),
        ("TST_HP_01", EnergySystems.s_process),
        ("TST_SRC_01", EnergySystems.s_process),
        ("TST_GRI_01", EnergySystems.s_process),
        ("TST_HP_01", EnergySystems.s_reset)
    ]

    expected_order_calculated = [
        ("TST_DEM_01", EnergySystems.s_reset),
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
        ("TST_GRI_01", EnergySystems.s_process)
    ]

    components = Resie.load_components(components_config)

    # from input file
    order_from_input_file = Resie.load_order_of_operations(order_of_operation, components)
    @test pwc_steps_astr(expected_order_from_input_file, order_from_input_file) == ""

    # from ooo calculation
    order_calculated = Resie.calculate_order_of_operations(components)
    @test pwc_steps_astr(expected_order_calculated, order_calculated) == ""

end


@testset "load_order_of_operation" begin
    load_order_of_operation()
end