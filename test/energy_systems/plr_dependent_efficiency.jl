using Debugger
using Test
using Resie
using Resie.EnergySystems
using Resie.Profiles

include("../test_util.jl")

function get_demand_energy_system_config()
    return Dict{String,Any}(
        "TST_DEM_01" => Dict{String,Any}(
            "type" => "Demand",
            "medium" => "m_h_w_ht1",
            "output_refs" => [],
            "constant_demand" => 4000,
            "scale" => 1,
        ),
        "TST_GRI_01" => Dict{String,Any}(
            "type" => "GridConnection",
            "medium" => "m_c_g_natgas",
            "output_refs" => ["TST_GB_01"],
            "is_source" => true,
        ),
        "TST_GB_01" => Dict{String,Any}(
            "type" => "FuelBoiler",
            "m_fuel_in" => "m_c_g_natgas",
            "output_refs" => ["TST_DEM_01"],
            "power_th" => 4000,
            "linear_interface" => "fuel_in",
            "efficiency_fuel_in" => "const:1.0",
            "efficiency_heat_out" => "poly:-0.9117,1.8795,0.0322",
            "nr_discretization_steps" => 20,
        ),
    )
end

function test_efficiency_parsing()
    efficiency = EnergySystems.parse_efficiency_function("const:0.314")
    @test efficiency(0.0) ≈ 0.314
    @test efficiency(0.5) ≈ 0.314
    @test efficiency(1.0) ≈ 0.314

    efficiency = EnergySystems.parse_efficiency_function("poly:0.5,0.314")
    @test efficiency(0.0) ≈ 0.314
    @test efficiency(0.5) ≈ 0.564
    @test efficiency(1.0) ≈ 0.814

    efficiency = EnergySystems.parse_efficiency_function("poly:-0.9117,1.8795,0.0322")
    @test efficiency(0.0) ≈ 0.0322
    @test efficiency(0.5) ≈ 0.744025
    @test efficiency(1.0) ≈ 1.0

    efficiency = EnergySystems.parse_efficiency_function("pwlin:0.5,0.9,1.0")
    @test efficiency(0.0) ≈ 0.5
    @test efficiency(0.25) ≈ 0.7
    @test efficiency(0.5) ≈ 0.9
    @test efficiency(0.8) ≈ 0.96
    @test efficiency(1.0) ≈ 1.0
end

@testset "test_efficiency_parsing" begin
    test_efficiency_parsing()
end

function test_cop_parsing()
    const_val, cop_func = EnergySystems.parse_cop_function("const:0.314")
    @test const_val == 0.314
    plr_func = cop_func(30.0, 50.0)
    @test plr_func(1.0) ≈ 0.314
    plr_func = cop_func(10.0, 70.0)
    @test plr_func(1.0) ≈ 0.314

    const_val, cop_func = EnergySystems.parse_cop_function("carnot:0.4")
    @test const_val === nothing
    plr_func = cop_func(30.0, 50.0)
    @test plr_func(1.0) ≈ 6.463
    plr_func = cop_func(10.0, 70.0)
    @test plr_func(1.0) ≈ 2.2876666666666

    field_def = "field:" *
                " 0, 0,10,20,30;" *
                " 0,30,20,10, 5;" *
                "10,30,30,20,10;" *
                "20,30,30,30,20"
    const_val, cop_func = EnergySystems.parse_cop_function(field_def)
    @test const_val === nothing
    @test cop_func(5.0, 5.0)(1.0) ≈ 30.0
    @test cop_func(10.0, 5.0)(1.0) ≈ 30.0
    @test cop_func(5.0, 10.0)(1.0) ≈ 25.0
    @test cop_func(12.5, 26.9)(1.0) ≈ 15.6

    error_occured = false
    try
        cop_func(9.0, 31.0)(1.0) # ≈ 5.0
    catch BoundsError
        error_occured = true
    end
    @test error_occured

    field_def = "field:" *
                " 0, 0,10,20,30;" *
                " 0,30,20,10, 5;" *
                "10,30,30,25,10;" *
                "20,30,30,30,20"
    const_val, cop_func = EnergySystems.parse_cop_function(field_def)
    @test const_val === nothing
    @test cop_func(12.5, 26.9)(1.0) ≈ 17.15
end

@testset "test_cop_parsing" begin
    test_cop_parsing()
end

function test_2dim_parsing()
    func_def = "poly-2:1.0,0.5,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0"
    func_2dim = EnergySystems.parse_2dim_function(func_def)
    @test func_2dim(12.0, 99.9) ≈ 7.0
    func_def = "poly-2:1.0,0.0,0.5,0.0,0.0,0.0,0.0,0.0,0.0,0.0"
    func_2dim = EnergySystems.parse_2dim_function(func_def)
    @test func_2dim(99.0, 16.0) ≈ 9.0
    func_def = "poly-2:2.0,0.5,0.5,0.0,0.1,0.0,0.0,0.0,0.0,0.0"
    func_2dim = EnergySystems.parse_2dim_function(func_def)
    @test func_2dim(10.0, 30.0) ≈ 52.0
    func_def = "poly-2:1.0,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1"
    func_2dim = EnergySystems.parse_2dim_function(func_def)
    @test func_2dim(2.0, 3.0) ≈ 1.0 + 0.2 + 0.3 + 0.4 + 0.9 + 0.6 + 0.8 + 1.2 + 1.8 + 2.7
end

@testset "test_2dim_parsing" begin
    test_2dim_parsing()
end

function test_bilinear_interpolate()
    #   x    1   2   3
    #  y   -----------
    #    1 | 0       1
    #    2 |     v
    #    3 | 1       3
    v = EnergySystems.bilinear_interpolate(1, 2, 3, 1, 2, 3, 0, 1, 1, 3)
    @test v ≈ 1.5

    #   x    1   2   3
    #  y   -----------
    #    1 | 0   v   1
    #    2 |
    #    3 | 1       3
    v = EnergySystems.bilinear_interpolate(1, 2, 3, 1, 1, 3, 0, 1, 1, 3)
    @test v ≈ 0.5

    #   x    1   2   3
    #  y   -----------
    #    1 | 0       3
    #    2 |         v
    #    3 | 1       3
    v = EnergySystems.bilinear_interpolate(1, 3, 3, 1, 2, 3, 0, 3, 1, 3)
    @test v ≈ 3.0

    #   x    1   2   3
    #  y   -----------
    #    1 | 3       1
    #    2 |         
    #    3 | 1       2
    #
    # five positions, four along each line of the box and one in the middle
    v = EnergySystems.bilinear_interpolate(1, 2, 3, 1, 1, 3, 3, 1, 1, 2)
    @test v ≈ 2.0
    v = EnergySystems.bilinear_interpolate(1, 3, 3, 1, 2, 3, 3, 1, 1, 2)
    @test v ≈ 1.5
    v = EnergySystems.bilinear_interpolate(1, 2, 3, 1, 3, 3, 3, 1, 1, 2)
    @test v ≈ 1.5
    v = EnergySystems.bilinear_interpolate(1, 1, 3, 1, 2, 3, 3, 1, 1, 2)
    @test v ≈ 2.0
    # the crease is along the axis from (1,1) to (3,3), thus the value between (1,3) and
    # (3,1) is not linear and can have a different value than expected from a linear inter-
    # polation between the values 1 and 1
    v = EnergySystems.bilinear_interpolate(1, 2, 3, 1, 2, 3, 3, 1, 1, 2)
    @test v ≈ 2.5

    # concave crease from (1,1) to (3,3)
    #
    #   x    1   2   3
    #  y   -----------
    #    1 | 1       3
    #    2 |     v   v
    #    3 | 2   v   1
    v = EnergySystems.bilinear_interpolate(1, 2, 3, 1, 3, 3, 1, 3, 2, 1)
    @test v ≈ 1.5
    v = EnergySystems.bilinear_interpolate(1, 2, 3, 1, 2, 3, 1, 3, 2, 1)
    @test v ≈ 1.0
    v = EnergySystems.bilinear_interpolate(1, 3, 3, 1, 2, 3, 1, 3, 2, 1)
    @test v ≈ 2.0
end

@testset "test_bilinear_interpolate" begin
    test_bilinear_interpolate()
end

function test_inverse_efficiency()
    components_config = get_demand_energy_system_config()
    components_config["TST_GB_01"]["nr_discretization_steps"] = 10
    eps = 1e-9
    simulation_parameters = get_default_sim_params()
    w2wh = simulation_parameters["watt_to_wh"]

    components = Resie.load_components(components_config, simulation_parameters)
    boiler = components["TST_GB_01"]

    @test abs(plr_from_energy(boiler, Symbol("fuel_in"), 0.0, w2wh)) < eps
    @test abs(plr_from_energy(boiler, Symbol("fuel_in"), 1000.0, w2wh) - 1.0) < eps
    plr = plr_from_energy(boiler, Symbol("fuel_in"), 450.0, w2wh)
    @test plr > 0.45 - eps && plr < 0.45 + eps

    @test abs(plr_from_energy(boiler, Symbol("heat_out"), 0.0, w2wh)) < eps
    @test abs(plr_from_energy(boiler, Symbol("heat_out"), 1000.0, w2wh) - 1.0) < eps
    plr = plr_from_energy(boiler, Symbol("heat_out"), 450.0, w2wh)
    @test plr > 0.5614073352582631 - eps && plr < 0.5614073352582631 + eps

    boiler.discretization_step = 1.0 / 20
    EnergySystems.initialise!(boiler, simulation_parameters)

    @test abs(plr_from_energy(boiler, Symbol("heat_out"), 0.0, w2wh)) < eps
    @test abs(plr_from_energy(boiler, Symbol("heat_out"), 1000.0, w2wh) - 1.0) < eps
    plr = plr_from_energy(boiler, Symbol("heat_out"), 450.0, w2wh)
    @test plr > 0.561969105640274 - eps && plr < 0.561969105640274 + eps
end

@testset "test_inverse_efficiency" begin
    test_inverse_efficiency()
end

function test_gas_boiler_demand_driven_plrd()
    components_config = get_demand_energy_system_config()
    simulation_parameters = get_default_sim_params()

    components = Resie.load_components(components_config, simulation_parameters)
    gasboiler = components["TST_GB_01"]
    grid = components["TST_GRI_01"]
    demand = components["TST_DEM_01"]

    # test if the thermal efficiency for a gas boiler in demand-driven mode is calculated
    # correctly depending on the part-load-ratio (PLR)

    # first time step: demand is exactly the max power of GasBoiler

    for unit in values(components)
        EnergySystems.reset(unit)
    end

    expected_efficiency = 1.0

    EnergySystems.control(demand, components, simulation_parameters)
    EnergySystems.control(gasboiler, components, simulation_parameters)
    EnergySystems.control(grid, components, simulation_parameters)

    EnergySystems.process(demand, simulation_parameters)
    @test demand.input_interfaces[demand.medium].balance ≈ -1000

    EnergySystems.process(gasboiler, simulation_parameters)
    @test gasboiler.output_interfaces[gasboiler.m_heat_out].balance ≈ 0
    @test gasboiler.input_interfaces[gasboiler.m_fuel_in].balance ≈ -1000 / expected_efficiency
    @test demand.input_interfaces[demand.medium].balance ≈ 0

    EnergySystems.process(grid, simulation_parameters)
    @test grid.output_interfaces[grid.medium].balance ≈ 0

    # second time step: demand is half the max power of GasBoiler

    for unit in values(components)
        EnergySystems.reset(unit)
    end

    demand.constant_demand = 2000
    expected_efficiency = 0.832290069922532027

    EnergySystems.control(demand, components, simulation_parameters)
    EnergySystems.control(gasboiler, components, simulation_parameters)
    EnergySystems.control(grid, components, simulation_parameters)

    EnergySystems.process(demand, simulation_parameters)
    @test demand.input_interfaces[demand.medium].balance ≈ -500

    EnergySystems.process(gasboiler, simulation_parameters)
    @test isapprox(gasboiler.output_interfaces[gasboiler.m_heat_out].balance, 0.0, atol=1e-9)
    @test gasboiler.input_interfaces[gasboiler.m_fuel_in].balance ≈ -500 / expected_efficiency
    @test isapprox(demand.input_interfaces[demand.medium].balance, 0.0, atol=1e-9)

    EnergySystems.process(grid, simulation_parameters)
    @test grid.output_interfaces[grid.medium].balance ≈ 0
end

@testset "test_gas_boiler_demand_driven_plrd" begin
    test_gas_boiler_demand_driven_plrd()
end

function test_gas_boiler_supply_driven_plrd()
    components_config = Dict{String,Any}(
        "TST_DEM_01" => Dict{String,Any}(
            "type" => "BoundedSink",
            "medium" => "m_h_w_ht1",
            "output_refs" => [],
            "max_power_profile_file_path" => "./profiles/tests/source_heat_max_power.prf",
            "scale" => 1.0,
        ),
        "TST_GRI_01" => Dict{String,Any}(
            "type" => "FixedSupply",
            "medium" => "m_c_g_natgas",
            "output_refs" => ["TST_GB_01"],
            "is_source" => true,
            "energy_profile_file_path" => "./profiles/tests/demand_heating_energy.prf",
            "supply" => 0,
            "scale" => 1.0,
        ),
        "TST_GB_01" => Dict{String,Any}(
            "type" => "FuelBoiler",
            "m_fuel_in" => "m_c_g_natgas",
            "output_refs" => ["TST_DEM_01"],
            "power_th" => 4000,
            "efficiency_fuel_in" => "poly:-0.9117,1.8795,0.0322",
            "nr_discretization_steps" => 20,
        ),
    )

    simulation_parameters = get_default_sim_params()

    components = Resie.load_components(components_config, simulation_parameters)
    gasboiler = components["TST_GB_01"]
    grid = components["TST_GRI_01"]
    demand = components["TST_DEM_01"]

    # test if the thermal efficiency for a gas boiler in supply-driven mode is calculated
    # correctly depending on the part-load-ratio (PLR)

    # first time step: demand is exactly the max power of GasBoiler

    for unit in values(components)
        EnergySystems.reset(unit)
    end

    expected_efficiency = 1.0

    EnergySystems.control(grid, components, simulation_parameters)
    grid.supply = 1000 / expected_efficiency
    grid.output_interfaces[grid.medium].max_energy.max_energy[1] = grid.supply
    EnergySystems.control(gasboiler, components, simulation_parameters)
    EnergySystems.control(demand, components, simulation_parameters)
    demand.input_interfaces[demand.medium].max_energy.max_energy[1] = 1000 * 10.0
    demand.max_energy = 1000 * 10

    EnergySystems.process(grid, simulation_parameters)
    @test grid.output_interfaces[grid.medium].balance ≈ 1000 / expected_efficiency

    EnergySystems.process(gasboiler, simulation_parameters)
    @test gasboiler.output_interfaces[gasboiler.m_heat_out].balance ≈ 1000
    @test gasboiler.input_interfaces[gasboiler.m_fuel_in].balance ≈ 0
    @test demand.input_interfaces[demand.medium].balance ≈ 1000

    EnergySystems.process(demand, simulation_parameters)
    @test demand.input_interfaces[demand.medium].balance ≈ 0

    # second time step: demand is half the max power of GasBoiler

    for unit in values(components)
        EnergySystems.reset(unit)
    end

    expected_efficiency = 0.68375260628542269248721

    EnergySystems.control(grid, components, simulation_parameters)
    grid.supply = 500 / expected_efficiency
    grid.output_interfaces[grid.medium].max_energy.max_energy[1] = grid.supply
    EnergySystems.control(gasboiler, components, simulation_parameters)
    EnergySystems.control(demand, components, simulation_parameters)
    demand.input_interfaces[demand.medium].max_energy.max_energy[1] = 500 * 10.0
    demand.max_energy = 500 * 10

    EnergySystems.process(grid, simulation_parameters)
    @test grid.output_interfaces[grid.medium].balance ≈ 500 / expected_efficiency
    @test gasboiler.output_interfaces[gasboiler.m_heat_out].balance ≈ 0

    EnergySystems.process(gasboiler, simulation_parameters)
    @test gasboiler.output_interfaces[gasboiler.m_heat_out].balance ≈ 776.5061852364258
    @test gasboiler.input_interfaces[gasboiler.m_fuel_in].balance ≈ 0.0
    @test demand.input_interfaces[demand.medium].balance ≈ 776.5061852364258

    EnergySystems.process(demand, simulation_parameters)
    @test gasboiler.output_interfaces[gasboiler.m_heat_out].balance ≈ 0
    @test demand.input_interfaces[demand.medium].balance ≈ 0
end

@testset "test_gas_boiler_supply_driven_plrd" begin
    test_gas_boiler_supply_driven_plrd()
end

function test_CHPP_el_eff_plrd()
    components_config = Dict{String,Any}(
        "TST_DEM_01" => Dict{String,Any}(
            "type" => "Demand",
            "medium" => "m_h_w_ht1",
            "output_refs" => [],
            "energy_profile_file_path" => "./profiles/tests/demand_heating_energy.prf",
            "scale" => 1000.0,
        ),
        "TST_GRI_01" => Dict{String,Any}(
            "type" => "GridConnection",
            "medium" => "m_c_g_natgas",
            "output_refs" => ["TST_CHP_01"],
            "is_source" => true,
        ),
        "TST_GRO_01" => Dict{String,Any}(
            "type" => "GridConnection",
            "medium" => "m_e_ac_230v",
            "output_refs" => [],
            "is_source" => false,
        ),
        "TST_CHP_01" => Dict{String,Any}(
            "type" => "CHPP",
            "m_fuel_in" => "m_c_g_natgas",
            "output_refs" => ["TST_DEM_01",
                              "TST_GRO_01"],
            "power_el" => 5000,
            "linear_interface" => "el_out",
            "min_power_fraction" => 0.1,
            "efficiency_fuel_in" => "const:2.5",
            "efficiency_heat_out" => "pwlin:0.8,0.9,1.0,0.8",
            "efficiency_el_out" => "const:1.0",
            "nr_discretization_steps" => 25,
        ),
    )

    simulation_parameters = get_default_sim_params()

    components = Resie.load_components(components_config, simulation_parameters)
    chpp = components["TST_CHP_01"]
    grid_in = components["TST_GRI_01"]
    grid_out = components["TST_GRO_01"]
    demand = components["TST_DEM_01"]

    # test if the thermal efficiency for a CHPP in electricity-defined efficiency mode is
    # calculated correctly depending on the part-load-ratio (PLR)

    # first time step: demand is exactly the thermal power of CHPP at PLR of 1.0

    for unit in values(components)
        EnergySystems.reset(unit)
    end

    EnergySystems.control(chpp, components, simulation_parameters)
    EnergySystems.control(demand, components, simulation_parameters)
    demand.input_interfaces[demand.medium].max_energy.max_energy[1] = 1000.0
    demand.demand = 1000
    EnergySystems.control(grid_in, components, simulation_parameters)
    EnergySystems.control(grid_out, components, simulation_parameters)

    EnergySystems.process(demand, simulation_parameters)
    @test demand.input_interfaces[demand.medium].balance ≈ -1000

    EnergySystems.process(chpp, simulation_parameters)
    @test chpp.output_interfaces[chpp.m_heat_out].balance ≈ 0
    @test chpp.output_interfaces[chpp.m_el_out].balance ≈ 1000 / 0.8
    @test chpp.input_interfaces[chpp.m_fuel_in].balance ≈ -1000 / 0.4 / 0.8
    @test demand.input_interfaces[demand.medium].balance ≈ 0

    EnergySystems.process(grid_in, simulation_parameters)
    @test grid_in.output_interfaces[grid_in.medium].balance ≈ 0
    EnergySystems.process(grid_out, simulation_parameters)
    @test grid_out.input_interfaces[grid_out.medium].balance ≈ 0

    # second time step: demand is at highest thermal efficiency at PLR of 2/3 (of linear
    # efficiency of electricity output). however due to the discretization step in the
    # piece-wise linear efficiency not being a multiple of that of the inverse calculations,
    # the peak at efficiency of 1.0 is "missed", leading to inexact calculations, which
    # makes a fairly significant difference in absolute energy values

    for unit in values(components)
        EnergySystems.reset(unit)
    end

    EnergySystems.control(chpp, components, simulation_parameters)
    EnergySystems.control(demand, components, simulation_parameters)
    demand.input_interfaces[demand.medium].max_energy.max_energy[1] = 1250 * 2.0 / 3.0
    demand.demand = 1250 * 2.0 / 3.0
    EnergySystems.control(grid_in, components, simulation_parameters)
    EnergySystems.control(grid_out, components, simulation_parameters)

    EnergySystems.process(demand, simulation_parameters)
    @test demand.input_interfaces[demand.medium].balance >= -1250 * 2.0 / 3.0 - 5.0
    @test demand.input_interfaces[demand.medium].balance <= -1250 * 2.0 / 3.0 + 5.0

    EnergySystems.process(chpp, simulation_parameters)
    @test chpp.output_interfaces[chpp.m_heat_out].balance >= -5.0
    @test chpp.output_interfaces[chpp.m_heat_out].balance <= 5.0
    @test chpp.output_interfaces[chpp.m_el_out].balance >= 1250 * 2.0 / 3.0 - 10.0
    @test chpp.output_interfaces[chpp.m_el_out].balance <= 1250 * 2.0 / 3.0 + 10.0
    @test chpp.input_interfaces[chpp.m_fuel_in].balance >= -1250 * 2.0 / 3.0 / 0.4 - 20.0
    @test chpp.input_interfaces[chpp.m_fuel_in].balance <= -1250 * 2.0 / 3.0 / 0.4 + 20.0
    @test demand.input_interfaces[demand.medium].balance >= -5.0
    @test demand.input_interfaces[demand.medium].balance <= 5.0

    EnergySystems.process(grid_in, simulation_parameters)
    @test grid_in.output_interfaces[grid_in.medium].balance >= -5.0
    @test grid_in.output_interfaces[grid_in.medium].balance <= 5.0
    EnergySystems.process(grid_out, simulation_parameters)
    @test grid_out.input_interfaces[grid_out.medium].balance >= -5.0
    @test grid_out.input_interfaces[grid_out.medium].balance <= 5.0
end

@testset "test_CHPP_el_eff_plrd" begin
    test_CHPP_el_eff_plrd()
end

function test_electrolyser_dispatch_units()
    components_config = Dict{String,Any}(
        "TST_DEM_01" => Dict{String,Any}(
            "type" => "Demand",
            "medium" => "m_h_w_ht1",
            "output_refs" => [],
            "energy_profile_file_path" => "./profiles/tests/demand_heating_energy.prf",
            "scale" => 500.0,
        ),
        "TST_GRI_01" => Dict{String,Any}(
            "type" => "GridConnection",
            "medium" => "m_e_ac_230v",
            "output_refs" => ["TST_ELY_01"],
            "is_source" => true,
        ),
        "TST_GRO_01" => Dict{String,Any}(
            "type" => "GridConnection",
            "medium" => "m_c_g_h2",
            "output_refs" => [],
            "is_source" => false,
        ),
        "TST_GRO_02" => Dict{String,Any}(
            "type" => "GridConnection",
            "medium" => "m_c_g_o2",
            "output_refs" => [],
            "is_source" => false,
        ),
        "TST_ELY_01" => Dict{String,Any}(
            "type" => "Electrolyser",
            "output_refs" => ["TST_DEM_01",
                              "TST_GRO_01",
                              "TST_GRO_02"],
            "power_el" => 4000,
            "heat_lt_is_usable" => false,
            "nr_switchable_units" => 4,
            "dispatch_strategy" => "equal_with_mpf",
            "min_power_fraction_total" => 0.3,
            "min_power_fraction" => 0.4,
        ),
    )

    simulation_parameters = get_default_sim_params()
    w2wh = simulation_parameters["watt_to_wh"]

    components = Resie.load_components(components_config, simulation_parameters)
    electrolyser = components["TST_ELY_01"]

    nr_units, plr = EnergySystems.dispatch_units(electrolyser, 0.5, Symbol("el_in"), 500.0, w2wh)
    @test nr_units == 4
    @test isapprox(plr, 0.5, atol=1e-9)
    nr_units, plr = EnergySystems.dispatch_units(electrolyser, 0.31, Symbol("el_in"), 310.0, w2wh)
    @test nr_units == 3
    @test isapprox(plr, 0.41333333333, atol=1e-9)

    electrolyser.dispatch_strategy = "all_equal"
    nr_units, plr = EnergySystems.dispatch_units(electrolyser, 0.5, Symbol("el_in"), 500.0, w2wh)
    @test nr_units == 4
    @test isapprox(plr, 0.5, atol=1e-9)
    nr_units, plr = EnergySystems.dispatch_units(electrolyser, 0.3, Symbol("el_in"), 300.0, w2wh)
    @test nr_units == 4
    @test isapprox(plr, 0.3, atol=1e-9)

    electrolyser.dispatch_strategy = "try_optimal"
    nr_units, plr = EnergySystems.dispatch_units(electrolyser, 0.5, Symbol("el_in"), 500.0, w2wh)
    @test nr_units == 3
    @test isapprox(plr, 0.66666666666, atol=1e-9)
    nr_units, plr = EnergySystems.dispatch_units(electrolyser, 0.3, Symbol("el_in"), 300.0, w2wh)
    @test nr_units == 2
    @test isapprox(plr, 0.6, atol=1e-9)
end

@testset "test_electrolyser_dispatch_units" begin
    test_electrolyser_dispatch_units()
end
