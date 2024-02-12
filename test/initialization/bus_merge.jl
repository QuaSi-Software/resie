using Test
using Resie
using Resie.EnergySystems

function energy_system_simple()::Dict{String,Any}
    return Dict{String,Any}(
        "TST_SRC_01" => Dict{String,Any}(
            "type" => "BoundedSupply",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "output_refs" => ["TST_BUS_01"],
            "constant_supply" => 1500,
            "constant_temperature" => 65
        ),
        "TST_SRC_02" => Dict{String,Any}(
            "type" => "BoundedSupply",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "output_refs" => ["TST_BUS_02"],
            "constant_supply" => 500,
            "constant_temperature" => 40
        ),
        "TST_BUS_01" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "connections" => Dict{String, Any}(
                "input_order" => [
                    "TST_SRC_01"
                ],
                "output_order" => [
                    "TST_DEM_01",
                    "TST_BUS_02"
                ],
                "energy_flow" => [
                    [1,1],
                ]
            )
        ),
        "TST_BUS_02" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "connections" => Dict{String, Any}(
                "input_order" => [
                    "TST_BUS_01",
                    "TST_SRC_02",
                ],
                "output_order" => [
                    "TST_DEM_02",
                ],
                "energy_flow" => [
                    [1],
                    [1],
                ]
            )
        ),
        "TST_DEM_01" => Dict{String,Any}(
            "type" => "Demand",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "output_refs" => [],
            "constant_demand" => 1000,
            "constant_temperature" => 65
        ),
        "TST_DEM_02" => Dict{String,Any}(
            "type" => "Demand",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "output_refs" => [],
            "constant_demand" => 1000,
            "constant_temperature" => 40
        ),
    )
end

function test_shallow_copy()
    components_config = energy_system_simple()

    simulation_parameters = Dict{String,Any}(
        "time_step_seconds" => 900,
        "time" => 0,
        "epsilon" => 1e-9
    )

    components = Resie.load_components(components_config, simulation_parameters)
    bus = components["TST_BUS_01"]
    new_bus = Bus(bus, true)

    bus.uac = "new_uac"
    @test new_bus.uac == "TST_BUS_01"

    bus.controller.parameter["foo"] = "bar"
    @test new_bus.controller.parameter["foo"] == "bar"
end

@testset "shallow_copy" begin
    test_shallow_copy()
end

function test_deep_copy()
    components_config = energy_system_simple()

    simulation_parameters = Dict{String,Any}(
        "time_step_seconds" => 900,
        "time" => 0,
        "epsilon" => 1e-9
    )

    components = Resie.load_components(components_config, simulation_parameters)
    bus = components["TST_BUS_01"]
    new_bus = EnergySystems.deepcopy(bus)

    bus.uac = "new_uac"
    @test new_bus.uac == "TST_BUS_01"

    bus.controller.parameter["foo"] = "bar"
    @test !haskey(new_bus.controller.parameter, "foo")

    bus.input_interfaces[1].source.uac = "new_uac"
    @test new_bus.input_interfaces[1].source.uac == "new_uac"
end

@testset "deep_copy" begin
    test_deep_copy()
end

function test_merge_busses()
    components_config = energy_system_simple()

    simulation_parameters = Dict{String,Any}(
        "time_step_seconds" => 900,
        "time" => 0,
        "epsilon" => 1e-9
    )

    components = Resie.load_components(components_config, simulation_parameters)
    new_bus = EnergySystems.merge(
        components["TST_BUS_01"],
        components["TST_BUS_02"],
        "TST_BUS_01"
    )

    expected = [
        "TST_SRC_01",
        "TST_SRC_02",
    ]
    inputs = [f.source.uac for f in new_bus.input_interfaces]
    @test inputs == expected
    @test new_bus.connectivity.input_order == expected

    expected = [
        "TST_DEM_01",
        "TST_DEM_02",
    ]
    outputs = [f.target.uac for f in new_bus.output_interfaces]
    @test outputs == expected
    @test new_bus.connectivity.output_order == expected

    expected = [
        ("TST_SRC_01", 1),
        ("TST_SRC_02", 2),
    ]
    inputs = [
        (row.source.uac, row.priority)
        for row in sort(collect(values(new_bus.balance_table_inputs)), by=x->x.priority)
    ]
    @test inputs == expected

    expected = [
        ("TST_DEM_01", 1),
        ("TST_DEM_02", 2),
    ]
    outputs = [
        (row.target.uac, row.priority)
        for row in sort(collect(values(new_bus.balance_table_outputs)), by=x->x.priority)
    ]
    @test outputs == expected

    expected = [
        [true,true],
        [false,true]
    ]
    @test new_bus.connectivity.energy_flow == expected
end

@testset "merge_busses" begin
    test_merge_busses()
end

function energy_system_complicated()::Dict{String,Any}
    return Dict{String,Any}(
        "TST_SRC_01" => Dict{String,Any}(
            "type" => "BoundedSupply",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "output_refs" => ["TST_BUS_01"],
            "constant_supply" => 1500,
            "constant_temperature" => 65
        ),
        "TST_SRC_02" => Dict{String,Any}(
            "type" => "BoundedSupply",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "output_refs" => ["TST_BUS_02"],
            "constant_supply" => 500,
            "constant_temperature" => 60
        ),
        "TST_SRC_03" => Dict{String,Any}(
            "type" => "BoundedSupply",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "output_refs" => ["TST_BUS_03"],
            "constant_supply" => 500,
            "constant_temperature" => 55
        ),
        "TST_SRC_04" => Dict{String,Any}(
            "type" => "BoundedSupply",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "output_refs" => ["TST_BUS_04"],
            "constant_supply" => 500,
            "constant_temperature" => 50
        ),
        "TST_SRC_05" => Dict{String,Any}(
            "type" => "BoundedSupply",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "output_refs" => ["TST_BUS_05"],
            "constant_supply" => 500,
            "constant_temperature" => 45
        ),
        "TST_BFT_01" => Dict{String,Any}(
            "type" => "Storage",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "output_refs" => ["TST_BUS_01"],
            "capacity" => 4000,
            "load" => 2000,
        ),
        "TST_BFT_02" => Dict{String,Any}(
            "type" => "Storage",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "output_refs" => ["TST_BUS_02"],
            "capacity" => 4000,
            "load" => 2000,
        ),
        "TST_BFT_03" => Dict{String,Any}(
            "type" => "Storage",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "output_refs" => ["TST_BUS_03"],
            "capacity" => 40000,
            "load" => 20000,
        ),
        "TST_BFT_04" => Dict{String,Any}(
            "type" => "Storage",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "output_refs" => ["TST_BUS_04"],
            "capacity" => 4000,
            "load" => 2000,
        ),
        "TST_BFT_05" => Dict{String,Any}(
            "type" => "Storage",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "output_refs" => ["TST_BUS_05"],
            "capacity" => 4000,
            "load" => 2000,
        ),
        "TST_BUS_01" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "connections" => Dict{String, Any}(
                "input_order" => [
                    "TST_SRC_01",
                    "TST_BFT_01"
                ],
                "output_order" => [
                    "TST_DEM_01",
                    "TST_BFT_01",
                    "TST_BUS_03"
                ],
                "energy_flow" => [
                    [1,1,1],
                    [1,0,0],
                ]
            )
        ),
        "TST_BUS_02" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "connections" => Dict{String, Any}(
                "input_order" => [
                    "TST_SRC_02",
                    "TST_BFT_02"
                ],
                "output_order" => [
                    "TST_DEM_02",
                    "TST_BFT_02",
                    "TST_BUS_03"
                ],
                "energy_flow" => [
                    [1,1,1],
                    [1,0,0],
                ]
            )
        ),
        "TST_BUS_03" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "connections" => Dict{String, Any}(
                "input_order" => [
                    "TST_SRC_03",
                    "TST_BUS_02",
                    "TST_BUS_01",
                    "TST_BFT_03",
                ],
                "output_order" => [
                    "TST_DEM_03",
                    "TST_BFT_03",
                    "TST_BUS_04",
                    "TST_BUS_05",
                ],
                "energy_flow" => [
                    [1,1,1,1],
                    [1,1,1,1],
                    [1,1,1,1],
                    [1,0,1,1],
                ]
            )
        ),
        "TST_BUS_04" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "connections" => Dict{String, Any}(
                "input_order" => [
                    "TST_SRC_04",
                    "TST_BUS_03",
                    "TST_BFT_04"
                ],
                "output_order" => [
                    "TST_DEM_04",
                    "TST_BFT_04",
                ],
                "energy_flow" => [
                    [1,1],
                    [1,1],
                    [1,0],
                ]
            )
        ),
        "TST_BUS_05" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "connections" => Dict{String, Any}(
                "input_order" => [
                    "TST_SRC_05",
                    "TST_BUS_03",
                    "TST_BFT_05"
                ],
                "output_order" => [
                    "TST_DEM_05",
                    "TST_BFT_05",
                ],
                "energy_flow" => [
                    [1,1],
                    [1,1],
                    [1,0],
                ]
            )
        ),
        "TST_DEM_01" => Dict{String,Any}(
            "type" => "Demand",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "output_refs" => [],
            "constant_demand" => 1000,
            "constant_temperature" => 65
        ),
        "TST_DEM_02" => Dict{String,Any}(
            "type" => "Demand",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "output_refs" => [],
            "constant_demand" => 1000,
            "constant_temperature" => 60
        ),
        "TST_DEM_03" => Dict{String,Any}(
            "type" => "Demand",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "output_refs" => [],
            "constant_demand" => 1000,
            "constant_temperature" => 55
        ),
        "TST_DEM_04" => Dict{String,Any}(
            "type" => "Demand",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "output_refs" => [],
            "constant_demand" => 1000,
            "constant_temperature" => 50
        ),
        "TST_DEM_05" => Dict{String,Any}(
            "type" => "Demand",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "output_refs" => [],
            "constant_demand" => 1000,
            "constant_temperature" => 45
        ),
    )
end

function test_merge_busses_cross_shape_half_manual()
    components_config = energy_system_complicated()

    simulation_parameters = Dict{String,Any}(
        "time_step_seconds" => 900,
        "time" => 0,
        "epsilon" => 1e-9
    )

    components = Resie.load_components(components_config, simulation_parameters)
    new_bus = EnergySystems.merge(
        components["TST_BUS_03"],
        components["TST_BUS_05"],
        "TST_BUS_03"
    )
    new_bus = EnergySystems.merge(
        new_bus,
        components["TST_BUS_04"],
        new_bus.uac
    )

    expected = [
        "TST_SRC_04",
        "TST_SRC_05",
        "TST_SRC_03",
        "TST_BUS_02",
        "TST_BUS_01",
        "TST_BFT_03",
        "TST_BFT_05",
        "TST_BFT_04",
    ]
    inputs = [f.source.uac for f in new_bus.input_interfaces]
    @test inputs == expected
    @test new_bus.connectivity.input_order == expected

    expected = [
        "TST_DEM_03",
        "TST_BFT_03",
        "TST_DEM_04",
        "TST_BFT_04",
        "TST_DEM_05",
        "TST_BFT_05",
    ]
    outputs = [f.target.uac for f in new_bus.output_interfaces]
    @test outputs == expected
    @test new_bus.connectivity.output_order == expected

    expected = [
        ("TST_SRC_04",1),
        ("TST_SRC_05",2),
        ("TST_SRC_03",3),
        ("TST_BUS_02",4),
        ("TST_BUS_01",5),
        ("TST_BFT_03",6),
        ("TST_BFT_05",7),
        ("TST_BFT_04",8),
    ]
    inputs = [
        (row.source.uac, row.priority)
        for row in sort(collect(values(new_bus.balance_table_inputs)), by=x->x.priority)
    ]
    @test inputs == expected

    expected = [
        ("TST_DEM_03",1),
        ("TST_BFT_03",2),
        ("TST_DEM_04",3),
        ("TST_BFT_04",4),
        ("TST_DEM_05",5),
        ("TST_BFT_05",6),
    ]
    outputs = [
        (row.target.uac, row.priority)
        for row in sort(collect(values(new_bus.balance_table_outputs)), by=x->x.priority)
    ]
    @test outputs == expected

    expected = [
        [0,0,1,1,0,0],
        [0,0,0,0,1,1],
        [1,1,1,1,1,1],
        [1,1,1,1,1,1],
        [1,1,1,1,1,1],
        [1,0,1,1,1,1],
        [0,0,0,0,1,0],
        [0,0,1,0,0,0],
    ]
    @test new_bus.connectivity.energy_flow == expected
end

@testset "merge_busses_cross_shape_half_manual" begin
    test_merge_busses_cross_shape_half_manual()
end

function test_merge_busses_cross_shape_full_manual()
    components_config = energy_system_complicated()
    simulation_parameters = Dict{String,Any}(
        "time_step_seconds" => 900,
        "time" => 0,
        "epsilon" => 1e-9
    )
    components = Resie.load_components(components_config, simulation_parameters)

    new_bus = EnergySystems.merge(
        components["TST_BUS_03"],
        components["TST_BUS_05"],
        "TST_BUS_03"
    )
    new_bus = EnergySystems.merge(
        new_bus,
        components["TST_BUS_04"],
        new_bus.uac
    )
    new_bus = EnergySystems.merge(
        components["TST_BUS_01"],
        new_bus,
        new_bus.uac
    )

    expected = [
        "TST_SRC_04",
        "TST_SRC_05",
        "TST_SRC_03",
        "TST_BUS_02",
        "TST_SRC_01",
        "TST_BFT_01",
        "TST_BFT_03",
        "TST_BFT_05",
        "TST_BFT_04",
    ]
    inputs = [f.source.uac for f in new_bus.input_interfaces]
    @test inputs == expected
    @test new_bus.connectivity.input_order == expected

    expected = [
        "TST_DEM_01",
        "TST_BFT_01",
        "TST_DEM_03",
        "TST_BFT_03",
        "TST_DEM_04",
        "TST_BFT_04",
        "TST_DEM_05",
        "TST_BFT_05",
    ]
    outputs = [f.target.uac for f in new_bus.output_interfaces]
    @test outputs == expected
    @test new_bus.connectivity.output_order == expected

    expected = [
        ("TST_SRC_04",1),
        ("TST_SRC_05",2),
        ("TST_SRC_03",3),
        ("TST_BUS_02",4),
        ("TST_SRC_01",5),
        ("TST_BFT_01",6),
        ("TST_BFT_03",7),
        ("TST_BFT_05",8),
        ("TST_BFT_04",9),
    ]
    inputs = [
        (row.source.uac, row.priority)
        for row in sort(collect(values(new_bus.balance_table_inputs)), by=x->x.priority)
    ]
    @test inputs == expected

    expected = [
        ("TST_DEM_01",1),
        ("TST_BFT_01",2),
        ("TST_DEM_03",3),
        ("TST_BFT_03",4),
        ("TST_DEM_04",5),
        ("TST_BFT_04",6),
        ("TST_DEM_05",7),
        ("TST_BFT_05",8),
    ]
    outputs = [
        (row.target.uac, row.priority)
        for row in sort(collect(values(new_bus.balance_table_outputs)), by=x->x.priority)
    ]
    @test outputs == expected

    expected = [
        [0,0,0,0,1,1,0,0],
        [0,0,0,0,0,0,1,1],
        [0,0,1,1,1,1,1,1],
        [0,0,1,1,1,1,1,1],
        [1,1,1,1,1,1,1,1],
        [1,0,0,0,0,0,0,0],
        [0,0,1,0,1,1,1,1],
        [0,0,0,0,0,0,1,0],
        [0,0,0,0,1,0,0,0],
    ]
    @test new_bus.connectivity.energy_flow == expected

    new_bus = EnergySystems.merge(
        components["TST_BUS_02"],
        new_bus,
        "TST_BUS_03"
    )

    expected = [
        "TST_SRC_04",
        "TST_SRC_05",
        "TST_SRC_03",
        "TST_SRC_02",
        "TST_BFT_02",
        "TST_SRC_01",
        "TST_BFT_01",
        "TST_BFT_03",
        "TST_BFT_05",
        "TST_BFT_04",
    ]
    inputs = [f.source.uac for f in new_bus.input_interfaces]
    @test inputs == expected
    @test new_bus.connectivity.input_order == expected

    expected = [
        "TST_DEM_02",
        "TST_BFT_02",
        "TST_DEM_01",
        "TST_BFT_01",
        "TST_DEM_03",
        "TST_BFT_03",
        "TST_DEM_04",
        "TST_BFT_04",
        "TST_DEM_05",
        "TST_BFT_05",
    ]
    outputs = [f.target.uac for f in new_bus.output_interfaces]
    @test outputs == expected
    @test new_bus.connectivity.output_order == expected

    expected = [
        ("TST_SRC_04",1),
        ("TST_SRC_05",2),
        ("TST_SRC_03",3),
        ("TST_SRC_02",4),
        ("TST_BFT_02",5),
        ("TST_SRC_01",6),
        ("TST_BFT_01",7),
        ("TST_BFT_03",8),
        ("TST_BFT_05",9),
        ("TST_BFT_04",10),
    ]
    inputs = [
        (row.source.uac, row.priority)
        for row in sort(collect(values(new_bus.balance_table_inputs)), by=x->x.priority)
    ]
    @test inputs == expected

    expected = [
        ("TST_DEM_02",1),
        ("TST_BFT_02",2),
        ("TST_DEM_01",3),
        ("TST_BFT_01",4),
        ("TST_DEM_03",5),
        ("TST_BFT_03",6),
        ("TST_DEM_04",7),
        ("TST_BFT_04",8),
        ("TST_DEM_05",9),
        ("TST_BFT_05",10),
    ]
    outputs = [
        (row.target.uac, row.priority)
        for row in sort(collect(values(new_bus.balance_table_outputs)), by=x->x.priority)
    ]
    @test outputs == expected

    expected = [
        [0,0,0,0,0,0,1,1,0,0],
        [0,0,0,0,0,0,0,0,1,1],
        [0,0,0,0,1,1,1,1,1,1],
        [1,1,0,0,1,1,1,1,1,1],
        [1,0,0,0,0,0,0,0,0,0],
        [0,0,1,1,1,1,1,1,1,1],
        [0,0,1,0,0,0,0,0,0,0],
        [0,0,0,0,1,0,1,1,1,1],
        [0,0,0,0,0,0,0,0,1,0],
        [0,0,0,0,0,0,1,0,0,0],
    ]
    @test new_bus.connectivity.energy_flow == expected

end

@testset "merge_busses_cross_shape_full_manual" begin
    test_merge_busses_cross_shape_full_manual()
end

function test_merge_busses_cross_shape_full_automatic()
    components_config = energy_system_complicated()
    simulation_parameters = Dict{String,Any}(
        "time_step_seconds" => 900,
        "time" => 0,
        "epsilon" => 1e-9
    )
    components = Resie.load_components(components_config, simulation_parameters)
    busses = Grouping()
    for (key,val) in pairs(components)
        if val.sys_function == EnergySystems.sf_bus
            busses[key] = val
        end
    end
    new_bus = EnergySystems.merge_busses(busses, components)

    expected = [
        "TST_SRC_04",
        "TST_SRC_05",
        "TST_SRC_03",
        "TST_SRC_02",
        "TST_BFT_02",
        "TST_SRC_01",
        "TST_BFT_01",
        "TST_BFT_03",
        "TST_BFT_05",
        "TST_BFT_04",
    ]
    inputs = [f.source.uac for f in new_bus.input_interfaces]
    @test inputs == expected
    @test new_bus.connectivity.input_order == expected

    expected = [
        "TST_DEM_02",
        "TST_BFT_02",
        "TST_DEM_01",
        "TST_BFT_01",
        "TST_DEM_03",
        "TST_BFT_03",
        "TST_DEM_04",
        "TST_BFT_04",
        "TST_DEM_05",
        "TST_BFT_05",
    ]
    outputs = [f.target.uac for f in new_bus.output_interfaces]
    @test outputs == expected
    @test new_bus.connectivity.output_order == expected

    expected = [
        ("TST_SRC_04",1),
        ("TST_SRC_05",2),
        ("TST_SRC_03",3),
        ("TST_SRC_02",4),
        ("TST_BFT_02",5),
        ("TST_SRC_01",6),
        ("TST_BFT_01",7),
        ("TST_BFT_03",8),
        ("TST_BFT_05",9),
        ("TST_BFT_04",10),
    ]
    inputs = [
        (row.source.uac, row.priority)
        for row in sort(collect(values(new_bus.balance_table_inputs)), by=x->x.priority)
    ]
    @test inputs == expected

    expected = [
        ("TST_DEM_02",1),
        ("TST_BFT_02",2),
        ("TST_DEM_01",3),
        ("TST_BFT_01",4),
        ("TST_DEM_03",5),
        ("TST_BFT_03",6),
        ("TST_DEM_04",7),
        ("TST_BFT_04",8),
        ("TST_DEM_05",9),
        ("TST_BFT_05",10),
    ]
    outputs = [
        (row.target.uac, row.priority)
        for row in sort(collect(values(new_bus.balance_table_outputs)), by=x->x.priority)
    ]
    @test outputs == expected

    expected = [
        [0,0,0,0,0,0,1,1,0,0],
        [0,0,0,0,0,0,0,0,1,1],
        [0,0,0,0,1,1,1,1,1,1],
        [1,1,0,0,1,1,1,1,1,1],
        [1,0,0,0,0,0,0,0,0,0],
        [0,0,1,1,1,1,1,1,1,1],
        [0,0,1,0,0,0,0,0,0,0],
        [0,0,0,0,1,0,1,1,1,1],
        [0,0,0,0,0,0,0,0,1,0],
        [0,0,0,0,0,0,1,0,0,0],
    ]
    @test new_bus.connectivity.energy_flow == expected
end

@testset "merge_busses_cross_shape_full_automatic" begin
    test_merge_busses_cross_shape_full_automatic()
end

function energy_system_busses_only()::Dict{String,Any}
    return Dict{String,Any}(
        "TST_BUS_01" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "connections" => Dict{String, Any}(
                "input_order" => [],
                "output_order" => [
                    "TST_BUS_03"
                ],
            )
        ),
        "TST_BUS_02" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "connections" => Dict{String, Any}(
                "input_order" => [],
                "output_order" => [
                    "TST_BUS_04"
                ],
            )
        ),
        "TST_BUS_03" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "connections" => Dict{String, Any}(
                "input_order" => [
                    "TST_BUS_01",
                ],
                "output_order" => [
                    "TST_BUS_04",
                ],
                "energy_flow" => [
                    [1],
                ]
            )
        ),
        "TST_BUS_04" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "connections" => Dict{String, Any}(
                "input_order" => [
                    "TST_BUS_02",
                    "TST_BUS_03",
                ],
                "output_order" => [
                    "TST_BUS_05",
                ],
                "energy_flow" => [
                    [1],
                    [1],
                ]
            )
        ),
        "TST_BUS_05" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "connections" => Dict{String, Any}(
                "input_order" => [
                    "TST_BUS_04",
                ],
                "output_order" => [
                    "TST_BUS_06",
                    "TST_BUS_07",
                ],
                "energy_flow" => [
                    [1,1],
                ]
            )
        ),
        "TST_BUS_06" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "connections" => Dict{String, Any}(
                "input_order" => [
                    "TST_BUS_05",
                ],
                "output_order" => [
                    "TST_BUS_08",
                ],
                "energy_flow" => [
                    [1],
                ]
            )
        ),
        "TST_BUS_07" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "connections" => Dict{String, Any}(
                "input_order" => [
                    "TST_BUS_05",
                ],
                "output_order" => [],
            )
        ),
        "TST_BUS_08" => Dict{String,Any}(
            "type" => "Bus",
            "medium" => "m_h_w_ht1",
            "control_refs" => [],
            "connections" => Dict{String, Any}(
                "input_order" => [
                    "TST_BUS_06",
                ],
                "output_order" => [],
            )
        ),
    )
end

function test_merge_busses_order_of_merges()
    components_config = energy_system_busses_only()
    simulation_parameters = Dict{String,Any}(
        "time_step_seconds" => 900,
        "time" => 0,
        "epsilon" => 1e-9
    )
    components = Resie.load_components(components_config, simulation_parameters)
    chains = Resie.find_chains(values(components), EnergySystems.sf_bus)
    chain = Grouping()
    for component in chains[1]
        chain[component.uac] = component
    end
    merged = EnergySystems.merge_busses(chain, components)
    @test merged isa Bus
end

@testset "merge_busses_order_of_merges" begin
    test_merge_busses_order_of_merges()
end