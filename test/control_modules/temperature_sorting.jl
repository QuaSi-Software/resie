using Test
using Resie
using Resie.EnergySystems

function test_coalesce_greater()
    @test EnergySystems.coalesce_greater(1.0, 2.0) == false
    @test EnergySystems.coalesce_greater(1.0, 1.0) == false
    @test EnergySystems.coalesce_greater(2.0, 1.0) == true
    @test EnergySystems.coalesce_greater(1.0, nothing) == true
    @test EnergySystems.coalesce_greater(nothing, 1.0) == false
    @test EnergySystems.coalesce_greater(nothing, nothing) == false
end

@testset "coalesce_greater" begin
    test_coalesce_greater()
end

function test_sorting_inputs_correctly()
    control_module = EnergySystems.CM_Temperature_Sorting(Dict{String,Any}(
                                                              "name" => "temperature_sorting",
                                                              "input_temps" => "max",
                                                              "input_order" => "desc",
                                                              "output_temps" => "min",
                                                              "output_order" => "none",
                                                          ), EnergySystems.Grouping(), Dict{String,Any}())
    inputs_min = [nothing, nothing, nothing]
    inputs_max = [20.0, 10.0, 30.0]
    @test EnergySystems.reorder_inputs(control_module, inputs_min, inputs_max) == [3, 1, 2]
    outputs_min = [30.0, 10.0, 20.0]
    outputs_max = [nothing, nothing, nothing]
    @test EnergySystems.reorder_outputs(control_module, outputs_min, outputs_max) == [1, 2, 3]
end

function test_sorting_outputs_correctly()
    control_module = EnergySystems.CM_Temperature_Sorting(Dict{String,Any}(
                                                              "name" => "temperature_sorting",
                                                              "input_temps" => "max",
                                                              "input_order" => "none",
                                                              "output_temps" => "min",
                                                              "output_order" => "asc",
                                                          ), EnergySystems.Grouping(), Dict{String,Any}())
    outputs_min = [30.0, 10.0, 20.0]
    outputs_max = [nothing, nothing, nothing]
    @test EnergySystems.reorder_outputs(control_module, outputs_min, outputs_max) == [2, 3, 1]
    inputs_min = [nothing, nothing, nothing]
    inputs_max = [20.0, 10.0, 30.0]
    @test EnergySystems.reorder_inputs(control_module, inputs_min, inputs_max) == [1, 2, 3]
end

function test_sorting_both_correctly()
    control_module = EnergySystems.CM_Temperature_Sorting(Dict{String,Any}(
                                                              "name" => "temperature_sorting",
                                                              "input_temps" => "min",
                                                              "input_order" => "asc",
                                                              "output_temps" => "max",
                                                              "output_order" => "desc",
                                                          ), EnergySystems.Grouping(), Dict{String,Any}())
    inputs_min = [30.0, 20.0, 10.0]
    inputs_max = [nothing, nothing, nothing]
    outputs_min = [nothing, nothing, nothing]
    outputs_max = [10.0, 20.0, 30.0]
    @test EnergySystems.reorder_inputs(control_module, inputs_min, inputs_max) == [3, 2, 1]
    @test EnergySystems.reorder_outputs(control_module, outputs_min, outputs_max) == [3, 2, 1]
end

function test_sorting_with_nothing_correctly()
    control_module = EnergySystems.CM_Temperature_Sorting(Dict{String,Any}(
                                                              "name" => "temperature_sorting",
                                                              "input_temps" => "max",
                                                              "input_order" => "desc",
                                                              "output_temps" => "min",
                                                              "output_order" => "none",
                                                          ), EnergySystems.Grouping(), Dict{String,Any}())
    inputs_min = [nothing, nothing, nothing]
    inputs_max = [nothing, 20.0, 10.0, nothing, 30.0, nothing]
    results = EnergySystems.reorder_inputs(control_module, inputs_min, inputs_max)
    # the nothings can be in any order, hence they are not tested
    @test results[1] == 5
    @test results[2] == 2
    @test results[3] == 3
end

@testset "sorting_correctly" begin
    test_sorting_inputs_correctly()
    test_sorting_outputs_correctly()
    test_sorting_both_correctly()
    test_sorting_with_nothing_correctly()
end
