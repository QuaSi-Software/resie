using Debugger
using Test
using Resie
using Resie.EnergySystems

@testset "project_loading_tests" begin
    include("initialization/loading.jl")
    include("initialization/bus.jl")
    include("initialization/bus_merge.jl")
    include("initialization/profiles.jl")
    include("initialization/parameters.jl")
    include("order_of_operations/load_order_of_operation.jl")
    include("order_of_operations/bus_output_priorities.jl")
    include("order_of_operations/bus_to_bus.jl")
    include("order_of_operations/one_bus_to_many_with_storage.jl")
    include("order_of_operations/heat_pumps_wrong.jl")
    include("order_of_operations/reorderings.jl")
    include("order_of_operations/storage_loading_switch.jl")
    include("order_of_operations/transformer_chains_ooo.jl")
    include("order_of_operations/transformer_chains_utils.jl")
end
