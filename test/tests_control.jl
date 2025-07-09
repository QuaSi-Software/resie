using Test

@testset "control_tests" begin
    include("control_modules/storage_driven.jl")
    include("control_modules/temperature_sorting.jl")
end
