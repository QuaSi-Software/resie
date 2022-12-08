using Test

@testset "tests_are_working" begin
    @test true
end

# there is an order to the includes in the sense that later tests might use functionality
# that is being tested in earlier test sets, in particular the project loading methods
include("tests_project_loading.jl")
include("tests_control.jl")
include("tests_energy_systems.jl")
