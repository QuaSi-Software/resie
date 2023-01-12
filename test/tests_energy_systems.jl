using Bran.EnergySystems
using Bran.Profiles

@testset "energy_systems_individual_tests" begin
    include("energy_systems/heat_pump_demand_driven.jl")
    include("energy_systems/demand_heating.jl")
end