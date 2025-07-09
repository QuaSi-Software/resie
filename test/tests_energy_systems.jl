using Test

@testset "energy_systems_individual_tests" begin
    include("energy_systems/balance_and_distribution/one_to_one.jl")
    include("energy_systems/balance_and_distribution/one_to_one_via_bus.jl")
    include("energy_systems/balance_and_distribution/many_to_one.jl")
    include("energy_systems/balance_and_distribution/one_to_many.jl")
    include("energy_systems/balance_and_distribution/many_plus_storage_to_one.jl")
    include("energy_systems/balance_and_distribution/one_plus_storage_to_many.jl")
    include("energy_systems/balance_and_distribution/many_to_many.jl")
    include("energy_systems/balance_and_distribution/one_bus_to_one_bus.jl")
    include("energy_systems/balance_and_distribution/one_bus_to_many_bus.jl")
    include("energy_systems/balance_and_distribution/bus_to_bus_distribution.jl")
    include("energy_systems/heat_pump.jl")
    include("energy_systems/energy_system_from_storage.jl")
    include("energy_systems/storage_loading_switch.jl")
    include("energy_systems/multiple_transformers.jl")
    include("energy_systems/gasboiler_demand_driven_with_bus.jl")
    include("energy_systems/plr_dependent_efficiency.jl")
end
