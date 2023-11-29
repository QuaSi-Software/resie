using Test

@testset "energy_systems_individual_tests" begin
    include("energy_systems/balance_and_distribution/one_to_one.jl")
    include("energy_systems/balance_and_distribution/one_to_one_via_bus.jl")
    include("energy_systems/balance_and_distribution/many_to_one.jl")
    include("energy_systems/balance_and_distribution/one_to_many.jl")
    include("energy_systems/heat_pump_demand_driven.jl")
    include("energy_systems/energy_system_from_storage.jl")
    include("energy_systems/bus_to_bus.jl")
    include("energy_systems/storage_loading_switch.jl")
    include("energy_systems/multiple_transformer_limited.jl")
    include("energy_systems/gasboiler_demand_driven_with_bus.jl")
    include("energy_systems/extended_storage_control_strategy.jl")
    include("energy_systems/gas_boiler_plr_efficiency.jl")
end