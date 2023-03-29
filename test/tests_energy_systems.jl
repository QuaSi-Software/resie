using Test

@testset "energy_systems_individual_tests" begin
    include("energy_systems/heat_pump_demand_driven.jl")
    include("energy_systems/demand_heating.jl")
    include("energy_systems/bus_to_bus.jl")
    include("energy_systems/storage_loading_switch.jl")
    include("energy_systems/gasboiler_demand_driven_with_bus.jl")
    include("energy_systems/extended_storage_control_strategy.jl")

end