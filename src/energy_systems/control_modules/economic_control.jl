using Dates
using ..Resie

"""
Control module for setting limits to the PLR of a component according to a price_profile and 
available storage capacity or demand.
"""
mutable struct CM_EconomicControl <: ControlModule
    name::String
    parameters::Dict{String,Any}
    price_profile::Profile
    new_connectivity::ConnectionMatrix
    original_connectivity::ConnectionMatrix
    new_OoO::OrderOfOperations

    function CM_EconomicControl(parameters::Dict{String,Any},
                                components::Grouping,
                                sim_params::Dict{String,Any},
                                unit_uac::String)
        default_parameters = Dict{String,Any}(
            "name" => "economic_control",
            "price_profile_path" => nothing,
            "limit_price" => 0.0,
            "new_connections" => nothing,
            "bus_uac" => unit_uac,
        )
        params = Base.merge(default_parameters, parameters)

        if params["price_profile_path"] === nothing
            @error "Required price profile path for control module economic_control not given"
        end
        price_profile = Profile(params["price_profile_path"], sim_params)

        if params["new_connections"] === nothing
            @error "new_connections must be given to be applied when the price is below " *
                   "the limit"
        end
        config = Dict{String,Any}(
            "connections" => params["new_connections"],
        )
        new_connectivity = ConnectionMatrix(config)
        original_connectivity = components[unit_uac].connectivity

        new_components = deepcopy(components)
        new_components[unit_uac].connectivity = new_connectivity
        # follow steps load_components() in project_loading to make sure new_OoO gets 
        # calculated correctly
        Resie.reorder_interfaces_of_bus!(new_components[unit_uac])
        initialise_components(new_components, sim_params)
        chains = Resie.find_chains(values(new_components), sf_bus)
        merge_bus_chains(chains, new_components, sim_params)
        # calculate new order_of_operations to have it available for later use
        new_OoO = Resie.calculate_order_of_operations(new_components)
        # remove new_components and chains from memory since they are not needed any more
        new_components = nothing
        chains = nothing

        return new("economic_control", params, price_profile, new_connectivity,
                   original_connectivity, new_OoO)
    end
end

# method for control module name on type-level
control_module_name(x::Type{CM_EconomicControl})::String = "economic_control"

function has_method_for(mod::CM_EconomicControl, func::ControlModuleFunction)::Bool
    return func == cmf_change_bus_priorities || func == cmf_reorder_operations
end

function update(mod::CM_EconomicControl)
    # nothing to do
end

function reorder_operations(mod::CM_EconomicControl,
                            order_of_operations::OrderOfOperations,
                            sim_params::Dict{String,Any})::OrderOfOperations
    if value_at_time(mod.price_profile, sim_params) <= mod.parameters["limit_price"]
        return mod.new_OoO
    end
    return order_of_operations
end

function change_bus_priorities!(mod::CM_EconomicControl,
                                components::Grouping,
                                sim_params::Dict{String,Any})
    bus = components[mod.parameters["bus_uac"]]
    if value_at_time(mod.price_profile, sim_params) <= mod.parameters["limit_price"]
        bus.connectivity = mod.new_connectivity
    else
        bus.connectivity = mod.original_connectivity
    end
    Resie.reorder_interfaces_of_bus!(bus)
    initialise!(bus, sim_params)
end
