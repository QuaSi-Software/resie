using Dates
include("../../order_of_operations.jl")
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
    new_OoO::StepInstructions

    function CM_EconomicControl(parameters::Dict{String,Any},
                                components::Grouping,
                                sim_params::Dict{String,Any},
                                unit_uac::String)
        default_parameters = Dict{String,Any}(
            "name" => "economic_control",
            "price_profile_path" => nothing,
            "limit_price" => 0.0,
            "new_connections" => nothing,
            "bus_uac" => unit_uac
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
        config = Dict{String, Any}(
            "connections" => params["new_connections"]
        )
        new_connectivity = ConnectionMatrix(config)
        original_connectivity = components[unit_uac].connectivity

        new_components = deepcopy(components)
        new_components[unit_uac].connectivity = new_connectivity
        # follow steps load_components() in project_loading to make sure new_OoO gets 
        # calculated correctly
        reorder_interfaces_of_bus!(new_components[unit_uac])
        initialise_components(new_components, sim_params)
        chains = find_chains(values(new_components), sf_bus)
        merge_bus_chains(chains, new_components, sim_params)
        # calculate new order_of_operations to have it available for later use
        new_OoO = calculate_order_of_operations(new_components)
        # remove new_components and chains from memory since they are not needed any more
        new_components = nothing
        chains = nothing

        return new("economic_control", params, price_profile, new_connectivity, 
                   original_connectivity, new_OoO)
    end

    function CM_EconomicControl()
        dummy_dates = DateTime[]
        push!(dummy_dates, DateTime(2024, 1, 1, 0, 0, 0))
        dummy_values = Float64[]
        push!(dummy_values, 0.0)
        dummy_profile = Profile("dummy",
                                Dict{String,Any}(
                                    "time_step_seconds" => 1,
                                    "start_date" => DateTime(2024, 1, 1, 0, 0, 0),
                                    "end_date" => DateTime(2024, 1, 1, 0, 0, 0),
                                    "force_profiles_to_repeat" => false,
                                );
                                given_profile_values=dummy_values,
                                given_timestamps=dummy_dates,
                                given_time_step=Second(1),
                                given_data_type="intensive")
        return new("economic_control", Dict{String,Any}(), dummy_profile, 
                   ConnectionMatrix([], [], nothing), ConnectionMatrix([], [], nothing))
    end
end

function has_method_for(mod::CM_EconomicControl, func::ControlModuleFunction)::Bool
    return func == cmf_change_priorities
end

function update(mod::CM_EconomicControl)
    # nothing to do
end

"""
reorder_interfaces_of_bus(components)

Reorder the input and output interfaces of a bus according to the input and output
priorities given in the connectivity matrix.
"""
function reorder_interfaces_of_bus!(bus::Bus)
    # get correct order according to connectivity matrix
    output_order = bus.connectivity.output_order
    input_order = bus.connectivity.input_order

    # Create a dictionary to map 'uac' to its correct position
    output_order_dict = Dict(uac => idx for (idx, uac) in enumerate(output_order))
    input_order_dict = Dict(uac => idx for (idx, uac) in enumerate(input_order))

    # Get the permutation indices that would sort the 'source'/'target' field by
    # 'uac' order
    output_perm_indices = sortperm([output_order_dict[bus.output_interfaces[i].target.uac]
                                    for i in 1:length(bus.output_interfaces)])
    input_perm_indices = sortperm([input_order_dict[bus.input_interfaces[i].source.uac]
                                    for i in 1:length(bus.input_interfaces)])

    # Reorder the input and output interfaces using the permutation indices
    bus.output_interfaces = bus.output_interfaces[output_perm_indices]
    bus.input_interfaces = bus.input_interfaces[input_perm_indices]
end

function change_priorities(mod::CM_EconomicControl, 
                           components::Grouping,
                           order_of_operations::StepInstructions, 
                           sim_params::Dict{String,Any})::StepInstructions
    bus = components[mod.parameters["bus_uac"]]
    if value_at_time(mod.price_profile, sim_params) <= mod.parameters["limit_price"]
        order_of_operations = mod.new_OoO
        bus.connectivity = mod.new_connectivity
    else
        bus.connectivity = mod.original_connectivity
    end
    reorder_interfaces_of_bus!(bus)
    initialise!(bus, sim_params)

    return order_of_operations
end
