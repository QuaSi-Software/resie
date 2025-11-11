using Dates
using ..Resie

"""
Control module for setting limits to the PLR of a component according to a price_profile and 
available storage capacity or demand.
"""
mutable struct CM_EconomicControl <: ControlModule
    name::String
    parameters::Dict{String,Any}
    state_machine::StateMachine
    price_profile::Profile
    connectivity_by_state::Dict{Int,ConnectionMatrix}
    ooo_by_state::Dict{Int,OrderOfOperations}

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

        # record original connectivity
        connectivities = Dict{Int,ConnectionMatrix}()
        connectivities[1] = components[unit_uac].connectivity

        # calculate new connectivity
        config = Dict{String,Any}(
            "connections" => params["new_connections"],
        )
        new_connectivity = ConnectionMatrix(config)
        new_components = deepcopy(components)
        new_components[unit_uac].connectivity = new_connectivity
        # follow steps load_components() in project_loading to make sure new order of
        # operations gets calculated correctly
        Resie.reorder_interfaces_of_bus!(new_components[unit_uac])
        initialise_components(new_components, sim_params)
        chains = Resie.find_chains(values(new_components), sf_bus)
        merge_bus_chains(chains, new_components, sim_params)

        # calculate new order_of_operations to have it available for later use
        ooo_by_state = Dict{Int,OrderOfOperations}()
        ooo_by_state[2] = Resie.calculate_order_of_operations(new_components)

        # set new connectivity and remove new_components and chains from memory since
        # they are not needed any more
        connectivities[2] = new_connectivity
        new_components = nothing
        chains = nothing

        # setup state machine for switching between states based on price
        # for now this is a simple two-state threshold switch, but it might implement a
        # more complicated calculation in the future
        table_state_off = TruthTable(;  # State: Off
                                     conditions=[function (state_machine)
                                                     return value_at_time(price_profile, sim_params) <=
                                                            params["limit_price"]
                                                 end],
                                     table_data=Dict{Tuple,UInt}(
                                         (false,) => 1,
                                         (true,) => 2,
                                     ))
        table_state_on = TruthTable(;  # State: On
                                    conditions=[function (state_machine)
                                                    return value_at_time(price_profile, sim_params) >
                                                           params["limit_price"]
                                                end],
                                    table_data=Dict{Tuple,UInt}(
                                        (true,) => 1,
                                        (false,) => 2,
                                    ))

        state_machine = StateMachine(UInt(1),               # state
                                     Dict{UInt,String}(     # state_names
                                         1 => "Off",
                                         2 => "On",
                                     ),
                                     Dict{UInt,TruthTable}( # transitions
                                         1 => table_state_off,
                                         2 => table_state_on,
                                     ))

        return new("economic_control", params, state_machine, price_profile,
                   connectivities, ooo_by_state)
    end
end

# method for control module name on type-level
control_module_name(x::Type{CM_EconomicControl})::String = "economic_control"

function has_method_for(mod::CM_EconomicControl, func::ControlModuleFunction)::Bool
    return func == cmf_change_bus_priorities || func == cmf_reorder_operations
end

function update(mod::CM_EconomicControl)
    # ordinarily we would update the state machine in the control modules's update step,
    # however this control module is special in that it's callbacks are used outside the
    # normal order of operations, so the update is performed by the callbacks instead
end

function reorder_operations(mod::CM_EconomicControl,
                            order_of_operations::OrderOfOperations,
                            sim_params::Dict{String,Any})::OrderOfOperations
    # record "original" ooo since it might be the first time the module has access to it
    # and also it might've changed due to other control modules
    if mod.state_machine.state == 1
        mod.ooo_by_state[1] = order_of_operations
    end
    return mod.ooo_by_state[mod.state_machine.state]
end

function change_bus_priorities!(mod::CM_EconomicControl,
                                components::Grouping,
                                sim_params::Dict{String,Any})
    # perform update outside of normal order of operations
    old_state = mod.state_machine.state
    move_state(mod.state_machine)

    # now reorder bus interfaces, but only if the state changed
    if old_state != mod.state_machine.state
        bus = components[mod.parameters["bus_uac"]]
        bus.connectivity = mod.connectivity_by_state[mod.state_machine.state]
        Resie.reorder_interfaces_of_bus!(bus)
        initialise!(bus, sim_params)
    end
end
