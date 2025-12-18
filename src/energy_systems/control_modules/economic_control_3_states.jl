using Dates
using ..Resie
using Infiltrator
include("../../file_output.jl")
using OrderedCollections: OrderedDict
"""
Control module for setting limits to the PLR of a component according to a price_profile and 
available storage capacity or demand.
"""
mutable struct CM_EconomicControl3States <: ControlModule
    name::String
    parameters::Dict{String,Any}
    state_machine::StateMachine
    price_profiles::Array{Profile}
    connectivity_by_state::Dict{Int,ConnectionMatrix}
    ooo_by_state::Dict{Int,Union{OrderOfOperations,Nothing}}

    function CM_EconomicControl3States(parameters::Dict{String,Any},
                                       components::Grouping,
                                       sim_params::Dict{String,Any},
                                       unit_uac::String)
        default_parameters = Dict{String,Any}(
            "name" => "economic_control_3_states",
            "price_profile_paths" => [nothing], # Strompreis, Regelreservevergütung
            "limit_prices" => [0.0], # Grenzpreis Strom, Kosten Eigenerzeugung, Mindestwert Regelreserve
            "new_connections_below_limits" => [nothing], # (Normalbetrieb), Sparbetrieb, Einnahmenbetrieb
            "bus_uac" => unit_uac,
            "storage_uac" => nothing,
            "demand_uac" => nothing,
            "hp_uac" => nothing,
            "boiler_uac" => nothing,
            "reserve_uac" => nothing 
        )
        params = Base.merge(default_parameters, parameters)

        if params["price_profile_paths"] === nothing
            @error "Required price profile paths for control module economic_control not" * 
                   "given"
        end

        price_profiles = Array{Profile}(undef, length(params["price_profile_paths"]))
        for idx in 1:length(params["price_profile_paths"])
            price_profiles[idx] = Profile(params["price_profile_paths"][idx], sim_params)
        end

        if params["new_connections_below_limits"] === nothing
            @error "new_connections_below_limits must be given to be applied when the " *
                   "price is below the limit"
        end

        # record original connectivity
        connectivities = Dict{Int,ConnectionMatrix}()
        connectivities[1] = components[unit_uac].connectivity

        N_connections = length(params["new_connections_below_limits"]) + 1
        ooo_by_state = Dict{Int,Union{OrderOfOperations,Nothing}}()
        # walk through all possible connections
        for idx in 1:N_connections
            if idx != 1
                # calculate new connectivity
                config = Dict{String,Any}("connections" => params["new_connections_below_limits"][idx-1])   
                # create and set new connectivity
                connectivities[idx] = ConnectionMatrix(config)
            end
            # check if ooo might need to be recalculated otherwise set to nothing and fill 
            # in reorder_operations() with original ooo
            new_components = deepcopy(components)
            new_components[unit_uac].connectivity = connectivities[idx]
            # follow steps load_components() in project_loading to make sure new order of
            # operations gets calculated correctly
            Resie.reorder_interfaces_of_bus!(new_components[unit_uac])
            initialise_components(new_components, sim_params)
            chains = Resie.find_chains(values(new_components), sf_bus)
            merge_bus_chains(chains, new_components, sim_params)

            # calculate new order_of_operations to have it available for later use
            ooo_by_state[idx] = Resie.calculate_order_of_operations(new_components)
        end
        
        # N_new_connections = length(params["new_connections_below_limits"])
        # for idx in 1:N_new_connections
        #     # calculate new connectivity
        #     config = Dict{String,Any}(
        #         "connections" => params["new_connections_below_limits"][idx],
        #     )
        #     # create and set new connectivity
        #     connectivities[idx+1] = ConnectionMatrix(config)
        #     # check if ooo might need to be recalculated otherwise set to nothing and fill 
        #     # in reorder_operations() with original ooo
        #     if connectivities[idx+1].input_order != connectivities[1].input_order || 
        #        connectivities[idx+1].output_order != connectivities[1].output_order
        #        # end of expression
        #         new_components = deepcopy(components)
        #         new_components[unit_uac].connectivity = connectivities[idx+1]
        #         # follow steps load_components() in project_loading to make sure new order of
        #         # operations gets calculated correctly
        #         Resie.reorder_interfaces_of_bus!(new_components[unit_uac])
        #         initialise_components(new_components, sim_params)
        #         chains = Resie.find_chains(values(new_components), sf_bus)
        #         merge_bus_chains(chains, new_components, sim_params)

        #         # calculate new order_of_operations to have it available for later use
        #         ooo_by_state[idx+1] = Resie.calculate_order_of_operations(new_components)
        #     else
        #         ooo_by_state[idx+1] = ooo_by_state[1]
        #     end
        # end

        # remove new_components and chains from memory sincethey are not needed any more
        new_components = nothing
        chains = nothing
        
        ooo_by_state[N_connections + 1] = ooo_by_state[N_connections]
        connectivities[N_connections + 1] = connectivities[N_connections]

        # setup state machine for switching between states based on price
        # for now this is a simple two-state threshold switch, but it might implement a
        # more complicated calculation in the future
        function time_is_4h(sim_params::Dict{String,Any})
            return Hour(sim_params["current_date"]).value % 4 == 0 && Minute(sim_params["current_date"]).value == 0
        end

        function smaller_than_price(pp_idx::Int, lp_idx::Int, sim_params::Dict{String,Any})
            return value_at_time(price_profiles[pp_idx], sim_params) <= params["limit_prices"][lp_idx]
        end

        function future_sim(mod_params, ooo_by_state, components::Grouping, sim_params::Dict{String,Any})
            # create deep copies of simulation relevant variables
            sp = deepcopy(sim_params)
            comps = deepcopy(components)
            ops = ooo_by_state[1]

            # define how long the future simulation will be; it's for now fixed for 4h
            end_date = sim_params["current_date"] + Hour(4) - Second(sim_params["time_step_seconds"])
            sim_range = sp["current_date"]:Second(sp["time_step_seconds"]):end_date
            
            # define relevant output_keys to be able to gather data after future simulation
            output_keys_dict = OrderedDict{String, Any}(
                mod_params["demand_uac"] => ["Demand"],
                mod_params["storage_uac"] => ["Load", "Capacity"],
                mod_params["hp_uac"] =>  ["m_power:IN", "m_heat:OUT"]
            )
            output_data_keys = output_keys(comps, output_keys_dict)
            output_data = zeros(Float64, length(sim_range), 1 + length(output_data_keys))
            output_data_keys_string = []
            for key in output_data_keys
                if key.medium === nothing
                    push!(output_data_keys_string, key.unit.uac * key.value_key)
                else
                    push!(output_data_keys_string, key.unit.uac * String(key.medium) * key.value_key)
                end
            end
            # output_data = OrderedDict()

            # excecute future simulation and save relevant data
            for (step, ts_ahead) in enumerate(sim_range)
                sp["current_date"] = ts_ahead
                perform_operations(comps, ops, sp)
                output_data[step, :] = gather_output_data(output_data_keys, sp["time_since_output"])
                sp["time_since_output"] += Int(sp["time_step_seconds"])
            end
            output_data = OrderedDict(zip(output_data_keys_string, eachcol(output_data)))
            @infiltrate
            # logic to evaluate if reserve control energy can be used
            storage = comps[mod_params["storage_uac"]]
            hp = comps[mod_params["storage_uac"]]
            available_storage_capacity = minimum(output_data[storage.uac * "Capacity"]) -
                                         maximum(output_data[storage.uac * "Load"]) 
            max_storage_charge = storage.max_load_rate * minimum(output_data[storage.uac * "Capacity"])
            #TODO caluclate on timestep basis and find min between storage_power and th_power on timestep basis
            available_storage_power = min(available_storage_capacity, max_storage_charge)
            available_th_power = maximum(hp.plr)hp.uac * "Effective_COP", hp.uac * "Avg_PLR", hp.uac * "Time_active"
            available_el_power = hp.uac * "Effective_COP"
            energy_indicator = true

            # cleanup sim relevant variables to save memory; TODO might not be significant
            # sp = nothing
            # comps = nothing
            # ops = nothing

            price_indicator = smaller_than_price(2, 3, sim_params)

            return energy_indicator && price_indicator
        end

        function requirements_balance(mod_params, 
                                      components::Grouping, 
                                      sim_params::Dict{String,Any})
            # get the future demand for the next 4h
            sp = deepcopy(sim_params)
            # sp["current_date"] = deepcopy(sim_params["current_date"])
            end_date = sim_params["current_date"] + Hour(4) - Second(sim_params["time_step_seconds"])
            demands_ahead_power = []
            demand = components[mod_params["demand_uac"]]
            for ts_ahead in sim_params["current_date"]:Second(sim_params["time_step_seconds"]):end_date
                sp["current_date"] = ts_ahead
                power = power_at_time(demand.energy_profile, sp) * demand.scaling_factor
                push!(demands_ahead_power, power)
            end
            max_power_demand = maximum(demands_ahead_power)

            # estimate the available storage capacity ignoring losses to make sure to be able
            # to take up the energy
            storage = components[mod_params["storage_uac"]]
            available_storage_capacity = storage.capacity - storage.load
            # check maximum continous power that can be taken by the storage for worst case 
            # of 1h continous reserve demand
            available_storage_power = min(available_storage_capacity, storage.max_load_rate * storage.capacity)
            
            hp = components[mod_params["hp_uac"]]
            boiler = components[mod_params["boiler_uac"]]
            available_th_power = hp.design_power_th + boiler.design_power_th - max_power_demand
            available_th_power = min(available_th_power, available_storage_power)
            # check if heat pump has capacity for providing control reserve
            if max_power_demand >= hp.design_power_th # only boiler used for control reserve
                available_el_power = available_th_power
            elseif available_th_power > 0 # boiler and heat pump used for control reserve
                available_power_hp = hp.design_power_th - max_power_demand
                cop = hp.dynamic_cop(30, demand.temperature) * hp.plf_function(1.0) #TODO variable input_temp for cop calculation
                available_el_power_hp = available_power_hp / cop
                available_el_power = available_el_power_hp + boiler.design_power_th
            else
                available_el_power = 0.0
            end

            #
            # available_el_power - pv_ertrag - wind_ertrag #TODO integrate 

            # check if minimum el power of 1 MW can be provided
            if available_el_power/10^6 >= 1.0
                energy_indicator = true
                # scale controlreserve profile with available_energy 
                components[mod_params["reserve_uac"]].scaling_factor = available_el_power/10^6 # TODO check basic scaling_factor
            else
                energy_indicator = false
            end
            price_indicator = smaller_than_price(2, 3, sim_params)

            return energy_indicator && price_indicator
        end
        base_conds = [function (state_machine)
                        return time_is_4h(sim_params)
                     end
                     function (state_machine)
                        return smaller_than_price(1, 1, sim_params) # Netzpreis < Grenzwert
                     end
                     function (state_machine)
                        return smaller_than_price(1, 2, sim_params) # Netzpreis < Eigenerzeugung
                     end
                     ]
        base_table = Dict{Tuple,UInt}(
                                    (true ,true ,true ) => 5,
                                    (true ,true ,false) => 5,
                                    (true ,false,true ) => 5,
                                    (true ,false,false) => 5,
                                    (false,true ,true ) => 3,
                                    (false,true ,false) => 1,
                                    (false,false,true ) => 2,
                                    (false,false,false) => 2
                                    )

        table_state_basic = TruthTable(;  # State: basic
                                conditions=base_conds,
                                table_data=base_table
                                )
        table_state_save = TruthTable(;  # State: Saving operation
                                conditions=base_conds,
                                table_data=base_table
                                )
        table_state_earn = TruthTable(;  # State: Earning operation
                                conditions=base_conds,
                                table_data=base_table
                                )
        table_state_earn_reserve = TruthTable(;  # State: Earning from reserve operation
                        conditions=[function (state_machine)
                                        return state_machine.time_in_state *
                                               sim_params["time_step_seconds"] >=
                                               240 # stay 4h in earn_reserve operation
                                    end
                                    function (state_machine)
                                        return smaller_than_price(1, 1, sim_params)
                                    end
                                    function (state_machine)
                                        return smaller_than_price(1, 2, sim_params)  
                                    end
                                    ],
                        table_data=Dict{Tuple,UInt}(
                                    (true ,true ,true ) => 3,
                                    (true ,true ,false) => 1,
                                    (true ,false,true ) => 2,
                                    (true ,false,false) => 2,
                                    (false,true ,true ) => 4,
                                    (false,true ,false) => 4,
                                    (false,false,true ) => 4,
                                    (false,false,false) => 4
                        ))
        table_state_4h_check = TruthTable(;  # State: Check at 4h mark
                        conditions=[# function (state_machine)
                                    #     return future_sim(params, ooo_by_state, components, sim_params)
                                    # end
                                    function (state_machine)
                                        return requirements_balance(params, components, sim_params)
                                    end
                                    function (state_machine)
                                        return smaller_than_price(1, 1, sim_params)
                                    end
                                    function (state_machine)
                                        return smaller_than_price(1, 2, sim_params)  
                                    end
                                    ],
                        table_data=Dict{Tuple,UInt}(
                                    (true ,true ,true ) => 4,
                                    (true ,true ,false) => 4,
                                    (true ,false,true ) => 4,
                                    (true ,false,false) => 4,
                                    (false,true ,true ) => 3,
                                    (false,true ,false) => 1,
                                    (false,false,true ) => 2,
                                    (false,false,false) => 2
                        ))

        state_machine = StateMachine(UInt(1),               # state
                                     Dict{UInt,String}(     # state_names
                                         1 => "Basic",
                                         2 => "Save",
                                         3 => "Earn",
                                         4 => "EarnReserve",
                                         5 => "4hCheck"
                                     ),
                                     Dict{UInt,TruthTable}( # transitions
                                         1 => table_state_basic,
                                         2 => table_state_save,
                                         3 => table_state_earn,
                                         4 => table_state_earn_reserve,
                                         5 => table_state_4h_check
                                     ))

        return new("economic_control_3_states", params, state_machine, price_profiles,
                   connectivities, ooo_by_state)
    end
end

# method for control module name on type-level
control_module_name(x::Type{CM_EconomicControl3States})::String = "economic_control_3_states"

function has_method_for(mod::CM_EconomicControl3States, func::ControlModuleFunction)::Bool
    return func == cmf_change_bus_priorities || func == cmf_reorder_operations
end

function update(mod::CM_EconomicControl3States)
    # ordinarily we would update the state machine in the control modules's update step,
    # however this control module is special in that it's callbacks are used outside the
    # normal order of operations, so the update is performed by the callbacks instead
end

function reorder_operations(mod::CM_EconomicControl3States,
                            order_of_operations::OrderOfOperations,
                            sim_params::Dict{String,Any})::OrderOfOperations
    # record "original" ooo since it might be the first time the module has access to it
    # and also it might've changed due to other control modules
    state = mod.state_machine.state
    if state == 1 || mod.ooo_by_state[state] === nothing
        mod.ooo_by_state[state] = order_of_operations
    end

    return mod.ooo_by_state[state]
end

function change_bus_priorities!(mod::CM_EconomicControl3States,
                                components::Grouping,
                                sim_params::Dict{String,Any})
    # perform update outside of normal order of operations
    old_state = mod.state_machine.state
    move_state(mod.state_machine)
    if mod.state_machine.state == 5
        move_state(mod.state_machine)
    end

    # now reorder bus interfaces, but only if the state changed
    if old_state != mod.state_machine.state
        bus = components[mod.parameters["bus_uac"]]
        bus.connectivity = mod.connectivity_by_state[mod.state_machine.state]
        Resie.reorder_interfaces_of_bus!(bus)
        initialise!(bus, sim_params)
    end
end
