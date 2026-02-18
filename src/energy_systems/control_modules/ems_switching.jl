using Dates
using ..Resie
using OrderedCollections: OrderedDict
using Infiltrator
# using Statistics
"""
Control module for setting limits to the PLR of a component according to a price_profile and 
available storage capacity or demand.
"""
mutable struct CM_EMS_Switching <: ControlModule
    name::String
    parameters::Dict{String,Any}
    state_machine::StateMachine
    connectivity_by_state::Dict{String,Dict{Int,ConnectionMatrix}}
    ooo_by_state::Dict{Int,Union{OrderOfOperations,Nothing}}

    function CM_EMS_Switching(parameters::Dict{String,Any},
                    components::Grouping,
                    sim_params::Dict{String,Any},
                    unit_uac::String)
        default_parameters = Dict{String,Any}(
            "name" => "EMS_Switching",
            "price_profile_paths" => [], # Strompreis, neg_Regelleistungsangebot, neg_Regelenergieangebot, pos_Regelleistungsangebot, pos_Regelenergieangebot
            "limit_prices" => [0.0], # Grenzpreis Strom, Kosten Eigenerzeugung, Mindestwert Regelleistung, Mindestwert Regelenergie
            "new_connections_below_limits" => Dict{String,Any}(), # (Normalbetrieb), Sparbetrieb, Einnahmenbetrieb
            "bus_uacs" => [],
            "storage_uac" => nothing,
            "hp_uac" => nothing,
            "boiler_uac" => nothing,
            "neg_reserve_uac" => nothing,
            "pos_reserve_uac" => nothing
        )
        params = Base.merge(default_parameters, parameters)

        if isempty(params["price_profile_paths"])
            @error "Required price profile paths for control module EMS_Switching not" * 
                   "given."
        end

        price_profiles = Array{Profile}(undef, length(params["price_profile_paths"]))
        for idx in 1:length(params["price_profile_paths"])
            price_profiles[idx] = Profile(params["price_profile_paths"][idx], sim_params)
        end

        if isempty(params["new_connections_below_limits"])
            @error "new_connections_below_limits must be given for the control module EMS_Switching" *
                   " to be applied when the price is below the limit."
        end

        if isempty(params["bus_uacs"])
            @error "At least one bus_uac in bus_uacs must be given for control module EMS_Switching."
        end
        bus_uacs = params["bus_uacs"]

        for (key, value) in params
            if endswith(key, "uac") && value === nothing
                @error "The parameter $key must be given for control module EMS_Switching."
            end
        end

        # for uac in [params["hp_uac"], params["boiler_uac"]]
        #     if components[uac].controller.modules[1].name != "profile_limited"
        #         @error "The component $uac must have the control module 'profile_limited' "*
        #                "for control module EMS_Switching to be able to limit the component power for "*
        #                "postive reserve power"
        #     end
        # end

        connectivity_by_state = Dict{String,Dict{Int,ConnectionMatrix}}()
        ooo_by_state = Dict{Int,Union{OrderOfOperations,Nothing}}()
        N_connections = length(params["new_connections_below_limits"][bus_uacs[1]]) + 1

        # walk through all possible connections
        for idx in 1:N_connections
            new_components = deepcopy(components)
            for uac in bus_uacs
                if idx == 1
                    connectivity_by_state[uac] = Dict{Int,ConnectionMatrix}()
                    # # record already calculated original connectivities for each bus
                    connectivity_by_state[uac][1] = components[uac].connectivity
                else
                    # calculate new connectivity
                    config = Dict{String,Any}("connections" => params["new_connections_below_limits"][uac][idx-1])   
                    # create and set new connectivity
                    connectivity_by_state[uac][idx] = ConnectionMatrix(config)
                    new_components[uac].connectivity = connectivity_by_state[uac][idx]
                end
                Resie.reorder_interfaces_of_bus!(new_components[uac])
            end
            # follow steps load_components() in project_loading to make sure new order of
            # operations gets calculated correctly
            initialise_components(new_components, sim_params)
            chains = Resie.find_chains(values(new_components), sf_bus)
            merge_bus_chains(chains, new_components, sim_params)

            # calculate new order_of_operations to have it available for later use
            ooo_by_state[idx] = Resie.calculate_order_of_operations(new_components)
        end

        # remove new_components and chains from memory since they are not needed any more
        new_components = nothing
        chains = nothing
        
        ooo_by_state[N_connections + 1] = ooo_by_state[N_connections]
        ooo_by_state[N_connections + 2] = ooo_by_state[N_connections] # TODO change to reflect different reserves
        for uac in bus_uacs
            connectivity_by_state[uac][N_connections + 1] = connectivity_by_state[uac][N_connections]
            connectivity_by_state[uac][N_connections + 2] = connectivity_by_state[uac][N_connections] # TODO change to reflect different reserves
        end

        # setup state machine for switching between states based on price and technical limitations
        function time_is_4h(sim_params::Dict{String,Any})
            return Hour(sim_params["current_date"]).value % 4 == 0 && Minute(sim_params["current_date"]).value == 0
        end

        function smaller_than_price(pp_idx::Int, lp_idx::Int, sim_params::Dict{String,Any})
            return value_at_time(price_profiles[pp_idx], sim_params) <= params["limit_prices"][lp_idx]
        end

        # TODO for testing
        time_array_1 = []
        time_array_2 = []
        time_array_3 = []

        function calc_reserve_power(sim_length::TimePeriod,
                                    mod_params::Dict{String,Any}, 
                                    ooo_by_state::Dict{Int,Union{OrderOfOperations,Nothing}}, 
                                    components::Grouping, 
                                    sim_params::Dict{String,Any})::Bool
            power_neg, power_pos, 
            baseline_el_hp, max_el_hp, 
            baseline_el_boiler, max_el_boiler = future_sim(sim_length, mod_params, 
                                                           ooo_by_state, components, 
                                                           sim_params)

            # decision logic for amount of reserve power
            minimal_reserve_power = 1e6
            power_neg_bool = power_neg < minimal_reserve_power
            power_pos_bool = power_pos < minimal_reserve_power
            if power_neg_bool || power_pos_bool
                # TODO FUZZY LOGIC RETURNS power_value 
                # power_value = fuzzy_control()
                power_value = -1e6
                power_values = fill(power_value, length(baseline_el_hp))
                if power_value > 0
                    components[mod_params["neg_reserve_uac"]].scaling_factor = floor(power_value, digits=-6) # TODO change floor digits to bidding size of virtual power plant
                    components[mod_params["pos_reserve_uac"]].scaling_factor = 0.0
                
                elseif power_value < 0
                    components[mod_params["pos_reserve_uac"]].scaling_factor = ceil(power_value, digits=-6) # TODO change floor digits to bidding size of virtual power plant
                    components[mod_params["neg_reserve_uac"]].scaling_factor = 0.0       

                    end_date = sim_params["current_date"] + sim_length - Second(sim_params["time_step_seconds"])  
                    el_power_boiler = max.(baseline_el_boiler .- power_values, 0)
                    left_power_values = max.(power_values .- baseline_el_boiler, 0)
                    el_power_hp = max.(baseline_el_hp .- left_power_values, 0)

                    for (idx, dt) in enumerate(sim_params["current_date"]:Second(sim_params["time_step_seconds"]):end_date)
                        components[mod_params["boiler_uac"]].controller.modules[1].profile.data[dt] = el_power_boiler[idx] / max_el_boiler[idx]
                        components[mod_params["hp_uac"]].controller.modules[1].profile.data[dt] = el_power_hp[idx] / max_el_hp[idx]
                    end
                end
            end

            #TODO 
            # alle 4h wird technische Machbarkeit geprüft
            #   ergibt maximale Leistung pos und neg über ganzen Zeitraum
            #   minimum 0.5 MW
            #   Wenn nur eine Richtung möglich dann wird angeboten in die Richtung
            #   Wenn beide Richtungen dann wird Preisvergleich mit Erlösprofilen * maximale Leistung
            #       Erlösprofil = (Leistungspreis + Arbeitspreis * Abrufdauer) * time_step_seconds/3600
            #       Auf 4h Aufsummieren für Gesamterlös über 4h
            #   Wenn keine der Richtungen möglich, dann wird in jedem Zeitschritt nach freebids gearbeitet
            #       pos und neg vergleichen: Arbeitspreis * Abrufdauer für jeden Zeitschritt
            # Parameter in Inputfile ob nur freebids oder beides

            # if available_el_power_neg >= 1e6
            #     energy_indicator = true
            #     # scale controlreserve profile with available_energy to full MW values
            #     components[mod_params["neg_reserve_uac"]].scaling_factor = floor(available_el_power_neg, digits=-6)
            # else
            #     energy_indicator = false
            #     components[mod_params["neg_reserve_uac"]].scaling_factor = 1e6
            # end
            # TODO Vergleich nicht nötig, da nicht ausschlaggebend für nächsten Zustand (Regelenergie kann trotzdem angeboten werden)
            # check if Regelleistungsangebot (2) > Mindestwert Regelleistung (3)
            # price_indicator = !smaller_than_price(2, 3, sim_params)

            return power_pos_bool || power_neg_bool
        end 

        function future_sim(sim_length::TimePeriod, 
                            mod_params::Dict{String,Any}, 
                            ooo_by_state::Dict{Int,Union{OrderOfOperations,Nothing}}, 
                            components::Grouping, 
                            sim_params::Dict{String,Any}
                           )::Tuple{Float64, Float64, Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}}
            # t1 = 0.0
            # t1 = @elapsed begin
            # create deep copies of simulation relevant variables
            comps = deepcopy(components)
            # end 
            # append!(time_array_1, t1*10^6)
            
            sp = deepcopy(sim_params)
            ops = ooo_by_state[1]

            # for bus in [unit for unit in values(components) if unit.sys_function === EnergySystems.sf_bus]
            #     t2 = 0.0
            #     t2 = @elapsed begin
            #         deepcopy(bus)
            #     end 
            #     append!(time_array_2, t2*10^6)
            # end      

            # t3 = 0.0
            # t3 = @elapsed begin
            # for (uac, unit) in pairs(components)
            #     # map(field -> setfield!(comps[uac], field, getfield(unit, field)), fieldnames(typeof(unit)))
            #     for field in fieldnames(typeof(unit))
            #         if typeof(field) != Profile
            #             setfield!(comps[uac], field, getfield(unit, field))
            #         end
            #     end
            # end    
            # end 
            # append!(time_array_3, t3*10^6)  

            hp_uac = mod_params["hp_uac"]
            boiler_uac = mod_params["boiler_uac"]
            storage_uac = mod_params["storage_uac"]

            # define how long the future simulation will be
            end_date = sim_params["current_date"] + sim_length - Second(sim_params["time_step_seconds"])
            sim_range = sp["current_date"]:Second(sp["time_step_seconds"]):end_date
            
            # define relevant output_keys to be able to gather data after future simulation
            output_keys_dict = OrderedDict{String, Any}(
                storage_uac => ["Load", "Capacity"],
                hp_uac => ["m_power:IN", "m_heat:OUT", "Avg_PLR", "COP", "secondary_m_heat:OUT", 
                           "MixingTemperature_Input", "MixingTemperature_Output"],
                boiler_uac => ["m_power:IN", "m_heat:OUT", "Avg_PLR", "COP", "secondary_m_heat:OUT"]
            )
            output_data_keys = Resie.output_keys(comps, output_keys_dict)
            output_data = zeros(Float64, length(sim_range), 1 + length(output_data_keys))
            output_data_keys_string = ["time"]
            for key in output_data_keys
                if key.medium === nothing
                    push!(output_data_keys_string, key.unit.uac * key.value_key)
                else
                    push!(output_data_keys_string, key.unit.uac * String(key.medium) * key.value_key)
                end
            end


            # excecute future simulation and save relevant data
            for (step, ts_ahead) in enumerate(sim_range)
                sp["current_date"] = ts_ahead
                try
                    perform_operations(comps, ops, sp)
                catch
                    @info "Undefined Error in future_sim()"
                end
                output_data[step, :] = Resie.gather_output_data(output_data_keys, sp["time_since_output"])
                sp["time_since_output"] += Int(sp["time_step_seconds"])
            end
            output_data = OrderedDict(zip(output_data_keys_string, eachcol(output_data)))

            # calculate available el. power for negative reserve control power
            # ---------------------------------------------------------------


            # calculate the availalble storage power to be put into the storage each timestep
            available_storage_capacity = sim_params["wh_to_watts"].(output_data[storage_uac * "Capacity"] .-
                                                                    output_data[storage_uac * "Load"])
            max_storage_charge = comps[storage_uac].max_load_rate .* output_data[storage_uac * "Capacity"]
            available_storage_power_neg = min.(available_storage_capacity, max_storage_charge)

            # calculate available power from heat pump
            hp = comps[hp_uac]
            used_th_power_hp = sim_params["wh_to_watts"].(output_data[hp_uac * "m_heat" * "OUT"]) .+ 
                               sim_params["wh_to_watts"].(output_data[hp_uac * "secondary_m_heat" * "OUT"])
            available_th_power_hp_neg = zeros(length(used_th_power_hp))
            available_el_power_hp_neg = zeros(length(used_th_power_hp))
            cops_hp = zeros(length(used_th_power_hp))
            for (idx, th_power) in enumerate(used_th_power_hp)
                if th_power == 0
                    src_temp = comps[hp.input_interfaces[hp.m_heat_in].source_uac].temperature_src_in
                    snk_temp = comps[storage_uac].high_temperature
                    available_th_power_hp_neg[idx] = hp.max_power_function(src_temp, snk_temp) * hp.design_power_th
                    available_th_power_hp_neg[idx] = min(available_storage_power_neg[idx], available_th_power_hp_neg[idx])
                    cops_hp[idx] = hp.dynamic_cop(src_temp, snk_temp)
                    available_el_power_hp_neg[idx] = available_th_power_hp_neg[idx] / cops_hp[idx]
                else
                    available_th_power_hp_neg[idx] = th_power / output_data[hp_uac * "Avg_PLR"][idx] - th_power
                    available_th_power_hp_neg[idx] = min(available_storage_power_neg[idx], available_th_power_hp_neg[idx])
                    cops_hp[idx] = output_data[hp_uac * "COP"]
                    available_el_power_hp_neg[idx] = available_th_power_hp_neg[idx] / cops_hp[idx]
                end
            end

            # calculate available power from boiler 
            used_th_power_boiler = sim_params["wh_to_watts"].(output_data[boiler_uac * "m_heat" * "OUT"]) .+ 
                                   sim_params["wh_to_watts"].(output_data[boiler_uac * "secondary_m_heat" * "OUT"])
            available_th_power_boiler_neg = zeros(length(used_th_power_boiler))
            for (idx, th_power) in enumerate(used_th_power_boiler)
                if th_power == 0
                    available_th_power_boiler_neg[idx] = comps[boiler_uac].design_power_th
                else
                    # how much power is left to use
                    available_th_power_boiler_neg[idx] = th_power / output_data[boiler_uac * "Avg_PLR"][idx] - th_power
                end
            end
            available_el_power_boiler_neg = min.(available_storage_power_neg .- available_th_power_hp_neg, available_th_power_boiler_neg)

            # add the feed-in power from renewables to the reserve control power since it 
            # can be marketed additionally TODO double check if that makes sense
            available_el_power_neg = minimum(available_el_power_hp_neg .+ available_el_power_boiler_neg)
            available_el_power_neg = min(available_storage_capacity[1] / length(available_storage_capacity), available_el_power_neg)


            # Calculate available el. power for positive reserve control power
            # ----------------------------------------------------------------

            # calculate available storage power to be discharged
            available_storage_load = sim_params["wh_to_watts"].(output_data[storage_uac * "Load"])
            max_storage_discharge = comps[storage_uac].max_unload_rate .* output_data[storage_uac * "Capacity"]
            available_storage_power_pos = min.(available_storage_load, max_storage_discharge)
            
            # calculate how much the power from the boiler can be reduced 
            available_th_power_boiler_pos = zeros(length(used_th_power_boiler))
            for (idx, th_power) in enumerate(used_th_power_boiler)
                if th_power > 0
                    available_th_power_boiler_pos[idx] = min(available_storage_power_pos[idx], th_power)
                end
            end
            available_el_power_boiler_pos = available_th_power_boiler_pos

            # calculate how much the power from the hp can be reduced additionally to the boiler
            available_th_power_hp_pos = zeros(length(used_th_power_hp))
            available_el_power_hp_pos = zeros(length(used_th_power_hp))
            for (idx, th_power) in enumerate(used_th_power_hp)
                if th_power > 0
                    available_th_power_hp_pos[idx] = min(available_storage_power_pos[idx] - available_th_power_boiler_pos[idx], th_power)
                    available_el_power_hp_pos[idx] = available_th_power_hp_pos[idx] / cops_hp[idx]
                end
            end

            available_el_power_pos = minimum(available_el_power_hp_pos .+ available_el_power_boiler_pos)
            available_el_power_pos = min(available_storage_load[1] / length(available_storage_load), available_el_power_pos)

            # calculate used and maximum el power for hp and boiler for calculation of max_plr for positve control reserve
            used_el_power_boiler = used_th_power_boiler
            max_el_power_boiler = fill(comps[boiler_uac].design_power_th, length(used_el_power_boiler))

            used_el_power_hp = used_th_power_hp ./ cops_hp
            max_el_power_hp = hp.design_power_th ./ cops_hp
            
            # if sim_params["time"] >= sim_params["number_of_time_steps"] * sim_params["time_step_seconds"] * 0.99
            #     time_array_1 = time_array_1[2:end]
            #     # @info time_array_1
            #     @info "Mean t1:$(round(mean(time_array_1), digits=3)) us; std: $(round(std(time_array_1), digits=3)) us; max:$(round(maximum(time_array_1), digits=3)) us"
            #     time_array_2 = time_array_2[2:end]
            #     # @info time_array_2
            #     @info "Mean t2:$(round(mean(time_array_2), digits=3)) us; std: $(round(std(time_array_2), digits=3)) us; max:$(round(maximum(time_array_2), digits=3)) us"
            #     time_array_3 = time_array_3[2:end]
            #     # @info time_array_3
            #     @info "Mean t3:$(round(mean(time_array_3), digits=3)) us; std: $(round(std(time_array_3), digits=3)) us; max:$(round(maximum(time_array_3), digits=3)) us"
            # end
            return available_el_power_neg, available_el_power_pos, used_el_power_hp, max_el_power_hp, used_el_power_boiler, max_el_power_boiler
        end

        base_conds = [function (state_machine)
                        return time_is_4h(sim_params)
                     end
                     function (state_machine)
                        # Stormpreis <= Grenzpreis Strom
                        # Strom sehr günstig
                        return smaller_than_price(1, 1, sim_params)
                     end
                     function (state_machine)
                        # Strompreis <= Kosten Eigenerzeugung
                        # Netzbezug günstiger als Eigenerzeugung TODO immer false, da Eigenerzeugung immer verbraucht werden soll
                        return smaller_than_price(1, 2, sim_params) # Netzpreis < Eigenerzeugung
                     end
                     ]
        base_table = Dict{Tuple,UInt}(
                                    (true ,true ,true ) => 6,
                                    (true ,true ,false) => 6,
                                    (true ,false,true ) => 6,
                                    (true ,false,false) => 6,
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
        table_state_earn_reserve_pos = TruthTable(;  # State: Earning from reserve operation
                conditions=[function (state_machine)
                                return state_machine.time_in_state *
                                        sim_params["time_step_seconds"] >=
                                        240 # stay 4h in earn_reserve operation
                            end
                            function (state_machine)
                                # Strompreis <= Grenzpreis Strom OR Regelenergieangebot > Mindestpreis Regelenergie
                                # Strom günstig OR Regelenergie lohnend
                                return smaller_than_price(1, 1, sim_params) || !smaller_than_price(3, 4, sim_params)
                            end
                            function (state_machine)
                                # Strompreis <= Kosten Eigenerzeugung
                                # Netzbezug günstiger als Eigenerzeugung TODO immer false, da Eigenerzeugung immer verbraucht werden soll
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
        table_state_earn_reserve_neg = TruthTable(;  # State: Earning from reserve operation
                conditions=[function (state_machine)
                                return state_machine.time_in_state *
                                        sim_params["time_step_seconds"] >=
                                        240 # stay 4h in earn_reserve operation
                            end
                            function (state_machine)
                                # Strompreis <= Grenzpreis Strom OR Regelenergieangebot > Mindestpreis Regelenergie
                                # Strom günstig OR Regelenergie lohnend
                                return smaller_than_price(1, 1, sim_params) || !smaller_than_price(3, 4, sim_params)
                            end
                            function (state_machine)
                                # Strompreis <= Kosten Eigenerzeugung
                                # Netzbezug günstiger als Eigenerzeugung TODO immer false, da Eigenerzeugung immer verbraucht werden soll
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
                conditions=[
                            function (state_machine)
                                return calc_reserve_power(Hour(4), params, ooo_by_state, components, sim_params)
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
                                         4 => "EarnReservePos",
                                         5 => "EarnReserveNeg",
                                         6 => "4hCheck"
                                     ),
                                     Dict{UInt,TruthTable}( # transitions
                                         1 => table_state_basic,
                                         2 => table_state_save,
                                         3 => table_state_earn,
                                         4 => table_state_earn_reserve_pos,
                                         5 => table_state_earn_reserve_neg,
                                         6 => table_state_4h_check
                                     ))

        return new("EMS_Switching", params, state_machine, connectivity_by_state, ooo_by_state)
    end
end

# method for control module name on type-level
control_module_name(x::Type{CM_EMS_Switching})::String = "EMS_Switching"

function has_method_for(mod::CM_EMS_Switching, func::ControlModuleFunction)::Bool
    return func == cmf_change_bus_priorities || func == cmf_reorder_operations
end

function update(mod::CM_EMS_Switching)
    # ordinarily we would update the state machine in the control modules's update step,
    # however this control module is special in that it's callbacks are used outside the
    # normal order of operations, so the update is performed by the callbacks instead
end

function reorder_operations(mod::CM_EMS_Switching,
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

function change_bus_priorities!(mod::CM_EMS_Switching,
                                components::Grouping,
                                sim_params::Dict{String,Any})
    # perform update outside of normal order of operations
    old_state = mod.state_machine.state
    if sim_params["time"] > 0
        move_state(mod.state_machine)
        if mod.state_machine.state == length(mod.state_machine.transitions)
            move_state(mod.state_machine)
        end
    end

    # now reorder bus interfaces, but only if the state changed
    if old_state != mod.state_machine.state
        for bus_uac in mod.parameters["bus_uacs"]
            bus = components[bus_uac]
            bus.connectivity = mod.connectivity_by_state[bus_uac][mod.state_machine.state]
            Resie.reorder_interfaces_of_bus!(bus)
            initialise!(bus, sim_params)
        end
    end
end
