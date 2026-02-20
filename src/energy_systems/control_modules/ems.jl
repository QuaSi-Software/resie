using Dates
using ..Resie
using OrderedCollections: OrderedDict

"""
Control module for setting limits to the PLR of a component according to a price_profiles and 
available storage capacity or demand.
"""
mutable struct CM_EMS <: ControlModule
    name::String
    parameters::Dict{String,Any}
    state_machine::StateMachine
    connectivity_by_state::Dict{String,Dict{UInt,ConnectionMatrix}}
    ooo_by_state::Dict{UInt,Union{OrderOfOperations,Nothing}}

    function CM_EMS(parameters::Dict{String,Any},
                    components::Grouping,
                    sim_params::Dict{String,Any},
                    unit_uac::String)
        default_parameters = Dict{String,Any}(
            "name" => "EMS",
            "grid_price_path" => nothing,               # Strompreis, 
            "reserve_power_price_neg_path" => nothing,  # neg_Regelleistungsangebot,
            "reserve_energy_price_neg_path" => nothing, # neg_Regelenergieangebot
            "reserve_call_neg_path" => nothing,         # neg_Abrufdauer
            "reserve_power_price_pos_path" => nothing,  # pos_Regelleistungsangebot,
            "reserve_energy_price_pos_path" => nothing, # pos_Regelenergieangebot
            "reserve_call_pos_path" => nothing,         # pos_Abrufdauer   
            "min_reserve_power" => 0.0,
            "offer_only_freebids" => false,
            "min_storage_load" => 0.0,
            "new_connections_below_limits" => Dict{String,Any}(), # (Normalbetrieb), Sparbetrieb, Einnahmenbetrieb
            "bus_uacs" => [],
            "storage_uac" => nothing,
            "hp_uac" => nothing,
            "boiler_uac" => nothing,
            "neg_reserve_uac" => nothing,
            "pos_reserve_uac" => nothing
        )
        params = Base.merge(default_parameters, parameters)

        price_profile_paths = [params["grid_price_path"], 
                               params["reserve_power_price_neg_path"], 
                               params["reserve_energy_price_neg_path"],
                               params["reserve_call_neg_path"], 
                               params["reserve_power_price_pos_path"],
                               params["reserve_energy_price_pos_path"], 
                               params["reserve_call_pos_path"]
                              ]

        if all(isnothing.(price_profile_paths))
            @error "Required price profile paths for control module EMS not" * 
                   "given."
        end

        params["price_profiles"] = Array{Profile}(undef, length(price_profile_paths))
        for idx in eachindex(price_profile_paths)
            params["price_profiles"][idx] = Profile(price_profile_paths[idx], sim_params)
        end

        if xor(isempty(params["new_connections_below_limits"]), isempty(params["bus_uacs"]))
            @error "bus_uacs and new_connections_below_limits must be given together for" * 
                   " the control module EMS."
        end

        for (key, value) in params
            if endswith(key, "uac") && value === nothing
                @error "The parameter $key must be given for control module EMS."
            end
        end

        # for uac in [params["hp_uac"], params["boiler_uac"]]
        #     if components[uac].controller.modules[1].name != "profile_limited"
        #         @error "The component $uac must have the control module 'profile_limited' "*
        #                "for control module EMS to be able to limit the component power for "*
        #                "postive and negative reserve power"
        #     end
        # end

        if isempty(params["bus_uacs"])
            bus_uacs = [unit_uac]
            params["new_connections_below_limits"][unit_uac] = []
        else
            bus_uacs = params["bus_uacs"]
        end

        connectivity_by_state = Dict{String,Dict{UInt,ConnectionMatrix}}()
        ooo_by_state = Dict{UInt,Union{OrderOfOperations,Nothing}}()
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

        # setup state machine for switching between states based on price and technical limitations
        function time_is_4h(sim_params::Dict{String,Any})
            return Hour(sim_params["current_date"]).value % 4 == 0 && Minute(sim_params["current_date"]).value == 0
        end

        # State: Base checks if time is at the 4h mark
        table_state_base =TruthTable(;
            conditions=[
                function (state_machine)
                    return time_is_4h(sim_params)
                end
                function (state_machine)
                    return params["offer_only_freebids"]
                end
            ],
            table_data=Dict{Tuple, UInt}(
                (true, true) => 2,
                (true, false) => 3,
                (false, true) => 2,
                (false, false) => 1
            )
        )

        # State: Offer reserve control as freebids in each time step
        table_state_freebids =TruthTable(;
            conditions=[
                function (state_machine)
                    return time_is_4h(sim_params)
                end
                function (state_machine)
                    return params["offer_only_freebids"]
                end
            ],
            table_data=Dict{Tuple, UInt}(
                (true, true) => 2,
                (true, false) => 3,
                (false, true) => 2,
                (false, false) => 2
            )
        )

        # State: Check possible reserve power at 4h mark
        table_state_4h_check = TruthTable(;  
            conditions=[
                function (state_machine)
                    return calc_reserve_power(Hour(4), params, ooo_by_state, components, 
                                              sim_params)
                end
            ],
            table_data=Dict{Tuple,UInt}(
                (true,) => 1,
                (false,) => 2
            )
        )

        state_machine = StateMachine(UInt(1),               # state
                                     Dict{UInt,String}(     # state_names
                                        1 => "Base",
                                        2 => "Freebids",
                                        3 => "4hCheck"
                                     ),
                                     Dict{UInt,TruthTable}( # transitions
                                        1 => table_state_base,
                                        2 => table_state_freebids,
                                        3 => table_state_4h_check
                                     ))
        # fill the connectivity and ooo to match the amount of states
        for key in keys(state_machine.transitions)
            if !haskey(ooo_by_state, key)
                ooo_by_state[key] = ooo_by_state[key-1]
                if N_connections > 1
                    for uac in bus_uacs
                        connectivity_by_state[uac][key] = connectivity_by_state[uac][key-1]
                    end
                end
            end
        end

        return new("EMS", params, state_machine, connectivity_by_state, ooo_by_state)
    end
end

function calc_reserve_power(sim_length::TimePeriod,
                            mod_params::Dict{String,Any}, 
                            ooo_by_state::Dict{UInt,Union{OrderOfOperations,Nothing}}, 
                            components::Grouping, 
                            sim_params::Dict{String,Any})::Bool

    power_neg, power_pos, 
    baseline_el_hp, max_el_hp, 
    baseline_el_boiler, max_el_boiler = future_sim(sim_length, mod_params, 
                                                    ooo_by_state, components, 
                                                    sim_params)

    # decision logic for amount of reserve power
    power_neg_bool = power_neg > mod_params["min_reserve_power"]
    power_pos_bool = power_pos > mod_params["min_reserve_power"]        
    if power_neg_bool || power_pos_bool
        end_date = min(add_ignoring_leap_days(sim_params["current_date"], 
                                              Second(sim_length) - Second(sim_params["time_step_seconds"])),
                       sim_params["end_date"])
        date_range = remove_leap_days(collect(range(sim_params["current_date"]; 
                                                    stop=end_date, 
                                                    step=Second(sim_params["time_step_seconds"])
                                                    )))
        price_profiles = mod_params["price_profiles"]

        if power_neg_bool && power_pos_bool
            # decide which direction to offer based on the earnings
            # calculate total revenue over length of future simulation
            rev_neg = 0
            rev_pos = 0
            if length(date_range) > 1
                for dt in date_range
                    rev_neg += (price_profiles[2].data[dt] + price_profiles[3].data[dt] * price_profiles[4].data[dt]) * 
                               sim_params["time_step_seconds"]/3600
                    rev_pos += (price_profiles[5].data[dt] + price_profiles[6].data[dt] * price_profiles[7].data[dt]) * 
                               sim_params["time_step_seconds"]/3600
                end
            else
                rev_neg = value_at_time(price_profiles[3], sim_params) * 
                          value_at_time(price_profiles[4], sim_params) * 
                          sim_params["time_step_seconds"]/3600
                rev_pos = value_at_time(price_profiles[3], sim_params) * 
                          value_at_time(price_profiles[4], sim_params) * 
                          sim_params["time_step_seconds"]/3600
            end
            power_neg_bool = rev_neg * power_neg_bool >= rev_pos * power_pos_bool
            power_pos_bool = !power_pos_bool
        end

        # control reserve realted values are always positive for ease of use and handle pos
        # and neg differntiation with power_neg_bool and power_pos_bool
        power_value = power_neg * power_neg_bool + power_pos * power_pos_bool
        energy_weighted_values = Array{Float64}(undef, length(date_range))
        for (idx, dt) in enumerate(date_range)
            if power_neg_bool
                energy_weighted_values[idx] = power_value * price_profiles[4].data[dt]
            elseif power_pos_bool
                energy_weighted_values[idx] = power_value * price_profiles[7].data[dt]
            end
        end

        # negative control reserve
        if power_neg_bool
            components[mod_params["neg_reserve_uac"]].scaling_factor = floor(power_value, digits=-6)
            components[mod_params["pos_reserve_uac"]].scaling_factor = 0.0
            
            #set plr_limit to baseline + offered negative control reserve
            el_power_hp = max.(baseline_el_hp .+ energy_weighted_values, max_el_hp)
            left_energy_weighted_values = energy_weighted_values .- (el_power_hp .- baseline_el_hp)
            el_power_boiler = max.(baseline_el_boiler .+ left_energy_weighted_values, max_el_boiler)

            for (idx, dt) in enumerate(date_range)
                components[mod_params["boiler_uac"]].controller.modules[1].profile.data[dt] = el_power_boiler[idx] / max_el_boiler[idx]
                components[mod_params["hp_uac"]].controller.modules[1].profile.data[dt] = el_power_hp[idx] / max_el_hp[idx]
            end
        # positive control reserve
        elseif power_pos_bool
            components[mod_params["pos_reserve_uac"]].scaling_factor = ceil(power_value, digits=-6)
            components[mod_params["neg_reserve_uac"]].scaling_factor = 0.0       
            
            #set plr_limit to baseline - offered positive control reserve
            el_power_boiler = max.(baseline_el_boiler .- energy_weighted_values, 0)
            left_energy_weighted_values = max.(energy_weighted_values .- baseline_el_boiler, 0)
            el_power_hp = max.(baseline_el_hp .- left_energy_weighted_values, 0)

            for (idx, dt) in enumerate(date_range)
                components[mod_params["boiler_uac"]].controller.modules[1].profile.data[dt] = el_power_boiler[idx] / max_el_boiler[idx]
                components[mod_params["hp_uac"]].controller.modules[1].profile.data[dt] = el_power_hp[idx] / max_el_hp[idx]
            end
        end
    end


    return power_pos_bool || power_neg_bool
end 

function future_sim(sim_length::TimePeriod, 
                    mod_params::Dict{String,Any}, 
                    ooo_by_state::Dict{UInt,Union{OrderOfOperations,Nothing}}, 
                    components::Grouping, 
                    sim_params::Dict{String,Any}
                    )::Tuple{Float64, Float64, Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}}

    comps = deepcopy(components)           
    sp = deepcopy(sim_params)
    ops = ooo_by_state[1]

    hp_uac = mod_params["hp_uac"]
    boiler_uac = mod_params["boiler_uac"]
    storage_uac = mod_params["storage_uac"]

    # define how long the future simulation will be
    end_date = min(add_ignoring_leap_days(sp["current_date"], Second(sim_length) - Second(sp["time_step_seconds"])),
                   sim_params["end_date"])
    sim_range = remove_leap_days(collect(range(sp["current_date"]; stop=end_date, step=Second(sp["time_step_seconds"]))))

    # disallow charging of storage for baseline simulation
    comps[storage_uac].controller.modules[1].parameters["charge_is_allowed"] = false
    
    # define relevant output_keys to be able to gather data after future simulation
    if comps[hp_uac].has_secondary_interface
        output_keys_dict = OrderedDict{String, Any}(
            storage_uac => ["Load", "Capacity"],
            hp_uac => ["m_power:IN", "m_heat:OUT", "Avg_PLR", "COP", "secondary_m_heat:OUT", 
                        "MixingTemperature_Input", "MixingTemperature_Output"],
            boiler_uac => ["m_power:IN", "m_heat:OUT", "Avg_PLR", "COP", "secondary_m_heat:OUT"],
        )
    else
        output_keys_dict = OrderedDict{String, Any}(
            storage_uac => ["Load", "Capacity"],
            hp_uac => ["m_power:IN", "m_heat:OUT", "Avg_PLR", "COP", 
                        "MixingTemperature_Input", "MixingTemperature_Output"],
            boiler_uac => ["m_power:IN", "m_heat:OUT", "Avg_PLR", "COP"],
        )
    end
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
    # shift BufferTank Load by one timestep since we need to know the load at the beginning
    # of each time_step; set first value to value before simulation
    storage_load_idx = findfirst(isequal(storage_uac * "Load"), output_data_keys_string)
    output_data[:, storage_load_idx] = circshift(output_data[:, storage_load_idx], 1)
    output_data[1, storage_load_idx] = components[storage_uac].load

    output_data = OrderedDict(zip(output_data_keys_string, eachcol(output_data)))   

    # calculate available el. power for negative reserve control power
    # ---------------------------------------------------------------


    # calculate the availalble storage power to be put into the storage each timestep
    available_storage_capacity = sim_params["wh_to_watts"].(output_data[storage_uac * "Capacity"] .-
                                                            output_data[storage_uac * "Load"])
    max_storage_charge = comps[storage_uac].max_load_rate .* output_data[storage_uac * "Capacity"]
    available_storage_power_neg = min.(available_storage_capacity[1] / length(available_storage_capacity), max_storage_charge)

    # calculate available power from heat pump
    hp = comps[hp_uac]
    if hp.has_secondary_interface
        used_th_power_hp = sim_params["wh_to_watts"].(output_data[hp_uac * "m_heat" * "OUT"]) .+ 
                            sim_params["wh_to_watts"].(output_data[hp_uac * "secondary_m_heat" * "OUT"])
    else
        used_th_power_hp = sim_params["wh_to_watts"].(output_data[hp_uac * "m_heat" * "OUT"])
    end
    available_th_power_hp_neg = zeros(length(used_th_power_hp))
    available_el_power_hp_neg = zeros(length(used_th_power_hp))
    cops_hp = zeros(length(used_th_power_hp))
    for (idx, th_power) in enumerate(used_th_power_hp)
        if th_power == 0
            src_temp = comps[hp.input_interfaces[hp.m_heat_in].source.uac].temperature_snk_out
            snk_temp = comps[storage_uac].high_temperature
            available_th_power_hp_neg[idx] = hp.max_power_function(src_temp, snk_temp) * hp.design_power_th
            available_th_power_hp_neg[idx] = min(available_storage_power_neg[idx], 
                                                    available_th_power_hp_neg[idx])
            cops_hp[idx] = hp.dynamic_cop(src_temp, snk_temp)
            available_el_power_hp_neg[idx] = available_th_power_hp_neg[idx] / cops_hp[idx]
        else
            available_th_power_hp_neg[idx] = th_power / output_data[hp_uac * "Avg_PLR"][idx] - th_power
            available_th_power_hp_neg[idx] = min(available_storage_power_neg[idx], 
                                                    available_th_power_hp_neg[idx])
            cops_hp[idx] = output_data[hp_uac * "COP"][idx]
            available_el_power_hp_neg[idx] = available_th_power_hp_neg[idx] / cops_hp[idx]
        end
    end

    # calculate available power from boiler 
    if comps[boiler_uac].has_secondary_interface
        used_th_power_boiler = sim_params["wh_to_watts"].(output_data[boiler_uac * "m_heat" * "OUT"]) .+ 
                                sim_params["wh_to_watts"].(output_data[boiler_uac * "secondary_m_heat" * "OUT"])
    else
        used_th_power_boiler = sim_params["wh_to_watts"].(output_data[boiler_uac * "m_heat" * "OUT"])
    end
    available_th_power_boiler_neg = zeros(length(used_th_power_boiler))
    for (idx, th_power) in enumerate(used_th_power_boiler)
        if th_power == 0
            available_th_power_boiler_neg[idx] = comps[boiler_uac].design_power_th
        else
            # how much power is left to use
            available_th_power_boiler_neg[idx] = th_power / output_data[boiler_uac * "Avg_PLR"][idx] - th_power
        end
    end
    available_el_power_boiler_neg = min.(available_storage_power_neg .- available_th_power_hp_neg, 
                                            available_th_power_boiler_neg)

    # calculate total marketable negative control reserve
    available_el_power_neg = minimum(available_el_power_hp_neg .+ available_el_power_boiler_neg)
    # available_el_power_neg = min(available_storage_capacity[1] / length(available_storage_capacity), 
    #                                 available_el_power_neg)


    # Calculate available el. power for positive reserve control power
    # ----------------------------------------------------------------

    # calculate available storage power to be discharged
    available_storage_load = sim_params["wh_to_watts"].(output_data[storage_uac * "Load"] .- 
                             mod_params["min_storage_load"] .* output_data[storage_uac * "Capacity"])
    max_storage_discharge = comps[storage_uac].max_unload_rate .* output_data[storage_uac * "Capacity"]
    available_storage_power_pos = min.(available_storage_load[1] / length(available_storage_load), max_storage_discharge)
    
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
            available_th_power_hp_pos[idx] = min(available_storage_power_pos[idx] - available_th_power_boiler_pos[idx], 
                                                    th_power)
            available_el_power_hp_pos[idx] = available_th_power_hp_pos[idx] / cops_hp[idx]
        end
    end

    available_el_power_pos = minimum(available_el_power_hp_pos .+ available_el_power_boiler_pos)
    # available_el_power_pos = min(available_storage_load[1] / length(available_storage_load), 
    #                                 available_el_power_pos)

    # calculate used and maximum el power for hp and boiler for calculation of 
    # max_plr for positve control reserve
    used_el_power_boiler = used_th_power_boiler
    max_el_power_boiler = fill(comps[boiler_uac].design_power_th, 
                                length(used_el_power_boiler))

    used_el_power_hp = used_th_power_hp ./ cops_hp
    max_el_power_hp = hp.design_power_th ./ cops_hp

    return available_el_power_neg, available_el_power_pos, used_el_power_hp, 
            max_el_power_hp, used_el_power_boiler, max_el_power_boiler
end

# method for control module name on type-level
control_module_name(x::Type{CM_EMS})::String = "EMS"

function has_method_for(mod::CM_EMS, func::ControlModuleFunction)::Bool
    return func == cmf_change_bus_priorities || func == cmf_reorder_operations
end

function update(mod::CM_EMS)
    # ordinarily we would update the state machine in the control modules's update step,
    # however this control module is special in that it's callbacks are used outside the
    # normal order of operations, so the update is performed by the callbacks instead
end

function reorder_operations(mod::CM_EMS,
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

function change_bus_priorities!(mod::CM_EMS,
                                components::Grouping,
                                sim_params::Dict{String,Any})
    # perform update outside of normal order of operations
    old_state = mod.state_machine.state
    
    move_state(mod.state_machine)
    if mod.state_machine.state == length(mod.state_machine.transitions)
        move_state(mod.state_machine)
    end

    if mod.state_machine.state == 2
        calc_reserve_power(Second(sim_params["time_step_seconds"]), mod.parameters, 
                           mod.ooo_by_state, components, sim_params)
    end

    # now reorder bus interfaces, but only if the state changed and the control module 
    # actually has multiple connectivities
    if old_state != mod.state_machine.state && length(mod.connectivity_by_state) > 1
        for bus_uac in mod.parameters["bus_uacs"]
            bus = components[bus_uac]
            bus.connectivity = mod.connectivity_by_state[bus_uac][mod.state_machine.state]
            Resie.reorder_interfaces_of_bus!(bus)
            initialise!(bus, sim_params)
        end
    end
end

function Base.deepcopy_internal(x::Profile, stackdict::IdDict)
    if haskey(stackdict, x)
        return stackdict[x]::typeof(x)
    end
    stackdict[x] = x
    return x
end
