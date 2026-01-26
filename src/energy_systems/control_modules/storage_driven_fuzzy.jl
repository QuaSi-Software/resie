"""
Control module for running a component depending on the state of a linked storage component.
In particular it switches to a state of allowing operation of the component when the load
of the linked storage falls below the lower threshold. The module stays in this state until
the load has reached the upper threshold and the minimum run time has passed.
"""

mutable struct CM_StorageDrivenFuzzy <: ControlModule
    name::String
    parameters::Dict{String,Any}
    state_machine::StateMachine

    function CM_StorageDrivenFuzzy(parameters::Dict{String,Any},
                                   components::Grouping,
                                   sim_params::Dict{String,Any},
                                   unit_uac::String)
        default_parameters = Dict{String,Any}(
            "name" => "storage_driven_fuzzy",
            "low_threshold" => 0.2,
            "high_threshold" => 0.95,
            "min_run_time" => 1800,
            "storage_uac" => nothing,
            "price_profile_path" => nothing
        )
        params = Base.merge(default_parameters, parameters)

        if !(params["storage_uac"] !== nothing
             && params["storage_uac"] in keys(components)
             && components[params["storage_uac"]] isa StorageComponent)
            @error "Required storage component for control module storage_driven_fuzzy not given"
        end
        params["storage"] = components[params["storage_uac"]]

        params["price_profile"] = Profile(params["price_profile_path"], sim_params)

        function run_fuzzy!(mod_params::Dict{String,Any}, sim_params::Dict{String,Any})    
            # get input values for Fuzzy control
            p_now = value_at_time(mod_params["price_profile"], sim_params)
            p_future = value_at_time(mod_params["price_profile"], sim_params)
            p_stability = value_at_time(mod_params["price_profile"], sim_params) / value_at_time(mod_params["price_profile"], sim_params)
            SOC = 100 * mod_params["storage"].load_end_of_last_timestep / mod_params["storage"].capacity

            # TODO FUZZY LOGIC RETURNS SOC_target 
            SOC_target = fuzzy_control(p_now, p_future, p_stability, SOC)
            if SOC_target <= SOC
                # change low_threshold of storage_driven_fuzzy controller
                mod_params["low_threshold"] = SOC_target
                mod_params["high_threshold"] = SOC
            else
                # change high_threshold of storage_driven_fuzzy controller
                mod_params["low_threshold"] = SOC
                mod_params["high_threshold"] = SOC_target
            end
        end 

        function fuzzy_control(p_now, p_future, p_stability, SOC_now)
            cheap = max(min((p_now - -100) / 50, (450 - p_now) / 500), 0)
            average = max(min((p_now - -50) / 500, (950 - p_now) / 500), 0)
            expensive = max(min((p_now - 450) / 500, (1500 - p_now) / 550), 0)
            cheap = max(min((p_future - -100) / 50, (450 - p_future) / 500), 0)
            average = max(min((p_future - -50) / 500, (950 - p_future) / 500), 0)
            expensive = max(min((p_future - 450) / 500, (1500 - p_future) / 550), 0)
            low = max(min((p_stability - -5) / 5, (5 - p_stability) / 5), 0)
            average = max(min((p_stability - 0) / 5, (10 - p_stability) / 5), 0)
            high = max(min((p_stability - 5) / 5, (15 - p_stability) / 5), 0)
            low = max(min((SOC_now - -1.0) / 1.0, (0.5 - SOC_now) / 0.5), 0)
            middle = max(min((SOC_now - 0.0) / 0.5, (1.0 - SOC_now) / 0.5), 0)
            high = max(min((SOC_now - 0.5) / 0.5, (2.0 - SOC_now) / 1.0), 0)
            ant1 = max(cheap, high)
            ant2 = max(expensive, high)
            ant3 = max(low, low)
            SOC_target_agg = collect(LinRange{Float64}(0.0, 1.0, 101))
            @inbounds for (i, x) = enumerate(SOC_target_agg)
                    low = max(min((x - -1.0) / 1.0, (0.5 - x) / 0.5), 0)
                    middle = max(min((x - 0.0) / 0.5, (1.0 - x) / 0.5), 0)
                    high = max(min((x - 0.5) / 0.5, (2.0 - x) / 1.0), 0)
                    SOC_target_agg[i] = max(max(min(ant1, low), min(ant2, high)), min(ant3, middle))
                end
            SOC_target = ((2 * sum((mfi * xi for (mfi, xi) = zip(SOC_target_agg, LinRange{Float64}(0.0, 1.0, 101)))) - first(SOC_target_agg) * 0) - last(SOC_target_agg) * 1) / ((2 * sum(SOC_target_agg) - first(SOC_target_agg)) - last(SOC_target_agg))
            return SOC_target
        end

        transition_off = TruthTable(;  # State: Off
                                    conditions=[function (state_machine)
                                                    run_fuzzy!(params, sim_params)
                                                    return params["storage"].load_end_of_last_timestep <
                                                        params["storage"].capacity *
                                                        params["low_threshold"]
                                                end],
                                    table_data=Dict{Tuple,UInt}(
                                        (false,) => 1,
                                        (true,) => 2,
                                    ))

        transition_load = TruthTable(;  # State: Load
                                     conditions=[function (state_machine)
                                                    run_fuzzy!(params, sim_params)
                                                    return params["storage"].load_end_of_last_timestep >=
                                                        params["storage"].capacity *
                                                        params["high_threshold"]
                                                end,
                                                function (state_machine)
                                                    return state_machine.time_in_state *
                                                        sim_params["time_step_seconds"] >=
                                                        params["min_run_time"]
                                                end],
                                     table_data=Dict{Tuple,UInt}(
                                        (true, true) => 1,
                                        (false, true) => 2,
                                        (true, false) => 2,
                                        (false, false) => 2,
                                     ))
        

        state_machine = StateMachine(UInt(1),           # state
                                     Dict{UInt,String}( # state_names
                                         1 => "Off",
                                         2 => "Load",
                                     ),
                                     Dict{UInt,TruthTable}( # transitions
                                         1 => transition_off,
                                         2 => transition_load,
                                     ))

        return new("storage_driven_fuzzy", params, state_machine)
    end
end

# method for control module name on type-level
control_module_name(::Type{CM_StorageDrivenFuzzy})::String = "storage_driven_fuzzy"

function has_method_for(mod::CM_StorageDrivenFuzzy, func::ControlModuleFunction)::Bool
    return func == cmf_upper_plr_limit
end

function update(mod::CM_StorageDrivenFuzzy)
    move_state(mod.state_machine)
end

function upper_plr_limit(mod::CM_StorageDrivenFuzzy, sim_params::Dict{String,Any})::Float64
    return mod.state_machine.state == 2 ? 1.0 : 0.0
end
