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

        function fuzzy_control(p_now, p_trend, p_volatility, SOC_now)
            low = max(min((p_now - -200) / 64, 1, (80 - p_now) / 25), 0)
            med = max(min((p_now - 55) / 25, (100 - p_now) / 20), 0)
            high = max(min((p_now - 80) / 20, 1, (1500 - p_now) / 500), 0)
            falling = max(min((p_trend - -270) / 10, 1, (-1 - p_trend) / 6), 0)
            stable = max(min((p_trend - -7) / 8, (6 - p_trend) / 5), 0)
            rising = max(min((p_trend - 1) / 5, 1, (270 - p_trend) / 20), 0)
            low = max(min((p_volatility - -10) / 10, 1, (8 - p_volatility) / 4), 0)
            med = max(min((p_volatility - 4) / 4, (15 - p_volatility) / 7), 0)
            high = max(min((p_volatility - 8) / 7, 1, (275 - p_volatility) / 10), 0)
            low = max(min((SOC_now - -0.5) / 0.5, (0.5 - SOC_now) / 0.5), 0)
            med = max(min((SOC_now - 0.0) / 0.5, (1.0 - SOC_now) / 0.5), 0)
            high = max(min((SOC_now - 0.5) / 0.5, (1.5 - SOC_now) / 0.5), 0)
            ant1 = min(low, min(falling, min(low, low)))
            ant2 = min(low, min(falling, min(low, med)))
            ant3 = min(low, min(falling, min(low, high)))
            ant4 = min(low, min(falling, min(med, low)))
            ant5 = min(low, min(falling, min(med, med)))
            ant6 = min(low, min(falling, min(med, high)))
            ant7 = min(low, min(falling, min(high, low)))
            ant8 = min(low, min(falling, min(high, med)))
            ant9 = min(low, min(falling, min(high, high)))
            ant10 = min(low, min(stable, min(low, low)))
            ant11 = min(low, min(stable, min(low, med)))
            ant12 = min(low, min(stable, min(low, high)))
            ant13 = min(low, min(stable, min(med, low)))
            ant14 = min(low, min(stable, min(med, med)))
            ant15 = min(low, min(stable, min(med, high)))
            ant16 = min(low, min(stable, min(high, low)))
            ant17 = min(low, min(stable, min(high, med)))
            ant18 = min(low, min(stable, min(high, high)))
            ant19 = min(low, min(rising, min(low, low)))
            ant20 = min(low, min(rising, min(low, med)))
            ant21 = min(low, min(rising, min(low, high)))
            ant22 = min(low, min(rising, min(med, low)))
            ant23 = min(low, min(rising, min(med, med)))
            ant24 = min(low, min(rising, min(med, high)))
            ant25 = min(low, min(rising, min(high, low)))
            ant26 = min(low, min(rising, min(high, med)))
            ant27 = min(low, min(rising, min(high, high)))
            ant28 = min(med, min(falling, min(low, low)))
            ant29 = min(med, min(falling, min(low, med)))
            ant30 = min(med, min(falling, min(low, high)))
            ant31 = min(med, min(falling, min(med, low)))
            ant32 = min(med, min(falling, min(med, med)))
            ant33 = min(med, min(falling, min(med, high)))
            ant34 = min(med, min(falling, min(high, low)))
            ant35 = min(med, min(falling, min(high, med)))
            ant36 = min(med, min(falling, min(high, high)))
            ant37 = min(med, min(stable, min(low, low)))
            ant38 = min(med, min(stable, min(low, med)))
            ant39 = min(med, min(stable, min(low, high)))
            ant40 = min(med, min(stable, min(med, low)))
            ant41 = min(med, min(stable, min(med, med)))
            ant42 = min(med, min(stable, min(med, high)))
            ant43 = min(med, min(stable, min(high, low)))
            ant44 = min(med, min(stable, min(high, med)))
            ant45 = min(med, min(stable, min(high, high)))
            ant46 = min(med, min(rising, min(low, low)))
            ant47 = min(med, min(rising, min(low, med)))
            ant48 = min(med, min(rising, min(low, high)))
            ant49 = min(med, min(rising, min(med, low)))
            ant50 = min(med, min(rising, min(med, med)))
            ant51 = min(med, min(rising, min(med, high)))
            ant52 = min(med, min(rising, min(high, low)))
            ant53 = min(med, min(rising, min(high, med)))
            ant54 = min(med, min(rising, min(high, high)))
            ant55 = min(high, min(falling, min(low, low)))
            ant56 = min(high, min(falling, min(low, med)))
            ant57 = min(high, min(falling, min(low, high)))
            ant58 = min(high, min(falling, min(med, low)))
            ant59 = min(high, min(falling, min(med, med)))
            ant60 = min(high, min(falling, min(med, high)))
            ant61 = min(high, min(falling, min(high, low)))
            ant62 = min(high, min(falling, min(high, med)))
            ant63 = min(high, min(falling, min(high, high)))
            ant64 = min(high, min(stable, min(low, low)))
            ant65 = min(high, min(stable, min(low, med)))
            ant66 = min(high, min(stable, min(low, high)))
            ant67 = min(high, min(stable, min(med, low)))
            ant68 = min(high, min(stable, min(med, med)))
            ant69 = min(high, min(stable, min(med, high)))
            ant70 = min(high, min(stable, min(high, low)))
            ant71 = min(high, min(stable, min(high, med)))
            ant72 = min(high, min(stable, min(high, high)))
            ant73 = min(high, min(rising, min(low, low)))
            ant74 = min(high, min(rising, min(low, med)))
            ant75 = min(high, min(rising, min(low, high)))
            ant76 = min(high, min(rising, min(med, low)))
            ant77 = min(high, min(rising, min(med, med)))
            ant78 = min(high, min(rising, min(med, high)))
            ant79 = min(high, min(rising, min(high, low)))
            ant80 = min(high, min(rising, min(high, med)))
            ant81 = min(high, min(rising, min(high, high)))
            P2H_agg = collect(LinRange{Float64}(-1.0, 1.0, 101))
            @inbounds for (i, x) = enumerate(P2H_agg)
                    lower = max(min((x - -1.5) / 0.5, (0.0 - x) / 1.0), 0)
                    keep = max(min((x - -1.0) / 1.0, (1.0 - x) / 1.0), 0)
                    rise = max(min((x - 0.0) / 1.0, (1.5 - x) / 0.5), 0)
                    P2H_agg[i] = max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(min(ant1, keep), min(ant2, lower)), min(ant3, lower)), min(ant4, rise)), min(ant5, rise)), min(ant6, rise)), min(ant7, rise)), min(ant8, rise)), min(ant9, rise)), min(ant10, rise)), min(ant11, rise)), min(ant12, rise)), min(ant13, rise)), min(ant14, rise)), min(ant15, rise)), min(ant16, rise)), min(ant17, 
            rise)), min(ant18, rise)), min(ant19, rise)), min(ant20, rise)), min(ant21, rise)), min(ant22, rise)), min(ant23, rise)), min(ant24, rise)), min(ant25, rise)), min(ant26, rise)), min(ant27, rise)), min(ant28, keep)), min(ant29, lower)), min(ant30, lower)), min(ant31, keep)), min(ant32, lower)), min(ant33, lower)), min(ant34, keep)), min(ant35, keep)), min(ant36, keep)), min(ant37, keep)), min(ant38, keep)), min(ant39, keep)), min(ant40, keep)), min(ant41, keep)), min(ant42, keep)), min(ant43, keep)), min(ant44, keep)), min(ant45, keep)), min(ant46, rise)), min(ant47, rise)), min(ant48, rise)), min(ant49, rise)), min(ant50, rise)), min(ant51, rise)), min(ant52, keep)), min(ant53, keep)), min(ant54, keep)), min(ant55, keep)), min(ant56, lower)), min(ant57, lower)), min(ant58, keep)), min(ant59, lower)), min(ant60, lower)), min(ant61, keep)), min(ant62, lower)), 
            min(ant63, lower)), min(ant64, keep)), min(ant65, keep)), min(ant66, keep)), min(ant67, keep)), min(ant68, keep)), min(ant69, keep)), min(ant70, keep)), min(ant71, keep)), min(ant72, keep)), min(ant73, rise)), min(ant74, rise)), min(ant75, rise)), min(ant76, rise)), min(ant77, rise)), min(ant78, keep)), min(ant79, keep)), min(ant80, keep)), min(ant81, keep))
            end
            P2H = ((2 * sum((mfi * xi for (mfi, xi) = zip(P2H_agg, LinRange{Float64}(-1.0, 1.0, 101)))) - first(P2H_agg) * -1) - last(P2H_agg) * 1) / ((2 * sum(P2H_agg) - first(P2H_agg)) - last(P2H_agg))
            SOC_target_agg = collect(LinRange{Float64}(0.0, 1.0, 101))
            @inbounds for (i, x) = enumerate(SOC_target_agg)
                low = max(min((x - -0.5) / 0.5, (0.5 - x) / 0.5), 0)
                med = max(min((x - 0.0) / 0.5, (1.0 - x) / 0.5), 0)
                high = max(min((x - 0.5) / 0.5, (1.5 - x) / 0.5), 0)
                SOC_target_agg[i] = max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(min(ant1, low), min(ant2, med)), min(ant3, high)), min(ant4, high)), min(ant5, high)), min(ant6, high)), min(ant7, high)), min(ant8, high)), min(ant9, high)), min(ant10, high)), min(ant11, high)), min(ant12, high)), min(ant13, high)), min(ant14, high)), min(ant15, high)), min(ant16, high)), min(ant17, high)), min(ant18, high)), min(ant19, high)), min(ant20, high)), min(ant21, high)), min(ant22, high)), min(ant23, high)), min(ant24, high)), min(ant25, high)), min(ant26, high)), min(ant27, high)), min(ant28, low)), min(ant29, med)), min(ant30, high)), min(ant31, low)), min(ant32, med)), min(ant33, high)), min(ant34, low)), min(ant35, med)), min(ant36, high)), min(ant37, low)), min(ant38, med)), min(ant39, high)), min(ant40, low)), min(ant41, med)), min(ant42, high)), min(ant43, low)), min(ant44, med)), min(ant45, high)), min(ant46, high)), min(ant47, high)), min(ant48, high)), min(ant49, high)), min(ant50, high)), min(ant51, high)), min(ant52, high)), min(ant53, high)), min(ant54, high)), min(ant55, low)), min(ant56, low)), min(ant57, low)), min(ant58, low)), min(ant59, low)), min(ant60, low)), min(ant61, low)), min(ant62, low)), min(ant63, low)), min(ant64, low)), min(ant65, med)), min(ant66, high)), min(ant67, low)), min(ant68, med)), min(ant69, high)), min(ant70, low)), min(ant71, med)), min(ant72, high)), min(ant73, med)), min(ant74, high)), min(ant75, high)), min(ant76, med)), min(ant77, high)), min(ant78, high)), min(ant79, low)), min(ant80, med)), min(ant81, high))
            end
            SOC_target = ((2 * sum((mfi * xi for (mfi, xi) = zip(SOC_target_agg, LinRange{Float64}(0.0, 1.0, 101)))) - first(SOC_target_agg) * 0) - last(SOC_target_agg) * 1) / ((2 * sum(SOC_target_agg) - first(SOC_target_agg)) - last(SOC_target_agg))
            $(Expr(:return, :P2H, :SOC_target))
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
