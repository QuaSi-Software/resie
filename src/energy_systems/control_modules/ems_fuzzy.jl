using Infiltrator
"""
Control module for running a component depending on the state of a linked storage component.
In particular it switches to a state of allowing operation of the component when the load
of the linked storage falls below the lower threshold. The module stays in this state until
the load has reached the upper threshold and the minimum run time has passed.
"""
mutable struct CM_EMSFuzzy <: ControlModule
    name::String
    parameters::Dict{String,Any}

    function CM_EMSFuzzy(parameters::Dict{String,Any},
                                   components::Grouping,
                                   sim_params::Dict{String,Any},
                                   unit_uac::String)
        default_parameters = Dict{String,Any}(
            "name" => "ems_fuzzy",
            "primary_source_uac" => nothing,
            "secondary_source_uac" => nothing,
            "storage_uac" => nothing,
            "demand_uac" => nothing,
            "price_profile_path" => nothing,
            "price_trend_profile_path" => nothing,
            "price_volatility_profile_path" => nothing,

        )
        params = Base.merge(default_parameters, parameters)

        if !(params["storage_uac"] !== nothing
             && params["storage_uac"] in keys(components)
             && components[params["storage_uac"]] isa StorageComponent)
            @error "Required storage component for control module ems_fuzzy not given"
            throw(InputError)
        end
        params["storage"] = components[params["storage_uac"]]
        params["price_profile"] = Profile(params["price_profile_path"], sim_params)
        params["price_trend_profile"] = Profile(params["price_trend_profile_path"], sim_params)
        params["price_volatility_profile"] = Profile(params["price_volatility_profile_path"], sim_params)
        params["plr_limit"] = 1.0
        params["primary_source"] = components[params["primary_source_uac"]]
        params["secondary_source"] = components[params["secondary_source_uac"]]
        params["demand"] = components[params["demand_uac"]]

        return new(params["name"], params)
    end
end

# method for control module name on type-level
control_module_name(::Type{CM_EMSFuzzy})::String = "ems_fuzzy"

function has_method_for(mod::CM_EMSFuzzy, func::ControlModuleFunction)::Bool
    return func == cmf_change_bus_priorities
end

function update(mod::CM_EMSFuzzy)
    # ordinarily we would update the state machine in the control modules's update step,
    # however this control module is special in that it's callbacks are used outside the
    # normal order of operations, so the update is performed by the callbacks instead
end

function change_bus_priorities!(mod::CM_EMSFuzzy,
                                components::Grouping,
                                sim_params::Dict{String,Any})
    # perform update outside of normal order of operations
    run_fuzzy_ems!(mod.parameters, sim_params)
    # set limits by using the profile limited control module thats attached to each source
    mod.parameters["primary_source"].controller.modules[1].profile.data[sim_params["current_date"]] = mod.parameters["plr_limit_primary"]
    mod.parameters["secondary_source"].controller.modules[1].profile.data[sim_params["current_date"]] = mod.parameters["plr_limit_secondary"]
end

function reorder_operations(mod::CM_EMSFuzzy,
                            order_of_operations::OrderOfOperations,
                            sim_params::Dict{String,Any})::OrderOfOperations
    # ooo doesn't get changed in this control module, so just return the original ooo
    return order_of_operations
end

function run_fuzzy_ems!(mod_params::Dict{String,Any}, sim_params::Dict{String,Any})    
    # get input values for Fuzzy control
    p_now = value_at_time(mod_params["price_profile"], sim_params)
    p_trend = value_at_time(mod_params["price_trend_profile"], sim_params)
    p_volatility = value_at_time(mod_params["price_volatility_profile"], sim_params)
    SOC_now = mod_params["storage"].load_end_of_last_timestep / mod_params["storage"].capacity
    demand_now = mod_params["demand"].demand

    # Fuzzy logic returns change in plr and SOC target value
    plr_diff, SOC_target = fuzzy_control_ems(p_now, p_trend, p_volatility, SOC_now)
    if isnan(plr_diff) || isnan(SOC_target)
        
        @error "Fuzzy Controller $(mod_params["name"]) couldn't be calculated. " *
               "Check if the parameters are inside their bounds. " *
               "p_now=$p_now, p_trend=$p_trend, p_volatility=$p_volatility, SOC_now=$SOC_now"
        throw(InputError)
    end

    plr_primary = mod_params["primary_source"].avg_plr
    plr_secondary = mod_params["secondary_source"].avg_plr
    primary_power = mod_params["primary_source"].design_power_th
    secondary_power = mod_params["secondary_source"].design_power_th
    total_power = primary_power + secondary_power

    if total_power <= 0
        mod_params["plr_limit_primary"] = 0.0
        mod_params["plr_limit_secondary"] = 0.0
        return
    end
    
    # control charging and discharging in the same way
    # plr_target_prim = clamp(plr_primary + rel_primary * plr_secondary + plr_diff * (1+rel_primary), 0.0, (1+rel_primary))

    if SOC_target <= SOC_now
        # reduce output power to empty the storage
        current_power = plr_primary * primary_power + plr_secondary * secondary_power
        target_power = clamp(current_power + plr_diff * total_power, 0.0, total_power)
    else
        # increase output power to fill the storage
        missing_power = sim_params["wh_to_watts"]((SOC_target - SOC_now) * mod_params["storage"].capacity)
        target_power = clamp(missing_power, 0.0, total_power)
    end

    # calculate the plrs of the sources with respect to their maximum power and priorities
    mod_params["plr_limit_primary"] = primary_power > 0 ? clamp(target_power / primary_power, 0.0, 1.0) : 0.0
    remaining_power = max(target_power - mod_params["plr_limit_primary"] * primary_power, 0.0)
    mod_params["plr_limit_secondary"] = secondary_power > 0 ? clamp(remaining_power / secondary_power, 0.0, 1.0) : 0.0

    # Check how much heat both heat sources and storage can provide
    storage_power_by_soc = SOC_now * sim_params["wh_to_watts"](mod_params["storage"].capacity)
    storage_power_limit = hasproperty(mod_params["storage"], :max_output_energy) ?
                          sim_params["wh_to_watts"](mod_params["storage"].max_output_energy) :
                          storage_power_by_soc
    storage_available_power = min(storage_power_by_soc, storage_power_limit)
    maximum_heat = mod_params["plr_limit_primary"] * mod_params["primary_source"].design_power_th +
                   mod_params["plr_limit_secondary"] * mod_params["secondary_source"].design_power_th +
                   storage_available_power
   
    # Fallback: Ensure plr_limit is never reduced if doing so would leave demand uncovered
    can_raise_primary = primary_power > 0 && mod_params["plr_limit_primary"] < 1.0
    can_raise_secondary = secondary_power > 0 && mod_params["plr_limit_secondary"] < 1.0
    if demand_now > maximum_heat && (can_raise_primary || can_raise_secondary)
        mod_params["plr_limit_primary"] = primary_power > 0 ? 1.0 : 0.0
        mod_params["plr_limit_secondary"] = secondary_power > 0 ? 1.0 : 0.0
        println("Fuzzy Controller set plr too small")
    end
end 

function fuzzy_control_ems(p_now, p_trend, p_volatility, SOC_now)
      cheap = max(min((p_now - -137) / 1, 1, (75 - p_now) / 65), 0)
      average = max(min((p_now - 55) / 25, (105 - p_now) / 25), 0)
      expensive = max(min((p_now - 100) / 20, 1, (1001 - p_now) / 1), 0)
      falling = max(min((p_trend - -258) / 1, 1, (1 - p_trend) / 8), 0)
      rising = max(min((p_trend - -1) / 7, 1, (247 - p_trend) / 1), 0)
      stable = max(min((p_volatility - -1) / 1, (10 - p_volatility) / 10), 0)
      volatile = max(min((p_volatility - 8) / 7, 1, (266 - p_volatility) / 1), 0)        
      empty = max(min((SOC_now - -0.5) / 0.5, 1, (0.3 - SOC_now) / 0.15), 0)
      middle = max(min((SOC_now - 0.2) / 0.3, (0.8 - SOC_now) / 0.30000000000000004), 0) 
      full = max(min((SOC_now - 0.7) / 0.15000000000000002, 1, (1.5 - SOC_now) / 0.5), 0)
      ant1 = min(cheap, min(falling, min(stable, empty)))
      ant2 = min(cheap, min(falling, min(stable, middle)))
      ant3 = min(cheap, min(falling, min(stable, full)))
      ant4 = min(cheap, min(falling, min(volatile, empty)))
      ant5 = min(cheap, min(falling, min(volatile, middle)))
      ant6 = min(cheap, min(falling, min(volatile, full)))
      ant7 = min(cheap, min(rising, min(stable, empty)))
      ant8 = min(cheap, min(rising, min(stable, middle)))
      ant9 = min(cheap, min(rising, min(stable, full)))
      ant10 = min(cheap, min(rising, min(volatile, empty)))
      ant11 = min(cheap, min(rising, min(volatile, middle)))
      ant12 = min(cheap, min(rising, min(volatile, full)))
      ant13 = min(average, min(falling, min(stable, empty)))
      ant14 = min(average, min(falling, min(stable, middle)))
      ant15 = min(average, min(falling, min(stable, full)))
      ant16 = min(average, min(falling, min(volatile, empty)))
      ant17 = min(average, min(falling, min(volatile, middle)))
      ant18 = min(average, min(falling, min(volatile, full)))
      ant19 = min(average, min(rising, min(stable, empty)))
      ant20 = min(average, min(rising, min(stable, middle)))
      ant21 = min(average, min(rising, min(stable, full)))
      ant22 = min(average, min(rising, min(volatile, empty)))
      ant23 = min(average, min(rising, min(volatile, middle)))
      ant24 = min(average, min(rising, min(volatile, full)))
      ant25 = min(expensive, min(falling, min(stable, empty)))
      ant26 = min(expensive, min(falling, min(stable, middle)))
      ant27 = min(expensive, min(falling, min(stable, full)))
      ant28 = min(expensive, min(falling, min(volatile, empty)))
      ant29 = min(expensive, min(falling, min(volatile, middle)))
      ant30 = min(expensive, min(falling, min(volatile, full)))
      ant31 = min(expensive, min(rising, min(stable, empty)))
      ant32 = min(expensive, min(rising, min(stable, middle)))
      ant33 = min(expensive, min(rising, min(stable, full)))
      ant34 = min(expensive, min(rising, min(volatile, empty)))
      ant35 = min(expensive, min(rising, min(volatile, middle)))
      ant36 = min(expensive, min(rising, min(volatile, full)))
      P2H_agg = collect(LinRange{Float64}(-2.0, 2.0, 101))
      @inbounds for (i, x) = enumerate(P2H_agg)
              lower = max(min((x - -2.0) / 1.0, (-0.2 - x) / 0.8), 0)
              keep = max(min((x - -0.3) / 0.3, (0.3 - x) / 0.3), 0)
              rise = max(min((x - 0.2) / 0.8, (2.0 - x) / 1.0), 0)
              P2H_agg[i] = max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(min(ant1, rise), min(ant2, keep)), min(ant3, keep)), min(ant4, rise)), min(ant5, rise)), min(ant6, rise)), min(ant7, rise)), min(ant8, rise)), min(ant9, rise)), min(ant10, rise)), min(ant11, rise)), min(ant12, rise)), min(ant13, rise)), min(ant14, lower)), min(ant15, lower)), min(ant16, rise)), min(ant17, keep)), min(ant18, keep)), min(ant19, rise)), min(ant20, rise)), min(ant21, rise)), min(ant22, rise)), min(ant23, rise)), min(ant24, rise)), min(ant25, rise)), min(ant26, lower)), min(ant27, lower)), min(ant28, keep)), min(ant29, keep)), min(ant30, keep)), min(ant31, rise)), min(ant32, rise)), min(ant33, keep)), min(ant34, rise)), min(ant35, keep)), min(ant36, keep))
          end
      P2H = ((2 * sum((mfi * xi for (mfi, xi) = zip(P2H_agg, LinRange{Float64}(-2.0, 2.0, 101)))) - first(P2H_agg) * -2.0) - last(P2H_agg) * 2.0) / ((2 * sum(P2H_agg) - first(P2H_agg)) - last(P2H_agg))
      SOC_target_agg = collect(LinRange{Float64}(-2.0, 2.0, 101))
      @inbounds for (i, x) = enumerate(SOC_target_agg)
              low = max(min((x - -1.0) / 0.5, (0.5 - x) / 1.0), 0)
              med = max(min((x - 0.0) / 0.5, (1.0 - x) / 0.5), 0)
              high = max(min((x - 0.5) / 0.5, (1.5 - x) / 0.5), 0)
              SOC_target_agg[i] = max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(min(ant1, med), min(ant2, med)), min(ant3, high)), min(ant4, high)), min(ant5, high)), min(ant6, high)), min(ant7, high)), min(ant8, high)), min(ant9, high)), min(ant10, high)), min(ant11, high)), min(ant12, high)), min(ant13, med)), min(ant14, low)), min(ant15, low)), min(ant16, med)), min(ant17, med)), min(ant18, high)), min(ant19, med)), min(ant20, high)), min(ant21, high)), min(ant22, med)), min(ant23, high)), min(ant24, high)), min(ant25, med)), min(ant26, low)), min(ant27, low)), min(ant28, med)), min(ant29, low)), min(ant30, low)), min(ant31, med)), min(ant32, med)), min(ant33, high)), min(ant34, high)), min(ant35, med)), min(ant36, med))
          end
      SOC_target = ((2 * sum((mfi * xi for (mfi, xi) = zip(SOC_target_agg, LinRange{Float64}(-2.0, 2.0, 101)))) - first(SOC_target_agg) * -2.0) - last(SOC_target_agg) * 2.0) / ((2 * sum(SOC_target_agg) - first(SOC_target_agg)) - last(SOC_target_agg))
      return P2H, SOC_target
  end
