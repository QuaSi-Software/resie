"""
Control module for running a component depending on the state of a linked storage component.
In particular it switches to a state of allowing operation of the component when the load
of the linked storage falls below the lower threshold. The module stays in this state until
the load has reached the upper threshold and the minimum run time has passed.
"""

mutable struct CM_StorageDrivenFuzzy <: ControlModule
    name::String
    parameters::Dict{String,Any}

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
            "price_profile_path" => nothing,
            "price_trend_profile_path" => nothing,
            "price_volatility_profile_path" => nothing,
            "demand_heat_profile_path" => nothing

        )
        params = Base.merge(default_parameters, parameters)

        if !(params["storage_uac"] !== nothing
             && params["storage_uac"] in keys(components)
             && components[params["storage_uac"]] isa StorageComponent)
            @error "Required storage component for control module storage_driven_fuzzy not given"
            throw(InputError)
        end
        params["storage"] = components[params["storage_uac"]]

        params["price_profile"] = Profile(params["price_profile_path"], sim_params)
        params["price_trend_profile"] = Profile(params["price_trend_profile_path"], sim_params)
        params["price_volatility_profile"] = Profile(params["price_volatility_profile_path"], sim_params)
        params["plr_limit"] = 1.0
        params["unit"] = components[unit_uac]
        params["demand_heat"] = Profile(params["demand_heat_profile_path"], sim_params)

        return new("storage_driven_fuzzy", params)
    end
end

# method for control module name on type-level
control_module_name(::Type{CM_StorageDrivenFuzzy})::String = "storage_driven_fuzzy"

function has_method_for(mod::CM_StorageDrivenFuzzy, func::ControlModuleFunction)::Bool
    return func == cmf_upper_plr_limit
end

function update(mod::CM_StorageDrivenFuzzy)
    # do nothing
end

function upper_plr_limit(mod::CM_StorageDrivenFuzzy, sim_params::Dict{String,Any})::Float64
    run_fuzzy!(mod.parameters, sim_params)
    return mod.parameters["plr_limit"]
end

function run_fuzzy!(mod_params::Dict{String,Any}, sim_params::Dict{String,Any})    
    # get input values for Fuzzy control
    p_now = value_at_time(mod_params["price_profile"], sim_params)
    p_trend = value_at_time(mod_params["price_trend_profile"], sim_params)
    p_volatility = value_at_time(mod_params["price_volatility_profile"], sim_params)
    SOC_now = mod_params["storage"].load_end_of_last_timestep / mod_params["storage"].capacity
    demand_now = mod_params["demand_heat"].value_at_time(sim_params)

    # Fuzzy logic returns change in plr and SOC target value
    plr_diff, SOC_target = fuzzy_control(p_now, p_trend, p_volatility, SOC_now)
    if isnan(plr_diff) || isnan(SOC_target)
        @error "Fuzzy Controller $(mod_params["name"]) couldn't be calculated. " *
               "Check if the parameters are inside their bounds."
        throw(InputError)
    end
    plr = mod_params["unit"].avg_plr

    # Check if demand is higher than what unit and storage can provide
    # TODO das nimmt an, dass WP und Boiler jeweils alleine (mit Speicher) Bedarf decken müssen?
    demand_coverable = plr * mod_params["unit"].design_heat_th + SOC_now * mod_params["storage"].capacity * sim_params["wh_to_watts"]
    
    if SOC_target <= SOC_now
        # change low_threshold of storage_driven_fuzzy controller
        mod_params["low_threshold"] = SOC_target
        mod_params["high_threshold"] = SOC_now
        mod_params["plr_limit"] = clamp(plr + plr_diff, 0.0, 1.0)
    else
        # change high_threshold of storage_driven_fuzzy controller
        mod_params["low_threshold"] = SOC_now
        mod_params["high_threshold"] = SOC_target
        missing_power = sim_params["wh_to_watts"]((SOC_target - SOC_now) * mod_params["storage"].capacity)
        mod_params["plr_limit"] = clamp(missing_power / mod_params["unit"].design_power_th, 0.0, 1.0)
    end
    
    # Ensure plr_limit is never reduced if doing so would leave demand uncovered
    if demand_now > demand_coverable
        mod_params["plr_limit"] = max(mod_params["plr_limit"], plr)
    end
end 

function fuzzy_control(p_now, p_trend, p_volatility, SOC_now)
      cheap = max(min((p_now - -137) / 1, 1, (80 - p_now) / 25), 0)
      average = max(min((p_now - 55) / 25, (100 - p_now) / 20), 0)
      expensive = max(min((p_now - 80) / 20, 1, (1001 - p_now) / 1), 0)
      falling = max(min((p_trend - -258) / 1, 1, (6 - p_trend) / 13), 0)
      rising = max(min((p_trend - -7) / 13, 1, (247 - p_trend) / 1), 0)
      stable = max(min((p_volatility - -1) / 1, 1, (15 - p_volatility) / 11), 0)  
      volatile = max(min((p_volatility - 4) / 11, 1, (266 - p_volatility) / 1), 0)
      empty = max(min((SOC_now - -0.5) / 0.5, (0.5 - SOC_now) / 0.5), 0)
      middle = max(min((SOC_now - 0.0) / 0.5, (1.0 - SOC_now) / 0.5), 0)
      full = max(min((SOC_now - 0.5) / 0.5, (1.5 - SOC_now) / 0.5), 0)
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
              lower = max(min((x - -2.0) / 1.0, (0.0 - x) / 1.0), 0)
              keep = max(min((x - -1.0) / 1.0, (1.0 - x) / 1.0), 0)
              rise = max(min((x - 0.0) / 1.0, (2.0 - x) / 1.0), 0)
              P2H_agg[i] = max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(min(ant1, rise), min(ant2, keep)), min(ant3, keep)), min(ant4, rise)), min(ant5, rise)), min(ant6, rise)), min(ant7, rise)), min(ant8, rise)), min(ant9, rise)), min(ant10, rise)), min(ant11, rise)), min(ant12, rise)), min(ant13, 
rise)), min(ant14, lower)), min(ant15, lower)), min(ant16, rise)), min(ant17, keep)), min(ant18, keep)), min(ant19, rise)), min(ant20, rise)), min(ant21, rise)), min(ant22, rise)), min(ant23, rise)), min(ant24, rise)), min(ant25, rise)), min(ant26, lower)), min(ant27, lower)), min(ant28, keep)), min(ant29, keep)), min(ant30, keep)), min(ant31, rise)), min(ant32, rise)), min(ant33, keep)), min(ant34, rise)), min(ant35, keep)), min(ant36, keep))
          end
      P2H = ((2 * sum((mfi * xi for (mfi, xi) = zip(P2H_agg, LinRange{Float64}(-2.0, 2.0, 101)))) - first(P2H_agg) * -2.0) - last(P2H_agg) * 2.0) / ((2 * sum(P2H_agg) - first(P2H_agg)) - last(P2H_agg))
      SOC_target_agg = collect(LinRange{Float64}(0.0, 1.0, 101))
      @inbounds for (i, x) = enumerate(SOC_target_agg)
              low = max(min((x - -0.5) / 0.5, (0.5 - x) / 0.5), 0)
              med = max(min((x - 0.0) / 0.5, (1.0 - x) / 0.5), 0)
              high = max(min((x - 0.5) / 0.5, (1.5 - x) / 0.5), 0)
              SOC_target_agg[i] = max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(min(ant1, med), min(ant2, med)), min(ant3, high)), min(ant4, high)), min(ant5, high)), min(ant6, high)), min(ant7, high)), min(ant8, high)), min(ant9, high)), min(ant10, high)), min(ant11, high)), min(ant12, high)), min(ant13, med)), min(ant14, low)), min(ant15, low)), min(ant16, med)), min(ant17, med)), min(ant18, high)), min(ant19, med)), min(ant20, high)), min(ant21, high)), min(ant22, med)), min(ant23, high)), min(ant24, high)), min(ant25, med)), min(ant26, low)), min(ant27, low)), min(ant28, med)), min(ant29, low)), min(ant30, low)), min(ant31, med)), min(ant32, med)), min(ant33, high)), min(ant34, high)), min(ant35, med)), min(ant36, med))
          end
      SOC_target = ((2 * sum((mfi * xi for (mfi, xi) = zip(SOC_target_agg, LinRange{Float64}(0.0, 1.0, 101)))) - first(SOC_target_agg) * 0) - last(SOC_target_agg) * 1) / ((2 * sum(SOC_target_agg) - 
first(SOC_target_agg)) - last(SOC_target_agg))
      return P2H, SOC_target
  end