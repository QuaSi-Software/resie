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
            "demand_uac" => nothing,
            "price_profile_path" => nothing,
            "price_trend_profile_path" => nothing,
            "price_volatility_profile_path" => nothing
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
        params["demand"] = components[params["demand_uac"]]

        return new(params["name"], params)
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
    demand_now = mod_params["demand"].demand

    # Fuzzy logic returns change in plr and SOC target value
    plr_diff, SOC_target = fuzzy_control(p_now, p_trend, p_volatility, SOC_now)
    if isnan(plr_diff) || isnan(SOC_target)
        @error "Fuzzy Controller $(mod_params["name"]) couldn't be calculated. " *
               "Check if the parameters are inside their bounds. " *
               "p_now=$p_now, p_trend=$p_trend, p_volatility=$p_volatility, SOC_now=$SOC_now"
        throw(InputError)
    end
    plr = mod_params["unit"].avg_plr
    
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

    # Check if demand is higher than what unit and storage can provide
    # TODO das nimmt an, dass WP und ElectrodeBoiler jeweils alleine (mit Speicher) Bedarf decken müssen?
    demand_coverable = mod_params["plr_limit"] * mod_params["unit"].design_power_th + SOC_now * sim_params["wh_to_watts"](mod_params["storage"].capacity) 
    
    # Ensure plr_limit is never reduced if doing so would leave demand uncovered
    if demand_now > demand_coverable && mod_params["plr_limit"] < 1.0
        mod_params["plr_limit"] = 1.0
        println("Fuzzy Controller set plr too small")
    end
end 

function fuzzy_control_ems(p_now, p_trend)
      cheap = max(min((p_now - -137) / 1, 1, (75 - p_now) / 65), 0)
      average = max(min((p_now - 55) / 25, (105 - p_now) / 25), 0)
      expensive = max(min((p_now - 100) / 20, 1, (1001 - p_now) / 1), 0)
      falling = max(min((p_trend - -258) / 1, 1, (1 - p_trend) / 8), 0)
      stable = max(min((p_trend - -7) / 7, (6 - p_trend) / 6), 0)
      rising = max(min((p_trend - -1) / 7, 1, (247 - p_trend) / 1), 0)
      ant1 = min(cheap, falling)
      ant2 = min(cheap, stable)
      ant3 = min(cheap, rising)
      ant4 = min(average, falling)
      ant5 = min(average, stable)
      ant6 = min(average, rising)
      ant7 = min(expensive, falling)
      ant8 = min(expensive, stable)
      ant9 = min(expensive, rising)
      plr_max_agg = collect(LinRange{Float64}(-1.0, 2.0, 101))
      @inbounds for (i, x) = enumerate(plr_max_agg)
              low = max(min((x - -1.0) / 1.0, (0.5 - x) / 0.5), 0)
              mid = max(min((x - 0.0) / 0.5, (1.0 - x) / 0.5), 0)
              high = max(min((x - 0.5) / 0.5, (2.0 - x) / 1.0), 0)
              plr_max_agg[i] = max(max(max(max(max(max(max(max(min(ant1, high), min(ant2, high)), min(ant3, high)), min(ant4, mid)), min(ant5, mid)), min(ant6, mid)), min(ant7, low)), min(ant8, low)), min(ant9, low))
          end
      plr_max = ((2 * sum((mfi * xi for (mfi, xi) = zip(plr_max_agg, LinRange{Float64}(-1.0, 2.0, 101)))) - first(plr_max_agg) * -1.0) - last(plr_max_agg) * 2.0) / ((2 * sum(plr_max_agg) - first(plr_max_agg)) - last(plr_max_agg))
      chargemode_agg = collect(LinRange{Float64}(-1.0, 101.0, 101))
      @inbounds for (i, x) = enumerate(chargemode_agg)
              off = max(min((x - -1.0) / 1.0, (51.0 - x) / 51.0), 0)
              on = max(min((x - 49.0) / 51.0, (101.0 - x) / 1.0), 0)
              chargemode_agg[i] = max(max(max(max(max(max(max(max(min(ant1, off), min(ant2, on)), min(ant3, on)), min(ant4, off)), min(ant5, on)), min(ant6, on)), min(ant7, off)), min(ant8, off)), min(ant9, on))
          end
      chargemode = ((2 * sum((mfi * xi for (mfi, xi) = zip(chargemode_agg, LinRange{Float64}(-1.0, 101.0, 101)))) - first(chargemode_agg) * -1.0) - last(chargemode_agg) * 101.0) / ((2 * 
sum(chargemode_agg) - first(chargemode_agg)) - last(chargemode_agg))
      return plr_max, chargemode
  end