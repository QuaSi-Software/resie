using Dates
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
    unit_uac::String

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
            "renewable_source_uacs" => [],
            "price_profile_path" => nothing,
            "price_trend_profile_path" => nothing,
            "price_volatility_profile_path" => nothing,
            "min_storage_load" => 0.0,
            "rolling_horizon_hours" => 8760.0

        )
        params = Base.merge(default_parameters, parameters)
        if haskey(parameters, "rolling_horizion_hours") && !haskey(parameters, "rolling_horizon_hours")
            params["rolling_horizon_hours"] = parameters["rolling_horizion_hours"]
        end

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
        params["renewable_sources"] = [components[uac] for uac in params["renewable_source_uacs"]]
        
        connectivity_by_state = Dict{String,Dict{Int,ConnectionMatrix}}()
        connectivity_by_state[unit_uac] = Dict{Int,ConnectionMatrix}()
        # record already calculated original connectivities for each bus
        connectivity_by_state[unit_uac][1] = components[unit_uac].connectivity

        new_connectivity = deepcopy(components[unit_uac].connectivity)

        storage_idx = findfirst(isequal(params["storage_uac"]), new_connectivity.output_order)
        if params["primary_source"].has_secondary_interface
            primary_idx = findfirst(isequal(params["primary_source_uac"] * "#secondary"), new_connectivity.input_order)
            new_connectivity.energy_flow[primary_idx][storage_idx] = 0.0
        end
        if params["secondary_source"].has_secondary_interface
            secondary_idx = findfirst(isequal(params["secondary_source_uac"] * "#secondary"), new_connectivity.input_order)
            new_connectivity.energy_flow[secondary_idx][storage_idx] = 0.0
        end

        connectivity_by_state[unit_uac][2] = new_connectivity
        params["connectivity_by_state"] = connectivity_by_state
        params["charging_state"] = 1
       
        return new(params["name"], params, unit_uac)
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
    # set connectivity to allow or disallow charging from secondary interfaces
    components[mod.unit_uac].connectivity = mod.parameters["connectivity_by_state"][mod.unit_uac][mod.parameters["charging_state"]]

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
    # p_volatility = value_at_time(mod_params["price_volatility_profile"], sim_params)
    # SOC_now = mod_params["storage"].load_end_of_last_timestep / mod_params["storage"].capacity

    # calculate price bounds for fuzzy control with rolling horizon.
    # Keep horizon length stable even when Feb 29 timestamps are skipped.
    horizon_hours = round(Int, mod_params["rolling_horizon_hours"])
    step = Dates.Second(round(Int, sim_params["time_step_seconds"]))
    desired_len = round(Int, horizon_hours * 3600 / sim_params["time_step_seconds"]) + 1
    rolling_range = DateTime[]
    dt = sim_params["current_date"]
    max_iterations = max(desired_len * 3, 1)
    iterations = 0

    while length(rolling_range) < desired_len && iterations < max_iterations
        if !(Dates.month(dt) == 2 && Dates.day(dt) == 29)
            push!(rolling_range, dt)
        end

        dt += step
        if dt > sim_params["end_date"]
            dt = sim_params["start_date"]
        end
        iterations += 1
    end

    if isempty(rolling_range)
        @error "Fuzzy Controller $(mod_params["name"]) rolling horizon has no valid timestamps."
        throw(InputError)
    end

    p_range = Array{Float64}(undef, length(rolling_range))
    for (idx, dt) in enumerate(rolling_range)
        p_range[idx] = mod_params["price_profile"].data[dt]
    end
    p_min_temp = minimum(p_range)
    p_max_temp = maximum(p_range)

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

    # Fuzzy logic returns absolute plr and chargemode
    plr, chargemode = fuzzy_control_chargemode(p_now, p_trend, p_min_temp, p_max_temp)
    if isnan(plr) || isnan(chargemode)
        
        @error "Fuzzy Controller $(mod_params["name"]) couldn't be calculated. " *
               "Check if the parameters are inside their bounds. " *
               "p_now=$p_now, p_trend=$p_trend"
        throw(InputError)
    end

    snk_temp = mod_params["storage"].high_temperature

    src_temp_primary = mod_params["primary_source"].input_interfaces[mod_params["primary_source"].m_heat_in].source.temperature_snk_out
    available_th_power_primary = mod_params["primary_source"].max_power_function(src_temp_primary, snk_temp) * primary_power
    cop_primary = mod_params["primary_source"].dynamic_cop(src_temp_primary, snk_temp)
    available_el_power_primary = available_th_power_primary / cop_primary

    src_temp_secondary = mod_params["secondary_source"].input_interfaces[mod_params["secondary_source"].m_heat_in].source.temperature_snk_out
    available_th_power_secondary = mod_params["secondary_source"].max_power_function(src_temp_secondary, snk_temp) * secondary_power
    cop_secondary = mod_params["secondary_source"].dynamic_cop(src_temp_secondary, snk_temp)
    available_el_power_secondary = available_th_power_secondary / cop_secondary

    renewable_el_energy = sum(source.supply for source in mod_params["renewable_sources"]; init=0.0)
    available_th_power_storage = sim_params["wh_to_watts"](min(mod_params["storage"].load - mod_params["min_storage_load"] * mod_params["storage"].capacity, mod_params["storage"].max_output_energy))
    th_power_demand = sim_params["wh_to_watts"](mod_params["demand"].scaling_factor * Profiles.work_at_time(mod_params["demand"].energy_profile, sim_params))

    total_th_power_timestep = available_th_power_primary + available_th_power_secondary
    total_el_power_timestep = available_el_power_primary + available_el_power_secondary
    target_power = clamp(plr * total_th_power_timestep, 0.0, total_th_power_timestep)


    if th_power_demand > (available_th_power_storage + target_power)
        # add security factor to account for calculation inaccuracies
        demand_cover_power = (th_power_demand * 1.1) - available_th_power_storage
    else
        min_plr = renewable_el_energy / total_el_power_timestep
        demand_cover_power = max(plr, min_plr) * total_th_power_timestep
    end

    storage_free_energy = max(mod_params["storage"].capacity - mod_params["storage"].load, 0.0)
    storage_charge_power_soc = sim_params["wh_to_watts"](storage_free_energy)
    storage_charge_power_rate = sim_params["wh_to_watts"](mod_params["storage"].max_load_rate * mod_params["storage"].capacity)
    max_storage_charge_power = min(storage_charge_power_soc, storage_charge_power_rate)
    charge_bonus_power = chargemode > 50 ? max_storage_charge_power : 0.0

    target_power = clamp(demand_cover_power + charge_bonus_power, 0.0, total_th_power_timestep)

    if chargemode < 50
        mod_params["charging_state"] = 2 # charging off
    else
        mod_params["charging_state"] = 1 # charging on
    end

    ### old Implementation with SOC_target got replaced by chargemode

    # Fuzzy logic returns change in plr and SOC target value
    # plr_diff, SOC_target = fuzzy_control_ems(p_now, p_trend, p_volatility, SOC_now)
    # if isnan(plr_diff) || isnan(SOC_target)
        
    #     @error "Fuzzy Controller $(mod_params["name"]) couldn't be calculated. " *
    #            "Check if the parameters are inside their bounds. " *
    #            "p_now=$p_now, p_trend=$p_trend, p_volatility=$p_volatility, SOC_now=$SOC_now"
    #     throw(InputError)
    # end

    ## control charging and discharging in the same way
    ## plr_target_prim = clamp(plr_primary + rel_primary * plr_secondary + plr_diff * (1+rel_primary), 0.0, (1+rel_primary))

    # if SOC_target <= SOC_now
    #     # reduce output power to empty the storage
    #     current_power = plr_primary * primary_power + plr_secondary * secondary_power
    #     target_power = clamp(current_power + plr_diff * total_power, 0.0, total_power)
    # else
    #     # increase output power to fill the storage
    #     missing_power = sim_params["wh_to_watts"]((SOC_target - SOC_now) * mod_params["storage"].capacity)
    #     target_power = clamp(missing_power, 0.0, total_power)
    # end

    # calculate the plrs of the sources with respect to their maximum power and priorities
    mod_params["plr_limit_primary"] = primary_power > 0 ? clamp(target_power / available_th_power_primary, 0.0, 1.0) : 0.0
    remaining_power = max(target_power - mod_params["plr_limit_primary"] * available_th_power_primary, 0.0)
    mod_params["plr_limit_secondary"] = secondary_power > 0 ? clamp(remaining_power / available_th_power_secondary, 0.0, 1.0) : 0.0

    if sim_params["time"] >= (11*30 + 12) * 24 * 3600 + 4 * 3600 @infiltrate end
    
    # # Check how much heat both heat sources and storage can provide
    # storage_power_by_soc = SOC_now * sim_params["wh_to_watts"](mod_params["storage"].capacity)
    # storage_power_limit = hasproperty(mod_params["storage"], :max_output_energy) ?
    #                       sim_params["wh_to_watts"](mod_params["storage"].max_output_energy) :
    #                       storage_power_by_soc
    # storage_available_power = min(storage_power_by_soc, storage_power_limit)
    # maximum_heat = mod_params["plr_limit_primary"] * mod_params["primary_source"].design_power_th +
    #                mod_params["plr_limit_secondary"] * mod_params["secondary_source"].design_power_th +
    #                storage_available_power
   
    # # Fallback: Ensure plr_limit is never reduced if doing so would leave demand uncovered
    # can_raise_primary = primary_power > 0 && mod_params["plr_limit_primary"] < 1.0
    # can_raise_secondary = secondary_power > 0 && mod_params["plr_limit_secondary"] < 1.0
    # if energy_demand > maximum_heat && (can_raise_primary || can_raise_secondary)
    #     mod_params["plr_limit_primary"] = primary_power > 0 ? 1.0 : 0.0
    #     mod_params["plr_limit_secondary"] = secondary_power > 0 ? 1.0 : 0.0
    #     println("Fuzzy Controller set plr too small")
    # end
end 


# TrapezoidalMF max(min((x - a) / (b - a), 1, (d - x) / (d - c)), 0)
# a:"left foot"
# b:"left shoulder" 
# c:"right shoulder"
# d:"right foot"
# TriangularMF max(min((x - a) / (b - a), (c - x) / (c - b)), 0)
    # a:"left foot"
    # b:"peak" 
    # c:"right foot"
function fuzzy_control_chargemode(p_now, p_trend, temp_min, temp_max)
   
    # definition of membership functions for price with rolling horizon bounds
      # helpers for p_now rolling horizon inputs
      # avoid degenerate membership functions when temp_min == temp_max
      min_temp = float(temp_min)
      max_temp = float(temp_max)
      if !isfinite(min_temp) || !isfinite(max_temp)
          return NaN, NaN
      end
      if max_temp <= min_temp
          delta = max(abs(min_temp) * 1e-6, 1e-6)
          min_temp -= delta
          max_temp += delta
      end

      a1 = min_temp - 1
      b1 = min_temp
      c1 = (max_temp + min_temp) / 2
      a2 = min_temp
      b2 = (max_temp + min_temp) / 2
      c2 = max_temp
      a3 = (max_temp + min_temp) / 2
      b3 = max_temp
      c3 = max_temp + 1
      
      cheap = max(min((p_now-a1) / (b1-a1), (c1-p_now) / (c1-b1)), 0)
      average = max(min((p_now-a2) / (b2-a2), (c2-p_now) / (c2-b2)), 0)
      expensive = max(min((p_now-a3) / (b3-a3), (c3-p_now) / (c3-b3)), 0) 
    
    # old version with fixed price bounds
        # cheap = max(min((p_now - -137) / 1, 1, (75 - p_now) / 65), 0)
        # average = max(min((p_now - 55) / 25, (105 - p_now) / 25), 0)
        # expensive = max(min((p_now - 100) / 20, 1, (1001 - p_now) / 1), 0)
    
    # definition of membership functions for price trend
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
      plr_num = (2 * sum((mfi * xi for (mfi, xi) = zip(plr_max_agg, LinRange{Float64}(-1.0, 2.0, 101)))) - first(plr_max_agg) * -1.0) - last(plr_max_agg) * 2.0
      plr_den = (2 * sum(plr_max_agg) - first(plr_max_agg)) - last(plr_max_agg)
      plr_max = plr_den == 0 ? NaN : plr_num / plr_den
      chargemode_agg = collect(LinRange{Float64}(-1.0, 101.0, 101))
      @inbounds for (i, x) = enumerate(chargemode_agg)
              off = max(min((x - -1.0) / 1.0, (51.0 - x) / 51.0), 0)
              on = max(min((x - 49.0) / 51.0, (101.0 - x) / 1.0), 0)
              chargemode_agg[i] = max(max(max(max(max(max(max(max(min(ant1, off), min(ant2, on)), min(ant3, on)), min(ant4, off)), min(ant5, on)), min(ant6, on)), min(ant7, off)), min(ant8, off)), min(ant9, on))
          end
      charge_num = (2 * sum((mfi * xi for (mfi, xi) = zip(chargemode_agg, LinRange{Float64}(-1.0, 101.0, 101)))) - first(chargemode_agg) * -1.0) - last(chargemode_agg) * 101.0
      charge_den = (2 * sum(chargemode_agg) - first(chargemode_agg)) - last(chargemode_agg)
      chargemode = charge_den == 0 ? NaN : charge_num / charge_den
      return plr_max, chargemode
  end
