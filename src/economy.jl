
export extend_profile

"""
Holds the energy profiles of each component that are relevant for economy calculation
"""
Base.@kwdef struct EconomyEmissionData{Component}
    component::Component
    energy_out::Union{Nothing,Vector{Float64}} = nothing
    energy_in::Union{Nothing,Vector{Float64}} = nothing
    energy_supply::Union{Nothing,Vector{Float64}} = nothing
    energy_demand::Union{Nothing,Vector{Float64}} = nothing
end

"""
Holds results of the economy calculation
"""
Base.@kwdef mutable struct EconomyResult
    annuity_capex::Float64 = 0.0
    annuity_opex::Float64 = 0.0
    annuity_energies::Float64 = 0.0
    breakdown::Dict{String,Any} = Dict{String,Any}()
end

function prepare_economy_emissions_data(components::Grouping,
                                        output_keys::Vector{EnergySystems.OutputKey},
                                        output_data::AbstractMatrix{<:Real})
    # Helper: Find index of required value_key in output_keys for a given component
    function find_key(component::EnergySystems.Component, output_keys::Vector{EnergySystems.OutputKey},
                      value_key::String)
        for (idx, key) in pairs(output_keys)
            if component == key.unit && value_key == key.value_key
                return idx
            end
        end
        @error "The output key `$value_key` could not be retrieved from component `$(component.uac)` but is " *
               "required for the calculation of economy and/or emissions."
    end

    data = EconomyEmissionData[]

    for (key, component) in pairs(components)
        sf = component.sys_function

        # ignore busses
        sf === EnergySystems.sf_bus && continue

        # add all other components without additional information
        entry = EconomyEmissionData(; component=component)

        # add fixed and flexible sources and sinks with actual input/output energies and for fixed ones also
        # with the demand or supply to calculate unmet demands/supplies.
        if sf === EnergySystems.sf_fixed_source
            entry = EconomyEmissionData(; component=component,
                                        energy_out=copy(@view output_data[:, find_key(component, output_keys, "OUT")]),
                                        energy_supply=copy(@view output_data[:,
                                                                             find_key(component, output_keys, "Supply")]))
        elseif sf === EnergySystems.sf_flexible_source
            entry = EconomyEmissionData(; component=component,
                                        energy_out=copy(@view output_data[:, find_key(component, output_keys, "OUT")]))

        elseif sf == EnergySystems.sf_fixed_sink
            entry = EconomyEmissionData(; component=component,
                                        energy_in=copy(@view output_data[:, find_key(component, output_keys, "IN")]),
                                        energy_demand=copy(@view output_data[:,
                                                                             find_key(component, output_keys, "Demand")]))
        elseif sf == EnergySystems.sf_flexible_sink
            entry = EconomyEmissionData(; component=component,
                                        energy_in=copy(@view output_data[:, find_key(component, output_keys, "IN")]))
        end

        push!(data, entry)
    end

    return data
end

function calculate_economy(shared_data::Vector{EconomyEmissionData},
                           sim_params::Dict{String,Any},
                           economy_parameter::AbstractDict{String,Any})
    # set initials
    result = EconomyResult()
    observation_period_in_years_economy = Year(economy_parameter["observation_period_in_years"])

    # get time stamp of extended simulation results
    economy_end_date = sim_params["start_date_output"] + Year(economy_parameter["observation_period_in_years"])
    number_of_timesteps = Int(Dates.value(Second(sub_ignoring_leap_days(economy_end_date,
                                                                        sim_params["start_date_output"]))) /
                              sim_params["time_step_seconds"])
    # simulation_result_timestamp = Vector{DateTime}(undef, number_of_timesteps)
    # for idx in 1:number_of_timesteps
    #     simulation_result_timestamp[idx] = add_ignoring_leap_days(sim_params["start_date_output"],
    #                                                               Second((idx - 1) * sim_params["time_step_seconds"]))
    #     # TODO not used currently
    # end

    # iterate over all components in the current energy system
    for item in shared_data
        component = item.component
        sf = component.sys_function

        if sf === EnergySystems.sf_bus
            continue
        end
        if sf in [EnergySystems.sf_fixed_source, EnergySystems.sf_flexible_source,
                  EnergySystems.sf_flexible_sink, EnergySystems.sf_fixed_sink]
            # get energy profiles
            if sf === EnergySystems.sf_fixed_source
                energy = extend_profile(.-item.energy_out, economy_parameter["observation_period_in_years"], sim_params) # source output
                energy_supply = extend_profile(.-item.energy_supply, economy_parameter["observation_period_in_years"], sim_params) # source supply
                energy_unmet = energy_supply .- energy # unmet supply of source
            elseif sf === EnergySystems.sf_flexible_source
                energy = extend_profile(.-item.energy_out, economy_parameter["observation_period_in_years"], sim_params) #  source output
                energy_unmet = nothing # unmet supply of source
            elseif sf === EnergySystems.sf_flexible_sink
                energy = extend_profile(item.energy_in, economy_parameter["observation_period_in_years"], sim_params) # sink input
                energy_demand = extend_profile(item.energy_demand, economy_parameter["observation_period_in_years"],
                                               sim_params) # sink demand
                energy_unmet = energy_demand .- energy # unmet demand of sink
            elseif sf === EnergySystems.sf_fixed_sink
                energy = extend_profile(item.energy_in, economy_parameter["observation_period_in_years"], sim_params) # sink input
                energy_unmet = nothing  # unmet demand of sink
            end

            # Note: (unmet) energies from sources are negative, (unmet) energies into sinks are positive at this point
            # TODO Unmet energies not considered yet!

            # calculate annuity for energies (demand-related costs)
            result.annuity_energies, result.breakdown = calculate_annuity_of_energies(energy, component,
                                                                                      sim_params, annuity_energies,
                                                                                      breakdown)

        elseif sf in [EnergySystems.sf_storage, EnergySystems.sf_transformer]
            # start and end energy of storage?! TODO
            # calculate capex, operation-related opex and optional additional component-specific opex
            result.annuity_capex, result.annuity_opex, result.breakdown = calculate_annuity_of_capex_and_opex(component,
                                                                                                              sim_params,
                                                                                                              result.annuity_capex,
                                                                                                              result.annuity_opex,
                                                                                                              result.breakdown)
        end
    end
    return result
end

function extend_profile(profile::Union{Nothing,Vector{Float64}}, observation_period_in_years_economy::Float64,
                        sim_params::Dict{String,Any})
    if profile === nothing
        return nothing
    end
    # currently only a very simple algorithm is used. The profile is taken as it is and it is repeated until the 
    # observation_period_in_years_economy is reached. If the profile does not cover a whole year, this is may not
    # the best way to extend the profile...
    # Add: Repeat only last year/month TODO
    economy_end_date = sim_params["start_date_output"] + Year(observation_period_in_years_economy)
    nr_to_repeat = Int(ceil((economy_end_date - sim_params["end_date"]) /
                            (sim_params["end_date"] - sim_params["start_date_output"])))
    profile_extended = repeat(profile, nr_to_repeat + 1)

    return profile_extended
end

function add_to_breakdown!(breakdown::Dict{String,Any}, uac::String, dict_to_add::Dict{String,Any})
    merge!(get!(breakdown, uac, Dict{String,Any}()), dict_to_add)
end

function calculate_annuity_of_energies(energy_profile::Vector{Float64}, component::EnergySystems.Component,
                                       sim_params::Dict{String,Any}, annuity_energies::Float64,
                                       breakdown::Dict{String,Any})
    # get price profile from component
    if isnothing(component.economy_parameter["constant_energy_price"])
        price_profile_energy = component.economy_parameter["energy_price_profile_scale"] .*
                               component.economy_parameter["energy_price_profile"]
    else
        price_profile_energy = fill(component.economy_parameter["constant_energy_price"],
                                    sim_params["number_of_time_steps_output"])
    end
    price_profile_energy = extend_profile(price_profile_energy, economy_parameter["observation_period_in_years"],
                                          sim_params)

    energy_price_change_rate_per_year = component.economy_parameter["energy_price_change_rate_per_year"]
    base_cost_per_year = component.economy_parameter["base_cost_per_year"]
    base_cost_change_rate_per_year = component.economy_parameter["base_cost_change_rate_per_year"]
    observation_period = sim_params["economy_parameter"]["observation_period_in_years"]

    # factors
    r_energy_costs = 1.0 + energy_price_change_rate_per_year  # price change factor of energy costs/revenues
    r_base_costs = 1.0 + base_cost_change_rate_per_year   # price change factor of base costs

    # calculate yearly costs / revenues of energies
    energy_costs_per_year = zeros(Float64, observation_period)
    energy_base_costs_per_year = zeros(Float64, observation_period)
    energy_costs_per_timestep = energy_profile .* price_profile_energy
    timesteps_per_year = Int(floor(length(energy_costs_per_timestep) / observation_period))
    for y in 1:observation_period
        energy_base_costs_per_year[year] = base_cost_per_year * r_base_costs^(y - 1)
        energy_costs_per_year[year] = sum(energy_costs_per_timestep[((y - 1) * timesteps_per_year):(y * timesteps_per_year)];
                                          init=1.0) * r_energy_costs^(y - 1)
    end

    # calculate annuity
    annuity_energies += get_annuity(energy_costs_per_year .+ energy_base_costs_per_year, sim_params)

    add_to_breakdown!(breakdown, component.uac,
                      Dict("annuity_energies" => annuity_energies,
                           "energy_costs_per_year" => energy_costs_per_year,
                           "energy_base_costs_per_year" => energy_base_costs_per_year))

    return annuity_energies, breakdown
end

function get_annuity(costs::Vector{Float64}, sim_params::Dict{String,Any})::Float64
    q = 1.0 + sim_params["economy_parameter"]["interest_rate"]
    T = length(costs)

    # Timing convention: costs[1] occurs at t=0, costs[k] at t=k-1
    # Present value of the given cost stream
    present_value = 0.0
    for k in 1:T
        present_value += costs[k] / q^(k - 1)
    end

    # Annuity factor a(q,T)
    a = if abs(sim_params["economy_parameter"]["interest_rate"]) < sim_params["epsilon"]
        1.0 / T
    else
        (q^T * (q - 1.0)) / (q^T - 1.0)
    end

    return present_value * a
end

function calculate_annuity_of_capex_and_opex(component::EnergySystems.Component, sim_params::Dict{String,Any},
                                             annuity_capex::Float64, annuity_opex::Float64, breakdown::Dict{String,Any})
    # calculate capex /(capital-related costs) including replacements, residual returns and subsidies
    investment_first_year,
    investments_per_year,
    residuals_per_year,
    subsidies_per_year = get_capex_from_component(component, sim_params)

    annuity_capex += get_annuity(investments_per_year .- residuals_per_year .- subsidies_per_year, sim_params)

    # calculate opex (operation-related costs) including inspections, servicing and repair 
    # (no energy flows here, as they are handled differently)
    maintenance_per_year,
    repair_per_year,
    labour_per_year = get_opex_from_component(component, investment_first_year, sim_params)

    annuity_opex += get_annuity(maintenance_per_year .+ repair_per_year .+ labour_per_year, sim_params)

    # get additional component-specific opex of material flows not represented by any interfaces
    # call component-specific function to get special opex, e.g. electrolysers water demand and oxygen production
    names, additional_opex_per_year = get_additional_opex_from_component(component, sim_params)
    for (idx, name) in enumerate(names)
        annuity_opex += get_annuity(additional_opex_per_year[idx], sim_params)
        add_to_breakdown!(breakdown, component.uac, Dict(name => additional_opex_per_year[idx]))
    end

    add_to_breakdown!(breakdown, component.uac,
                      Dict("annuity_capex" => annuity_capex,
                           "annuity_opex" => annuity_opex,
                           "investments_cost_per_year" => investments_per_year,
                           "residuals_cost_per_year" => residuals_per_year,
                           "subsidies_cost_per_year" => subsidies_per_year,
                           "maintenance_cost_per_year" => maintenance_per_year,
                           "repairs_cost_per_year" => repair_per_year,
                           "labour_cost_per_year" => labour_per_year
                           ))

    return annuity_capex, annuity_opex, breakdown
end

function get_capex_from_component(component::EnergySystems.Component, sim_params::Dict{String,Any})
    # capex including replacements, residuals and subsidies
    # this is probably a component-specific implementation with a default implementation? No, I guess its only general.

    # parameter component
    lifetime_years = component.economy_parameter["lifetime_years"]
    capex_specific = component.economy_parameter["capex_specific"]
    capex_price_change_rate_per_year = component.economy_parameter["capex_price_change_rate_per_year"]
    subsidy_rate_of_capex = component.economy_parameter["subsidy_rate_of_capex"]
    subsidy_max = component.economy_parameter["subsidy_max"]

    capex_reference = get_capex_reference(component)  # e.g. installed power, area, etc.

    # parameter general
    observation_period = sim_params["economy_parameter"]["observation_period_in_years"]

    # factors
    r = 1.0 + capex_price_change_rate_per_year  # price change factor of capex

    # results
    investments_per_year = zeros(Float64, observation_period)
    residuals_per_year = zeros(Float64, observation_period)
    subsidies_per_year = zeros(Float64, observation_period)

    # calculate initial investment (A0)
    investment_first_year = capex_specific * capex_reference   # € at time 0
    investments_per_year[1] = investment_first_year

    # subsidies as one-time event at the beginning of the observation period
    cap = (subsidy_max > 0.0) ? subsidy_max : Inf
    subsidies_first_year = min(investment_first_year * subsidy_rate_of_capex, cap)
    subsidies_per_year[1] += subsidies_first_year

    # calculate capex of replacements
    # number of replacements n within observation period
    n = max(0, Int(floor((observation_period - 1e-9) / lifetime_years)))
    for j in 1:n
        t = j * lifetime_years      # years since start
        idx = Int(round(t)) + 1     # index for t=0 -> 1, t=1 -> 2, ...

        if 1 <= idx <= observation_period
            investments_per_year[idx] += investment_first_year * r^t
        end
    end

    # calculate residual value
    t_last = n * lifetime_years
    remaining_fraction = max(0.0, ((n + 1) * lifetime_years - observation_period) / lifetime_years)
    residuals_per_year[observation_period] = (investment_first_year * r^t_last) * remaining_fraction

    return investment_first_year, investments_per_year, residuals_per_year, subsidies_per_year
end

function get_opex_from_component(component::EnergySystems.Component, investment_first_year::Float64,
                                 sim_params::Dict{String,Any})
    # opex, here only for component-specific opex, energy are calculated separately.
    maintenance_inspection_rate_per_year = component.economy_parameter["maintenance_inspection_rate_per_year"]
    maintenance_inspection_price_change_rate_per_year = component.economy_parameter["maintenance_inspection_price_change_rate_per_year"]
    repair_rate_per_year = component.economy_parameter["repair_rate_per_year"]
    repair_price_change_rate_per_year = component.economy_parameter["repair_price_change_rate_per_year"]
    operational_labour_hours_per_year = component.economy_parameter["operational_labour_hours_per_year"]
    labour_costs_per_hour = sim_params["economy_parameter"]["labour_costs_per_hour"]
    labour_costs_price_change_rate_per_year = sim_params["economy_parameter"]["labour_costs_price_change_rate_per_year"]
    observation_period = sim_params["economy_parameter"]["observation_period_in_years"]

    # calculate factors
    r_maintenance = 1.0 + maintenance_inspection_price_change_rate_per_year
    r_repair = 1.0 + repair_price_change_rate_per_year
    r_labour = 1.0 + labour_costs_price_change_rate_per_year

    # results (nominal yearly cashflows)
    # convention: index 1 corresponds to the first modeled year.
    maintenance_per_year = zeros(Float64, observation_period)
    repair_per_year = zeros(Float64, observation_period)
    labour_per_year = zeros(Float64, observation_period)

    # build yearly series
    # year y uses exponent (y-1): year1 -> r^0, year2 -> r^1, ...
    # Note: annual maintenance/repair rates are applied to the initial investment.
    for y in 1:observation_period
        maintenance_per_year[y] = investment_first_year * maintenance_inspection_rate_per_year * r_maintenance^(y - 1)
        repair_per_year[y] = investment_first_year * repair_rate_per_year * r_repair^(y - 1)
        labour_per_year[y] = operational_labour_hours_per_year * labour_costs_per_hour * r_labour^(y - 1)
    end

    return maintenance_per_year, repair_per_year, labour_per_year
end

# ----------------------------------------------------------------------------------------------
# Emissions

Base.@kwdef mutable struct EmissionsResult
    total::Float64 = 0.0
    breakdown::Dict{Any,Any} = Dict{Any,Any}()
end

function calculate_emissions(shared_data::Vector{EconomyEmissionData},
                             sim_params::Dict{String,Any},
                             emissions_parameter::AbstractDict{String,Any})
    result = EmissionsResult()

    for item in shared_data
        component = item.component
        sf = item.sys_function

        if sf === EnergySystems.sf_bus
            continue
        end
        if sf in [EnergySystems.sf_fixed_source, EnergySystems.sf_flexible_source]
            energy_output = extend_profile(item.energy_out, emissions_parameter["observation_period_in_years"],
                                           sim_params)
            # emission calculation here

        elseif sf in [EnergySystems.sf_flexible_sink, EnergySystems.sf_fixed_sink]
            energy_input = extend_profile(item.energy_in, emissions_parameter["observation_period_in_years"],
                                          sim_params)
            # emission calculation here

        elseif sf in [EnergySystems.sf_storage, EnergySystems.sf_transformer]
            # if no emissions, do nothing
        end
    end
end

# # ----------------------------------------------------------------------------------------------
# # Von Adrian

# ############################################################
# #  VDI 2067 MATHEMATICAL FUNCTIONS
# ############################################################

# # Annuity factor 
# # a = q^T * (q − 1) / (q^T − 1)
# function annuity_factor(i, T)
#     q = 1 + i
#     T == 0 && return 0.0
#     isapprox(q, 1.0) && return 1 / T
#     return (q^T * (q - 1)) / (q^T - 1)
# end

# # Price change factor 
# # b = (1 − (qv / q)^T) / (q − qv)
# function price_change_factor(r, i, T)
#     q = 1 + i
#     qv = 1 + r
#     isapprox(q, qv) && return T / q
#     return (1 - (qv / q)^T) / (q - qv)
# end

# # Accept scalar or vector price input
# vecize_price(p, N) = isa(p, Number) ? fill(p, N) :
#                      (length(p) == N ? p : error("Price vector length mismatch."))

# ############################################################
# #  CAPITAL COSTS — REPLACEMENTS AND RESIDUAL VALUE
# ############################################################

# function npv_replacements(A0, TN, T, i, r)      # net present value of component replacements
#     npv = 0.0
#     t = TN
#     while t < T
#         cost_t = A0 * (1 + r)^t
#         npv += cost_t / (1 + i)^t
#         t += TN
#     end
#     return npv
# end

# function residual_value(A0, TN, T, i, r)    # residual value at end of evaluation period
#     q = 1 + i
#     qv = 1 + r

#     n = div(T, TN)                  # amount of component replacements in evaluation period
#     t_last_repl = n * TN            # year of the last component replacement

#     if t_last_repl == T
#         return 0.0                       # replacement at the end of the evaluation period → no residual value
#     end

#     A_last = A0 * qv^t_last_repl    # cost of last replacement in the evaluation period

#     write_off = ((n + 1) * TN - T) / TN     # function for linear write off / depreciation

#     RW_T = A_last * write_off       # residual value of last replacement at the end of evaluation period

#     return RW_T / q^T   # Discounting residual value to t=0
# end

# function capital_annuity(comp::VDIComponent, p::VDIParams)      # capital cost-related annuity
#     return (comp.A0 +
#             npv_replacements(comp.A0, comp.TN, p.T, p.i_cap, p.r_cap) -
#             residual_value(comp.A0, comp.TN, p.T, p.i_cap, p.r_cap)) *
#            annuity_factor(p.i_cap, p.T)
# end

# function capital_annuity_subsidy(comp::VDIComponent, p::VDIParams)      # capital cost-related annuity with subsidy considered
#     cap = (comp.subsidy_max <= 0) ? Inf : comp.subsidy_max
#     S0 = min(comp.A0 * comp.subsidy_p, cap)   # subsidy amount at t=0
#     return ((comp.A0 - S0) +
#             npv_replacements(comp.A0, comp.TN, p.T, p.i_cap, p.r_cap) -
#             residual_value(comp.A0, comp.TN, p.T, p.i_cap, p.r_cap)) *
#            annuity_factor(p.i_cap, p.T)
# end

# ############################################################
# #  OPERATING COSTS: OPERATION (B) + MAINTENANCE (IN)
# ############################################################

# function op_annuity(components::Vector{VDIComponent}, p::VDIParams)     # operation cost-related annuity

#     # first year operating costs (f_bedien in h/a, hourly labor cost -> 30 EUR/h)
#     A_B1 = sum((c.f_bedien * 30.0) for c in components if c.A0 > 0; init=0.0)  # only use f_bedien if component is actually used -> p_th > 0

#     # First-year inspection, maintenance, repair and reinstatement cost
#     A_IN1 = sum(c.A0 * (c.f_instand + c.f_wartung) for c in components if c.A0 > 0; init=0.0)

#     a = annuity_factor(p.i_cap, p.T)
#     b_B = price_change_factor(p.r_op, p.i_cap, p.T)
#     b_IN = price_change_factor(p.r_inst, p.i_cap, p.T)

#     return A_B1 * a * b_B + A_IN1 * a * b_IN
# end

# ############################################################
# #  MISCELLANEOUS COSTS (INSURANCE, ADMIN,…)
# ############################################################

# function misc_annuity(components::Vector{VDIComponent}, p::VDIParams)
#     total_cap = sum(c.A0 for c in components if c.A0 > 0; init=0.0)
#     A1 = 0.02 * total_cap        # rule-of-thumb: 2% of investment

#     a = annuity_factor(p.i_cap, p.T)
#     b = price_change_factor(p.r_misc, p.i_cap, p.T)

#     return A1 * a * b
# end

# ############################################################
# #  ENERGY COSTS — TIME-RESOLVED → ANNUITY
# ############################################################

# function energy_annuity(sim::Dict, p::VDIParams)
#     IN = sim["Grid_IN"] .* 1e-6        # convert Wh time series in MWh
#     base_price = vecize_price(sim["Grid_price"], length(IN))   # €/MWh (market price)   # fix price = 214.0                          

#     # A_V1: energy costs of first year [EUR]
#     # MWh * EUR/MWh → EUR
#     A1 = sum(IN .* base_price)

#     a = annuity_factor(p.i_cap, p.T)
#     b = price_change_factor(p.r_energy, p.i_cap, p.T)

#     return A1 * a * b
# end

# ############################################################
# #  REVENUE — CONTROL ENERGY + CONTROL CAPACITY
# ############################################################

# #TODO implement positive and negatice control power
# #TODO implement costs for aggregator /virtua power plant operator (percentage of revenues)

# function revenue_control(sim::Dict, p::VDIParams)

#     # Fixed model time step
#     Δt = 0.25  # hours (15 minutes)

#     # 1) CONTROL ENERGY (WORK) REVENUE
#     A_E_energy = 0.0

#     if haskey(sim, "Control_energy") && haskey(sim, "Control_energy_price")
#         E = sim["Control_energy"] .* 1e-6     # convert Wh in MWh
#         price_E = vecize_price(sim["Control_energy_price"], length(E))

#         # MWh * EUR/MWh → EUR
#         A_E_energy = sum(E .* price_E)
#     end

#     # 2) CONTROL CAPACITY (POWER) REVENUE
#     A_E_capacity = 0.0

#     if haskey(sim, "Control_energy") && haskey(sim, "Control_power_price")
#         E = sim["Control_energy"]     # Wh
#         price_P = vecize_price(sim["Control_power_price"], length(E))

#         # Offered capacity inferred from reserve energy:
#         # P_t = E_t / Δt  → Wh / h = W → MW
#         P_MW = (E ./ Δt) .* 1e-6

#         # MW * EUR/MW → EUR
#         A_E_capacity = sum(P_MW .* price_P)
#     end

#     # 3) TOTAL FIRST-YEAR CONTROL RESERVE REVENUE
#     A_E1 = A_E_energy + A_E_capacity

#     # 4) ANNUITIZATION
#     a = annuity_factor(p.i_cap, p.T)
#     b = price_change_factor(p.r_rev, p.i_cap, p.T)

#     return A_E1 * a * b
# end

# ############################################################
# #  REVENUES — FEED-IN
# ############################################################

# function revenue_feedin(sim::Dict, p::VDIParams)

#     # 1) PHOTOVOLTAIC
#     A_E_PV = 0.0

#     if haskey(sim, "Grid_Out_PV")
#         E_PV = sim["Grid_Out_PV"] .* 1e-6      # Wh -> MWh                                  
#         market_value_PV = vecize_price(sim["Market_Value_PV"], length(E_PV))

#         # Effective remuneration per timestep: max(MW, AW)
#         # the bigger number (market value or AW) per timestep is used
#         price_eff_PV = max.(market_value_PV, p.AW_PV)

#         # MWh * €/MWh → EUR
#         A_E_PV = sum(E_PV .* price_eff_PV)
#     end

#     # 2) WIND
#     A_E_Wind = 0.0

#     if haskey(sim, "Grid_Out_Wind")
#         E_Wind = sim["Grid_Out_Wind"] .* 1e-6     # Wh -> MWh  
#         market_value_Wind = vecize_price(sim["Market_Value_Wind"], length(E_Wind))

#         # Effective remuneration per timestep: max(MW, AW)
#         # the bigger number (market value or AW) per timestep is used
#         price_eff_Wind = max.(market_value_Wind, p.AW_Wind)

#         # MWh * €/MWh → EUR
#         A_E_Wind = sum(E_Wind .* price_eff_Wind)
#     end

#     # 3) TOTAL FIRST-YEAR FEED-IN REVENUE
#     A_E2 = A_E_PV + A_E_Wind

#     # 4) VDI 2067 ANNUITIZATION
#     a = annuity_factor(p.i_cap, p.T)
#     b = price_change_factor(p.r_rev, p.i_cap, p.T)

#     return A_E2 * a * b
# end

# ############################################################
# #  MAIN FUNCTION
# ############################################################

# function vdi2067_annuity(sim::Union{Dict,OrderedDict}, components::Vector{VDIComponent}, p::VDIParams)
#     # sim_new = Dict()
#     # sim_new["Grid_IN"] = sim["m_power EnergyFlow Grid_IN->HeatPump"] .+
#     #                      sim["m_power EnergyFlow Grid_IN->ElectrodeBoiler"] .+
#     #                      sim["m_power EnergyFlow Grid_IN->Demand_Power"] .+ sim["m_power EnergyFlow Grid_IN->Battery"]
#     # sim_new["Power_Demand_P2H"] = sim["m_power EnergyFlow Grid_IN->HeatPump"] .+
#     #                               sim["m_power EnergyFlow Grid_IN->ElectrodeBoiler"] .+
#     #                               sim["m_power EnergyFlow Photovoltaic->HeatPump"] .+
#     #                               sim["m_power EnergyFlow Photovoltaic->ElectrodeBoiler"] .+
#     #                               sim["m_power EnergyFlow WindFarm->HeatPump"] .+
#     #                               sim["m_power EnergyFlow WindFarm->ElectrodeBoiler"]
#     # sim_new["Power_Demand_Demand"] = sim["m_power EnergyFlow Grid_IN->Demand_Power"] .+
#     #                                  sim["m_power EnergyFlow Photovoltaic->Demand_Power"] .+
#     #                                  sim["m_power EnergyFlow WindFarm->Demand_Power"] .+
#     #                                  sim["m_power EnergyFlow Battery->Demand_Power"]
#     # sim_new["Grid_Out_PV"] = sim["m_power EnergyFlow Photovoltaic->Grid_OUT"]
#     # sim_new["Grid_Out_Wind"] = sim["m_power EnergyFlow WindFarm->Grid_OUT"]
#     # sim_new["Control_energy"] = sim["NegControlReserve m_power OUT"]
#     # sim_new["Control_power_price"] = sim["Reserve_Power_Price_Neg"]
#     # sim_new["Control_energy_price"] = sim["Reserve_Energy_Price_Neg"]
#     # sim_new["Grid_price"] = sim["Stock_Price"]
#     # sim_new["Market_Value_PV"] = sim["Market_Price_PV"]
#     # sim_new["Market_Value_Wind"] = sim["Market_Price_Wind"]
#     # sim_new["CO2_Grid"] = sim["CO2_Grid"]

#     # sim = OrderedDict()

#     A_cap = sum(capital_annuity(c, p) for c in components)
#     A_op = op_annuity(components, p)
#     A_misc = misc_annuity(components, p)
#     A_energy = energy_annuity(sim_new, p)
#     A_cap_subsidy = sum(capital_annuity_subsidy(c, p) for c in components)
#     A_rev_control = revenue_control(sim_new, p)
#     A_rev_feed = revenue_feedin(sim_new, p)

#     A_total = A_cap + A_op + A_misc + A_energy -
#               (A_rev_control + A_rev_feed)

#     A_total_subsidy = A_cap_subsidy + A_op + A_misc + A_energy - (A_rev_control + A_rev_feed)

#     # CO2_yearly = co2_yearly(sim_new) #TODO implement

#     return OrderedDict("A_cap" => A_cap,
#                        "A_cap_subsidy" => A_cap_subsidy,
#                        "A_op" => A_op,
#                        "A_misc" => A_misc,
#                        "A_energy" => A_energy,
#                        "A_rev_control" => A_rev_control,
#                        "A_rev_feed" => A_rev_feed,
#                        "A_total" => A_total,
#                        "A_total_subsidy" => A_total_subsidy
#                        # "CO2_yearly" => CO2_yearly  #TODO implement
#                        )
# end
