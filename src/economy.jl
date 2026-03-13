
Base.@kwdef struct EconomyEmissionData{Component}
    component::Component
    energy_out::Union{Nothing,Vector{Float64}} = nothing
    energy_in::Union{Nothing,Vector{Float64}} = nothing
    energy_supply::Union{Nothing,Vector{Float64}} = nothing
    energy_demand::Union{Nothing,Vector{Float64}} = nothing
end

Base.@kwdef mutable struct EconomyResult
    annuity_capex::Float64 = 0.0
    annuity_opex::Float64 = 0.0
    annuity_energies::Float64 = 0.0
    breakdown::Dict{Any,Any} = Dict{Any,Any}()
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
    economy_end_date = sim_params["start_date_output"] + observation_period_in_years_economy
    number_of_timesteps = Int(Dates.value(Second(sub_ignoring_leap_days(economy_end_date,
                                                                        sim_params["start_date_output"]))) /
                              sim_params["time_step_seconds"])
    simulation_result_timestamp = Vector{DateTime}(undef, number_of_timesteps)
    for idx in 1:number_of_timesteps
        simulation_result_timestamp[idx] = add_ignoring_leap_days(sim_params["start_date_output"],
                                                                  Second((idx - 1) * sim_params["time_step_seconds"]))
    end

    # iterate over all components in the current energy system
    for item in shared_data
        component = item.component
        sf = component.sys_function

        if sf === EnergySystems.sf_bus
            continue
        end
        if sf in [EnergySystems.sf_fixed_source, EnergySystems.sf_flexible_source,
                  EnergySystems.sf_flexible_sink, EnergySystems.sf_fixed_sink]
            if sf === EnergySystems.sf_fixed_source
                energy = extend_profile(.-item.energy_out, observation_period_in_years_economy, sim_params) # source output
                energy_supply = extend_profile(.-item.energy_supply, observation_period_in_years_economy, sim_params) # source supply
                energy_unmet = energy_supply .- energy # unmet supply of source
            elseif sf === EnergySystems.sf_flexible_source
                energy = extend_profile(.-item.energy_out, observation_period_in_years_economy, sim_params) #  source output
                energy_unmet = nothing # unmet supply of source
            elseif sf === EnergySystems.sf_flexible_sink
                energy = extend_profile(item.energy_in, observation_period_in_years_economy, sim_params) # sink input
                energy_demand = extend_profile(item.energy_demand, observation_period_in_years_economy, sim_params) # sink demand
                energy_unmet = energy_demand .- energy # unmet demand of sink
            elseif sf === EnergySystems.sf_fixed_sink
                energy = extend_profile(item.energy_in, observation_period_in_years_economy, sim_params) # sink input
                energy_unmet = nothing  # unmet demand of sink
            end

            # Note: (unmet) energies from sources are negative, (unmet) energies into sinks are positive at this point

            # calculate economy using correct price(profile)
            price_profile_energy = component.price_profile_energy  # TODO 
            price_profile_energy = extend_profile(price_profile_energy, observation_period_in_years_economy, sim_params)

            # price for unmet energy
            price_unmet = 0 # TODO

            # calculate annuity for opex
            result.annuity_energies, result.breakdown = calculate_annuity_of_energies(energy, price_profile_energy,
                                                                                      sim_params, annuity_energies,
                                                                                      breakdown)

        elseif sf in [EnergySystems.sf_storage, EnergySystems.sf_transformer]
            # start and end energy of storage?! TODO
            # calculate capex and optional additional component-specific opex
            result.annuity_capex, result.annuity_opex, result.breakdown = calculate_annuity_of_capex_and_opex(component::EnergySystems.Component,
                                                                                                              sim_params::Dict{String,
                                                                                                                               Any},
                                                                                                              annuity_capex::Float64,
                                                                                                              annuity_opex::Float64,
                                                                                                              breakdown_economy::Dict{String,
                                                                                                                                      Any})
        end
    end
    return result
end

function extend_profile(profile::Union{Nothing,Vector{Float64}}, observation_period_in_years_economy::DateTime,
                        sim_params::Dict{String,Any})
    if profile === nothing
        return nothing
    end
    # currently only a very simple algorithm is used. The profile is taken as it is and it is repeated until the 
    # observation_period_in_years_economy is reached. If the profile does not cover a whole year, this is may not
    # the best way to extend the profile...
    # Add: Repeat only last year/month TODO
    economy_end_date = sim_params["start_date_output"] + observation_period_in_years_economy
    nr_to_repeat = Int(ceil((economy_end_date - sim_params["end_date"]) /
                            (sim_params["end_date"] - sim_params["start_date_output"])))
    profile_extended = repeat(profile, nr_to_repeat + 1)

    return profile_extended
end

function add_to_breakdown!(breakdown::Dict{String,Any}, uac::String, dict_to_add::Dict{String,Any})
    merge!(get!(breakdown, uac, Dict{String,Any}()), dict_to_add)
end

function calculate_annuity_of_energies(energy::Vector{Float64}, price_profile::Vector{Float64},
                                       sim_params::Dict{String,Any}, annuity_energies::Float64,
                                       breakdown::Dict{String,Any})
    # calculate opex of energies
    annuity_energies_per_year = calculate_annuity(energy, price_profile, sim_params) # TODO
    annuity_energies += get_annuity(annuity_energies_per_year, sim_params)

    add_to_breakdown!(breakdown, uac, Dict("annuity_energies" => annuity_energies))

    return annuity_energies, breakdown
end

function calculate_annuity_of_capex_and_opex(component::EnergySystems.Component, sim_params::Dict{String,Any},
                                             annuity_capex::Float64,
                                             annuity_opex::Float64, breakdown::Dict{String,Any})
    # calculate capex including replacements and residual returns and subsidies
    capex_per_year, capex_first_year, residuals_per_year, subsidies_per_year = get_capex_from_component(component,
                                                                                                        sim_params)
    annuity_capex += get_annuity(capex_per_year, sim_params)
    annuity_capex += get_annuity(residuals_per_year, sim_params)
    annuity_capex += get_annuity(subsidies_per_year, sim_params)

    # calculate opex (no energy flows here, as they are handled differently)
    opex_per_year = get_opex_from_component(component, capex_first_year, sim_params)
    annuity_opex += get_annuity(opex_per_year, sim_params)

    add_to_breakdown!(breakdown, uac,
                      Dict("annuity_capex" => annuity_capex,
                           "annuity_opex" => annuity_opex,
                           "capex" => capex_per_year,
                           "residuals" => residuals_per_year,
                           "subsidies" => subsidies_per_year,
                           "opex" => opex_per_year
                           ))

    return annuity_capex, annuity_opex, breakdown
end

function get_capex_from_component(component::EnergySystems.Component, sim_params::Dict{String,Any})
    # capex including replacements, residuals and subsidies

    # this is probably a component-specific implementation with a default implementation
    # check if parameters are provided by the input:
    component.capex_specific
    component.lifetime_year
    component.capex_price_change_rate_per_year
    component.subsidy_rate_of_capex
    component.subsidy_max

    # get installed power and calculate A0

    # calculate capex of replacements and convert

    # calculate residual value

    # calculate total cash value of all capex under consideration of subsidies

end

function get_opex_from_component(component::EnergySystems.Component, capex_first_year::Float64,
                                 sim_params::Dict{String,Any})
    # opex, here only for component-specific opex, energy are calculated separately.
    # So here only for e.g. electrolysers water demand, labour, maintenance?
    component.maintenance_repair_rate_per_year
    component.operational_labour_hour_per_year
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

    # get period sto cover the whole economy observation_period_in_years
    observation_period_in_years_emissions = Years(emissions_parameter["observation_period_in_years"])

    for item in shared_data
        component = item.component
        sf = item.sys_function

        if sf === EnergySystems.sf_bus
            continue
        end
        if sf in [EnergySystems.sf_fixed_source, EnergySystems.sf_flexible_source]
            energy_output = extend_profile(item.energy_out, observation_period_in_years_emissions, sim_params)
            # emission calculation here

        elseif sf in [EnergySystems.sf_flexible_sink, EnergySystems.sf_fixed_sink]
            energy_input = extend_profile(item.energy_in, observation_period_in_years_emissions, sim_params)
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
