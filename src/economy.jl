using PlotlyJS, Printf

export extend_profile

"""
Holds the energy profiles of each component that are relevant for economy calculation
"""
Base.@kwdef struct EconomyEmissionsData{Component}
    component::Component
    energy_out::Union{Nothing,Vector{Float64}} = nothing
    energy_in::Union{Nothing,Vector{Float64}} = nothing
    energy_supply::Union{Nothing,Vector{Float64}} = nothing
    energy_demand::Union{Nothing,Vector{Float64}} = nothing
end

"""
Holds results of the economy calculation
"""
Base.@kwdef mutable struct EconomicResult
    total_annuity::Float64 = 0.0
    annuity_capex::Float64 = 0.0
    annuity_opex::Float64 = 0.0
    annuity_energies::Float64 = 0.0
    breakdown::Dict{String,Any} = Dict{String,Any}()
end

# Components that are only modeled by their energy flow and not with capex
ConnectionComponent = Union{GridInput,GridOutput,PVPlant,FixedSink,FixedSupply,
                            FlexibleSink,FlexibleSupply,GenericHeatSource}

function prepare_economic_emissions_data(components::Grouping,
                                         output_keys::Vector{EnergySystems.OutputKey},
                                         output_data::AbstractMatrix{<:Real})
    # Helper: Find index of required value_key in output_keys for a given component
    function find_key(component::EnergySystems.Component, output_keys::Vector{EnergySystems.OutputKey},
                      value_key::String)
        for (idx, key) in pairs(output_keys)
            if component == key.unit && value_key == key.value_key
                return idx + 1
            end
        end
        @error "The output key `$value_key` could not be retrieved from component `$(component.uac)` but is " *
               "required for the calculation of economy and/or emissions."
    end

    data = EconomyEmissionsData[]

    for (key, component) in pairs(components)
        sf = component.sys_function

        # ignore busses
        sf === EnergySystems.sf_bus && continue

        # add all other components without additional information
        entry = EconomyEmissionsData(; component=component)

        if isa(component, ConnectionComponent)
            # add fixed and flexible sources and sinks with actual input/output energies and for fixed ones also
            # with the demand or supply to calculate unmet demands/supplies.
            if sf === EnergySystems.sf_fixed_source
                entry = EconomyEmissionsData(; component=component,
                                             energy_out=copy(@view output_data[:,
                                                                               find_key(component, output_keys, "OUT")]),
                                             energy_supply=copy(@view output_data[:,
                                                                                  find_key(component, output_keys,
                                                                                           "Supply")]))
            elseif sf === EnergySystems.sf_flexible_source
                entry = EconomyEmissionsData(; component=component,
                                             energy_out=copy(@view output_data[:,
                                                                               find_key(component, output_keys, "OUT")]))

            elseif sf == EnergySystems.sf_fixed_sink
                entry = EconomyEmissionsData(; component=component,
                                             energy_in=copy(@view output_data[:,
                                                                              find_key(component, output_keys, "IN")]),
                                             energy_demand=copy(@view output_data[:,
                                                                                  find_key(component, output_keys,
                                                                                           "Demand")]))
            elseif sf == EnergySystems.sf_flexible_sink
                entry = EconomyEmissionsData(; component=component,
                                             energy_in=copy(@view output_data[:,
                                                                              find_key(component, output_keys, "IN")]))
            end
        end

        push!(data, entry)
    end

    return data
end

function calculate_economy(shared_data::Vector{EconomyEmissionsData}, sim_params::Dict{String,Any})
    # set initials
    result = EconomicResult()
    economic_parameters = sim_params["economic_parameters"]

    # iterate over all components in the current energy system
    for item in shared_data
        component = item.component
        sf = component.sys_function

        if sf === EnergySystems.sf_bus
            continue
        end
        if isa(component, ConnectionComponent)
            # get energy profiles
            if sf === EnergySystems.sf_fixed_source
                energy = extend_profile(item.energy_out, economic_parameters["observation_period_in_years"], sim_params) # source output
                energy_supply = extend_profile(item.energy_supply, economic_parameters["observation_period_in_years"],
                                               sim_params) # source supply
                energy_unmet = energy_supply .- energy # unmet supply of source
                type = :source
            elseif sf === EnergySystems.sf_flexible_source
                energy = extend_profile(item.energy_out, economic_parameters["observation_period_in_years"], sim_params) #  source output
                energy_unmet = nothing # unmet supply of source
                type = :source
            elseif sf === EnergySystems.sf_flexible_sink
                energy = extend_profile(item.energy_in, economic_parameters["observation_period_in_years"], sim_params) # sink input
                energy_unmet = nothing  # unmet demand of sink
                type = :sink
            elseif sf === EnergySystems.sf_fixed_sink
                energy = extend_profile(item.energy_in, economic_parameters["observation_period_in_years"], sim_params) # sink input
                energy_demand = extend_profile(item.energy_demand, economic_parameters["observation_period_in_years"],
                                               sim_params) # sink demand
                energy_unmet = energy_demand .- energy # unmet demand of sink
                type = :sink
            end
            # Note: Both (unmet) energies from sources AND (unmet) energies into sinks are positive at this point

            # calculate annuity for energies (demand-related costs)
            result.annuity_energies, result.breakdown = calculate_annuity_of_energies(energy,
                                                                                      component,
                                                                                      sim_params,
                                                                                      result.annuity_energies,
                                                                                      result.breakdown,
                                                                                      type)

            # calculate annuity for unmet energies
            result.annuity_energies, result.breakdown = calculate_annuity_of_unmet_energies(energy_unmet,
                                                                                            component,
                                                                                            sim_params,
                                                                                            result.annuity_energies,
                                                                                            result.breakdown,
                                                                                            type)
        end

        # calculate capex, operation-related opex and optional additional component-specific opex for all components
        result.annuity_capex, result.annuity_opex, result.breakdown = calculate_annuity_of_capex_and_opex(component,
                                                                                                          sim_params,
                                                                                                          result.annuity_capex,
                                                                                                          result.annuity_opex,
                                                                                                          result.breakdown)

        # start and end energy of storage?! TODO
    end
    result.total_annuity = result.annuity_capex + result.annuity_opex + result.annuity_energies
    return result
end

function extend_profile(profile::Union{Nothing,Vector{Float64}}, observation_period_in_years::Union{Float64,Int64},
                        sim_params::Dict{String,Any})
    if profile === nothing || profile == []
        return nothing
    end
    if sim_params["start_date_output"] + Year(observation_period_in_years) <= sim_params["end_date"]
        # no repeat required, whole period simulated
        return profile
    end

    economic_end_date = sim_params["start_date_output"] + Year(observation_period_in_years)
    if sim_params["start_date_output"] == sim_params["end_date"]
        # special case where only 1 time step is simulated
        nr_to_repeat = Int(ceil(Dates.value(Dates.Second(sub_ignoring_leap_days(economic_end_date,
                                                                                sim_params["end_date"]))) /
                                sim_params["time_step_seconds"]))
        data_to_repeat = profile
    else
        repeat_period = sim_params["economic_parameters"]["repeat_period"] # Days
        nr_to_repeat = Int(ceil(Dates.value(Dates.Second(sub_ignoring_leap_days(economic_end_date,
                                                                                sim_params["end_date"]))) /
                                Dates.value(Dates.Second(repeat_period))))
        start_idx = length(profile) + 1 -
                    Int(ceil(Dates.value(Dates.Second((repeat_period))) / sim_params["time_step_seconds"]))
        end_idx = length(profile)
        data_to_repeat = profile[start_idx:end_idx]
    end
    profile_extended = vcat(profile, repeat(data_to_repeat, nr_to_repeat))

    # cut profile to observation_period_in_years * 365 days
    number_of_timesteps = Int(ceil(observation_period_in_years * 365 * 24 * 60 * 60 / sim_params["time_step_seconds"]))
    profile_extended = profile_extended[1:number_of_timesteps]

    return profile_extended
end

function add_to_breakdown!(breakdown::Dict{String,Any}, uac::String, dict_to_add::Dict{String,<:Any})
    merge!(get!(breakdown, uac, Dict{String,Any}()), dict_to_add)
end

function calculate_annuity_of_energies(energy_profile::Vector{Float64}, component::EnergySystems.Component,
                                       sim_params::Dict{String,Any}, annuity_energies::Float64,
                                       breakdown::Dict{String,Any}, type::Symbol)
    energy_price_change_rate_per_year = component.economic_parameters["energy_price_change_rate_per_year"]
    base_cost_per_year = component.economic_parameters["base_cost_per_year"]
    base_cost_change_rate_per_year = component.economic_parameters["base_cost_change_rate_per_year"]
    observation_period = Int(sim_params["economic_parameters"]["observation_period_in_years"])

    if isnothing(component.economic_parameters["constant_energy_price"])
        # get price profile from component (one of them is given)
        step = Millisecond(Second(sim_params["time_step_seconds"]))
        times = collect(sim_params["start_date_output"]:step:sim_params["end_date"])
        times = filter(t -> !(month(t) == 2 && day(t) == 29), times)  # skip leap days
        price_profile_energy = component.economic_parameters["energy_price_profile_scale"] .*
                               [component.economic_parameters["energy_price_profile"].data[t] for t in times]
    else
        # create profile from constant energy price
        price_profile_energy = fill(component.economic_parameters["constant_energy_price"],
                                    sim_params["number_of_time_steps_output"])
    end
    price_profile_energy = extend_profile(price_profile_energy, observation_period, sim_params)

    # factors
    r_energy_costs = 1.0 + energy_price_change_rate_per_year  # price change factor of energy costs/revenues
    r_base_costs = 1.0 + base_cost_change_rate_per_year   # price change factor of base costs

    # calculate yearly costs / revenues of energies
    energy_costs_per_year = zeros(Float64, observation_period)
    energy_base_costs_per_year = zeros(Float64, observation_period)
    energy_costs_per_timestep = energy_profile .* price_profile_energy
    timesteps_per_year = Int(floor(length(energy_costs_per_timestep) / observation_period))
    for y in 1:observation_period
        energy_base_costs_per_year[y] = base_cost_per_year * r_base_costs^(y - 1)
        energy_costs_per_year[y] = sum(energy_costs_per_timestep[((y - 1) * timesteps_per_year + 1):(y * timesteps_per_year)];
                                       init=0.0) * r_energy_costs^(y - 1)
    end

    if type == :source
        # costs remain positive
        energy_costs_per_year_result = energy_costs_per_year
        energy_base_costs_per_year_result = energy_base_costs_per_year
        energy_name = "energy_costs_per_year"
        base_name = "energy_base_costs_per_year"
    elseif type == :sink
        # revenues will be negative
        energy_costs_per_year_result = .-energy_costs_per_year
        energy_base_costs_per_year_result = .-energy_base_costs_per_year
        energy_name = "energy_revenues_per_year"
        base_name = "energy_base_revenues_per_year"
    end

    # calculate annuity
    component_energies_annuity = get_annuity(energy_costs_per_year_result .+ energy_base_costs_per_year_result,
                                             sim_params)
    annuity_energies += component_energies_annuity

    add_to_breakdown!(breakdown, component.uac,
                      Dict("annuity_energies" => component_energies_annuity,
                           energy_name => energy_costs_per_year_result,
                           base_name => energy_base_costs_per_year_result))

    return annuity_energies, breakdown
end

function calculate_annuity_of_unmet_energies(unmet_energy_profile::Union{Nothing,Vector{Float64}},
                                             component::EnergySystems.Component,
                                             sim_params::Dict{String,Any}, annuity_energies::Float64,
                                             breakdown::Dict{String,Any}, type::Symbol)
    if unmet_energy_profile === nothing
        return annuity_energies, breakdown
    end

    unmet_energy_price_change_rate_per_year = component.economic_parameters["unmet_energy_price_change_rate_per_year"]
    observation_period = Int(sim_params["economic_parameters"]["observation_period_in_years"])

    if isnothing(component.economic_parameters["constant_unmet_energy_price"])
        # get price profile from component (one of them is given)
        step = Millisecond(round(Int, sim_params["time_step_seconds"] * 1000))
        times = collect(sim_params["start_date_output"]:step:sim_params["end_date"])
        times = filter(t -> !(month(t) == 2 && day(t) == 29), times)  # skip leap days
        price_profile_unmet_energy = component.economic_parameters["unmet_energy_price_profile_scale"] .*
                                     [component.economic_parameters["unmet_energy_price_profile"].data[t]
                                      for t in times]
    else
        # create profile from constant energy price
        price_profile_unmet_energy = fill(component.economic_parameters["constant_unmet_energy_price"],
                                          sim_params["number_of_time_steps_output"])
    end
    price_profile_unmet_energy = extend_profile(price_profile_unmet_energy, observation_period, sim_params)

    # factors
    r_unmet_energy_costs = 1.0 + unmet_energy_price_change_rate_per_year  # price change factor of unmet energy costs/revenues

    # calculate yearly costs / revenues of unmet energies
    unmet_energy_costs_per_year = zeros(Float64, observation_period)
    unmet_energy_costs_per_timestep = unmet_energy_profile .* price_profile_unmet_energy
    timesteps_per_year = Int(floor(length(unmet_energy_costs_per_timestep) / observation_period))
    for y in 1:observation_period
        unmet_energy_costs_per_year[y] = sum(unmet_energy_costs_per_timestep[((y - 1) * timesteps_per_year + 1):(y * timesteps_per_year)];
                                             init=0.0) * r_unmet_energy_costs^(y - 1)
    end

    if type == :source
        # fixed source: unmet energies from sink are defined as costs, so they are positive
        unmet_energy_costs_per_year_result = unmet_energy_costs_per_year
        energy_name = "unmet_energy_revenues_per_year"
    elseif type == :sink
        # fixed sink: unmet energies from demands are also costs and therefore positive
        unmet_energy_costs_per_year_result = unmet_energy_costs_per_year
        energy_name = "unmet_energy_costs_per_year"
    end

    # calculate annuity
    component_unmet_energies_annuity = get_annuity(unmet_energy_costs_per_year_result, sim_params)
    annuity_energies += component_unmet_energies_annuity

    add_to_breakdown!(breakdown, component.uac,
                      Dict("annuity_unmet_energies" => component_unmet_energies_annuity,
                           energy_name => unmet_energy_costs_per_year_result))

    return annuity_energies, breakdown
end

function get_annuity(costs::Vector{Float64}, sim_params::Dict{String,Any})::Float64
    q = 1.0 + sim_params["economic_parameters"]["interest_rate"]
    years = length(costs)

    # Timing convention: costs[1] occurs at t=0, costs[k] at t=k-1
    # Present value of the given cost stream
    present_value = 0.0
    for k in 1:years
        present_value += costs[k] / q^(k - 1)
    end

    # Annuity factor a(q,years)
    a = if abs(sim_params["economic_parameters"]["interest_rate"]) < sim_params["epsilon"]
        1.0 / years
    else
        (q^years * (q - 1.0)) / (q^years - 1.0)
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

    component_capex_annuity = get_annuity(investments_per_year .- residuals_per_year .- subsidies_per_year, sim_params)
    annuity_capex += component_capex_annuity

    # calculate opex (operation-related costs) including inspections, servicing and repair 
    # (no energy flows here, as they are handled differently)
    maintenance_per_year,
    repair_per_year,
    labour_per_year = get_opex_from_component(component, investment_first_year, sim_params)

    component_opex_annuity = get_annuity(maintenance_per_year .+ repair_per_year .+ labour_per_year, sim_params)
    annuity_opex += component_opex_annuity

    # get additional component-specific opex of material flows not represented by any interfaces
    # call component-specific function to get special opex, e.g. electrolysers water demand and oxygen production
    names, additional_opex_per_year = get_additional_opex_from_component(component, sim_params)
    for (idx, name) in enumerate(names)
        component_additional_opex_annuity = get_annuity(additional_opex_per_year[idx], sim_params)
        component_opex_annuity += component_additional_opex_annuity
        annuity_opex += component_additional_opex_annuity
        add_to_breakdown!(breakdown, component.uac, Dict(name => additional_opex_per_year[idx]))
    end

    add_to_breakdown!(breakdown, component.uac,
                      Dict("annuity_capex" => component_capex_annuity,
                           "annuity_opex" => component_opex_annuity,
                           "investments_costs_per_year" => investments_per_year,
                           "residual_revenues_per_year" => .-residuals_per_year,
                           "subsidy_revenues_per_year" => .-subsidies_per_year,
                           "maintenance_costs_per_year" => maintenance_per_year,
                           "repairs_costs_per_year" => repair_per_year,
                           "labour_costs_per_year" => labour_per_year
                           ))

    return annuity_capex, annuity_opex, breakdown
end

function get_capex_from_component(component::EnergySystems.Component, sim_params::Dict{String,Any})
    # capex including replacements, residuals and subsidies
    # parameter of the component
    lifetime_years = component.economic_parameters["lifetime_years"]
    subsidy_rate_of_capex = component.economic_parameters["subsidy_rate_of_capex"]
    subsidy_max = component.economic_parameters["subsidy_max"]
    capex_reference = get_reference_for_capex_and_embodied_emissions(component)  # e.g. installed power, area, etc.

    # parameter general
    observation_period = Int(sim_params["economic_parameters"]["observation_period_in_years"])

    # factors
    r = 1.0 + component.economic_parameters["capex_price_change_rate_per_year"]  # price change factor of capex

    # results
    investments_per_year = zeros(Float64, observation_period)
    residuals_per_year = zeros(Float64, observation_period)
    subsidies_per_year = zeros(Float64, observation_period)

    # calculate initial investment (A0)
    investment_first_year = component.economic_parameters["capex_specific"](capex_reference)   # € at time 0
    investments_per_year[1] = investment_first_year

    # subsidies as one-time event at the beginning of the observation period
    if !isnothing(subsidy_rate_of_capex)
        cap = (isnothing(subsidy_max) || subsidy_max ≤ 0) ? Inf : subsidy_max
        subsidies_first_year = min(investment_first_year * subsidy_rate_of_capex, cap)
        subsidies_per_year[1] += subsidies_first_year
    end

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
    maintenance_inspection_rate_per_year = component.economic_parameters["maintenance_inspection_rate_per_year"]
    maintenance_inspection_price_change_rate_per_year = component.economic_parameters["maintenance_inspection_price_change_rate_per_year"]
    repair_rate_per_year = component.economic_parameters["repair_rate_per_year"]
    repair_price_change_rate_per_year = component.economic_parameters["repair_price_change_rate_per_year"]
    operational_labour_hours_per_year = component.economic_parameters["operational_labour_hours_per_year"]
    labour_costs_per_hour = sim_params["economic_parameters"]["labour_costs_per_hour"]
    labour_costs_price_change_rate_per_year = sim_params["economic_parameters"]["labour_costs_price_change_rate_per_year"]
    observation_period = Int(sim_params["economic_parameters"]["observation_period_in_years"])

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

function plot_economic_results(result::EconomicResult, output_file_path::String, sim_params::Dict{String,Any},
                               cost_type::String)
    suffix = "_per_year"
    # Note: costs are positive and revenues are negative at this point! 
    # We will keep this in the figure for now...

    # collect yearly series from result struct
    series = Dict{String,Vector{Float64}}()
    for (component, component_results) in result.breakdown
        for (key, values) in component_results
            values isa AbstractVector{<:Real} || continue   # consider only vectors, no floats
            endswith(key, suffix) || continue               # consider entries with a key ending with suffix
            all(iszero, values) && continue                 # consider only entries that contain values

            name = "$(component) | $(key)"
            series[name] = copy(values)
        end
    end
    isempty(series) && return false

    if cost_type == "cashflows"
        # do nothing, values are already cashflows
        cost_name = "yearly cashflows"
    elseif cost_type == "present_values"
        # calculate net present values
        cost_name = "present values"
        q = 1.0 + sim_params["economic_parameters"]["interest_rate"]

        for (key, values) in series
            present_values = similar(values)
            for t in eachindex(values)
                present_values[t] = values[t] / q^(t - 1)
            end
            series[key] = present_values
        end
    end

    # parameter 
    observation_period_in_years = Int(sim_params["economic_parameters"]["observation_period_in_years"])
    start_year = year(sim_params["start_date_output"])

    # x-axis label
    years = collect(start_year:(start_year + observation_period_in_years - 1))

    # Net and cumulative costs
    net = zeros(observation_period_in_years)
    for v in values(series)
        net .+= v
    end
    cum = cumsum(net)

    # calculate overlays
    total_costs_over_period = 0.0
    for v in values(series)
        total_costs_over_period += sum(v)  # total costs - revenues
    end

    # Breakeven year (first sign change / zero crossing of cumulative)
    breakeven = nothing
    for t in 2:observation_period_in_years
        if cum[t - 1] == 0.0
            breakeven = start_year + t - 2
            break
        elseif cum[t - 1] * cum[t] < 0
            breakeven = start_year + t - 1
            break
        end
    end
    breakeven === nothing && cum[1] == 0.0 && (breakeven = start_year)

    # build traces
    traces = PlotlyJS.AbstractTrace[]

    # stacked bars
    for (name, v) in sort(collect(series); by=first)
        push!(traces, bar(; x=years, y=v, name=name))
    end

    # net line
    push!(traces, scatter(; x=years, y=net, mode="lines+markers",
                          name="Net per year", line=attr(; width=3)))

    # cumulative line on secondary axis
    push!(traces,
          scatter(; x=years, y=cum, mode="lines",
                  name="Cumulative net", yaxis="y2",
                  line=attr(; width=3, dash="dot")))

    # layout with dual axis and optional breakeven marker
    shapes = Any[]
    ann = Any[]
    if breakeven !== nothing
        push!(shapes,
              attr(; type="line", xref="x", yref="paper",
                   x0=breakeven, x1=breakeven, y0=0, y1=1,
                   line=attr(; width=1, dash="dot")))
        push!(ann,
              attr(; x=breakeven, y=1.02, xref="x", yref="paper",
                   text="Breakeven ≈ year $(breakeven) (in $(breakeven - start_year + 1)th year)", showarrow=false))
    end

    layout = Layout(;
                    title=attr(;
                               text="Economic results ($cost_name). Costs are positive, revenues are negative." *
                                    "<br><sup>Total yearly annuity: $(Int(round(result.total_annuity))) €/a, " *
                                    "Total costs ($cost_name) over period: $(round(total_costs_over_period; digits=0)) €</sup>"),
                    xaxis_title_text="Year",
                    yaxis_title_text="$cost_name [€/year]",
                    barmode="relative",
                    xaxis=attr(; dtick=1),
                    yaxis2=attr(; title="Cumulative $cost_name [€]", overlaying="y", side="right"),
                    legend=attr(; x=1.05, y=1.0, xanchor="left", yanchor="top"),
                    shapes=shapes,
                    annotations=ann)

    p = plot(traces, layout)
    savefig(p, output_file_path)

    return true
end

function write_economic_results_to_CSV(economic_result::EconomicResult, filepath::String, sim_params::Dict{String,Any})
    observation_period_in_years = Int(sim_params["economic_parameters"]["observation_period_in_years"])
    start_year = year(sim_params["start_date_output"])

    # write data to file
    open(filepath, "w") do file_handle
        write(file_handle, "\ufeff")  # UTF-8 BOM for Excel
        write(file_handle, "Economic results\n")
        write(file_handle, "Note: Costs are positive, revenues are negative.\n")

        write(file_handle, replace(@sprintf("Total Annuity:;%.2f;€\n", economic_result.total_annuity), "." => ","))
        write(file_handle, replace(@sprintf("Annuity of Capex:;%.2f;€\n", economic_result.annuity_capex), "." => ","))
        write(file_handle,
              replace(@sprintf("Annuity of Opex (without energies):;%.2f;€\n", economic_result.annuity_opex),
                      "." => ","))
        write(file_handle,
              replace(@sprintf("Annuity of Energies:;%.2f;€\n", economic_result.annuity_energies), "." => ","))

        write(file_handle, "\n")
        write(file_handle, "Yearly cashflows in [€/year] (without discounting):\n")

        # write years
        for year in 1:observation_period_in_years
            current_year = Int(year - 1 + start_year)
            write(file_handle, ";$(current_year)")
        end
        write(file_handle, "\n")

        # write yearly component-specific cashflows
        for (uac, component_results) in sort(collect(economic_result.breakdown); by=first)
            for (variable_name, entry) in sort(collect(component_results); by=first)
                if entry isa AbstractVector{<:Real}
                    write(file_handle, "$(uac) - $(variable_name):;")
                    write(file_handle, replace(join((@sprintf("%.2f", x) for x in entry), ";") * "\n", "." => ","))
                end
            end
        end

        # write annuities per category and component
        write(file_handle, "\n")
        write(file_handle, "Annuity breakdown per component:\n")
        for (uac, component_results) in sort(collect(economic_result.breakdown); by=first)
            for (variable_name, entry) in sort(collect(component_results); by=first)
                if entry isa Float64
                    write(file_handle, "$(uac) - $(variable_name):;")
                    write(file_handle, replace(@sprintf("%.2f;€\n", entry), "." => ","))
                end
            end
        end
    end
    return true
end

# ----------------------------------------------------------------------------------------------
# Emissions

Base.@kwdef mutable struct EmissionsResult
    total_emissions::Float64 = 0.0
    emissions_energies::Float64 = 0.0
    embodied_emissions::Float64 = 0.0
    breakdown::Dict{String,Any} = Dict{String,Any}()
end

function calculate_emissions(shared_data::Vector{EconomyEmissionsData}, sim_params::Dict{String,Any})
    # set initials
    result = EmissionsResult()
    emissions_parameters = sim_params["emissions_parameters"]

    # iterate over all components in the current energy system
    for item in shared_data
        component = item.component
        sf = component.sys_function

        if sf === EnergySystems.sf_bus
            continue
        end
        if isa(component, ConnectionComponent)
            # get energy profiles
            if sf in [EnergySystems.sf_fixed_source, EnergySystems.sf_flexible_source]
                energy = extend_profile(item.energy_out, emissions_parameters["observation_period_in_years"],
                                        sim_params) # source output
                type = :source
            elseif sf in [EnergySystems.sf_flexible_sink, EnergySystems.sf_fixed_sink]
                energy = extend_profile(item.energy_in, emissions_parameters["observation_period_in_years"], sim_params) # sink input
                type = :sink
            end
            # Note: Both energies from sources AND energies into sinks are positive at this point

            # calculate emissions and emissions credits of energies
            result.emissions_energies, result.breakdown = calculate_emissions_of_energies(energy,
                                                                                          component,
                                                                                          sim_params,
                                                                                          result.emissions_energies,
                                                                                          result.breakdown,
                                                                                          type)
        elseif emissions_parameters["include_embodied_emissions"]
            # calculate embodied emissions of transformers and storages
            result.embodied_emissions, result.breakdown = calculate_embodied_emissions(component,
                                                                                       sim_params,
                                                                                       result.embodied_emissions,
                                                                                       result.breakdown)
        end
    end
    result.total_emissions = result.emissions_energies + result.embodied_emissions
    return result
end

function calculate_emissions_of_energies(energy_profile::Vector{Float64}, component::EnergySystems.Component,
                                         sim_params::Dict{String,Any}, emissions_energies::Float64,
                                         breakdown::Dict{String,Any}, type::Symbol)
    observation_period = Int(sim_params["emissions_parameters"]["observation_period_in_years"])

    if isnothing(component.emissions_parameters["constant_energy_emissions"])
        # get price profile from component (one of them is given)
        step = Millisecond(Second(sim_params["time_step_seconds"]))
        times = collect(sim_params["start_date_output"]:step:sim_params["end_date"])
        times = filter(t -> !(month(t) == 2 && day(t) == 29), times)  # skip leap days
        emissions_profile_energy = component.emissions_parameters["energy_emissions_profile_scale"] .*
                                   [component.emissions_parameters["energy_emissions_profile"].data[t] for t in times]
    else
        # create profile from constant energy price
        emissions_profile_energy = fill(component.emissions_parameters["constant_energy_emissions"],
                                        sim_params["number_of_time_steps_output"])
    end
    emissions_profile_energy = extend_profile(emissions_profile_energy, observation_period, sim_params)

    # factors
    r_energy_emissions = 1.0 + component.emissions_parameters["energy_emissions_change_rate_per_year"]  # change factor of energy emissions

    # calculate yearly emissions of energies
    energy_emissions_per_year = zeros(Float64, observation_period)
    energy_emissions_per_timestep = energy_profile .* emissions_profile_energy
    timesteps_per_year = Int(floor(length(energy_emissions_per_timestep) / observation_period))
    for y in 1:observation_period
        energy_emissions_per_year[y] = sum(energy_emissions_per_timestep[((y - 1) * timesteps_per_year + 1):(y * timesteps_per_year)];
                                           init=0.0) * r_energy_emissions^(y - 1)
    end

    if type == :source
        # emissions remain positive
        energy_emissions_per_year_result = energy_emissions_per_year
        emissions_name = "energy_emissions_per_year"
    elseif type == :sink
        # emissions will be negative
        energy_emissions_per_year_result = .-energy_emissions_per_year
        emissions_name = "energy_emissions_credits_per_year"
    end

    # calculate annuity
    component_energies_emissions = sum(energy_emissions_per_year_result; init=0.0)
    emissions_energies += component_energies_emissions

    add_to_breakdown!(breakdown, component.uac,
                      Dict("emissions_energy" => component_energies_emissions,
                           emissions_name => energy_emissions_per_year_result))

    return emissions_energies, breakdown
end

function calculate_embodied_emissions(component::EnergySystems.Component, sim_params::Dict{String,Any},
                                      embodied_emissions::Float64, breakdown::Dict{String,Any})
    # calculate capex /(capital-related costs) including replacements, residual returns and subsidies
    embodied_emissions_per_year = get_embodied_emissions_from_component(component, sim_params)

    component_embodied_emissions = sum(embodied_emissions_per_year; init=0.0)
    embodied_emissions += component_embodied_emissions

    # no additional emissions per component for now

    add_to_breakdown!(breakdown, component.uac,
                      Dict("embodied_emissions" => component_embodied_emissions,
                           "embodied_emissions_per_year" => embodied_emissions_per_year
                           ))

    return embodied_emissions, breakdown
end

function get_embodied_emissions_from_component(component::EnergySystems.Component, sim_params::Dict{String,Any})
    # embodied emissions including replacements and residuals
    # parameter of the component
    lifetime_years = component.emissions_parameters["lifetime_years"]
    emissions_reference = get_reference_for_capex_and_embodied_emissions(component)  # e.g. installed power, area, etc.

    # parameter general
    observation_period = Int(sim_params["emissions_parameters"]["observation_period_in_years"])

    # factors
    r = 1.0 + component.emissions_parameters["embodied_emissions_change_rate_per_year"]  # embodied emissions change factor per year

    # results
    emissions_per_year = zeros(Float64, observation_period)

    # calculate initial emissions
    emissions_first_year = component.emissions_parameters["embodied_emissions_specific"] * emissions_reference # € at time 0
    emissions_per_year[1] = emissions_first_year

    # calculate emissions of replacements
    # number of replacements n within observation period
    n = max(0, Int(floor((observation_period - 1e-9) / lifetime_years)))
    for j in 1:n
        t = j * lifetime_years      # years since start
        idx = Int(round(t)) + 1     # index for t=0 -> 1, t=1 -> 2, ...

        if 1 <= idx <= observation_period
            emissions_per_year[idx] += emissions_first_year * r^t
        end
    end

    # calculate residual value
    t_last = n * lifetime_years
    remaining_fraction = max(0.0, ((n + 1) * lifetime_years - observation_period) / lifetime_years)
    emissions_per_year[observation_period] = (emissions_first_year * r^t_last) * remaining_fraction

    return emissions_per_year
end

function plot_emissions_results(result::EmissionsResult,
                                output_file_path::String,
                                sim_params::Dict{String,Any})
    suffix = "_per_year"
    emissions_unit = "kgCO₂e"
    emissions_factor = 1e-3
    # Convention:
    #   positive values  -> emissions
    #   negative values  -> credits / avoided emissions

    # parameters
    observation_period_in_years = Int(sim_params["emissions_parameters"]["observation_period_in_years"])
    start_year = year(sim_params["start_date_output"])
    years = collect(start_year:(start_year + observation_period_in_years - 1))

    # collect yearly series from result struct
    series = Dict{String,Vector{Float64}}()

    for (component, component_results) in result.breakdown
        component_results isa AbstractDict || continue

        for (key, values) in component_results
            values isa AbstractVector{<:Real} || continue
            endswith(String(key), suffix) || continue
            all(iszero, values) && continue

            v = Float64.(values) .* emissions_factor
            name = "$(component) | $(key)"
            series[name] = v
        end
    end

    isempty(series) && return false

    # net and cumulative emissions
    net = zeros(observation_period_in_years)
    for v in values(series)
        net .+= v
    end
    cum = cumsum(net)

    total_emissions_over_period = sum(net)

    # build traces
    traces = PlotlyJS.AbstractTrace[]

    # stacked bars
    for (name, v) in sort(collect(series); by=first)
        push!(traces, bar(; x=years, y=v, name=name))
    end

    # net line
    push!(traces,
          scatter(; x=years, y=net, mode="lines+markers",
                  name="Net emissions per year",
                  line=attr(; width=3)))

    # cumulative line on secondary axis
    push!(traces,
          scatter(; x=years, y=cum, mode="lines",
                  name="Cumulative net emissions", yaxis="y2",
                  line=attr(; width=3, dash="dot")))

    layout = Layout(;
                    title=attr(;
                               text="Emissions results. Positive values are emissions, negative values are credits." *
                                    "<br><sup>Total emissions: $(round(emissions_factor*result.total_emissions; digits=2)) $(emissions_unit), " *
                                    "Energy emissions: $(round(emissions_factor*result.emissions_energies; digits=2)) $(emissions_unit), " *
                                    "Embodied emissions: $(round(emissions_factor*result.embodied_emissions; digits=2)) $(emissions_unit), " *
                                    "Net sum of plotted yearly values: $(round(total_emissions_over_period; digits=2)) $(emissions_unit)</sup>"),
                    xaxis_title_text="Year",
                    yaxis_title_text="Yearly emissions [$emissions_unit/year]",
                    barmode="relative",
                    xaxis=attr(; dtick=1),
                    yaxis2=attr(;
                                title="Cumulative emissions [$emissions_unit]",
                                overlaying="y",
                                side="right"),
                    legend=attr(; x=1.05, y=1.0, xanchor="left", yanchor="top"))

    p = plot(traces, layout)
    savefig(p, output_file_path)

    return true
end

function write_emissions_results_to_CSV(emissions_result::EmissionsResult,
                                        filepath::String,
                                        sim_params::Dict{String,Any})
    emissions_unit = "kgCO₂e"
    emissions_factor = 1e-3

    observation_period_in_years = Int(sim_params["emissions_parameters"]["observation_period_in_years"])
    start_year = year(sim_params["start_date_output"])

    open(filepath, "w") do file_handle
        write(file_handle, "\ufeff")  # UTF-8 BOM for Excel
        write(file_handle, "Emissions results\n")
        write(file_handle, "Note: Positive values are emissions, negative values are credits / avoided emissions.\n")

        write(file_handle,
              replace(@sprintf("Total emissions:;%.2f;%s\n",
                               emissions_factor * emissions_result.total_emissions, emissions_unit), "." => ","))
        write(file_handle,
              replace(@sprintf("Energy emissions:;%.2f;%s\n",
                               emissions_factor * emissions_result.emissions_energies, emissions_unit), "." => ","))
        write(file_handle,
              replace(@sprintf("Embodied emissions:;%.2f;%s\n",
                               emissions_factor * emissions_result.embodied_emissions, emissions_unit), "." => ","))

        write(file_handle, "\n")
        write(file_handle, "Yearly emissions in [$(emissions_unit)/year]:\n")

        # write years
        for year in 1:observation_period_in_years
            current_year = Int(year - 1 + start_year)
            write(file_handle, ";$(current_year)")
        end
        write(file_handle, "\n")

        # write yearly component-specific emissions series
        for (component, component_results) in sort(collect(emissions_result.breakdown); by=first)
            component_results isa AbstractDict || continue

            for (variable_name, entry) in sort(collect(component_results); by=first)
                if entry isa AbstractVector{<:Real}
                    values = emissions_factor .* Float64.(entry)

                    # make vector length consistent with observation period
                    if length(values) < observation_period_in_years
                        values = vcat(values, zeros(observation_period_in_years - length(values)))
                    elseif length(values) > observation_period_in_years
                        values = values[1:observation_period_in_years]
                    end

                    write(file_handle, "$(component) - $(variable_name):;")
                    write(file_handle,
                          replace(join((@sprintf("%.2f", x) for x in values), ";") * "\n", "." => ","))
                end
            end
        end

        # write scalar breakdown per category and component
        write(file_handle, "\n")
        write(file_handle, "Emissions breakdown per component:\n")

        for (component, component_results) in sort(collect(emissions_result.breakdown); by=first)
            component_results isa AbstractDict || continue

            for (variable_name, entry) in sort(collect(component_results); by=first)
                if entry isa Real && !(entry isa AbstractVector)
                    write(file_handle, "$(component) - $(variable_name):;")
                    write(file_handle,
                          replace(@sprintf("%.2f;%s\n", emissions_factor * Float64(entry), emissions_unit), "." => ","))
                end
            end
        end
    end

    return true
end
