using PlotlyJS, Printf

# Calculate emissions balance, including emissions and credits of energies and embodied emission.
# Replacement investments of components induce new embodied emissions, while 
# possible remaining usable time of a components creates a credit of emissions. Note that no
# explicit end-of-life emissions are included, they need to be part of the total embodied emissions.
#
# Convention in output:
# emissions to the environment: positive values
# credits / avoided emissions: negative values

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
            # get energy profiles of connection components
            if sf in [EnergySystems.sf_fixed_source, EnergySystems.sf_flexible_source]
                energy = extend_profile(item.energy_out, emissions_parameters["observation_period_in_years"],
                                        emissions_parameters["repeat_period"], sim_params) # source output
                type = :source
            elseif sf in [EnergySystems.sf_flexible_sink, EnergySystems.sf_fixed_sink]
                energy = extend_profile(item.energy_in, emissions_parameters["observation_period_in_years"],
                                        emissions_parameters["repeat_period"], sim_params) # sink input
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
        end

        if emissions_parameters["include_embodied_emissions"]
            # calculate embodied emissions of all components
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

    if type == :source
        param_name = "energy_emissions"
    elseif type == :sink
        param_name = "energy_emissions_credits"
    end

    if isnothing(component.emissions_parameters["constant_" * param_name])
        # get emission profile from component (one of them is given)
        step = Millisecond(Second(sim_params["time_step_seconds"]))
        profile_end_date = maximum(keys(component.emissions_parameters[param_name * "_profile"].data))
        times = collect(sim_params["start_date_output"]:step:profile_end_date)
        times = filter(t -> !(month(t) == 2 && day(t) == 29), times)  # skip leap days
        emissions_profile_energy = component.emissions_parameters[param_name * "_profile_scale"] .*
                                   [component.emissions_parameters[param_name * "_profile"].data[t] for t in times]

        # Deal with cases where the price profile is longer than the simulation period.
        if (sub_ignoring_leap_days(profile_end_date, sim_params["start_date"]) >
            sim_params["emissions_parameters"]["repeat_period"]) &&
           sim_params["emissions_parameters"]["repeat_method"] === "all"
            emissions_repeat_period = sub_ignoring_leap_days(profile_end_date, sim_params["start_date_output"])
        else
            emissions_repeat_period = sim_params["emissions_parameters"]["repeat_period"]
        end
    else
        # create profile from constant emissions
        emissions_profile_energy = fill(component.emissions_parameters["constant_" * param_name],
                                        sim_params["number_of_time_steps_output"])
        emissions_repeat_period = sim_params["emissions_parameters"]["repeat_period"]
    end
    emissions_profile_energy = extend_profile(emissions_profile_energy, observation_period, emissions_repeat_period,
                                              sim_params)

    emissions_profile_energy_effective = apply_yearly_emissions_change_to_profile(emissions_profile_energy,
                                                                                  observation_period,
                                                                                  component.emissions_parameters[param_name * "_change_rate_per_year"])

    # calculate yearly emissions of energies
    energy_emissions_per_year = zeros(Float64, observation_period)
    energy_emissions_per_timestep = energy_profile .* emissions_profile_energy_effective
    timesteps_per_year = Int(floor(length(energy_emissions_per_timestep) / observation_period))
    for y in 1:observation_period
        energy_emissions_per_year[y] = sum(energy_emissions_per_timestep[((y - 1) * timesteps_per_year + 1):(y * timesteps_per_year)];
                                           init=0.0)
    end

    if type == :source
        # emissions remain positive
        energy_emissions_per_year_result = energy_emissions_per_year
        emissions_name = "emissions_energy"
        emissions_name_per_year = "energy_emissions_per_year"

        # save profile for later use
        store_extended_emissions_profile!(breakdown,
                                          component,
                                          "energy_emissions_profile_effective",
                                          emissions_profile_energy_effective)
    elseif type == :sink
        # emissions will be negative
        energy_emissions_per_year_result = .-energy_emissions_per_year
        emissions_name = "emissions_credits_energy"
        emissions_name_per_year = "energy_emissions_credits_per_year"

        # save profile for later use
        store_extended_emissions_profile!(breakdown,
                                          component,
                                          "energy_emissions_credits_profile_effective",
                                          .-emissions_profile_energy_effective)
    end

    # calculate annuity
    component_energies_emissions = sum(energy_emissions_per_year_result; init=0.0)
    emissions_energies += component_energies_emissions

    add_to_breakdown!(breakdown, component.uac,
                      Dict(emissions_name => component_energies_emissions,
                           emissions_name_per_year => energy_emissions_per_year_result))

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
    emissions_first_year = component.emissions_parameters["embodied_emissions_specific"](emissions_reference) # € at time 0
    emissions_per_year[1] = emissions_first_year

    # calculate emissions of replacements and treat them as credits
    # number of replacements n within observation period
    n = max(0, Int(floor((observation_period - 1e-9) / lifetime_years)))
    for j in 1:n
        t = j * lifetime_years      # years since start
        idx = Int(round(t)) + 1     # index for t=0 -> 1, t=1 -> 2, ...

        if 1 <= idx <= observation_period
            emissions_per_year[idx] += emissions_first_year * r^t
        end
    end

    # calculate residual emissions
    t_last = n * lifetime_years
    remaining_fraction = max(0.0, ((n + 1) * lifetime_years - observation_period) / lifetime_years)
    emissions_per_year[observation_period] -= (emissions_first_year * r^t_last) * remaining_fraction

    return emissions_per_year
end

function store_extended_emissions_profile!(breakdown::Dict{String,Any},
                                           component::EnergySystems.Component,
                                           profile_name::String,
                                           profile::AbstractVector{<:Real})
    component_results = get!(breakdown, component.uac, Dict{String,Any}())

    profiles_any = get!(component_results,
                        "extended_emissions_profiles",
                        Dict{String,Vector{Float64}}())

    profiles = profiles_any::Dict{String,Vector{Float64}}
    profiles[profile_name] = Float64.(profile)

    return nothing
end

function apply_yearly_emissions_change_to_profile(profile::AbstractVector{<:Real},
                                                  observation_period::Int,
                                                  change_rate_per_year::Real)
    effective_profile = Float64.(profile)

    observation_period <= 0 && return effective_profile

    r = 1.0 + Float64(change_rate_per_year)
    timesteps_per_year = Int(floor(length(effective_profile) / observation_period))

    timesteps_per_year <= 0 && return effective_profile

    for y in 1:observation_period
        idx_start = (y - 1) * timesteps_per_year + 1
        idx_end = y * timesteps_per_year

        idx_start > length(effective_profile) && break
        idx_end = min(idx_end, length(effective_profile))

        effective_profile[idx_start:idx_end] .*= r^(y - 1)
    end

    return effective_profile
end

function plot_emissions_results(result::EmissionsResult,
                                output_file_path::String,
                                sim_params::Dict{String,Any},
                                fixed_output_precision::Int)
    suffix = "_per_year"
    # unit and factor for output. Note that internally, everything is based on [g/W] or [g/Wh]
    emissions_unit = "kg CO₂e"
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

    # Optional fixed output precision for plotted values.
    round_for_plot(v) = fixed_output_precision > 0 ? round.(v; digits=fixed_output_precision) : v

    # net and cumulative emissions
    net = zeros(observation_period_in_years)
    for v in values(series)
        net .+= v
    end
    cum = cumsum(net)

    # build traces
    traces = PlotlyJS.AbstractTrace[]

    # stacked bars
    for (name, v) in sort(collect(series); by=first)
        push!(traces, bar(; x=years, y=round_for_plot(v), name=name))
    end

    # net line
    push!(traces,
          scatter(; x=years, y=round_for_plot(net), mode="lines+markers",
                  name="Net emissions per year",
                  line=attr(; width=3)))

    # cumulative line on secondary axis
    push!(traces,
          scatter(; x=years, y=round_for_plot(cum), mode="lines",
                  name="Cumulative net emissions", yaxis="y2",
                  line=attr(; width=3, dash="dot")))

    layout = Layout(;
                    title=attr(;
                               text="Emissions results. Positive values are emissions, negative values are credits." *
                                    "<br><sup>Total emissions: $(round(emissions_factor*result.total_emissions; digits=2)) $(emissions_unit), " *
                                    "Energy emissions: $(round(emissions_factor*result.emissions_energies; digits=2)) $(emissions_unit), " *
                                    "Embodied emissions: $(round(emissions_factor*result.embodied_emissions; digits=2)) $(emissions_unit)</sup>"),
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
    # unit and factor for output. Note that internally, everything is based on [g/W] or [g/Wh]
    emissions_unit = "kg CO₂e"
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
