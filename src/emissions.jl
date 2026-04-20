using PlotlyJS, Printf

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
