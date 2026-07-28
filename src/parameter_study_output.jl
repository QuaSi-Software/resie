# this file contains functionality for writing output and creating plots of the parameter-study results
# including parameter variation, optimisation and sensitivity analysis.

const PARAMETER_STUDY_PLOT_INSERTIONS_PATH = joinpath(@__DIR__, "parameter_study_plots_insertions")

"""
    render_parameter_study_plot_insertion(file_name, replacements)

Takes a separate html insertion file for interactive plots and processes it.
"""
function render_parameter_study_plot_insertion(file_name::String,
                                               replacements::Pair{String,String}...)::String
    insertion_path = joinpath(PARAMETER_STUDY_PLOT_INSERTIONS_PATH, file_name)
    insertion = read(insertion_path, String)

    for (name, value) in replacements
        placeholder = "{{$name}}"
        if !occursin(placeholder, insertion)
            @error "Missing placeholder $placeholder in parameter-study plot insertion $file_name."
        end
        insertion = replace(insertion, placeholder => value)
    end

    unresolved_placeholder = match(r"\{\{[A-Z0-9_]+\}\}", insertion)
    if !isnothing(unresolved_placeholder)
        @error "Unresolved placeholder $(unresolved_placeholder.match) in parameter-study plot insertion $file_name."
    end
    return insertion
end

"""
    parameter_study_csv_path(sim_params, io_settings, suffix)

Create a CSV path derived from `parameter_study_csv_path`.
"""
function parameter_study_csv_path(sim_params::Dict{String,Any},
                                  io_settings::Dict{String,Any},
                                  suffix::String)::String
    base_path = sim_params["run_path"](io_settings["parameter_study_csv_path"])
    directory, filename = splitdir(base_path)
    root, _ = splitext(filename)
    return joinpath(directory, "$(root)_$(suffix).csv")
end

function sensitivity_csv_string(value)::String
    if value isa Real
        return replace(string(Float64(value)), "." => ",")
    end
    return string(value)
end

function sensitivity_csv_escape(value)::String
    text = sensitivity_csv_string(value)
    if occursin(';', text) || occursin('"', text) || occursin('\n', text) || occursin('\r', text)
        return "\"" * replace(text, "\"" => "\"\"") * "\""
    end
    return text
end

function write_sensitivity_csv(file_path::String,
                               header::Vector{String},
                               rows::Vector{<:AbstractVector})::String
    directory = dirname(file_path)
    isempty(directory) || mkpath(directory)

    open(file_path, "w") do io
        println(io, join(sensitivity_csv_escape.(header), ';'))
        for row in rows
            println(io, join(sensitivity_csv_escape.(row), ';'))
        end
    end

    return file_path
end

"""
    write_global_sensitivity_csv(file_path, parameter_keys, bounds, S_total, S_first, rel_rmse, r2)

Write one row per parameter. `interaction_effect` is the part of the total-order Sobol index
that is not explained by the first-order effect.
"""
function write_global_sensitivity_csv(file_path::String,
                                      parameter_keys::Vector{String},
                                      bounds::AbstractMatrix{<:Real},
                                      S_total::Vector{Float64},
                                      S_first::Vector{Float64},
                                      rel_rmse::Float64,
                                      r2::Float64)::String
    header = ["parameter",
              "bounds_min",
              "bounds_max",
              "S_first",
              "S_total",
              "interaction_effect",
              "surrogate_relative_rmse",
              "surrogate_r2"]

    rows = Vector{Vector{Any}}()
    for i in eachindex(parameter_keys)
        push!(rows,
              Any[parameter_keys[i],
                  bounds[i, 1],
                  bounds[i, 2],
                  S_first[i],
                  S_total[i],
                  S_total[i] - S_first[i],
                  rel_rmse,
                  r2])
    end

    return write_sensitivity_csv(file_path, header, rows)
end

"""
    write_local_sensitivity_csv(file_path, sensitivity_results)

Write one row per varied parameter. Relative parameter and objective changes are written as
percent values to make the CSV directly comparable with the local-sensitivity plots.
"""
function write_local_sensitivity_csv(file_path::String,
                                     sensitivity_results::Vector{Dict{String,Any}})::String
    header = ["parameter",
              "lower_value",
              "reference_value",
              "upper_value",
              "lower_parameter_change_percent",
              "upper_parameter_change_percent",
              "lower_objective",
              "reference_objective",
              "upper_objective",
              "lower_objective_change",
              "upper_objective_change",
              "lower_objective_change_percent",
              "upper_objective_change_percent",
              "gradient",
              "elasticity"]

    rows = Vector{Vector{Any}}()
    for result in sensitivity_results
        reference_value = Float64(result["reference_value"])
        lower_parameter_change_percent = reference_value == 0.0 ? NaN :
                                         100.0 * (Float64(result["lower_value"]) - reference_value) /
                                         abs(reference_value)
        upper_parameter_change_percent = reference_value == 0.0 ? NaN :
                                         100.0 * (Float64(result["upper_value"]) - reference_value) /
                                         abs(reference_value)

        push!(rows,
              Any[result["parameter"],
                  result["lower_value"],
                  reference_value,
                  result["upper_value"],
                  lower_parameter_change_percent,
                  upper_parameter_change_percent,
                  result["lower_objective"],
                  result["reference_objective"],
                  result["upper_objective"],
                  result["lower_change"],
                  result["upper_change"],
                  100.0 * Float64(result["lower_change_relative"]),
                  100.0 * Float64(result["upper_change_relative"]),
                  result["gradient"],
                  result["elasticity"]])
    end

    return write_sensitivity_csv(file_path, header, rows)
end

# Format percentages consistently in plot labels. A decimal comma is used to keep
# compact labels such as `+2,34 %` easy to read.
function format_plot_percent(value::Real; signed::Bool=false)::String
    display_value = abs(Float64(value)) < 0.005 ? 0.0 : Float64(value)
    formatted = if signed
        @sprintf("%+.2f %%", display_value)
    else
        @sprintf("%.2f %%", display_value)
    end
    return replace(formatted, "." => ",")
end

"""
    create_global_sensitivity_plot(parameter_keys, S_total, S_first, rel_rmse, r2,
                                   io_settings, sim_params)

Create a vertical comparison of first-order and total-order Sobol indices. Parameters are
ordered by their total-order index. The gap between `S_first` and `S_total` indicates
interaction effects with other parameters.
"""
function create_global_sensitivity_plot(parameter_keys::Vector{String},
                                        S_total::Vector{Float64},
                                        S_first::Vector{Float64},
                                        rel_rmse::Float64,
                                        r2::Float64,
                                        io_settings::Dict{String,Any},
                                        sim_params::Dict{String,Any})::String
    if isempty(parameter_keys) ||
       length(parameter_keys) != length(S_total) ||
       length(parameter_keys) != length(S_first)
        @error "Cannot create global-sensitivity plot: result dimensions do not match."
        return ""
    end

    order = sortperm(S_total; rev=true)
    names = parameter_keys[order]
    first_order = S_first[order]
    total_order = S_total[order]

    first_labels = [format_plot_percent(100.0 * value) for value in first_order]
    total_labels = [format_plot_percent(100.0 * value) for value in total_order]

    first_trace = bar(; x=names,
                      y=first_order,
                      name="First-order effect (S_first)",
                      text=first_labels,
                      textposition="outside",
                      cliponaxis=false,
                      hovertemplate="%{x}<br>S_first = %{y:.5f}<extra></extra>")

    total_trace = bar(; x=names,
                      y=total_order,
                      name="Total effect (S_total)",
                      text=total_labels,
                      textposition="outside",
                      cliponaxis=false,
                      hovertemplate="%{x}<br>S_total = %{y:.5f}<extra></extra>")

    quality_text = "Surrogate relative RMSE=$(round(rel_rmse; digits=3)), " *
                   "R²=$(round(r2; digits=3))"
    maximum_index = maximum(vcat(first_order, total_order))
    y_upper = max(0.05, 1.18 * maximum_index)

    layout = Layout(;
                    title=attr(;
                               text="Global sensitivity: Sobol indices" *
                                    "<br><sup>The gap between S_first and S_total indicates interactions; $quality_text</sup>",
                               x=0.5,
                               xanchor="center"),
                    barmode="group",
                    xaxis=attr(; title="Parameter",
                               tickangle=-30,
                               automargin=true),
                    yaxis=attr(; title="Sobol sensitivity index",
                               zeroline=true,
                               range=[0.0, y_upper]),
                    legend=attr(; orientation="h",
                                x=0.5,
                                xanchor="center",
                                y=1.02,
                                yanchor="bottom"),
                    margin=attr(; t=155,
                                b=190,
                                l=90,
                                r=45),
                    height=720,
                    autosize=true)

    plot_object = plot(GenericTrace[first_trace, total_trace], layout)
    file_path = parameter_study_plot_path(sim_params,
                                          io_settings,
                                          "global_sensitivity")
    mkpath(dirname(file_path))
    savefig(plot_object, file_path)
    return file_path
end

"""
    create_local_sensitivity_response_overview_plot(sensitivity_results, io_settings, sim_params)

Create a normalized dumbbell plot of the local-sensitivity response. Each parameter is shown
on one row. The lower and upper markers show the corresponding relative objective changes,
and the connecting line makes the response magnitude and asymmetry visible. A vertical line at
zero marks the reference objective.

The absolute three-point response panels are additionally written by
`create_local_sensitivity_absolute_response_plot`.
"""
function create_local_sensitivity_response_overview_plot(sensitivity_results::Vector{Dict{String,Any}},
                                                         io_settings::Dict{String,Any},
                                                         sim_params::Dict{String,Any})::String
    valid_results = filter(sensitivity_results) do result
        values = (result["lower_value"],
                  result["reference_value"],
                  result["upper_value"],
                  result["lower_objective"],
                  result["reference_objective"],
                  result["upper_objective"],
                  result["lower_change_relative"],
                  result["upper_change_relative"])
        all(value -> value isa Real && isfinite(value), values)
    end

    if isempty(valid_results)
        @error "Cannot create local-sensitivity response overview plot: no finite response points."
        return ""
    end

    importance = [max(abs(Float64(result["lower_change_relative"])),
                      abs(Float64(result["upper_change_relative"])))
                  for result in valid_results]
    order = sortperm(importance; rev=true)
    ordered_results = valid_results[order]

    names = String[result["parameter"] for result in ordered_results]
    lower_changes = Float64[100.0 * result["lower_change_relative"]
                            for result in ordered_results]
    upper_changes = Float64[100.0 * result["upper_change_relative"]
                            for result in ordered_results]
    lower_values = Float64[result["lower_value"] for result in ordered_results]
    upper_values = Float64[result["upper_value"] for result in ordered_results]
    lower_objectives = Float64[result["lower_objective"] for result in ordered_results]
    upper_objectives = Float64[result["upper_objective"] for result in ordered_results]

    traces = GenericTrace[]

    # One connector per parameter keeps the relationship between lower and upper
    # responses explicit without introducing another scale.
    for i in eachindex(names)
        push!(traces,
              scatter(; x=[lower_changes[i], upper_changes[i]],
                      y=[names[i], names[i]],
                      mode="lines",
                      line=attr(; width=3,
                                color="rgba(120,120,120,0.55)"),
                      hoverinfo="skip",
                      showlegend=false))
    end

    lower_labels = [format_plot_percent(value; signed=true) for value in lower_changes]
    upper_labels = [format_plot_percent(value; signed=true) for value in upper_changes]

    # Position labels relative to the actual left-to-right order of the two
    # response points. This keeps labels outside the connector even when the
    # lower and upper responses switch sides. Points that are close on the
    # common x-axis scale are separated vertically to prevent label overlap.
    all_changes = vcat(lower_changes, upper_changes, [0.0])
    change_minimum = minimum(all_changes)
    change_maximum = maximum(all_changes)
    change_span = max(change_maximum - change_minimum,
                      maximum(abs.(all_changes)),
                      1.0)

    lower_label_positions = String[]
    upper_label_positions = String[]

    for i in eachindex(names)
        if lower_changes[i] <= upper_changes[i]
            push!(lower_label_positions, "middle left")
            push!(upper_label_positions, "middle right")
        else
            push!(lower_label_positions, "middle right")
            push!(upper_label_positions, "middle left")
        end
    end

    lower_hover = ["$(names[i])<br>" *
                   "Point: Lower<br>" *
                   "Parameter value: $(lower_values[i])<br>" *
                   "Objective: $(lower_objectives[i])<br>" *
                   "Objective change: $(round(lower_changes[i]; digits=4)) %"
                   for i in eachindex(names)]

    upper_hover = ["$(names[i])<br>" *
                   "Point: Upper<br>" *
                   "Parameter value: $(upper_values[i])<br>" *
                   "Objective: $(upper_objectives[i])<br>" *
                   "Objective change: $(round(upper_changes[i]; digits=4)) %"
                   for i in eachindex(names)]

    push!(traces,
          scatter(; x=lower_changes,
                  y=names,
                  mode="markers+text",
                  name="Lower parameter value",
                  marker=attr(; size=11,
                              symbol="circle"),
                  text=lower_labels,
                  textposition=lower_label_positions,
                  cliponaxis=false,
                  hovertext=lower_hover,
                  hovertemplate="%{hovertext}<extra></extra>"))

    push!(traces,
          scatter(; x=upper_changes,
                  y=names,
                  mode="markers+text",
                  name="Upper parameter value",
                  marker=attr(; size=11,
                              symbol="diamond"),
                  text=upper_labels,
                  textposition=upper_label_positions,
                  cliponaxis=false,
                  hovertext=upper_hover,
                  hovertemplate="%{hovertext}<extra></extra>"))

    x_padding = 0.25 * change_span

    reference_objective = Float64(first(ordered_results)["reference_objective"])
    layout = Layout(;
                    title=attr(;
                               text="Local sensitivity responses overview" *
                                    "<br><sup>Relative objective change from the reference " *
                                    "($(round(reference_objective; sigdigits=8))).",
                               x=0.5,
                               xanchor="center"),
                    xaxis=attr(; title="Objective change [%]",
                               zeroline=true,
                               zerolinewidth=2,
                               range=[change_minimum - x_padding,
                                      change_maximum + x_padding]),
                    yaxis=attr(; title="Parameter",
                               autorange="reversed",
                               automargin=true),
                    legend=attr(; orientation="h",
                                x=0.5,
                                xanchor="center",
                                y=1.02,
                                yanchor="bottom"),
                    margin=attr(; t=150,
                                b=90,
                                l=230,
                                r=80),
                    height=max(560, 72 * length(names) + 230),
                    autosize=true,
                    hovermode="closest")

    plot_object = plot(traces, layout)
    file_path = parameter_study_plot_path(sim_params,
                                          io_settings,
                                          "local_sensitivity_response_overview")
    mkpath(dirname(file_path))
    savefig(plot_object, file_path)

    return file_path
end

"""
    create_local_sensitivity_response_trends_plot(sensitivity_results, io_settings, sim_params)

Create one vertically stacked three-point absolute response plot for every parameter. Each
panel uses the parameter value on the x-axis and the absolute objective value on the y-axis.
Separate y-axis ranges prevent large responses from hiding smaller effects.
"""
function create_local_sensitivity_response_trends_plot(sensitivity_results::Vector{Dict{String,Any}},
                                                       io_settings::Dict{String,Any},
                                                       sim_params::Dict{String,Any})::String
    valid_results = filter(sensitivity_results) do result
        values = (result["lower_value"],
                  result["reference_value"],
                  result["upper_value"],
                  result["lower_objective"],
                  result["reference_objective"],
                  result["upper_objective"])
        all(value -> value isa Real && isfinite(value), values)
    end

    if isempty(valid_results)
        @error "Cannot create absolute local-sensitivity response trend plot: no finite response points."
        return ""
    end

    n_parameters = length(valid_results)
    vertical_spacing = n_parameters > 1 ? min(0.045, 0.16 / (n_parameters - 1)) : 0.0
    panel_height = (1.0 - (n_parameters - 1) * vertical_spacing) / n_parameters

    traces = GenericTrace[]
    shapes = Any[]
    layout_values = Dict{Symbol,Any}()

    for (i, result) in enumerate(valid_results)
        parameter = String(result["parameter"])
        parameter_values = Float64[result["lower_value"],
                                   result["reference_value"],
                                   result["upper_value"]]
        objective_values = Float64[result["lower_objective"],
                                   result["reference_objective"],
                                   result["upper_objective"]]
        point_names = ["Lower", "Reference", "Upper"]
        objective_changes = Float64[100.0 * result["lower_change_relative"],
                                    0.0,
                                    100.0 * result["upper_change_relative"]]
        hover_text = ["$parameter<br>" *
                      "Point: $(point_names[j])<br>" *
                      "Parameter value: $(parameter_values[j])<br>" *
                      "Objective: $(objective_values[j])<br>" *
                      "Objective change: $(round(objective_changes[j]; digits=4)) %"
                      for j in eachindex(point_names)]
        point_labels = [format_plot_percent(objective_changes[1]; signed=true),
                        "",
                        format_plot_percent(objective_changes[3]; signed=true)]
        label_positions = [objective_changes[1] >= 0.0 ? "top center" : "bottom center",
                           "top center",
                           objective_changes[3] >= 0.0 ? "top center" : "bottom center"]

        objective_minimum = minimum(objective_values)
        objective_maximum = maximum(objective_values)
        objective_span = objective_maximum - objective_minimum
        objective_padding = objective_span > eps(Float64) ?
                            0.18 * objective_span :
                            max(0.01 * abs(objective_values[2]), 1.0)

        x_reference = i == 1 ? "x" : "x$i"
        y_reference = i == 1 ? "y" : "y$i"
        x_layout_key = i == 1 ? :xaxis : Symbol("xaxis$i")
        y_layout_key = i == 1 ? :yaxis : Symbol("yaxis$i")

        domain_end = 1.0 - (i - 1) * (panel_height + vertical_spacing)
        domain_start = domain_end - panel_height

        push!(traces,
              scatter(; x=parameter_values,
                      y=objective_values,
                      mode="lines+markers+text",
                      line=attr(; width=3),
                      marker=attr(; size=10),
                      text=point_labels,
                      textposition=label_positions,
                      textfont=attr(; size=12),
                      cliponaxis=false,
                      hovertext=hover_text,
                      hovertemplate="%{hovertext}<extra></extra>",
                      xaxis=x_reference,
                      yaxis=y_reference,
                      showlegend=false))

        push!(shapes,
              attr(; type="line",
                   xref=x_reference,
                   yref=y_reference,
                   x0=minimum(parameter_values),
                   x1=maximum(parameter_values),
                   y0=objective_values[2],
                   y1=objective_values[2],
                   line=attr(; width=1,
                             dash="dot")))

        layout_values[x_layout_key] = attr(; domain=[0.0, 1.0],
                                           anchor=y_reference,
                                           title=parameter,
                                           automargin=true,
                                           tickformat=".7g")
        layout_values[y_layout_key] = attr(; domain=[domain_start, domain_end],
                                           anchor=x_reference,
                                           title="Objective",
                                           automargin=true,
                                           tickformat=".7g",
                                           range=[objective_minimum - objective_padding,
                                                  objective_maximum + objective_padding])
    end

    layout_values[:title] = attr(;
                                 text="Trend of local sensitivity responses" *
                                      "<br><sup>The dotted line marks the reference objective.</sup>",
                                 x=0.5,
                                 xanchor="center")
    layout_values[:shapes] = shapes
    layout_values[:showlegend] = false
    layout_values[:hovermode] = "closest"
    layout_values[:margin] = attr(; t=105,
                                  b=80,
                                  l=115,
                                  r=55)
    layout_values[:height] = max(650, 340 * n_parameters)
    layout_values[:autosize] = true

    plot_object = plot(traces, Layout(; layout_values...))
    file_path = parameter_study_plot_path(sim_params,
                                          io_settings,
                                          "local_sensitivity_response_trends")
    mkpath(dirname(file_path))
    savefig(plot_object, file_path)
    return file_path
end

"""
    create_matrix_plot(results, io_settings, sim_params; ...)

Create a lower-triangular parameter-study matrix.

The diagonal contains histograms of evaluated parameter values. Pairwise scatter plots are
shown below the diagonal; the redundant upper triangle remains empty. Scatter colours use the
selected parameter-study parameter, objective, objective component or scalar result quantity. For a single objective, the best point
is highlighted. For multiple objectives, all Pareto-optimal points are highlighted.

All panels that display the same parameter are linked with Plotly's `matches` axes. Zooming or
panning one panel therefore updates the corresponding horizontal or vertical axes in the other
visible panels.
"""
function create_matrix_plot(results::Vector{Any},
                            io_settings::Dict{String,Any},
                            sim_params::Dict{String,Any};
                            objective_keys=nothing,
                            objective_senses=nothing,
                            color_key=nothing,
                            histogram_bins::Union{Nothing,Int}=nothing)
    param_names = String.(sim_params["parameter_study"]["runtime"]["parameter_keys"])

    if isempty(param_names)
        @error("Cannot create matrix plot: no parameter-study parameters were found.")
        return ""
    end

    spec = parameter_study_objective_spec(results,
                                          sim_params;
                                          objective_keys=objective_keys)
    results_dict = spec.results_dict
    isempty(spec.objective_keys) && return ""

    missing_params = filter(param -> !haskey(results_dict, param), param_names)
    if !isempty(missing_params)
        @error("Cannot create matrix plot. Missing parameters: " *
               join(missing_params, ", "))
        return ""
    end

    configured_color_key = color_key === nothing ? nothing : String(color_key)

    if configured_color_key === nothing ||
       !haskey(results_dict, configured_color_key) ||
       !any(is_finite_number, results_dict[configured_color_key])
        configured_color_key = first(spec.objective_keys)
    end

    senses = parameter_study_objective_senses(spec.objective_keys,
                                              sim_params;
                                              objective_senses=objective_senses)
    isempty(senses) && return ""

    # Use the complete scalar parameter-study result set for the colour selector.
    # This includes the configured objective, named objective components and
    # other finite scalar objective/KPI result columns. Variable parameter-study
    # parameters remain in their own dropdown group.
    result_color_keys = parameter_study_result_axis_keys(results,
                                                         results_dict,
                                                         sim_params,
                                                         spec.objective_keys)
    objective_color_keys = unique(vcat(result_color_keys,
                                       [configured_color_key]))
    selectable_color_keys = unique(vcat(param_names,
                                        objective_color_keys))

    required_keys = selectable_color_keys

    valid_idx = [idx
                 for idx in eachindex(results)
                 if all(key -> is_finite_number(results_dict[key][idx]),
                        required_keys)]

    if isempty(valid_idx)
        @error("Cannot create matrix plot: no complete numeric result rows.")
        return ""
    end

    parameter_values = Dict(
        key => Float64[Float64(results_dict[key][idx]) for idx in valid_idx]
        for key in param_names
    )

    objective_matrix = hcat([Float64[Float64(results_dict[key][idx]) for idx in valid_idx]
                             for key in spec.objective_keys]...)

    color_values_by_key = Dict(
        key => Float64[Float64(results_dict[key][idx]) for idx in valid_idx]
        for key in selectable_color_keys
    )

    color_values = color_values_by_key[configured_color_key]

    highlight_mask = if spec.is_multiobjective
        pareto_front_mask(objective_matrix,
                          [senses[key] for key in spec.objective_keys])
    else
        mask = falses(length(valid_idx))
        sense = get(senses, first(spec.objective_keys), :min)
        best_index = sense == :min ?
                     argmin(objective_matrix[:, 1]) :
                     argmax(objective_matrix[:, 1])
        mask[best_index] = true
        mask
    end

    color_sense = get(senses, configured_color_key, :min)

    # Preserve the matrix plot's established colour grading: approximately the
    # best quarter of the selected quantity receives the complete colour scale.
    sorted_color_values = sort(color_values)
    number_of_values = length(sorted_color_values)
    number_of_best_values = max(1, cld(number_of_values, 4))

    if color_sense == :max
        first_best_index = number_of_values - number_of_best_values + 1
        cmin = sorted_color_values[first_best_index]
        cmax = last(sorted_color_values)
    else
        cmin = first(sorted_color_values)
        cmax = sorted_color_values[number_of_best_values]
    end

    if !isfinite(cmin) || !isfinite(cmax) || cmax <= cmin
        cmin = first(sorted_color_values)
        cmax = last(sorted_color_values)
    end

    if cmax <= cmin
        delta = max(abs(cmin), 1.0) * 1.0e-9
        cmin -= delta
        cmax += delta
    end

    hover_keys = unique(vcat(param_names, spec.objective_keys))
    hover_text = ["run $(valid_idx[local_index])" *
                  join(("<br>$key = $(results_dict[key][valid_idx[local_index]])"
                        for key in hover_keys))
                  for local_index in eachindex(valid_idx)]

    number_of_bins = if histogram_bins === nothing
        clamp(ceil(Int, sqrt(length(valid_idx))), 6, 20)
    else
        histogram_bins > 0 ||
            throw(ArgumentError("histogram_bins must be greater than zero."))
        histogram_bins
    end

    function padded_parameter_range(values::Vector{Float64})::Tuple{Float64,Float64}
        minimum_value, maximum_value = extrema(values)
        scale = max(abs(minimum_value), abs(maximum_value), 1.0)

        padding = if maximum_value <= minimum_value
            scale * 0.05
        else
            (maximum_value - minimum_value) * 0.03
        end

        return minimum_value - padding, maximum_value + padding
    end

    function histogram_counts(values::Vector{Float64},
                              lower::Float64,
                              upper::Float64,
                              n_bins::Int)::Tuple{Vector{Float64},
                                                  Vector{Float64},
                                                  Vector{Int},
                                                  Vector{Float64},
                                                  Vector{Float64}}
        bin_width = (upper - lower) / n_bins
        bin_width > 0.0 ||
            throw(ArgumentError("Histogram range must have positive width."))

        lower_edges = [lower + (index - 1) * bin_width for index in 1:n_bins]
        upper_edges = [lower + index * bin_width for index in 1:n_bins]
        centers = (lower_edges .+ upper_edges) ./ 2.0
        widths = fill(bin_width, n_bins)
        counts = zeros(Int, n_bins)

        for value in values
            bin_index = if value >= upper
                n_bins
            else
                floor(Int, (value - lower) / bin_width) + 1
            end

            counts[clamp(bin_index, 1, n_bins)] += 1
        end

        return centers, widths, counts, lower_edges, upper_edges
    end

    function axis_reference(prefix::String, axis_index::Int)::String
        return axis_index == 1 ? prefix : prefix * string(axis_index)
    end

    function axis_layout_key(prefix::String, axis_index::Int)::Symbol
        return axis_index == 1 ?
               Symbol(prefix * "axis") :
               Symbol(prefix * "axis" * string(axis_index))
    end

    parameter_ranges = Dict(
        parameter => padded_parameter_range(parameter_values[parameter])
        for parameter in param_names
    )

    n_parameters = length(param_names)

    # Assign one Plotly x/y-axis pair to every visible cell in the lower triangle.
    cell_axis_index = Dict{Tuple{Int,Int},Int}()
    next_axis_index = 0

    for row in 1:n_parameters
        for column in 1:row
            next_axis_index += 1
            cell_axis_index[(row, column)] = next_axis_index
        end
    end

    # Every column shares one x parameter. The diagonal histogram is the master x axis.
    x_master_reference = Dict(
        column => axis_reference("x", cell_axis_index[(column, column)])
        for column in 1:n_parameters
    )

    # Every scatter row shares one y parameter. Its left-most scatter cell is the master.
    y_master_reference = Dict(
        row => axis_reference("y", cell_axis_index[(row, 1)])
        for row in 2:n_parameters
    )

    horizontal_spacing = n_parameters > 1 ?
                         min(0.025, 0.10 / (n_parameters - 1)) :
                         0.0

    vertical_spacing = n_parameters > 1 ?
                       min(0.030, 0.12 / (n_parameters - 1)) :
                       0.0

    # Reserve space on the right for a clearly labelled colour bar.
    plot_x_max = 0.88

    cell_width = (plot_x_max -
                  (n_parameters - 1) * horizontal_spacing) /
                 n_parameters

    cell_height = (1.0 -
                   (n_parameters - 1) * vertical_spacing) /
                  n_parameters

    traces = GenericTrace[]
    annotations = Any[]
    panel_shapes = Any[]
    layout_values = Dict{Symbol,Any}()
    highlight_legend_added = false
    scatter_trace_added = false

    for row in 1:n_parameters
        for column in 1:row
            axis_index = cell_axis_index[(row, column)]
            x_reference = axis_reference("x", axis_index)
            y_reference = axis_reference("y", axis_index)
            x_layout_key = axis_layout_key("x", axis_index)
            y_layout_key = axis_layout_key("y", axis_index)

            x_domain_start = (column - 1) * (cell_width + horizontal_spacing)
            x_domain_end = x_domain_start + cell_width

            y_domain_end = 1.0 -
                           (row - 1) * (cell_height + vertical_spacing)
            y_domain_start = y_domain_end - cell_height

            x_parameter = param_names[column]
            y_parameter = param_names[row]
            x_range = parameter_ranges[x_parameter]
            y_range = parameter_ranges[y_parameter]

            panel_fill = row == column ?
                         "rgba(245,247,250,0.95)" :
                         "rgba(230,236,246,0.72)"

            push!(panel_shapes,
                  attr(; type="rect",
                       xref="paper",
                       yref="paper",
                       x0=x_domain_start,
                       x1=x_domain_end,
                       y0=y_domain_start,
                       y1=y_domain_end,
                       fillcolor=panel_fill,
                       line=attr(; color="rgba(90,110,140,0.18)",
                                 width=1),
                       layer="below"))

            xaxis_values = Dict{Symbol,Any}(
                :domain => [x_domain_start, x_domain_end],
                :anchor => y_reference,
                :range => [x_range[1], x_range[2]],
                :showgrid => true,
                :gridcolor => "rgba(110,130,160,0.20)",
                :zeroline => false,
                :showline => false,
                :ticks => "outside",
                :tickfont => attr(; size=9),
                :showticklabels => row == n_parameters || row == column,
                :automargin => true,
                :fixedrange => false,
            )

            if x_reference != x_master_reference[column]
                xaxis_values[:matches] = x_master_reference[column]
            end

            if row == n_parameters
                xaxis_values[:title] = attr(; text=x_parameter,
                                            font=attr(; size=11))
            end

            layout_values[x_layout_key] = attr(; xaxis_values...)

            yaxis_values = Dict{Symbol,Any}(
                :domain => [y_domain_start, y_domain_end],
                :anchor => x_reference,
                :zeroline => false,
                :showline => false,
                :ticks => "outside",
                :tickfont => attr(; size=9),
                :automargin => true,
                :fixedrange => false,
            )

            if row == column
                # The diagonal y axis is a histogram count axis and must not be
                # linked to the parameter-valued y axes in the scatter row.
                yaxis_values[:showgrid] = false
                yaxis_values[:showticklabels] = column == 1

                if column == 1
                    yaxis_values[:title] = attr(; text="Count",
                                                font=attr(; size=11))
                end
            else
                yaxis_values[:range] = [y_range[1], y_range[2]]
                yaxis_values[:showgrid] = true
                yaxis_values[:gridcolor] = "rgba(110,130,160,0.20)"
                yaxis_values[:showticklabels] = column == 1

                if y_reference != y_master_reference[row]
                    yaxis_values[:matches] = y_master_reference[row]
                end

                if column == 1
                    yaxis_values[:title] = attr(; text=y_parameter,
                                                font=attr(; size=11))
                end
            end

            layout_values[y_layout_key] = attr(; yaxis_values...)

            if row == column
                centers, bar_widths, counts, lower_edges, upper_edges = histogram_counts(parameter_values[x_parameter],
                                                                                         x_range[1],
                                                                                         x_range[2],
                                                                                         number_of_bins)

                bin_customdata = hcat(lower_edges, upper_edges)

                push!(traces,
                      bar(; x=centers,
                          y=counts,
                          width=bar_widths,
                          marker=attr(; color="rgba(70,130,180,0.72)",
                                      line=attr(; color="rgba(55,90,120,0.85)",
                                                width=0.7)),
                          customdata=bin_customdata,
                          xaxis=x_reference,
                          yaxis=y_reference,
                          showlegend=false,
                          hovertemplate=("$x_parameter<br>" *
                                         "Range: %{customdata[0]:.6g} - " *
                                         "%{customdata[1]:.6g}<br>" *
                                         "Runs: %{y}<extra></extra>")))

                # Draw red histogram outlines only for bins that actually contain
                # highlighted runs. Omitting zero-height bins prevents the red
                # baseline that Plotly otherwise draws at y = 0.
                if any(highlight_mask)
                    highlighted_values = parameter_values[x_parameter][highlight_mask]
                    _, _, highlighted_counts, _, _ = histogram_counts(highlighted_values,
                                                                      x_range[1],
                                                                      x_range[2],
                                                                      number_of_bins)

                    nonzero_bins = findall(count -> count > 0, highlighted_counts)

                    if !isempty(nonzero_bins)
                        push!(traces,
                              bar(; x=centers[nonzero_bins],
                                  y=highlighted_counts[nonzero_bins],
                                  width=bar_widths[nonzero_bins],
                                  marker=attr(; color="rgba(255,255,255,0)",
                                              line=attr(; color="red",
                                                        width=2)),
                                  customdata=bin_customdata[nonzero_bins, :],
                                  xaxis=x_reference,
                                  yaxis=y_reference,
                                  name=spec.is_multiobjective ?
                                       "Pareto solutions" :
                                       "Best solution",
                                  showlegend=(!highlight_legend_added),
                                  hovertemplate=("$x_parameter - highlighted<br>" *
                                                 "Range: %{customdata[0]:.6g} - " *
                                                 "%{customdata[1]:.6g}<br>" *
                                                 "Highlighted runs: %{y}<extra></extra>")))

                        highlight_legend_added = true
                    end
                end

                push!(annotations,
                      attr(; text="<b>$x_parameter</b>",
                           x=(x_domain_start + x_domain_end) / 2,
                           y=y_domain_end - 0.012,
                           xref="paper",
                           yref="paper",
                           showarrow=false,
                           xanchor="center",
                           yanchor="top",
                           bgcolor="rgba(255,255,255,0.80)",
                           borderpad=2,
                           font=attr(; size=11,
                                     color="rgb(45,65,90)")))

                continue
            end

            push!(traces,
                  scatter(; x=parameter_values[x_parameter],
                          y=parameter_values[y_parameter],
                          mode="markers",
                          marker=attr(; size=5,
                                      opacity=0.72,
                                      color=color_values,
                                      coloraxis="coloraxis",
                                      line=attr(; width=0)),
                          text=hover_text,
                          hovertemplate="%{text}<extra></extra>",
                          xaxis=x_reference,
                          yaxis=y_reference,
                          meta="matrix-colour-runs",
                          showlegend=false))

            scatter_trace_added = true

            if any(highlight_mask)
                push!(traces,
                      scatter(; x=parameter_values[x_parameter][highlight_mask],
                              y=parameter_values[y_parameter][highlight_mask],
                              mode="markers",
                              marker=attr(; size=9,
                                          opacity=1.0,
                                          color=color_values[highlight_mask],
                                          coloraxis="coloraxis",
                                          line=attr(; color="red",
                                                    width=2)),
                              text=hover_text[highlight_mask],
                              hovertemplate="%{text}<extra></extra>",
                              xaxis=x_reference,
                              yaxis=y_reference,
                              meta="matrix-colour-highlight",
                              name=spec.is_multiobjective ?
                                   "Pareto solutions" :
                                   "Best solution",
                              showlegend=(!highlight_legend_added)))

                highlight_legend_added = true
            end
        end
    end

    # With one parameter-study parameter there is no scatter cell. Add a fully
    # transparent marker so Plotly still renders the shared colour bar.
    if !scatter_trace_added
        first_parameter = first(param_names)

        push!(traces,
              scatter(; x=[first(parameter_values[first_parameter])],
                      y=[0.0],
                      mode="markers",
                      marker=attr(; size=0,
                                  opacity=0.0,
                                  color=[first(color_values)],
                                  coloraxis="coloraxis"),
                      hoverinfo="skip",
                      xaxis="x",
                      yaxis="y",
                      meta="matrix-colour-dummy",
                      showlegend=false))
    end

    number_of_highlighted = count(identity, highlight_mask)
    highlight_description = if spec.is_multiobjective
        number_of_highlighted == 1 ?
        "1 Pareto solution" :
        "$number_of_highlighted Pareto solutions"
    else
        "best solution"
    end

    layout_values[:title] = attr(;
                                 text=("Parameter-study matrix - " *
                                       "$(length(valid_idx)) valid runs; " *
                                       "red outline: $highlight_description"),
                                 x=0.01,
                                 xanchor="left",
                                 font=attr(; size=17))

    layout_values[:annotations] = annotations
    layout_values[:shapes] = panel_shapes

    # Shared colour scale and explicit title/legend for all scatter panels.
    layout_values[:coloraxis] = attr(; colorscale="Viridis",
                                     reversescale=color_sense != :max,
                                     cmin=cmin,
                                     cmax=cmax,
                                     showscale=true,
                                     colorbar=attr(; title=attr(; text=configured_color_key,
                                                                side="right"),
                                                   x=0.915,
                                                   xanchor="left",
                                                   y=0.5,
                                                   len=0.88,
                                                   thickness=18,
                                                   ticks="",
                                                   ticklen=0,
                                                   outlinecolor="rgba(70,70,70,0.65)",
                                                   outlinewidth=1))

    # Let JavaScript size the figure to the current browser viewport.
    # No fixed width or height is stored in the Plotly layout, so browser
    # zoom and small screens can use all available screen space.
    layout_values[:autosize] = true
    layout_values[:margin] = attr(; t=75,
                                  b=70,
                                  l=90,
                                  r=150)
    layout_values[:paper_bgcolor] = "white"
    layout_values[:plot_bgcolor] = "rgba(0,0,0,0)"
    layout_values[:hovermode] = "closest"
    layout_values[:dragmode] = "zoom"
    layout_values[:showlegend] = any(highlight_mask)
    layout_values[:legend] = attr(; x=0.90,
                                  xanchor="left",
                                  y=1.0,
                                  yanchor="top",
                                  bgcolor="rgba(255,255,255,0.85)")
    layout_values[:barmode] = "overlay"
    layout_values[:bargap] = 0.06

    p = plot(traces, Layout(; layout_values...))
    file_path = parameter_study_plot_path(sim_params,
                                          io_settings,
                                          "matrix_plot")
    savefig(p, file_path)

    # Add the same toolbar style used by the other interactive parameter-study
    # plots. The graph occupies the viewport space remaining below the toolbar.
    if endswith(lowercase(file_path), ".html")
        html = read(file_path, String)

        json_for_html(value) = replace(JSON.json(value),
                                       "</" => "<\\/")

        color_data_json = json_for_html(color_values_by_key)
        parameter_keys_json = json_for_html(param_names)
        objective_keys_json = json_for_html(objective_color_keys)
        initial_color_json = json_for_html(configured_color_key)
        highlight_mask_json = json_for_html(collect(highlight_mask))
        color_senses_json = json_for_html(Dict(
                                              key => string(get(senses, key, :min))
                                              for key in selectable_color_keys
                                          ))

        controls_injection = render_parameter_study_plot_insertion("matrix_controls.html",
                                                                   "COLOR_DATA_JSON" => color_data_json,
                                                                   "PARAMETER_KEYS_JSON" => parameter_keys_json,
                                                                   "OBJECTIVE_KEYS_JSON" => objective_keys_json,
                                                                   "HIGHLIGHT_MASK_JSON" => highlight_mask_json,
                                                                   "COLOR_SENSES_JSON" => color_senses_json,
                                                                   "INITIAL_COLOR_JSON" => initial_color_json)

        if occursin("</head>", html)
            html = replace(html,
                           "</head>" => controls_injection * "\n</head>";
                           count=1)
        else
            html *= controls_injection
        end

        write(file_path, html)
    end

    return file_path
end

"""
    result_value(result, key)

Read a result entry using either a string or symbol key. Return `missing` when the key is absent.
"""
function result_value(result, key::String)
    if haskey(result, key)
        return result[key]
    end

    symbol_key = Symbol(key)
    if haskey(result, symbol_key)
        return result[symbol_key]
    end

    return missing
end

"""
    parameter_study_results_dict(results)

Collect every result column in a dictionary of vectors. Missing entries are represented by
`missing`. String and symbol result keys are handled consistently.
"""
function parameter_study_results_dict(results::Vector{Any})::Dict{String,Vector{Any}}
    if isempty(results)
        @error("Cannot create parameter-study plots: results are empty.")
        return Dict{String,Vector{Any}}()
    end

    all_keys = unique([String(key) for result in results for key in keys(result)])

    return Dict(key => [result_value(result, key) for result in results]
                for key in all_keys)
end

"""
    is_finite_number(value)

Return `true` for finite scalar real values.
"""
is_finite_number(value)::Bool = value isa Real && isfinite(Float64(value))

"""
    parameter_study_plot_path(sim_params, io_settings, suffix)

Create a plot path derived from `parameter_study_plots_path`.
"""
function parameter_study_plot_path(sim_params::Dict{String,Any},
                                   io_settings::Dict{String,Any},
                                   suffix::String)::String
    base_path = sim_params["run_path"](io_settings["parameter_study_plots_path"])
    dir, filename = splitdir(base_path)
    root, _ = splitext(filename)
    ext = ".html"
    return joinpath(dir, "$(root)_$(suffix)$(ext)")
end

"""
    parameter_study_color_bounds(values)

Return colour limits spanning the complete finite value range.
"""
function parameter_study_color_bounds(values::Vector{Float64})::Tuple{Float64,Float64}
    if isempty(values)
        @error("Cannot determine colour bounds from an empty vector.")
        return 0.0, 1.0
    end

    cmin = minimum(values)
    cmax = maximum(values)

    if cmax <= cmin
        delta = max(abs(cmin), 1.0) * 1.0e-9
        cmin -= delta
        cmax += delta
    end

    return cmin, cmax
end

"""
    objective_colorscale(sense)

Use bright colours for favourable values: low values for minimisation and high values for
maximisation.
"""
function objective_colorscale(sense::Symbol)
    sense == :max ? ColorSchemes.viridis : reverse(ColorSchemes.viridis)
end

"""
    parameter_study_objective_spec(results, sim_params; objective_keys=nothing)

Resolve the scalar objective columns used by the plots.

Single-objective results keep the scalar `"objective"` column, including objectives derived
from sum, mean, economic or emissions calculations.

For vector-valued multi-objective results, the objective names are taken from
`objective_keys` when supplied, otherwise from
`sim_params["parameter_study"]["runtime"]["objective_params_keys"]`. The named scalar columns already stored
in every parameter-study result are used directly; no numerical value matching or predefined
objective categories are required.
"""
function parameter_study_objective_spec(results::Vector{Any},
                                        sim_params::Dict{String,Any};
                                        objective_keys=nothing)
    results_dict = parameter_study_results_dict(results)
    isempty(results_dict) &&
        return (; results_dict,
                objective_keys=String[],
                is_multiobjective=false,
                vector_objective=false)

    objective_column = get(results_dict,
                           "objective",
                           Any[missing for _ in eachindex(results)])

    vector_lengths = Int[length(value)
                         for value in objective_column
                         if value isa AbstractVector || value isa Tuple]

    vector_objective = !isempty(vector_lengths)

    if !vector_objective
        if haskey(results_dict, "objective") &&
           any(is_finite_number, results_dict["objective"])
            return (; results_dict,
                    objective_keys=["objective"],
                    is_multiobjective=false,
                    vector_objective=false)
        end

        @error("Cannot create parameter-study plots: no finite scalar objective was found.")
        return (; results_dict,
                objective_keys=String[],
                is_multiobjective=false,
                vector_objective=false)
    end

    n_objectives = first(vector_lengths)

    if any(length_value -> length_value != n_objectives,
           vector_lengths)
        @error("Cannot create parameter-study plots: vector-valued objective entries have inconsistent lengths.")
        return (; results_dict,
                objective_keys=String[],
                is_multiobjective=n_objectives > 1,
                vector_objective=true)
    end

    selected_keys = if objective_keys !== nothing
        String.(objective_keys)
    elseif haskey(sim_params, "parameter_study") &&
           haskey(sim_params["parameter_study"]["runtime"], "objective_params_keys")
        String.(sim_params["parameter_study"]["runtime"]["objective_params_keys"])
    else
        String[]
    end

    if isempty(selected_keys)
        @error("Multi-objective plotting requires the ordered objective parameter keys in " *
               "sim_params[\"parameter_study\"][\"runtime\"][\"objective_params_keys\"] or via objective_keys.")
        return (; results_dict,
                objective_keys=String[],
                is_multiobjective=n_objectives > 1,
                vector_objective=true)
    end

    if length(selected_keys) != n_objectives
        @error("The number of objective parameter keys ($(length(selected_keys))) does not match " *
               "the vector objective length ($n_objectives).")
        return (; results_dict,
                objective_keys=String[],
                is_multiobjective=n_objectives > 1,
                vector_objective=true)
    end

    duplicate_keys = [key
                      for key in unique(selected_keys)
                      if count(==(key), selected_keys) > 1]

    if !isempty(duplicate_keys)
        @error("Objective parameter keys must be unique: " *
               join(duplicate_keys, ", "))
        return (; results_dict,
                objective_keys=String[],
                is_multiobjective=n_objectives > 1,
                vector_objective=true)
    end

    missing_keys = [key
                    for key in selected_keys
                    if !haskey(results_dict, key)]

    if !isempty(missing_keys)
        @error("Cannot create multi-objective plots. Objective result columns are missing: " *
               join(missing_keys, ", "))
        return (; results_dict,
                objective_keys=String[],
                is_multiobjective=n_objectives > 1,
                vector_objective=true)
    end

    invalid_keys = [key
                    for key in selected_keys
                    if !any(is_finite_number, results_dict[key])]

    if !isempty(invalid_keys)
        @error("Cannot create multi-objective plots. Objective result columns contain no finite " *
               "scalar values: " * join(invalid_keys, ", "))
        return (; results_dict,
                objective_keys=String[],
                is_multiobjective=n_objectives > 1,
                vector_objective=true)
    end

    return (; results_dict,
            objective_keys=selected_keys,
            is_multiobjective=n_objectives > 1,
            vector_objective=true)
end

"""
    parameter_study_result_axis_keys(results, results_dict, sim_params, objective_keys)

Return selectable scalar result axes. True objective columns are listed first, followed by other
finite scalar result or KPI columns. Parameters and metadata are excluded.
"""
function parameter_study_result_axis_keys(results::Vector{Any},
                                          results_dict::Dict{String,Vector{Any}},
                                          sim_params::Dict{String,Any},
                                          objective_keys::Vector{String})::Vector{String}
    param_names = String.(sim_params["parameter_study"]["runtime"]["parameter_keys"])
    excluded_keys = Set(vcat(param_names,
                             ["error",
                              "run",
                              "run_id",
                              "sample_id"]))

    ordered_result_keys = unique([String(key) for result in results for key in keys(result)])
    ordered_keys = unique(vcat(objective_keys,
                               ordered_result_keys,
                               collect(keys(results_dict))))

    return [key
            for key in ordered_keys
            if !(key in excluded_keys) &&
                   haskey(results_dict, key) &&
                   any(is_finite_number, results_dict[key])]
end

"""
    parameter_study_objective_senses(objective_keys, sim_params; objective_senses=nothing)

Resolve one `:min` or `:max` sense for each objective. The keyword may be a single symbol, a
vector in objective order or a dictionary keyed by objective name. If omitted, the function
first uses `sim_params["parameter_study"]["runtime"]["objective_senses"]` when available. Otherwise it infers
the direction from the sign of `objective_factors` (positive means minimise, negative means
maximise) and defaults to `:min` when no factor is available.
"""
function parameter_study_objective_senses(objective_keys::Vector{String},
                                          sim_params::Dict{String,Any};
                                          objective_senses=nothing)::Dict{String,Symbol}
    parameter_study = get(sim_params, "parameter_study", Dict{String,Any}())
    runtime = get(parameter_study, "runtime", Dict{String,Any}())
    configured = objective_senses

    # An explicit sense configuration has priority.
    if configured === nothing && haskey(runtime, "objective_senses")
        configured = runtime["objective_senses"]
    end

    resolved = Dict{String,Symbol}()

    if configured === nothing
        # When no explicit senses are supplied, infer the direction from
        # objective_factors if they are available. Optimisation backends minimise the
        # factor-weighted objective values:
        #
        #   positive factor -> minimise the original objective
        #   negative factor -> maximise the original objective
        #
        # The scalar aggregate "objective" remains a minimisation quantity.
        factors = get(parameter_study, "objective_factors", nothing)

        for key in objective_keys
            if key == "objective"
                resolved[key] = :min
                continue
            end

            factor = if factors isa AbstractDict && haskey(factors, key)
                factors[key]
            elseif factors isa AbstractDict && haskey(factors, Symbol(key))
                factors[Symbol(key)]
            else
                nothing
            end

            if factor isa Real && isfinite(Float64(factor))
                numeric_factor = Float64(factor)

                if numeric_factor > 0.0
                    resolved[key] = :min
                elseif numeric_factor < 0.0
                    resolved[key] = :max
                else
                    @warn("Objective factor for \"$key\" is zero. " *
                          "The objective has no effective direction; using :min for plotting.")
                    resolved[key] = :min
                end
            else
                resolved[key] = :min
            end
        end
    elseif configured isa Symbol
        for key in objective_keys
            resolved[key] = configured
        end
    elseif configured isa AbstractVector
        if length(configured) != length(objective_keys)
            @error("The number of objective senses ($(length(configured))) does not match the number of objectives ($(length(objective_keys))).")
            return Dict{String,Symbol}()
        end

        for (key, sense) in zip(objective_keys, configured)
            resolved[key] = sense isa Symbol ? sense : Symbol(sense)
        end
    elseif configured isa AbstractDict
        for key in objective_keys
            if haskey(configured, key)
                resolved[key] = configured[key] isa Symbol ? configured[key] : Symbol(configured[key])
            elseif haskey(configured, Symbol(key))
                resolved[key] = configured[Symbol(key)] isa Symbol ? configured[Symbol(key)] :
                                Symbol(configured[Symbol(key)])
            else
                resolved[key] = :min
            end
        end
    else
        @error("objective_senses must be nothing, :min, :max, a vector or a dictionary.")
        return Dict{String,Symbol}()
    end

    invalid = [key for key in objective_keys if !(resolved[key] in (:min, :max))]
    if !isempty(invalid)
        @error("Objective senses must be :min or :max. Invalid entries: " *
               join(["$key => $(resolved[key])" for key in invalid], ", "))
        return Dict{String,Symbol}()
    end

    return resolved
end

"""
    pareto_front_mask(objective_values, senses)

Return a mask identifying nondominated rows. Every objective may independently be minimised or
maximised.
"""
function pareto_front_mask(objective_values::Matrix{Float64},
                           senses::Vector{Symbol})::BitVector
    n_points, n_objectives = size(objective_values)

    if length(senses) != n_objectives
        throw(ArgumentError("The number of objective senses must match the objective columns."))
    end

    losses = copy(objective_values)

    for objective_index in 1:n_objectives
        if senses[objective_index] == :max
            losses[:, objective_index] .*= -1.0
        elseif senses[objective_index] != :min
            throw(ArgumentError("Objective sense must be :min or :max."))
        end
    end

    is_pareto = trues(n_points)

    for candidate_index in 1:n_points
        for other_index in 1:n_points
            candidate_index == other_index && continue

            no_worse = all(losses[other_index, :] .<= losses[candidate_index, :])
            strictly_better = any(losses[other_index, :] .< losses[candidate_index, :])

            if no_worse && strictly_better
                is_pareto[candidate_index] = false
                break
            end
        end
    end

    return is_pareto
end

"""
    create_objective_convergence_plot(results, io_settings, sim_params; ...)

Create one convergence plot.

For a single-objective run, the existing convergence plot is retained. For a multi-objective
run, one HTML file is created with an objective dropdown. The selected objective determines the
evaluation points, objective-specific best-so-far line, y-axis label and linear/log axis mode.
"""
function create_objective_convergence_plot(results::Vector{Any},
                                           io_settings::Dict{String,Any},
                                           sim_params::Dict{String,Any};
                                           objective_keys=nothing,
                                           objective_senses=nothing)
    spec = parameter_study_objective_spec(results,
                                          sim_params;
                                          objective_keys=objective_keys)
    results_dict = spec.results_dict
    isempty(spec.objective_keys) && return ""

    senses = parameter_study_objective_senses(spec.objective_keys,
                                              sim_params;
                                              objective_senses=objective_senses)
    isempty(senses) && return ""

    objective_data = Dict{String,Vector{Union{Nothing,Float64}}}()
    best_data = Dict{String,Vector{Union{Nothing,Float64}}}()
    use_log_axis = Dict{String,Bool}()
    available_objective_keys = String[]

    for objective_key in spec.objective_keys
        objective = Union{Nothing,Float64}[is_finite_number(value) ?
                                           Float64(value) :
                                           nothing
                                           for value in results_dict[objective_key]]

        valid_objective = Float64[value
                                  for value in objective
                                  if value !== nothing]

        if isempty(valid_objective)
            @warn "Skipping convergence objective: no finite values." objective_key
            continue
        end

        sense = senses[objective_key]
        best_so_far = Vector{Union{Nothing,Float64}}(undef,
                                                     length(objective))
        fill!(best_so_far, nothing)

        current_best = sense == :min ? Inf : -Inf

        for idx in eachindex(objective)
            value = objective[idx]

            if value !== nothing
                current_best = sense == :min ?
                               min(current_best, value) :
                               max(current_best, value)
            end

            if isfinite(current_best)
                best_so_far[idx] = current_best
            end
        end

        objective_data[objective_key] = objective
        best_data[objective_key] = best_so_far
        use_log_axis[objective_key] = all(value -> value > 0.0, valid_objective)

        push!(available_objective_keys,
              objective_key)
    end

    if isempty(available_objective_keys)
        @error("Cannot create convergence plot: no finite objective values.")
        return ""
    end

    initial_objective_key = first(available_objective_keys)

    initial_objective = objective_data[initial_objective_key]
    initial_best = best_data[initial_objective_key]

    objective_hover_text = ["run $idx<br>$initial_objective_key = " *
                            format_convergence_value(initial_objective[idx])
                            for idx in eachindex(initial_objective)]

    best_hover_text = ["run $idx<br>best $initial_objective_key = " *
                       format_convergence_value(initial_best[idx])
                       for idx in eachindex(initial_best)]

    objective_trace = scatter(; x=collect(eachindex(initial_objective)),
                              y=initial_objective,
                              mode="markers",
                              name=initial_objective_key,
                              text=objective_hover_text,
                              hovertemplate="%{text}<extra></extra>")

    best_trace = scatter(; x=collect(eachindex(initial_best)),
                         y=initial_best,
                         mode="lines",
                         name="Best so far",
                         text=best_hover_text,
                         hovertemplate="%{text}<extra></extra>")

    title = spec.is_multiobjective ?
            "Objective convergence: $initial_objective_key" :
            "Objective convergence"

    layout = Layout(; title=title,
                    xaxis=attr(; title="Run"),
                    yaxis=attr(; title=initial_objective_key,
                               type=use_log_axis[initial_objective_key] ?
                                    "log" :
                                    "linear"))

    plot_object = plot([objective_trace, best_trace], layout)

    # Use the same output suffix for single- and multi-objective runs.
    file_path = parameter_study_plot_path(sim_params,
                                          io_settings,
                                          "convergence")

    savefig(plot_object, file_path)

    if spec.is_multiobjective &&
       length(available_objective_keys) > 1
        inject_convergence_objective_controls!(file_path,
                                               objective_data,
                                               best_data,
                                               available_objective_keys,
                                               use_log_axis,
                                               collect(eachindex(results)),
                                               initial_objective_key)
    end

    return file_path
end

"""
    format_convergence_value(value)

Format a convergence value for initial Plotly hover text.
"""
function format_convergence_value(value::Union{Nothing,Float64})::String
    value === nothing && return "missing"
    return string(value)
end

"""
    inject_convergence_objective_controls!(file_path, objective_data, best_data,
                                           objective_keys, use_log_axis, run_ids,
                                           initial_objective_key)

Inject an objective dropdown into a multi-objective convergence HTML plot. Changing the selected
objective updates the evaluation points, objective-specific best-so-far line, plot title, hover
text, y-axis label and linear/log scale.
"""
function inject_convergence_objective_controls!(file_path::String,
                                                objective_data::Dict{String,Vector{Union{Nothing,Float64}}},
                                                best_data::Dict{String,Vector{Union{Nothing,Float64}}},
                                                objective_keys::Vector{String},
                                                use_log_axis::Dict{String,Bool},
                                                run_ids::Vector{Int},
                                                initial_objective_key::String)
    html = read(file_path, String)

    json_for_html(value) = replace(JSON.json(value),
                                   "</" => "<\\/")

    objective_data_json = json_for_html(objective_data)
    best_data_json = json_for_html(best_data)
    objective_keys_json = json_for_html(objective_keys)
    use_log_axis_json = json_for_html(use_log_axis)
    run_ids_json = json_for_html(run_ids)
    initial_objective_json = json_for_html(initial_objective_key)

    injection = render_parameter_study_plot_insertion("convergence_controls.html",
                                                      "OBJECTIVE_DATA_JSON" => objective_data_json,
                                                      "BEST_DATA_JSON" => best_data_json,
                                                      "OBJECTIVE_KEYS_JSON" => objective_keys_json,
                                                      "USE_LOG_AXIS_JSON" => use_log_axis_json,
                                                      "RUN_IDS_JSON" => run_ids_json,
                                                      "INITIAL_OBJECTIVE_JSON" => initial_objective_json)

    if occursin("</body>", html)
        html = replace(html,
                       "</body>" =>
                           injection *
                           "\n</body>")
    else
        html *= injection
    end

    write(file_path, html)

    return file_path
end

"""
    create_objective_parameter_plots(results, io_settings, sim_params; ...)

Create a single interactive 2D objective/parameter explorer. Both axes and the marker colour
can be selected independently from all finite scalar parameter-study parameters, objectives and
result/KPI columns.

Single-objective runs highlight the best solution. Multi-objective runs highlight the global
Pareto set.
"""
function create_objective_parameter_plots(results::Vector{Any},
                                          io_settings::Dict{String,Any},
                                          sim_params::Dict{String,Any};
                                          objective_keys=nothing,
                                          objective_senses=nothing,
                                          color_key=nothing)
    param_names = String.(sim_params["parameter_study"]["runtime"]["parameter_keys"])
    configured_objective_params = if haskey(sim_params["parameter_study"]["runtime"],
                                            "objective_params_keys")
        unique(String.(sim_params["parameter_study"]["runtime"]["objective_params_keys"]))
    else
        String[]
    end

    spec = parameter_study_objective_spec(results,
                                          sim_params;
                                          objective_keys=objective_keys)
    results_dict = spec.results_dict
    isempty(spec.objective_keys) && return ""

    senses = parameter_study_objective_senses(spec.objective_keys,
                                              sim_params;
                                              objective_senses=objective_senses)
    isempty(senses) && return ""

    configured_senses = if objective_senses !== nothing
        objective_senses
    elseif haskey(sim_params["parameter_study"]["runtime"], "objective_senses")
        sim_params["parameter_study"]["runtime"]["objective_senses"]
    else
        nothing
    end

    function explorer_objective_sense(key::String)::Symbol
        if haskey(senses, key)
            return senses[key]
        end

        raw_sense = if configured_senses isa Symbol
            configured_senses
        elseif configured_senses isa AbstractDict && haskey(configured_senses, key)
            configured_senses[key]
        elseif configured_senses isa AbstractDict && haskey(configured_senses, Symbol(key))
            configured_senses[Symbol(key)]
        else
            :min
        end

        sense = raw_sense isa Symbol ?
                raw_sense :
                Symbol(lowercase(String(raw_sense)))

        if !(sense in (:min, :max))
            @warn "Invalid objective sense for explorer colour; defaulting to :min." key raw_sense
            return :min
        end

        return sense
    end

    available_param_names = [key
                             for key in param_names
                             if haskey(results_dict, key) && any(is_finite_number, results_dict[key])]

    if isempty(available_param_names)
        @error("Cannot create objective/parameter explorer: no finite parameter-study parameters were found.")
        return ""
    end

    complete_objective_idx = [idx
                              for idx in eachindex(results)
                              if all(key -> is_finite_number(results_dict[key][idx]),
                                     spec.objective_keys)]

    highlight_global_indices = Set{Int}()

    if spec.is_multiobjective && !isempty(complete_objective_idx)
        objective_matrix = hcat([Float64[Float64(results_dict[key][idx])
                                         for idx in complete_objective_idx]
                                 for key in spec.objective_keys]...)

        pareto_mask = pareto_front_mask(objective_matrix,
                                        [senses[key]
                                         for key in spec.objective_keys])

        for local_index in findall(pareto_mask)
            push!(highlight_global_indices,
                  complete_objective_idx[local_index])
        end
    elseif !spec.is_multiobjective
        primary_objective = first(spec.objective_keys)
        valid_objective_idx = [idx
                               for idx in eachindex(results)
                               if is_finite_number(results_dict[primary_objective][idx])]

        if !isempty(valid_objective_idx)
            objective_values = Float64[Float64(results_dict[primary_objective][idx])
                                       for idx in valid_objective_idx]

            best_local_index = senses[primary_objective] == :min ?
                               argmin(objective_values) :
                               argmax(objective_values)

            push!(highlight_global_indices,
                  valid_objective_idx[best_local_index])
        end
    end

    available_result_keys = parameter_study_result_axis_keys(results,
                                                             results_dict,
                                                             sim_params,
                                                             spec.objective_keys)

    interactive_keys = unique(vcat(available_param_names,
                                   available_result_keys))

    plot_data = Dict{String,Vector{Union{Nothing,Float64}}}(
        key => Union{Nothing,Float64}[is_finite_number(value) ? Float64(value) : nothing
                                      for value in results_dict[key]]
        for key in interactive_keys
        if haskey(results_dict, key)
    )

    function candidate_color_keys()::Vector{String}
        return [key
                for key in interactive_keys
                if haskey(plot_data, key) && any(!isnothing, plot_data[key])]
    end

    initial_x = first(available_param_names)
    initial_y = first(spec.objective_keys)
    initial_color_candidates = candidate_color_keys()

    preferred_initial_color_candidates = if spec.is_multiobjective
        [key for key in spec.objective_keys if key != initial_y]
    elseif length(configured_objective_params) > 1
        [key for key in configured_objective_params
         if key != initial_y && key != "objective"]
    else
        String[]
    end

    filter!(key -> key in initial_color_candidates,
            preferred_initial_color_candidates)

    initial_color = if color_key !== nothing &&
                       String(color_key) in initial_color_candidates
        String(color_key)
    elseif isempty(preferred_initial_color_candidates)
        nothing
    else
        first(preferred_initial_color_candidates)
    end

    function valid_indices(x_key::String,
                           y_key::String,
                           selected_color::Union{Nothing,String})::Vector{Int}
        return [idx
                for idx in eachindex(results)
                if plot_data[x_key][idx] !== nothing &&
                       plot_data[y_key][idx] !== nothing &&
                       (selected_color === nothing || plot_data[selected_color][idx] !== nothing)]
    end

    initial_valid_idx = valid_indices(initial_x,
                                      initial_y,
                                      initial_color)

    if isempty(initial_valid_idx)
        @error("Cannot create objective/parameter explorer: no complete numeric rows were found for the initial selection.")
        return ""
    end

    x_values = Float64[plot_data[initial_x][idx] for idx in initial_valid_idx]
    y_values = Float64[plot_data[initial_y][idx] for idx in initial_valid_idx]

    function base_hover_text(x_key::String,
                             y_key::String,
                             valid_idx::Vector{Int})::Vector{String}
        return ["run $(valid_idx[local_index])" *
                "<br>$x_key = $(plot_data[x_key][valid_idx[local_index]])" *
                "<br>$y_key = $(plot_data[y_key][valid_idx[local_index]])"
                for local_index in eachindex(valid_idx)]
    end

    base_text = base_hover_text(initial_x,
                                initial_y,
                                initial_valid_idx)

    initial_hover_text = if initial_color === nothing
        base_text
    else
        [base_text[local_index] *
         "<br>$initial_color = $(plot_data[initial_color][initial_valid_idx[local_index]])"
         for local_index in eachindex(initial_valid_idx)]
    end

    marker_settings = if initial_color === nothing
        attr(; color="rgb(31,119,180)",
             size=7,
             opacity=0.78,
             showscale=false)
    else
        color_values = Float64[plot_data[initial_color][idx] for idx in initial_valid_idx]
        cmin, cmax = parameter_study_color_bounds(color_values)

        initial_color_sense = explorer_objective_sense(initial_color)

        attr(; color=color_values,
             colorscale=objective_colorscale(initial_color_sense),
             cmin=cmin,
             cmax=cmax,
             size=7,
             opacity=0.78,
             showscale=true,
             colorbar=attr(; title=initial_color,
                           thickness=16,
                           x=1.02,
                           xanchor="left",
                           y=0.50,
                           len=0.76))
    end

    runs_trace = scatter(; x=x_values,
                         y=y_values,
                         mode="markers",
                         name="Runs",
                         marker=marker_settings,
                         text=initial_hover_text,
                         hovertemplate="%{text}<extra></extra>")

    highlight_local_idx = [local_index
                           for (local_index, global_index) in pairs(initial_valid_idx)
                           if global_index in highlight_global_indices]

    highlight_name = spec.is_multiobjective ?
                     "Pareto solutions" :
                     "Best solution"

    highlight_trace = scatter(; x=x_values[highlight_local_idx],
                              y=y_values[highlight_local_idx],
                              mode="markers",
                              name=highlight_name,
                              marker=attr(; symbol=spec.is_multiobjective ? "diamond" : "x",
                                          size=spec.is_multiobjective ? 10 : 14,
                                          color="black",
                                          line=attr(; color="white", width=1.2)),
                              text=[base_text[idx] for idx in highlight_local_idx],
                              hovertemplate="%{text}<extra></extra>")

    layout = Layout(; title=attr(; text="Interactive Objective/Parameter Explorer",
                                 x=0.5,
                                 xanchor="center"),
                    xaxis=attr(; title=initial_x,
                               autorange=true),
                    yaxis=attr(; title=initial_y,
                               type=all(value -> value > 0.0, y_values) ? "log" : "linear",
                               autorange=true),
                    legend=attr(; x=1.18,
                                xanchor="left",
                                y=1.0,
                                yanchor="top"),
                    margin=attr(; t=70,
                                b=70,
                                l=85,
                                r=300),
                    autosize=true)

    p = plot(GenericTrace[runs_trace, highlight_trace], layout)
    file_path = parameter_study_plot_path(sim_params,
                                          io_settings,
                                          "objective_parameter_explorer")
    if lowercase(splitext(file_path)[2]) != ".html"
        @error("The objective/parameter explorer requires an HTML output path.")
        return ""
    end

    mkpath(dirname(file_path))
    savefig(p, file_path)

    explorer_sense_keys = unique(vcat(spec.objective_keys,
                                      configured_objective_params))
    sense_strings = Dict(key => String(explorer_objective_sense(key))
                         for key in explorer_sense_keys)

    inject_objective_parameter_controls!(file_path,
                                         plot_data,
                                         interactive_keys,
                                         available_param_names,
                                         collect(eachindex(results)),
                                         collect(highlight_global_indices),
                                         sense_strings,
                                         spec.is_multiobjective,
                                         initial_x,
                                         initial_y,
                                         initial_color)

    return file_path
end

"""
    inject_objective_parameter_controls!(file_path, plot_data, variable_keys,
                                         parameter_keys, run_ids, highlight_run_ids,
                                         objective_senses, is_multiobjective, initial_x,
                                         initial_y, initial_color)

Inject grouped x-axis, y-axis and colour selectors into the objective/parameter explorer.
Every finite scalar parameter-study parameter, objective and result/KPI column remains available.
Variables are shown under the two dropdown groups "Parameter-study parameters" and "Objectives".
"""
function inject_objective_parameter_controls!(file_path::String,
                                              plot_data::Dict{String,Vector{Union{Nothing,Float64}}},
                                              variable_keys::Vector{String},
                                              parameter_keys::Vector{String},
                                              run_ids::Vector{Int},
                                              highlight_run_ids::Vector{Int},
                                              objective_senses::Dict{String,String},
                                              is_multiobjective::Bool,
                                              initial_x::String,
                                              initial_y::String,
                                              initial_color::Union{Nothing,String})
    html = read(file_path, String)

    json_for_html(value) = replace(JSON.json(value),
                                   "</" => "<\\/")

    plot_data_json = json_for_html(plot_data)
    variable_keys_json = json_for_html(variable_keys)
    parameter_keys_json = json_for_html(parameter_keys)
    run_ids_json = json_for_html(run_ids)
    highlight_run_ids_json = json_for_html(highlight_run_ids)
    objective_senses_json = json_for_html(objective_senses)
    initial_x_json = json_for_html(initial_x)
    initial_y_json = json_for_html(initial_y)
    initial_color_json = json_for_html(initial_color)
    is_multiobjective_json = json_for_html(is_multiobjective)

    injection = render_parameter_study_plot_insertion("objective_parameter_controls.html",
                                                      "PLOT_DATA_JSON" => plot_data_json,
                                                      "VARIABLE_KEYS_JSON" => variable_keys_json,
                                                      "PARAMETER_KEYS_JSON" => parameter_keys_json,
                                                      "RUN_IDS_JSON" => run_ids_json,
                                                      "HIGHLIGHT_RUN_IDS_JSON" => highlight_run_ids_json,
                                                      "OBJECTIVE_SENSES_JSON" => objective_senses_json,
                                                      "IS_MULTIOBJECTIVE_JSON" => is_multiobjective_json,
                                                      "INITIAL_X_JSON" => initial_x_json,
                                                      "INITIAL_Y_JSON" => initial_y_json,
                                                      "INITIAL_COLOR_JSON" => initial_color_json)

    if occursin("</body>", html)
        html = replace(html,
                       "</body>" => injection * "\n</body>")
    else
        html *= injection
    end

    write(file_path, html)

    return file_path
end

"""
    create_parallel_coordinates_plot(results, io_settings, sim_params; ...)

Create a filterable parallel-coordinates plot. Parameters and all finite scalar result axes are
shown. Multi-objective vector entries are represented by the configured named scalar objective columns.
Line colour defaults to the first objective and may be changed with `color_key`.
"""
function create_parallel_coordinates_plot(results::Vector{Any},
                                          io_settings::Dict{String,Any},
                                          sim_params::Dict{String,Any};
                                          objective_keys=nothing,
                                          objective_senses=nothing,
                                          color_key=nothing)
    param_names = String.(sim_params["parameter_study"]["runtime"]["parameter_keys"])
    spec = parameter_study_objective_spec(results,
                                          sim_params;
                                          objective_keys=objective_keys)
    results_dict = spec.results_dict
    isempty(spec.objective_keys) && return ""

    senses = parameter_study_objective_senses(spec.objective_keys,
                                              sim_params;
                                              objective_senses=objective_senses)
    isempty(senses) && return ""

    result_axis_keys = parameter_study_result_axis_keys(results,
                                                        results_dict,
                                                        sim_params,
                                                        spec.objective_keys)

    color_source_keys = unique(vcat(param_names, result_axis_keys))

    selected_color_key = color_key === nothing ?
                         first(spec.objective_keys) :
                         String(color_key)

    if !(selected_color_key in color_source_keys)
        @warn "Requested parallel-coordinates color key is unavailable; using the first objective." selected_color_key
        selected_color_key = first(spec.objective_keys)
    end

    required_keys = unique(vcat(param_names, result_axis_keys))
    missing_keys = filter(key -> !haskey(results_dict, key), required_keys)

    if !isempty(missing_keys)
        @error("Cannot create parallel-coordinates plot. Missing result keys: " *
               join(missing_keys, ", "))
        return ""
    end

    valid_idx = [idx
                 for idx in eachindex(results)
                 if all(key -> is_finite_number(results_dict[key][idx]), required_keys)]

    if isempty(valid_idx)
        @error("Cannot create parallel-coordinates plot: no complete numeric result rows.")
        return ""
    end

    dimensions = Any[]

    for param_name in param_names
        values = Float64[Float64(results_dict[param_name][idx]) for idx in valid_idx]
        push!(dimensions, attr(; label=param_name, values=values))
    end

    for key in result_axis_keys
        values = Float64[Float64(results_dict[key][idx]) for idx in valid_idx]
        push!(dimensions, attr(; label=key, values=values))
    end

    color_values = Float64[Float64(results_dict[selected_color_key][idx])
                           for idx in valid_idx]

    color_sense = get(senses, selected_color_key, :min)
    cmin, cmax = parameter_study_color_bounds(color_values)
    n_parameters = length(param_names)
    n_dimensions = length(dimensions)
    n_results = n_dimensions - n_parameters

    plot_x_min = 0.035
    plot_x_max = 0.955
    plot_y_max = 0.92
    group_boundary = if n_dimensions > 1
        plot_x_min +
        ((n_parameters - 0.5) / (n_dimensions - 1)) *
        (plot_x_max - plot_x_min)
    else
        plot_x_max
    end

    trace = parcoords(; ids=string.(valid_idx),
                      domain=attr(; x=[plot_x_min, plot_x_max],
                                  y=[0.0, plot_y_max]),
                      labelfont=attr(; size=11),
                      line=attr(; color=color_values,
                                colorscale=objective_colorscale(color_sense),
                                reversescale=false,
                                cmin=cmin,
                                cmax=cmax,
                                showscale=true,
                                colorbar=attr(; title=selected_color_key,
                                              len=plot_y_max,
                                              y=plot_y_max / 2,
                                              x=0.985,
                                              xanchor="left",
                                              xpad=2,
                                              thickness=13)),
                      unselected=attr(; line=attr(; color="rgb(145,145,145)",
                                                  opacity=0.10)),
                      dimensions=dimensions)

    group_shapes, group_annotations = parallel_group_decorations(n_parameters,
                                                                 n_results,
                                                                 plot_x_min,
                                                                 group_boundary,
                                                                 plot_x_max,
                                                                 plot_y_max)

    layout = Layout(;
                    title=attr(; text="Interactive Parameter Study Design Space",
                               x=0.5,
                               xanchor="center",
                               y=0.995,
                               yanchor="top",
                               font=attr(; size=18)),
                    margin=attr(; t=26, b=24, l=28, r=58),
                    shapes=group_shapes,
                    annotations=group_annotations)

    p = plot(trace, layout)
    file_path = parameter_study_plot_path(sim_params,
                                          io_settings,
                                          "parallel_coordinates")
    savefig(p, file_path)
    inject_parallel_axis_zoom_controls!(file_path,
                                        param_names,
                                        spec.objective_keys,
                                        selected_color_key,
                                        senses)

    return file_path
end

"""
    parallel_group_decorations(n_parameters, n_results, plot_x_min, group_boundary,
                               plot_x_max, plot_y_max)

Create background regions, group labels and the separator for a parallel-coordinates plot.
"""
function parallel_group_decorations(n_parameters::Int,
                                    n_results::Int,
                                    plot_x_min::Float64,
                                    group_boundary::Float64,
                                    plot_x_max::Float64,
                                    plot_y_max::Float64)
    shapes = Any[]
    annotations = Any[]
    groups = [(enabled=n_parameters > 0,
               label="Variable parameters",
               x0=plot_x_min,
               x1=group_boundary,
               fill="rgba(70,130,180,0.07)",
               text_color="rgb(55,90,120)"),
              (enabled=n_results > 0,
               label="Results",
               x0=group_boundary,
               x1=plot_x_max,
               fill="rgba(220,140,50,0.07)",
               text_color="rgb(140,85,30)")]

    for group in groups
        group.enabled || continue

        push!(shapes,
              attr(; type="rect",
                   xref="paper",
                   yref="paper",
                   x0=group.x0,
                   x1=group.x1,
                   y0=0.0,
                   y1=1.0,
                   fillcolor=group.fill,
                   line=attr(; width=0),
                   layer="below"))

        push!(annotations,
              attr(; text="<b>$(group.label)</b>",
                   x=(group.x0 + group.x1) / 2,
                   y=0.99,
                   xref="paper",
                   yref="paper",
                   showarrow=false,
                   xanchor="center",
                   yanchor="top",
                   align="center",
                   font=attr(; size=11, color=group.text_color)))
    end

    if n_parameters > 0 && n_results > 0
        push!(shapes,
              attr(; type="line",
                   xref="paper",
                   yref="paper",
                   x0=group_boundary,
                   x1=group_boundary,
                   y0=0.0,
                   y1=1.0,
                   line=attr(; color="rgba(80,80,80,0.65)",
                             width=2,
                             dash="dot"),
                   layer="above"))
    end

    return shapes, annotations
end

"""
    parameter_study_objective_axis_keys(results, sim_params, primary_obj_key="objective";
                                     objective_keys=nothing)

Backward-compatible helper returning all scalar objective or KPI axes. The vector-valued
`"objective"` column itself is excluded; the configured named scalar objective columns are included.
"""
function parameter_study_objective_axis_keys(results::Vector{Any},
                                             sim_params::Dict{String,Any},
                                             primary_obj_key::String="objective";
                                             objective_keys=nothing)::Vector{String}
    spec = parameter_study_objective_spec(results,
                                          sim_params;
                                          objective_keys=objective_keys)

    isempty(spec.objective_keys) && return String[]

    return parameter_study_result_axis_keys(results,
                                            spec.results_dict,
                                            sim_params,
                                            spec.objective_keys)
end

function inject_parallel_axis_zoom_controls!(file_path::String,
                                             parameter_keys::Vector{String},
                                             objective_keys::Vector{String},
                                             initial_color_key::String,
                                             objective_senses::Dict{String,Symbol})
    html = read(file_path, String)
    parameter_keys_json = replace(JSON.json(parameter_keys),
                                  "</" => "<\\/")
    objective_keys_json = replace(JSON.json(objective_keys),
                                  "</" => "<\\/")
    initial_color_key_json = replace(JSON.json(initial_color_key),
                                     "</" => "<\\/")
    objective_senses_json = replace(JSON.json(Dict(key => String(value)
                                                   for (key, value) in objective_senses)),
                                    "</" => "<\\/")

    injection = render_parameter_study_plot_insertion("parallel_coordinates_controls.html",
                                                      "PARAMETER_KEYS_JSON" => parameter_keys_json,
                                                      "OBJECTIVE_KEYS_JSON" => objective_keys_json,
                                                      "OBJECTIVE_SENSES_JSON" => objective_senses_json,
                                                      "INITIAL_COLOR_KEY_JSON" => initial_color_key_json)

    if occursin("</body>", html)
        html = replace(html,
                       "</body>" => injection * "\n</body>")
    else
        html *= injection
    end

    write(file_path, html)

    return file_path
end

"""
    create_interactive_3d_parameter_study_plot(results, io_settings, sim_params;
                                            objective_keys=nothing, x_key=nothing,
                                            y_key=nothing, z_key=nothing, color_key=nothing,
                                            objective_sense=:min)

Create an interactive 3D parameter-study explorer with selectable axes, colour variable and
independent axis zoom controls.
"""

"""
    create_interactive_3d_parameter_study_plot(results, io_settings, sim_params; ...)

Create an interactive 3D parameter-study explorer with selectable parameter and result axes,
selectable colour, two fixed control rows and independent axis zoom controls. Single-objective
runs highlight the best solution. Multi-objective runs highlight all Pareto-optimal solutions.
"""
function create_interactive_3d_parameter_study_plot(results::Vector{Any},
                                                    io_settings::Dict{String,Any},
                                                    sim_params::Dict{String,Any};
                                                    objective_keys=nothing,
                                                    objective_senses=nothing,
                                                    x_key=nothing,
                                                    y_key=nothing,
                                                    z_key=nothing,
                                                    color_key=nothing,
                                                    objective_sense::Symbol=:min)
    if !(objective_sense in (:min, :max))
        @error("objective_sense must be :min or :max.")
        return ""
    end

    param_names = String.(sim_params["parameter_study"]["runtime"]["parameter_keys"])
    spec = parameter_study_objective_spec(results,
                                          sim_params;
                                          objective_keys=objective_keys)
    results_dict = spec.results_dict
    isempty(spec.objective_keys) && return ""

    has_configured_senses = haskey(sim_params, "parameter_study") &&
                            haskey(sim_params["parameter_study"]["runtime"], "objective_senses")

    senses_input = objective_senses === nothing && !has_configured_senses ?
                   objective_sense :
                   objective_senses

    senses = parameter_study_objective_senses(spec.objective_keys,
                                              sim_params;
                                              objective_senses=senses_input)
    isempty(senses) && return ""

    result_axis_keys = parameter_study_result_axis_keys(results,
                                                        results_dict,
                                                        sim_params,
                                                        spec.objective_keys)
    axis_keys = unique(vcat(param_names, result_axis_keys))

    if isempty(axis_keys)
        @error("The 3D plot requires at least one selectable numeric quantity.")
        return ""
    end

    required_keys = axis_keys
    missing_keys = filter(key -> !haskey(results_dict, key), required_keys)

    if !isempty(missing_keys)
        @error("Cannot create 3D plot. Missing result keys: $(join(missing_keys, ", ")).")
        return ""
    end

    valid_idx = [idx
                 for idx in eachindex(results)
                 if all(key -> is_finite_number(results_dict[key][idx]), required_keys)]

    if isempty(valid_idx)
        @error("Cannot create 3D plot: no complete numeric result rows.")
        return ""
    end

    axis_data = Dict(
        key => Float64[Float64(results_dict[key][idx]) for idx in valid_idx]
        for key in axis_keys
    )

    function select_axis_key(explicit_key,
                             preferred_keys::Vector{String},
                             used_keys::Vector{String})::String
        if explicit_key !== nothing
            selected = String(explicit_key)

            if !(selected in axis_keys)
                @error("Axis key \"$selected\" is not available.")
                return ""
            end

            return selected
        end

        candidates = unique(vcat(preferred_keys, axis_keys))

        for candidate in candidates
            if candidate in axis_keys && !(candidate in used_keys)
                return candidate
            end
        end

        for candidate in candidates
            if candidate in axis_keys
                return candidate
            end
        end

        @error("Cannot find an available axis key.")
        return ""
    end

    default_objective_key = first(spec.objective_keys)
    x_selected = select_axis_key(x_key, param_names, String[])
    y_preferred = length(param_names) >= 2 ? param_names[2:end] : axis_keys
    y_selected = select_axis_key(y_key, y_preferred, [x_selected])
    z_selected = select_axis_key(z_key,
                                 [default_objective_key],
                                 [x_selected, y_selected])

    if any(isempty, [x_selected, y_selected, z_selected])
        return ""
    end

    color_selected = color_key === nothing ?
                     default_objective_key :
                     String(color_key)

    if !(color_selected in axis_keys)
        @error("Colour key \"$color_selected\" is not available.")
        return ""
    end

    objective_matrix = hcat([Float64[Float64(results_dict[key][idx]) for idx in valid_idx]
                             for key in spec.objective_keys]...)

    highlight_indices = if spec.is_multiobjective
        findall(pareto_front_mask(objective_matrix,
                                  [senses[key]
                                   for key in spec.objective_keys]))
    else
        sense = senses[first(spec.objective_keys)]
        [sense == :min ?
         argmin(objective_matrix[:, 1]) :
         argmax(objective_matrix[:, 1])]
    end

    highlight_name = spec.is_multiobjective ?
                     "Pareto solutions" :
                     "Best objective"

    color_values = axis_data[color_selected]
    cmin, cmax = parameter_study_color_bounds(color_values)
    color_sense = get(senses, color_selected, :min)

    hover_keys = unique(vcat(param_names,
                             [x_selected,
                              y_selected,
                              z_selected,
                              color_selected]))

    hover_text = ["run $(valid_idx[idx])" *
                  join(("<br>$key = $(axis_data[key][idx])" for key in hover_keys))
                  for idx in eachindex(valid_idx)]

    runs_trace = PlotlyJS.scatter(; type="scatter3d",
                                  x=axis_data[x_selected],
                                  y=axis_data[y_selected],
                                  z=axis_data[z_selected],
                                  mode="markers",
                                  name="Parameter-study runs",
                                  marker=attr(; size=5,
                                              opacity=0.75,
                                              color=color_values,
                                              colorscale=objective_colorscale(color_sense),
                                              cmin=cmin,
                                              cmax=cmax,
                                              showscale=true,
                                              colorbar=attr(; title=color_selected,
                                                            thickness=16,
                                                            x=1.02,
                                                            xanchor="left",
                                                            y=0.50,
                                                            len=0.72),
                                              line=attr(; width=0.3,
                                                        color="rgba(50,50,50,0.35)")),
                                  text=hover_text,
                                  hovertemplate="%{text}<extra></extra>")

    highlight_trace = PlotlyJS.scatter(; type="scatter3d",
                                       x=axis_data[x_selected][highlight_indices],
                                       y=axis_data[y_selected][highlight_indices],
                                       z=axis_data[z_selected][highlight_indices],
                                       mode="markers",
                                       name=highlight_name,
                                       marker=attr(; size=9,
                                                   symbol="diamond",
                                                   color="black",
                                                   line=attr(; color="white", width=1.5)),
                                       text=hover_text[highlight_indices],
                                       hovertemplate="%{text}<extra></extra>")

    layout = Layout(; title=attr(; text="Interactive 3D Parameter Study Explorer",
                                 x=0.5,
                                 xanchor="center"),
                    scene=attr(; xaxis=attr(; title=x_selected, autorange=true),
                               yaxis=attr(; title=y_selected, autorange=true),
                               zaxis=attr(; title=z_selected, autorange=true),
                               aspectmode="cube",
                               camera=attr(; eye=attr(; x=1.35,
                                                      y=1.35,
                                                      z=1.10))),
                    legend=attr(; x=1.16,
                                xanchor="left",
                                y=1.0,
                                yanchor="top"),
                    margin=attr(; t=60, b=20, l=20, r=250),
                    height=720)

    p = plot([runs_trace, highlight_trace], layout)
    file_path = parameter_study_plot_path(sim_params,
                                          io_settings,
                                          "interactive_3d")

    if lowercase(splitext(file_path)[2]) != ".html"
        @error("The interactive 3D plot requires an HTML output path.")
        return ""
    end

    mkpath(dirname(file_path))
    savefig(p, file_path)

    inject_3d_axis_selection_controls!(file_path,
                                       axis_data,
                                       axis_keys,
                                       param_names,
                                       valid_idx,
                                       highlight_indices,
                                       highlight_name,
                                       x_selected,
                                       y_selected,
                                       z_selected,
                                       color_selected)

    return file_path
end

"""
    inject_3d_axis_selection_controls!(file_path, axis_data, axis_keys,
                                       parameter_keys, run_ids, highlight_indices,
                                       highlight_name, initial_x, initial_y,
                                       initial_z, initial_color)

Inject axis selection, colour selection, camera reset and axis-zoom controls into a saved
3D parameter-study HTML plot.
"""
function inject_3d_axis_selection_controls!(file_path::String,
                                            axis_data::Dict{String,Vector{Float64}},
                                            axis_keys::Vector{String},
                                            parameter_keys::Vector{String},
                                            run_ids::Vector{Int},
                                            highlight_indices::Vector{Int},
                                            highlight_name::String,
                                            initial_x::String,
                                            initial_y::String,
                                            initial_z::String,
                                            initial_color::String)
    html = read(file_path, String)

    json_for_html(value) = replace(JSON.json(value),
                                   "</" => "<\\/")

    data_json = json_for_html(axis_data)
    keys_json = json_for_html(axis_keys)
    parameter_keys_json = json_for_html(parameter_keys)
    runs_json = json_for_html(run_ids)

    x_json = json_for_html(initial_x)
    y_json = json_for_html(initial_y)
    z_json = json_for_html(initial_z)
    color_json = json_for_html(initial_color)
    highlight_indices_json = json_for_html(highlight_indices .- 1)
    highlight_name_json = json_for_html(highlight_name)

    injection = render_parameter_study_plot_insertion("parameter_3d_controls.html",
                                                      "DATA_JSON" => data_json,
                                                      "KEYS_JSON" => keys_json,
                                                      "PARAMETER_KEYS_JSON" => parameter_keys_json,
                                                      "RUNS_JSON" => runs_json,
                                                      "HIGHLIGHT_INDICES_JSON" => highlight_indices_json,
                                                      "HIGHLIGHT_NAME_JSON" => highlight_name_json,
                                                      "X_JSON" => x_json,
                                                      "Y_JSON" => y_json,
                                                      "Z_JSON" => z_json,
                                                      "COLOR_JSON" => color_json)

    if occursin("</body>", html)
        html = replace(html,
                       "</body>" => injection * "\n</body>")
    else
        html *= injection
    end

    write(file_path, html)

    return file_path
end

"""
    create_parameter_study_diagnostic_plots(results, io_settings, sim_params; ...)

Create all parameter-study diagnostic plots and return their output paths. Objective columns are
resolved from the scalar `"objective"` column for single-objective runs and from
`sim_params["parameter_study"]["runtime"]["objective_params_keys"]` for multi-objective runs. Optional
objective names, senses and colour selection are forwarded consistently to all figures.
"""
function create_parameter_study_diagnostic_plots(results::Vector{Any},
                                                 io_settings::Dict{String,Any},
                                                 sim_params::Dict{String,Any};
                                                 objective_keys=nothing,
                                                 objective_senses=nothing,
                                                 color_key=nothing)
    if !io_settings["output_parameter_study_plots"]
        @globalInfo("Generation of parameter-study plots is deactivated. Set `output_parameter_study_plots` to true.")
        return
    end

    @globalInfo("Preparing parameter-study figures...")

    matrix = create_matrix_plot(results,
                                io_settings,
                                sim_params;
                                objective_keys=objective_keys,
                                objective_senses=objective_senses,
                                color_key=color_key)
    @globalInfo "Parameter-study matrix plot created and saved to $matrix"

    convergence = create_objective_convergence_plot(results,
                                                    io_settings,
                                                    sim_params;
                                                    objective_keys=objective_keys,
                                                    objective_senses=objective_senses)
    @globalInfo "Parameter-study convergence plot created and saved to $convergence"

    parameter_plots = create_objective_parameter_plots(results,
                                                       io_settings,
                                                       sim_params;
                                                       objective_keys=objective_keys,
                                                       objective_senses=objective_senses,
                                                       color_key=color_key)
    @globalInfo "Parameter-study parameter plot created and saved to $parameter_plots"

    parallel_coordinates = create_parallel_coordinates_plot(results,
                                                            io_settings,
                                                            sim_params;
                                                            objective_keys=objective_keys,
                                                            objective_senses=objective_senses,
                                                            color_key=color_key)
    @globalInfo "Parameter-study parallel plot created and saved to $parallel_coordinates"

    interactive_3d = create_interactive_3d_parameter_study_plot(results,
                                                                io_settings,
                                                                sim_params;
                                                                objective_keys=objective_keys,
                                                                objective_senses=objective_senses,
                                                                color_key=color_key)

    @globalInfo "Parameter-study 3D plot created and saved to $interactive_3d"
end
