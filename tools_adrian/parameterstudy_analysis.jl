using CSV
using DataFrames
using Plots
using Tables
using XLSX
using Dates
using JSON
import PlotlyJS

const CSV_PATH = "C:/Users/jenter/Documents/resie/output/parameterstudy/results_125runs_260308_154521.csv"
const OUTDIR = "c:/Users/jenter/Documents/resie/output/parameterstudy/plots/"
const OUT_CSV_PATH = "C:/Users/jenter/Documents/resie/output/out.csv"
const INPUT_JSON_PATH = "C:/Users/jenter/Documents/resie/inputfiles/inputfile_base_fuzzy_ems.json"

const XCOLS = [
    "HeatPump_Power / W",
    "ElectrodeBoiler_Power / W",
    "BufferTank_Capacity / Wh",
    "Battery_Capacity / Wh",
]
const YCOL_NO = "annuity_no A_total / €"
const CO2_COL = "CO2_yearly / t/a"
const BAL_COLS = ["balance_power", "balance_heat"]
const BAL_THRESHOLD = 1.0
const ANNUITY_COST_COLS = [
    "annuity_no A_cap / €",
    "annuity_no A_misc / €",
    "annuity_no A_op / €",
    "annuity_no A_energy / €",
]
const ANNUITY_REV_COLS = [
    "annuity_no A_rev_control / €",
    "annuity_no A_rev_feed / €",
]

to_million(x::Real) = x / 1e6

to_float(x) = x === missing ? missing : (
    x isa Real ? Float64(x) : (
        x isa AbstractString ? something(tryparse(Float64, replace(strip(x), "," => ".")), missing) : missing
    )
)

function coerce_columns!(df::DataFrame, cols::Vector{String})
    for c in cols
        c in names(df) && (df[!, c] = to_float.(df[!, c]))
    end
    return df
end

function filter_balance(df::DataFrame; bal_cols::Vector{String}=BAL_COLS, threshold::Float64=BAL_THRESHOLD)
    dfv = df
    for bcol in bal_cols
        @assert bcol in names(dfv) "Spalte '$bcol' fehlt in der CSV."
        dfv = dfv[.!ismissing.(dfv[!, bcol]) .&& (abs.(dfv[!, bcol]) .<= threshold), :]
    end
    return dfv
end

function topn_annuity(df::DataFrame, n::Int; xcols::Vector{String}=XCOLS, ycol::String=YCOL_NO)
    dfv = filter_balance(df)
    dfv = dfv[.!ismissing.(dfv[!, ycol]), :]
    sort!(dfv, ycol)

    nshow = min(n, nrow(dfv))
    keep_cols = [xcols; ANNUITY_COST_COLS; ANNUITY_REV_COLS; [ycol]]
    bestn = copy(dfv[1:nshow, keep_cols])

    for c in [xcols; ANNUITY_COST_COLS; ANNUITY_REV_COLS; [ycol]]
        bestn[!, c] = to_million.(bestn[!, c])
    end

    return bestn
end

function topn_annuity_rows(df::DataFrame, n::Int; ycol::String=YCOL_NO)
    dfv = filter_balance(df)
    dfv = dfv[.!ismissing.(dfv[!, ycol]), :]
    sort!(dfv, ycol)
    nshow = min(n, nrow(dfv))
    return copy(dfv[1:nshow, :])
end

function topn_co2(df::DataFrame, n::Int; xcols::Vector{String}=XCOLS, co2col::String=CO2_COL, ycol::String=YCOL_NO)
    dfv = filter_balance(df)
    dfv = dfv[.!ismissing.(dfv[!, co2col]), :]
    sort!(dfv, co2col)

    nshow = min(n, nrow(dfv))
    keep_cols = [xcols; [co2col, ycol]]
    bestn = copy(dfv[1:nshow, keep_cols])

    for c in xcols
        bestn[!, c] = to_million.(bestn[!, c])
    end
    bestn[!, co2col] = to_million.(bestn[!, co2col])
    bestn[!, ycol] = to_million.(bestn[!, ycol])

    return bestn
end

function save_topn_xlsx(topn::DataFrame; outdir::String=OUTDIR, filename::String="top100_annuity_no.xlsx")
    isdir(outdir) || mkpath(outdir)
    path = joinpath(outdir, filename)
    XLSX.writetable(path, Tables.columntable(topn); sheetname="Top100")
    return path
end

function save_topn_csv(topn::DataFrame; outdir::String=OUTDIR, filename::String)
    isdir(outdir) || mkpath(outdir)
    path = joinpath(outdir, filename)
    CSV.write(path, topn)
    return path
end

function plot_top10_annuity(best10::DataFrame; ycol::String=YCOL_NO, outdir::String=OUTDIR)
    isdir(outdir) || mkpath(outdir)
    x = collect(1:nrow(best10))
    x_cost = x .- 0.18
    x_rev = x .+ 0.18
    bw = 0.28
    total = Vector{Float64}(best10[!, ycol])
    cap = Vector{Float64}(best10[!, ANNUITY_COST_COLS[1]])
    misc = Vector{Float64}(best10[!, ANNUITY_COST_COLS[2]])
    op = Vector{Float64}(best10[!, ANNUITY_COST_COLS[3]])
    energy = Vector{Float64}(best10[!, ANNUITY_COST_COLS[4]])
    rev_control = Vector{Float64}(best10[!, ANNUITY_REV_COLS[1]])
    rev_feed = Vector{Float64}(best10[!, ANNUITY_REV_COLS[2]])

    p = plot(
        xlabel="Top-Rang",
        ylabel="Annuitaet [Mio. EUR]",
        legend=:outerright,
        grid=true,
        gridalpha=0.3,
    )

    stack_up = zeros(Float64, length(x))
    bar!(p, x_cost, stack_up .+ cap; fillrange=stack_up, bar_width=bw, label="A_cap")
    stack_up .+= cap
    bar!(p, x_cost, stack_up .+ misc; fillrange=stack_up, bar_width=bw, label="A_misc")
    stack_up .+= misc
    bar!(p, x_cost, stack_up .+ op; fillrange=stack_up, bar_width=bw, label="A_op")
    stack_up .+= op
    bar!(p, x_cost, stack_up .+ energy; fillrange=stack_up, bar_width=bw, label="A_energy")
    stack_up .+= energy

    top_down = copy(stack_up)
    bar!(p, x_rev, top_down; fillrange=top_down .- rev_control, bar_width=bw, label="A_rev_control")
    top_down .-= rev_control
    bar!(p, x_rev, top_down; fillrange=top_down .- rev_feed, bar_width=bw, label="A_rev_feed")

    for i in eachindex(x)
        plot!(
            p,
            [x_cost[i], x_rev[i]],
            [total[i], total[i]],
            linewidth=2.5,
            color=:black,
            label=(i == 1 ? "A_total" : false),
        )
    end
    xticks!(p, x, string.(x))

    savefig(p, joinpath(outdir, "top10_annuity_no.png"))
    return p
end

function plot_3d_parameter_space(df::DataFrame; outdir::String=OUTDIR)
    dfv = filter_balance(df)

    col_hp = XCOLS[1]
    col_eb = XCOLS[2]
    col_buf = XCOLS[3]
    col_obj = YCOL_NO

    dfv = dfv[
        .!ismissing.(dfv[!, col_hp]) .&&
        .!ismissing.(dfv[!, col_eb]) .&&
        .!ismissing.(dfv[!, col_buf]) .&&
        .!ismissing.(dfv[!, col_obj]),
        :
    ]

    x = Vector{Float64}(dfv[!, col_hp]) ./ 1e6
    y = Vector{Float64}(dfv[!, col_eb]) ./ 1e6
    z = Vector{Float64}(dfv[!, col_buf]) ./ 1e6
    objective = Vector{Float64}(dfv[!, col_obj]) ./ 1e6

    trace = PlotlyJS.scatter3d(
        x=x,
        y=y,
        z=z,
        mode="markers",
        marker=PlotlyJS.attr(
            size=5,
            color=objective,
            colorscale="Viridis",
            showscale=true,
            colorbar=PlotlyJS.attr(title="A_total [Mio. EUR]"),
        ),
        hovertemplate="HP=%{x:.2f} MW<br>EB=%{y:.2f} MW<br>Buffer=%{z:.2f} MWh<br>A_total=%{marker.color:.2f} Mio. EUR<extra></extra>",
    )

    layout = PlotlyJS.Layout(
        title="3D Parameterraum",
        scene=PlotlyJS.attr(
            xaxis=PlotlyJS.attr(title="HeatPump Power [MW]"),
            yaxis=PlotlyJS.attr(title="ElectrodeBoiler Power [MW]"),
            zaxis=PlotlyJS.attr(title="Buffer Capacity [MWh]"),
        ),
    )

    p = PlotlyJS.Plot(trace, layout)
    isdir(outdir) || mkpath(outdir)
    PlotlyJS.savefig(p, joinpath(outdir, "3d_parameter_space.html"))
    return p
end

function create_matrix_plot(
    df::DataFrame;
    xcols::Vector{String}=XCOLS,
    objective_col::String=YCOL_NO,
    outdir::String=OUTDIR,
    filename::String="matrix_plot.html",
)
    dfv = filter_balance(df)

    needed = [xcols; [objective_col]]
    mask = trues(nrow(dfv))
    for c in needed
        mask .&= .!ismissing.(dfv[!, c])
    end
    d = dfv[mask, :]

    values = Dict(c => Vector{Float64}(d[!, c]) for c in needed)

    dims = [
        PlotlyJS.attr(
            label=c,
            values=(contains(c, "/W") || contains(c, "/Wh") ? values[c] ./ 1e6 : values[c]),
        ) for c in xcols
    ]

    objective = values[objective_col] ./ 1e6

    trace = PlotlyJS.splom(
        dimensions=dims,
        marker=PlotlyJS.attr(
            color=objective,
            colorscale="Viridis",
            size=10,
            showscale=true,
            colorbar=PlotlyJS.attr(title="A_total [Mio. EUR]"),
        ),
    )

    p = PlotlyJS.plot(trace, PlotlyJS.Layout(title="Parameter-Matrix (SPLOM)"))
    isdir(outdir) || mkpath(outdir)
    PlotlyJS.savefig(p, joinpath(outdir, filename))
    return p
end

function read_prf_values(path::String)
    vals = Float64[]
    for line in eachline(path)
        s = strip(line)
        isempty(s) && continue
        startswith(s, "#") && continue
        parts = split(s, ';')
        length(parts) < 2 && continue
        push!(vals, something(tryparse(Float64, replace(strip(parts[2]), "," => ".")), 0.0))
    end
    return vals
end

function fit_length(v::Vector{Float64}, n::Int)
    if length(v) == n
        return v
    elseif isempty(v)
        return zeros(n)
    elseif length(v) > n
        return v[1:n]
    else
        out = Vector{Float64}(undef, n)
        for i in 1:n
            out[i] = v[mod1(i, length(v))]
        end
        return out
    end
end

function clean_col(df::DataFrame, col::String)
    @assert col in names(df) "Spalte '$col' fehlt in out.csv."
    return [ismissing(v) || isnan(v) ? 0.0 : abs(v) for v in to_float.(df[!, col])]
end

function detect_dt_hours(df::DataFrame)
    tcol = "Time [dd.mm.yyyy HH:MM:SS]"
    if !(tcol in names(df)) || nrow(df) < 2
        return 0.25
    end
    fmt = dateformat"dd.mm.yyyy HH:MM:SS"
    t1 = DateTime(df[1, tcol], fmt)
    t2 = DateTime(df[2, tcol], fmt)
    return Dates.value(t2 - t1) / 3_600_000
end

function remunerated_power_4h(power_w::Vector{Float64}, dt_h::Float64)
    steps = max(1, Int(round(4 / dt_h)))
    held = zeros(Float64, length(power_w))
    for idx in 1:steps:length(power_w)
        idx_end = min(idx + steps - 1, length(power_w))
        block = power_w[idx:idx_end]
        if !isempty(block) && all(isapprox.(block, block[1]; atol=1e-3, rtol=0.0))
            held[idx:idx_end] .= block[1]
        end
    end
    return held, steps
end

function count_held_blocks_4h(held_w::Vector{Float64}, steps_per_block::Int)
    c = 0
    for idx in 1:steps_per_block:length(held_w)
        idx_end = min(idx + steps_per_block - 1, length(held_w))
        if any(>(1e-9), held_w[idx:idx_end])
            c += 1
        end
    end
    return c
end


function main()
    df = CSV.read(CSV_PATH, DataFrame; delim=';', decimal=',')

    required_cols = vcat(XCOLS, ANNUITY_COST_COLS, ANNUITY_REV_COLS, [YCOL_NO, CO2_COL], BAL_COLS)
    missing_required = [c for c in required_cols if !(c in names(df))]
    @assert isempty(missing_required) "Fehlende Pflichtspalten: $(join(missing_required, ", "))"

    coerce_columns!(df, required_cols)

    best100 = topn_annuity(df, 100)
    xlsx_path = save_topn_xlsx(best100)
    top10_co2 = topn_co2(df, 10)
    co2_csv_path = save_topn_csv(top10_co2; filename="top10_lowest_co2_yearly.csv")

    best10 = topn_annuity(df, 10)
    best25_rows = topn_annuity_rows(df, 25)
    plot_top10_annuity(best10)
    plot_3d_parameter_space(df)
    create_matrix_plot(best25_rows; xcols=XCOLS[1:3], filename="matrix_plot_top25_no_battery.html")


    println("\nFertig. Dateien in: $OUTDIR")
    println("- $(basename(xlsx_path))")
    println("- $(basename(co2_csv_path))")
    println("- top10_annuity_no.png")
    println("- 3d_parameter_space.html")
    println("- matrix_plot_top25_no_battery.html")
    
end

main()
