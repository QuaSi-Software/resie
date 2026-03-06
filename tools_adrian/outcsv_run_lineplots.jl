using CSV
using DataFrames
using Dates
using JSON
import PlotlyJS

const DEFAULT_INFILE = "output/out.csv"
const DEFAULT_OUTDIR = "output/analysis_outcsv/run_lineplots"
const TIME_COLUMN = "Time [dd.mm.yyyy HH:MM:SS]"
const TIME_FORMAT = dateformat"dd.mm.yyyy HH:MM:SS"

# --------------------------- user configuration ---------------------------
# Set to `nothing` to disable start/end filtering.
# Input format: dd.mm.yyyy HH:MM:SS
const START_TIME = "01.01.2024 00:00:00"
const END_TIME = "31.12.2024 23:59:59"

# If your CSV has a run column, set it here (e.g. "run_id"), else keep `nothing`.
const RUN_COLUMN = nothing

const LEFT_AXIS_LABEL = "Leistung [kW]"
const RIGHT_AXIS_LABEL = "Preis [EUR/MWh] / SOC (%)"

# `col`: exact CSV column name
# `name`: legend display name (free text)
# `scale`: multiply raw values for plotting
const LEFT_SERIES = [
    #Demands
    (col="Demand_Heat Demand", name="Wärmebedarf", scale=1 / 1e3 * 4),
    (col="Demand_Power Demand", name="Strombedarf", scale=1 / 1e3 * 4),
    #Regelreserve
    (col="NegControlReserve_in Scaling_Factor", name="Negative Kapazitätsvorhaltung", scale=1 / 1e3),
    (col="PosControlReserve_out Scaling_Factor", name="Positive Kapazitätsvorhaltung", scale=1 / 1e3),
    (col="NegControlReserve_in m_power OUT", name="Negative Regelreserve", scale=1 / 1e3 * 4),
    (col="PosControlReserve_out m_power IN", name="Positive Regelreserve", scale=1 / 1e3 * 4),
    #Grid
    (col="Grid_IN m_power OUT", name="Strombezug aus Netz", scale=1 / 1e3 * 4),
    (col="Grid_OUT m_power IN", name="Stromeinspeisung ins Netz", scale=1 / 1e3 * 4),
    #PV
    (col="Photovoltaic m_power OUT", name="PV Erzeugung", scale=1 / 1e3 * 4),
    (col="m_power EnergyFlow Photovoltaic->Battery", name="PV -> Batterie", scale=1 / 1e3 * 4),
    (col="m_power EnergyFlow Photovoltaic->Demand_Power", name="PV -> Strombedarf", scale=1 / 1e3 * 4),
    (col="m_power EnergyFlow Photovoltaic->HeatPump", name="PV -> Wärmepumpe", scale=1 / 1e3 * 4),
    (col="m_power EnergyFlow Photovoltaic->ElectrodeBoiler", name="PV -> Elektrodenkessel", scale=1 / 1e3 * 4),
    (col="m_power EnergyFlow Photovoltaic->Grid_OUT", name="PV -> Stromnetz", scale=1 / 1e3 * 4),
    #Wind
    (col="WindFarm m_power OUT", name="Windenergie Erzeugung", scale=1 / 1e3 * 4),
    (col="m_power EnergyFlow WindFarm->Battery", name="Windenergie -> Batterie", scale=1 / 1e3 * 4),
    (col="m_power EnergyFlow WindFarm->Demand_Power", name="Windenergie -> Strombedarf", scale=1 / 1e3 * 4),
    (col="m_power EnergyFlow WindFarm->HeatPump", name="Windenergie -> Wärmepumpe", scale=1 / 1e3 * 4),
    (col="m_power EnergyFlow WindFarm->ElectrodeBoiler", name="Windenergie -> Elektrodenkessel", scale=1 / 1e3 * 4),
    (col="m_power EnergyFlow WindFarm->Grid_OUT", name="Windenergie -> Stromnetz", scale=1 / 1e3 * 4),
    #HeatPump
    (col="HeatPump_baseline_in m_power OUT", name="Wärmepumpe Strom-Input Baseline", scale=1 / 1e3 * 4),
    (col="HeatPump m_power IN", name="Wärmepumpe Strom-Input", scale=1 / 1e3 * 4),
    (col="HeatPump m_heat OUT", name="Wärmepumpe Wärme-Output #1", scale=1 / 1e3 * 4),
    (col="HeatPump secondary_m_heat OUT", name="Wärmepumpe Wärme-Output #2", scale=1 / 1e3 * 4),
    (col="m_heat EnergyFlow HeatPump->BufferTank", name="Wärmepumpe #1 -> Wärmespeicher", scale=1 / 1e3 * 4),
    (col="secondary_m_heat EnergyFlow secondary_HeatPump->BufferTank", name="Wärmepumpe #2 -> Wärmespeicher", scale=1 / 1e3 * 4),
    (col="m_heat EnergyFlow HeatPump->Demand_Heat", name="Wärmepumpe #1 -> Wärmebedarf", scale=1 / 1e3 * 4),
    (col="secondary_m_heat EnergyFlow secondary_HeatPump->Demand_Heat", name="Wärmepumpe #2 -> Wärmebedarf", scale=1 / 1e3 * 4),
    #ElectrodeBoiler
    (col="ElectrodeBoiler_baseline_in m_power OUT", name="Elektrodenkessel Strom-Input Baseline", scale=1 / 1e3 * 4),
    (col="ElectrodeBoiler m_power IN", name="Elektrodenkessel Strom-Input", scale=1 / 1e3 * 4),
    (col="ElectrodeBoiler m_heat OUT", name="Elektrodenkessel Wärme-Output #1", scale=1 / 1e3 * 4),
    (col="ElectrodeBoiler secondary_m_heat OUT", name="Elektrodenkessel Wärme-Output #2", scale=1 / 1e3 * 4),
    (col="m_heat EnergyFlow ElectrodeBoiler->BufferTank", name="Elektrodenkessel #1 -> Wärmespeicher", scale=1 / 1e3 * 4),
    (col="secondary_m_heat EnergyFlow secondary_ElectrodeBoiler->BufferTank", name="Elektrodenkessel #2 -> Wärmespeicher", scale=1 / 1e3 * 4),
    (col="m_heat EnergyFlow ElectrodeBoiler->Demand_Heat", name="Elektrodenkessel #1 -> Wärmebedarf", scale=1 / 1e3 * 4),
    (col="secondary_m_heat EnergyFlow secondary_ElectrodeBoiler->Demand_Heat", name="Elektrodenkessel #2 -> Wärmebedarf", scale=1 / 1e3 * 4),
    #BufferTank
    (col="m_heat EnergyFlow BufferTank->Demand_Heat", name="Wärmespeicher -> Bedarf", scale=1 / 1e3 * 4),
    #Battery
    (col="m_power EnergyFlow Battery->Demand_Power", name="Batterie -> Strombedarf", scale=1 / 1e3 * 4),
    (col="m_power EnergyFlow Battery->Grid_OUT", name="Batterie -> Stromnetz", scale=1 / 1e3 * 4),
    (col="m_power EnergyFlow Grid_IN->Battery", name="Stromnetz -> Batterie", scale=1 / 1e3 * 4),
]

const RIGHT_SERIES = [
    #Börsenstrompreis
    (col="StockMarket_IN m_money OUT", name="Börsenstrompreis", scale=1.0),
    #Speicherstände
    (col="BufferTank Load%", name="SOC Wärmespeicher", scale=1.0),
    (col="Battery Load%", name="SOC Batterie", scale=1.0),   
]

# Available out.csv columns (copy/paste helper):

# - "Battery Load%"
# - "Battery m_power IN"
# - "Battery m_power OUT"
# - "Battery SOC"
# - "m_heat EnergyFlow BufferTank->Demand_Heat"

# - "BufferTank Load"
# - "BufferTank Load%"
# - "BufferTank m_heat IN"
# - "BufferTank m_heat OUT"
# - "m_power EnergyFlow Battery->Demand_Power"

# - "Demand_Heat Demand"
# - "Demand_Power Demand"

# - "Grid_IN m_power OUT"
# - "Grid_OUT m_power IN"
# - "m_power EnergyFlow Grid_IN->Demand_Power"
# - "m_power EnergyFlow Grid_IN->ElectrodeBoiler"
# - "m_power EnergyFlow Grid_IN->HeatPump"
# - "StockMarket_IN m_money OUT"

# - "NegControlReserve_in m_power OUT"
# - "NegControlReserve_in Scaling_Factor"

# - "PosControlReserve_out m_power IN"
# - "PosControlReserve_out Scaling_Factor"

# - "Photovoltaic m_power OUT"
# - "m_power EnergyFlow Photovoltaic->Battery"
# - "m_power EnergyFlow Photovoltaic->Demand_Power"
# - "m_power EnergyFlow Photovoltaic->ElectrodeBoiler"
# - "m_power EnergyFlow Photovoltaic->Grid_OUT"
# - "m_power EnergyFlow Photovoltaic->HeatPump"

# - "WindFarm m_power OUT"
# - "m_power EnergyFlow WindFarm->Battery"
# - "m_power EnergyFlow WindFarm->Demand_Power"
# - "m_power EnergyFlow WindFarm->ElectrodeBoiler"
# - "m_power EnergyFlow WindFarm->Grid_OUT"
# - "m_power EnergyFlow WindFarm->HeatPump"

# - "HeatPump_baseline_in m_power OUT"
# - "HeatPump m_power IN"
# - "HeatPump m_heat OUT"
# - "HeatPump secondary_m_heat OUT"
# - "m_heat EnergyFlow HeatPump->BufferTank"
# - "m_heat EnergyFlow HeatPump->Demand_Heat"
# - "secondary_m_heat EnergyFlow secondary_HeatPump->BufferTank"
# - "secondary_m_heat EnergyFlow secondary_HeatPump->Demand_Heat"

# - "ElectrodeBoiler_baseline_in m_power OUT"
# - "ElectrodeBoiler m_power IN"
# - "ElectrodeBoiler m_heat OUT"
# - "ElectrodeBoiler secondary_m_heat OUT"
# - "m_heat EnergyFlow ElectrodeBoiler->BufferTank"
# - "m_heat EnergyFlow ElectrodeBoiler->Demand_Heat"
# - "secondary_m_heat EnergyFlow secondary_ElectrodeBoiler->BufferTank"
# - "secondary_m_heat EnergyFlow secondary_ElectrodeBoiler->Demand_Heat"


# Optional PNG export requires Plotly image backend setup.
const EXPORT_PNG = false
# ------------------------------------------------------------------------

to_float(v) = v === missing ? missing : (
    v isa Real ? Float64(v) : (
        v isa AbstractString ? something(tryparse(Float64, replace(strip(v), "," => ".")), missing) : missing
    )
)

function parse_input_datetime(s::Union{Nothing,AbstractString})
    s === nothing && return nothing
    x = tryparse(DateTime, strip(s), dateformat"dd.mm.yyyy HH:MM:SS")
    x === nothing && error("Ungueltiges Datumsformat '$s'. Erwartet: dd.mm.yyyy HH:MM:SS")
    return x
end

function sanitize_filename(s::AbstractString)
    replace(strip(s), r"[^A-Za-z0-9_\\-]+" => "_")
end

function read_out_csv(path::String)
    df = CSV.read(path, DataFrame; delim=';', decimal=',')
    @assert TIME_COLUMN in names(df) "Zeitspalte fehlt: $TIME_COLUMN"

    ts = [DateTime(x, TIME_FORMAT) for x in df[!, TIME_COLUMN]]
    for c in names(df)
        c == TIME_COLUMN && continue
        df[!, c] = to_float.(df[!, c])
    end
    return df, ts
end

function infer_run_index(ts::Vector{DateTime})
    ridx = Vector{Int}(undef, length(ts))
    run = 1
    ridx[1] = run
    for i in 2:length(ts)
        ts[i] <= ts[i - 1] && (run += 1)
        ridx[i] = run
    end
    return ridx
end

function get_run_groups(df::DataFrame, ts::Vector{DateTime})
    if RUN_COLUMN !== nothing && string(RUN_COLUMN) in names(df)
        rc = df[!, string(RUN_COLUMN)]
        ids = unique(rc)
        out = Vector{Tuple{String,Vector{Int}}}()
        for id in ids
            idx = findall(==(id), rc)
            push!(out, (string(id), idx))
        end
        return out
    end

    ridx = infer_run_index(ts)
    out = Vector{Tuple{String,Vector{Int}}}()
    for r in 1:maximum(ridx)
        idx = findall(==(r), ridx)
        push!(out, ("run_" * lpad(string(r), 3, '0'), idx))
    end
    return out
end

function filter_by_time(ts::Vector{DateTime}, start_dt::Union{Nothing,DateTime}, end_dt::Union{Nothing,DateTime})
    mask = trues(length(ts))
    start_dt !== nothing && (mask .&= ts .>= start_dt)
    end_dt !== nothing && (mask .&= ts .<= end_dt)
    findall(mask)
end

function clean_series(df::DataFrame, col::String, idx::Vector{Int}, scale::Real)
    col in names(df) || return nothing
    raw = df[idx, col]
    y = Vector{Union{Missing,Float64}}(undef, length(raw))
    for i in eachindex(raw)
        v = raw[i]
        if ismissing(v) || (v isa AbstractFloat && isnan(v))
            y[i] = missing
        else
            y[i] = Float64(v) * scale
        end
    end
    return y
end

function leapday_rangebreaks(xdt::Vector{DateTime})
    first_dt = first(xdt)
    last_dt = last(xdt)
    years = year(first_dt):year(last_dt)
    breaks = Any[]

    for y in years
        isleapyear(y) || continue
        leap_date = Date(y, 2, 29)
        leap_start = DateTime(leap_date)
        leap_end = leap_start + Day(1)
        in_window = (leap_start >= first_dt) && (leap_start <= last_dt)
        in_window || continue

        has_data_on_leapday = any(t -> (t >= leap_start && t < leap_end), xdt)
        if !has_data_on_leapday
            push!(breaks, PlotlyJS.attr(values=[Dates.format(leap_date, dateformat"yyyy-mm-dd")]))
        end
    end

    return breaks
end

function write_fullpage_html(fig::PlotlyJS.Plot, html_path::String)
    data_json = JSON.json(fig.data)
    layout_json = JSON.json(fig.layout)
    config_json = JSON.json(Dict(
        "responsive" => true,
        "displaylogo" => false,
        "scrollZoom" => true,
    ))

    html = """
<!doctype html>
<html lang="de">
<head>
  <meta charset="utf-8" />
  <meta name="viewport" content="width=device-width, initial-scale=1" />
  <title>Out.csv Linienplot</title>
  <script src="https://cdn.plot.ly/plotly-2.35.2.min.js"></script>
  <style>
    html, body {
      margin: 0;
      padding: 0;
      width: 100%;
      height: 100%;
      overflow: hidden;
      font-family: Arial, sans-serif;
    }
    #plot {
      width: 100vw;
      height: 100vh;
    }
  </style>
</head>
<body>
  <div id="plot"></div>
  <script>
    const data = $data_json;
    const layout = $layout_json;
    const config = $config_json;
    Plotly.newPlot('plot', data, layout, config);
  </script>
</body>
</html>
"""
    open(html_path, "w") do io
        write(io, html)
    end
end

function make_run_plot(
    df::DataFrame,
    ts::Vector{DateTime},
    run_label::String,
    run_idx::Vector{Int},
    outdir::String,
    start_dt::Union{Nothing,DateTime},
    end_dt::Union{Nothing,DateTime},
)
    ts_run = ts[run_idx]
    loc = filter_by_time(ts_run, start_dt, end_dt)
    isempty(loc) && return nothing

    idx = run_idx[loc]
    xdt = ts[idx]
    x_iso = Dates.format.(xdt, dateformat"yyyy-mm-ddTHH:MM:SS")
    x_start = first(x_iso)
    x_end = last(x_iso)
    rangebreaks = leapday_rangebreaks(xdt)

    left_colors = ["#E3A300", "#E54A00", "#1C8C5F", "#2A5DB0", "#8F3DB8"]
    right_colors = ["#F57C00", "#C62828", "#039BE5", "#6A1B9A", "#6D4C41"]
    traces = PlotlyJS.GenericTrace[]

    plotted = false
    for (i, s) in enumerate(LEFT_SERIES)
        y = clean_series(df, s.col, idx, s.scale)
        y === nothing && continue
        push!(
            traces,
            PlotlyJS.scatter(
                x=x_iso,
                y=y,
                mode="lines",
                name=s.name * " (left)",
                line=PlotlyJS.attr(width=2, color=left_colors[mod1(i, length(left_colors))]),
                yaxis="y",
                hovertemplate="%{x|%d.%m.%Y %H:%M:%S}<br>%{y:.3f}<extra>%{fullData.name}</extra>",
            ),
        )
        plotted = true
    end

    for (i, s) in enumerate(RIGHT_SERIES)
        y = clean_series(df, s.col, idx, s.scale)
        y === nothing && continue
        push!(
            traces,
            PlotlyJS.scatter(
                x=x_iso,
                y=y,
                mode="lines",
                name=s.name * " (right)",
                line=PlotlyJS.attr(width=2, color=right_colors[mod1(i, length(right_colors))]),
                yaxis="y2",
                hovertemplate="%{x|%d.%m.%Y %H:%M:%S}<br>%{y:.3f}<extra>%{fullData.name}</extra>",
            ),
        )
        plotted = true
    end

    plotted || return nothing

    layout = PlotlyJS.Layout(
        title="Out.csv Linienplot - " * run_label,
        template="plotly_white",
        xaxis=PlotlyJS.attr(
            title="",
            type="date",
            range=[x_start, x_end],
            domain=[0.0, 0.72],
            rangebreaks=rangebreaks,
            tickformat="%d.%m.%Y",
            showgrid=true,
            gridcolor="rgba(110,110,110,0.28)",
            gridwidth=1.0,
            griddash="solid",
            minor=PlotlyJS.attr(
                showgrid=true,
                gridcolor="rgba(110,110,110,0.14)",
                gridwidth=0.7,
                griddash="dot",
            ),
            showline=true,
            linewidth=1.4,
            linecolor="rgba(40,40,40,0.9)",
            mirror=true,
            ticks="outside",
            ticklen=6,
            tickwidth=1.2,
            tickcolor="rgba(40,40,40,0.9)",
            automargin=true,
            zeroline=false,
        ),
        yaxis=PlotlyJS.attr(
            title=LEFT_AXIS_LABEL,
            showgrid=true,
            gridcolor="rgba(110,110,110,0.28)",
            gridwidth=1.0,
            griddash="solid",
            minor=PlotlyJS.attr(
                showgrid=true,
                gridcolor="rgba(110,110,110,0.14)",
                gridwidth=0.7,
                griddash="dot",
            ),
            showline=true,
            linewidth=1.4,
            linecolor="rgba(40,40,40,0.9)",
            mirror=true,
            ticks="outside",
            ticklen=6,
            tickwidth=1.2,
            tickcolor="rgba(40,40,40,0.9)",
            automargin=true,
            zeroline=false,
        ),
        yaxis2=PlotlyJS.attr(
            title=RIGHT_AXIS_LABEL,
            overlaying="y",
            anchor="x",
            side="right",
            showgrid=false,
            zeroline=false,
            showline=true,
            linewidth=1.4,
            linecolor="rgba(40,40,40,0.9)",
            ticks="outside",
            ticklen=6,
            tickwidth=1.2,
            tickcolor="rgba(40,40,40,0.9)",
            automargin=true,
        ),
        legend=PlotlyJS.attr(
            x=0.75,
            y=1.0,
            xanchor="left",
            yanchor="top",
            bgcolor="rgba(255,255,255,0.9)",
            bordercolor="rgba(0,0,0,0.2)",
            borderwidth=1,
            traceorder="normal",
            itemsizing="trace",
        ),
        margin=PlotlyJS.attr(l=90, r=40, t=70, b=70),
        autosize=true,
    )

    fig = PlotlyJS.Plot(traces, layout)

    start_txt = Dates.format(first(xdt), "yyyymmdd_HHMMSS")
    end_txt = Dates.format(last(xdt), "yyyymmdd_HHMMSS")
    base = "lineplot_" * sanitize_filename(run_label) * "_" * start_txt * "__" * end_txt

    html_path = joinpath(outdir, base * ".html")
    write_fullpage_html(fig, html_path)

    if EXPORT_PNG
        try
            PlotlyJS.savefig(fig, joinpath(outdir, base * ".png"))
        catch err
            @warn "PNG-Export fehlgeschlagen. HTML wurde erstellt." exception=(err, catch_backtrace())
        end
    end

    return html_path
end

function main()
    inpath = length(ARGS) >= 1 ? ARGS[1] : DEFAULT_INFILE
    outdir = length(ARGS) >= 2 ? ARGS[2] : DEFAULT_OUTDIR
    mkpath(outdir)

    start_dt = parse_input_datetime(START_TIME)
    end_dt = parse_input_datetime(END_TIME)
    if start_dt !== nothing && end_dt !== nothing
        start_dt <= end_dt || error("START_TIME liegt nach END_TIME.")
    end

    df, ts = read_out_csv(inpath)
    @assert !isempty(ts) "Keine Daten in out.csv."

    groups = get_run_groups(df, ts)
    created = String[]
    skipped = String[]

    for (label, idx) in groups
        out = make_run_plot(df, ts, label, idx, outdir, start_dt, end_dt)
        out === nothing ? push!(skipped, label) : push!(created, out)
    end

    println("Lineplot-Export (Plotly) abgeschlossen.")
    println("Input: $inpath")
    println("Output-Ordner: $outdir")
    println("Datumsfenster: $(start_dt === nothing ? "-inf" : string(start_dt)) bis $(end_dt === nothing ? "+inf" : string(end_dt))")
    println("Runs erkannt: $(length(groups))")
    println("Plots erstellt: $(length(created))")
    !isempty(skipped) && println("Runs ohne Daten im Zeitfenster: " * join(skipped, ", "))
    for p in created
        println("- " * p)
    end
end

main()
