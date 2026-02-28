using CSV
using DataFrames
using Dates
using Statistics
using Printf
using Plots

const DEFAULT_INFILE = "output/out.csv"
const DEFAULT_OUTDIR = "output/analysis_outcsv"
const EPS = 1e-9

to_float(v) = v === missing ? missing : (
    v isa Real ? Float64(v) : (
        v isa AbstractString ? something(tryparse(Float64, replace(strip(v), "," => ".")), missing) : missing
    )
)

function parse_time_column(col)
    fmt = dateformat"dd.mm.yyyy HH:MM:SS"
    return [DateTime(x, fmt) for x in col]
end

function detect_dt_hours(ts::Vector{DateTime})
    if length(ts) < 2
        return 0.25
    end
    diffs = [Dates.value(ts[i + 1] - ts[i]) for i in 1:(length(ts) - 1)] # milliseconds
    return median(diffs) / 3_600_000
end

function read_out_csv(path::String)
    df = CSV.read(path, DataFrame; delim=';', decimal=',')
    @assert "Time [dd.mm.yyyy HH:MM:SS]" in names(df) "Zeitspalte fehlt: Time [dd.mm.yyyy HH:MM:SS]"

    time_col = "Time [dd.mm.yyyy HH:MM:SS]"
    ts = parse_time_column(df[!, time_col])

    for c in names(df)
        c == time_col && continue
        df[!, c] = to_float.(df[!, c])
    end

    dt_h = detect_dt_hours(ts)
    return df, ts, dt_h
end

colsum(df::DataFrame, col::String) = col in names(df) ? sum(skipmissing(replace(df[!, col], NaN => missing)); init=0.0) : 0.0
hascol(df::DataFrame, col::String) = col in names(df)

function positive_hours(series, dt_h::Real)
    vals = collect(skipmissing(replace(series, NaN => missing)))
    return count(>(EPS), vals) * dt_h
end

function peak_power_w(series, dt_h::Real)
    vals = collect(skipmissing(replace(series, NaN => missing)))
    isempty(vals) && return 0.0
    return maximum(vals) / dt_h
end

function make_duration_curve(df::DataFrame, dt_h::Real, outdir::String)
    generators = Dict(
        "HeatPump" => "HeatPump m_heat OUT",
        "ElectrodeBoiler" => "ElectrodeBoiler m_heat OUT",
    )
    available = filter(p -> hascol(df, p.second), collect(generators))
    isempty(available) && return nothing

    p = plot(
        xlabel = "Dauer [h]",
        ylabel = "Waermeleistung [MW]",
        title = "Jahresdauerlinien Waermeerzeuger",
        legend = :topright,
        linewidth = 2,
    )

    for (name, col) in available
        vals = [v for v in replace(df[!, col], NaN => missing) if !ismissing(v) && v >= 0]
        isempty(vals) && continue
        sorted_vals = sort(vals; rev=true)
        x = (1:length(sorted_vals)) .* dt_h
        y = sorted_vals ./ dt_h ./ 1e6
        plot!(p, x, y; label=name)
    end

    png_path = joinpath(outdir, "jahresdauerlinien_waermeerzeuger.png")
    savefig(p, png_path)
    return png_path
end

function unit_statistics(df::DataFrame, dt_h::Real)
    units = [
        ("HeatPump", "Waerme", "HeatPump m_heat OUT"),
        ("ElectrodeBoiler", "Waerme", "ElectrodeBoiler m_heat OUT"),
        ("Photovoltaic", "Strom", "Photovoltaic m_power OUT"),
        ("WindFarm", "Strom", "WindFarm m_power OUT"),
        ("Battery (Discharge)", "Strom", "Battery m_power OUT"),
        ("BufferTank (Discharge)", "Waerme", "BufferTank m_heat OUT"),
    ]

    rows = NamedTuple[]
    for (unit, medium, col) in units
        hascol(df, col) || continue
        e_wh = colsum(df, col)
        h = positive_hours(df[!, col], dt_h)
        p_peak_w = peak_power_w(df[!, col], dt_h)
        flh = p_peak_w > EPS ? e_wh / p_peak_w : missing
        push!(rows, (
            anlage = unit,
            medium = medium,
            energie_mwh = e_wh / 1e6,
            betriebsstunden_h = h,
            spitzenleistung_mw = p_peak_w / 1e6,
            volllaststunden_h = flh,
        ))
    end
    return DataFrame(rows)
end

function partial_load_statistics(df::DataFrame, dt_h::Real)
    units = [
        ("HeatPump", "HeatPump Avg_PLR", "HeatPump m_heat OUT"),
        ("ElectrodeBoiler", "ElectrodeBoiler Avg_PLR", "ElectrodeBoiler m_heat OUT"),
    ]
    bins = collect(0.0:0.1:1.0)
    rows = NamedTuple[]

    for (unit, plr_col, fallback_energy_col) in units
        plr = nothing
        if hascol(df, plr_col)
            vals = replace(df[!, plr_col], NaN => missing)
            plr = [v for v in vals if !ismissing(v) && v >= 0]
        elseif hascol(df, fallback_energy_col)
            p_peak = peak_power_w(df[!, fallback_energy_col], dt_h)
            if p_peak > EPS
                vals = replace(df[!, fallback_energy_col], NaN => missing)
                plr = [min(max((v / dt_h) / p_peak, 0.0), 1.0) for v in vals if !ismissing(v)]
            end
        end

        (plr === nothing || isempty(plr)) && continue
        active = [v for v in plr if v > EPS]
        isempty(active) && continue

        for i in 1:(length(bins) - 1)
            lo, hi = bins[i], bins[i + 1]
            inbin = i < length(bins) - 1 ? count(v -> (v >= lo && v < hi), active) : count(v -> (v >= lo && v <= hi), active)
            hours = inbin * dt_h
            share = inbin / length(active)
            push!(rows, (
                anlage = unit,
                teillastklasse = @sprintf("%.0f-%.0f%%", lo * 100, hi * 100),
                stunden_h = hours,
                anteil_aktive_zeit = share,
            ))
        end
    end
    return DataFrame(rows)
end

function storage_cycles(df::DataFrame)
    rows = NamedTuple[]

    battery_cap = hascol(df, "Battery Capacity") ? median(skipmissing(replace(df[!, "Battery Capacity"], NaN => missing))) : 0.0
    battery_in = colsum(df, "Battery m_power IN")
    battery_out = colsum(df, "Battery m_power OUT")
    push!(rows, (
        speicher = "Battery",
        kapazitaet_mwh = battery_cap / 1e6,
        zyklen_entladung = battery_cap > EPS ? battery_out / battery_cap : 0.0,
        zyklen_durchsatz = battery_cap > EPS ? (battery_in + battery_out) / (2 * battery_cap) : 0.0,
    ))

    buff_cap = hascol(df, "BufferTank Capacity") ? median(skipmissing(replace(df[!, "BufferTank Capacity"], NaN => missing))) : 0.0
    buff_in = colsum(df, "BufferTank m_heat IN")
    buff_out = colsum(df, "BufferTank m_heat OUT")
    push!(rows, (
        speicher = "BufferTank",
        kapazitaet_mwh = buff_cap / 1e6,
        zyklen_entladung = buff_cap > EPS ? buff_out / buff_cap : 0.0,
        zyklen_durchsatz = buff_cap > EPS ? (buff_in + buff_out) / (2 * buff_cap) : 0.0,
    ))

    return DataFrame(rows)
end

function power_key_figures(df::DataFrame)
    pv_gen = colsum(df, "Photovoltaic m_power OUT")
    wind_gen = colsum(df, "WindFarm m_power OUT")

    pv_self = colsum(df, "m_power EnergyFlow Photovoltaic->Demand_Power") +
              colsum(df, "m_power EnergyFlow Photovoltaic->Battery") +
              colsum(df, "m_power EnergyFlow Photovoltaic->HeatPump") +
              colsum(df, "m_power EnergyFlow Photovoltaic->ElectrodeBoiler")

    wind_self = colsum(df, "m_power EnergyFlow WindFarm->Demand_Power") +
                colsum(df, "m_power EnergyFlow WindFarm->Battery") +
                colsum(df, "m_power EnergyFlow WindFarm->HeatPump") +
                colsum(df, "m_power EnergyFlow WindFarm->ElectrodeBoiler")

    grid_total = colsum(df, "m_power EnergyFlow Grid_IN->Demand_Power") +
                 colsum(df, "m_power EnergyFlow Grid_IN->HeatPump") +
                 colsum(df, "m_power EnergyFlow Grid_IN->ElectrodeBoiler") +
                 colsum(df, "m_power EnergyFlow Grid_IN->Battery")

    own_total = pv_self + wind_self
    autarky_total = (own_total + grid_total) > EPS ? own_total / (own_total + grid_total) : missing

    own_to_demand = colsum(df, "m_power EnergyFlow Photovoltaic->Demand_Power") +
                    colsum(df, "m_power EnergyFlow WindFarm->Demand_Power") +
                    colsum(df, "m_power EnergyFlow Battery->Demand_Power")
    grid_to_demand = colsum(df, "m_power EnergyFlow Grid_IN->Demand_Power")
    autarky_demand = (own_to_demand + grid_to_demand) > EPS ? own_to_demand / (own_to_demand + grid_to_demand) : missing

    return DataFrame([(
        pv_zeugung_mwh = pv_gen / 1e6,
        wind_zeugung_mwh = wind_gen / 1e6,
        pv_eigenverbrauchsanteil = pv_gen > EPS ? pv_self / pv_gen : missing,
        wind_eigenverbrauchsanteil = wind_gen > EPS ? wind_self / wind_gen : missing,
        strombezug_netz_mwh = grid_total / 1e6,
        strombezug_eigenerzeugung_mwh = own_total / 1e6,
        autarkiegrad_gesamtstrom = autarky_total,
        autarkiegrad_demand_power = autarky_demand,
    )])
end

function reserve_metrics(df::DataFrame, dt_h::Real)
    req = ["PosControlReserve_in m_power OUT", "NegControlReserve_in m_power OUT",
           "PosControlReserve_in Scaling_Factor", "NegControlReserve_in Scaling_Factor"]
    if !all(c -> hascol(df, c), req)
        return DataFrame()
    end

    function vec_clean(col::String)
        [x for x in replace(df[!, col], NaN => missing) if !ismissing(x)]
    end

    e_pos = sum(abs, vec_clean("PosControlReserve_in m_power OUT"))
    e_neg = sum(abs, vec_clean("NegControlReserve_in m_power OUT"))
    ratio_pos_neg = e_neg > EPS ? e_pos / e_neg : missing

    p_pos = vec_clean("PosControlReserve_in Scaling_Factor")
    p_neg = vec_clean("NegControlReserve_in Scaling_Factor")
    n = min(length(p_pos), length(p_neg))
    n == 0 && return DataFrame()
    p_pos = abs.(p_pos[1:n])
    p_neg = abs.(p_neg[1:n])

    function remunerated_power_4h(power_w::Vector{Float64}, dt_h::Real)
        steps_per_block = max(1, Int(round(4 / dt_h)))
        held = zeros(Float64, length(power_w))
        for idx in 1:steps_per_block:length(power_w)
            idx_end = min(idx + steps_per_block - 1, length(power_w))
            block = power_w[idx:idx_end]
            if !isempty(block) && all(isapprox.(block, block[1]; atol=1e-3, rtol=0.0))
                held[idx:idx_end] .= block[1]
            end
        end
        return held
    end

    hold_pos = remunerated_power_4h(p_pos, dt_h)
    hold_neg = remunerated_power_4h(p_neg, dt_h)

    offer_mwh = sum((p_pos .+ p_neg) .* dt_h) / 1e6
    hold_mwh = sum((hold_pos .+ hold_neg) .* dt_h) / 1e6
    freebid_mwh = max(offer_mwh - hold_mwh, 0.0)
    ratio_freebid_hold = hold_mwh > EPS ? freebid_mwh / hold_mwh : missing

    return DataFrame([(
        regelenergie_pos_mwh = e_pos / 1e6,
        regelenergie_neg_mwh = e_neg / 1e6,
        verhaeltnis_pos_zu_neg_regelenergie = ratio_pos_neg,
        angebotene_reserve_mwh = offer_mwh,
        kapazitaetsvorhaltung_mwh = hold_mwh,
        freebids_mwh = freebid_mwh,
        verhaeltnis_freebids_zu_kapazitaetsvorhaltung = ratio_freebid_hold,
    )])
end

function main()
    inpath = length(ARGS) >= 1 ? ARGS[1] : DEFAULT_INFILE
    outdir = length(ARGS) >= 2 ? ARGS[2] : DEFAULT_OUTDIR
    mkpath(outdir)

    df, ts, dt_h = read_out_csv(inpath)
    @assert !isempty(ts) "Keine Zeitreihen in out.csv gefunden."

    duration_path = make_duration_curve(df, dt_h, outdir)
    df_units = unit_statistics(df, dt_h)
    df_partload = partial_load_statistics(df, dt_h)
    df_cycles = storage_cycles(df)
    df_power = power_key_figures(df)
    df_reserve = reserve_metrics(df, dt_h)

    CSV.write(joinpath(outdir, "anlagen_kennzahlen.csv"), df_units; delim=';')
    CSV.write(joinpath(outdir, "teillast_statistik.csv"), df_partload; delim=';')
    CSV.write(joinpath(outdir, "ladezyklen.csv"), df_cycles; delim=';')
    CSV.write(joinpath(outdir, "strom_kennzahlen.csv"), df_power; delim=';')
    if !isempty(df_reserve)
        CSV.write(joinpath(outdir, "regelreserve_kennzahlen.csv"), df_reserve; delim=';')
    end

    println("Analyse abgeschlossen.")
    println("Input: $inpath")
    println("Zeitauflosung: $(round(dt_h; digits=4)) h")
    println("Ausgabeordner: $outdir")
    duration_path !== nothing && println("Plot: $(basename(duration_path))")
    println("CSVs:")
    println("- anlagen_kennzahlen.csv")
    println("- teillast_statistik.csv")
    println("- ladezyklen.csv")
    println("- strom_kennzahlen.csv")
    !isempty(df_reserve) && println("- regelreserve_kennzahlen.csv")
end

main()
