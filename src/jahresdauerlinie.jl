using CSV
using DataFrames
using Dates
using Plots

const DATE_FMT = dateformat"dd.mm.yyyy HH:MM:SS"

function resolve_column_name(df::DataFrame, candidates::Vector{String})
    names_str = String.(names(df))
    for c in candidates
        idx = findfirst(==(c), names_str)
        idx !== nothing && return names(df)[idx]
    end
    error("Keine passende Spalte gefunden. Erwartet eine von: $(join(candidates, ", "))")
end

function finite_or_zero(x)
    if x isa Missing
        return 0.0
    end
    v = Float64(x)
    return isfinite(v) ? v : 0.0
end

function infer_timestep_hours(df::DataFrame)
    time_col = Symbol("Time [dd.mm.yyyy HH:MM:SS]")
    (time_col in names(df)) || return 1.0
    nrow(df) < 2 && return 1.0
    t1_raw = df[1, time_col]
    t2_raw = df[2, time_col]
    t1 = t1_raw isa DateTime ? t1_raw : DateTime(String(t1_raw), DATE_FMT)
    t2 = t2_raw isa DateTime ? t2_raw : DateTime(String(t2_raw), DATE_FMT)
    dt_h = Dates.value(t2 - t1) / 3_600_000
    return dt_h > 0 ? dt_h : 1.0
end

function sorted_desc(vals::Vector{Float64})
    out = copy(vals)
    sort!(out; rev = true)
    return out
end

function energy_wh_to_power_w(vals_wh::Vector{Float64}, factor::Float64 = 4.0)
    return vals_wh .* factor
end

function main()
    input_path = length(ARGS) >= 1 ? ARGS[1] : "output/out.csv"
    output_old_csv = length(ARGS) >= 2 ? ARGS[2] : "output/jahresdauerlinie_heat_getrennt.csv"
    output_old_png = length(ARGS) >= 3 ? ARGS[3] : "output/jahresdauerlinie_heat_getrennt.png"
    output_new_csv = length(ARGS) >= 4 ? ARGS[4] : "output/jahresdauerlinie_heat_gestapelt.csv"
    output_new_png = length(ARGS) >= 5 ? ARGS[5] : "output/jahresdauerlinie_heat_gestapelt.png"

    df = DataFrame(CSV.File(input_path; delim = ';', decimal = ','))

    col_eb = resolve_column_name(df, [
        "m_heat EnergyFlow ElectrodeBoiler->Demand",
        "m_heat EnergyFlow ElectrodeBoiler->Demand_Heat",
    ])
    col_hp = resolve_column_name(df, [
        "m_heat EnergyFlow Heatpump->Demand",
        "m_heat EnergyFlow HeatPump->Demand",
        "m_heat EnergyFlow Heatpump->Demand_Heat",
        "m_heat EnergyFlow HeatPump->Demand_Heat",
    ])
    col_bt = resolve_column_name(df, [
        "m_heat EnergyFlow BufferTank->Demand",
        "m_heat EnergyFlow BufferTank->Demand_Heat",
    ])

    eb_raw = [finite_or_zero(x) for x in df[!, Symbol(col_eb)]]
    hp_raw = [finite_or_zero(x) for x in df[!, Symbol(col_hp)]]
    bt_raw = [finite_or_zero(x) for x in df[!, Symbol(col_bt)]]

    n = nrow(df)
    rank = collect(1:n)
    duration_h = collect(range(0.0, stop = 8760.0, length = n))

    # Alte Variante: jede Reihe separat als Jahresdauerlinie sortieren.
    eb_old = sorted_desc(eb_raw)
    hp_old = sorted_desc(hp_raw)
    bt_old = sorted_desc(bt_raw)
    eb_old_w = energy_wh_to_power_w(eb_old)
    hp_old_w = energy_wh_to_power_w(hp_old)
    bt_old_w = energy_wh_to_power_w(bt_old)
    old = DataFrame(
        rank = rank,
        duration_h = duration_h,
        heatpump_to_demand_heat_W = hp_old_w,
        electrodeboiler_to_demand_heat_W = eb_old_w,
        buffertank_to_demand_heat_W = bt_old_w,
    )
    CSV.write(output_old_csv, old)

    p_old = plot(
        duration_h,
        old.heatpump_to_demand_heat_W ./ 1e6,
        label = "HeatPump -> Demand",
        xlabel = "Dauer [h]",
        ylabel = "Leistung [MW]",
        linewidth = 2,
    )
    plot!(p_old, duration_h, old.electrodeboiler_to_demand_heat_W ./ 1e6, label = "ElectrodeBoiler -> Demand", linewidth = 2)
    plot!(p_old, duration_h, old.buffertank_to_demand_heat_W ./ 1e6, label = "BufferTank -> Demand", linewidth = 2)
    savefig(p_old, output_old_png)

    # Neue Variante: pro Zeitschritt zusammenfassen, nach Summe sortieren, gestapelt plotten.
    total_raw = eb_raw .+ hp_raw .+ bt_raw
    order = sortperm(total_raw; rev = true)

    eb = eb_raw[order]
    hp = hp_raw[order]
    bt = bt_raw[order]
    total = total_raw[order]
    eb_w = energy_wh_to_power_w(eb)
    hp_w = energy_wh_to_power_w(hp)
    bt_w = energy_wh_to_power_w(bt)
    total_w = energy_wh_to_power_w(total)

    out = DataFrame(
        rank = rank,
        duration_h = duration_h,
        heatpump_to_demand_heat_W = hp_w,
        electrodeboiler_to_demand_heat_W = eb_w,
        buffertank_to_demand_heat_W = bt_w,
        total_W = total_w,
    )
    CSV.write(output_new_csv, out)

    hp_mw = hp_w ./ 1e6
    eb_mw = eb_w ./ 1e6
    bt_mw = bt_w ./ 1e6
    cum1 = hp_mw
    cum2 = hp_mw .+ eb_mw
    cum3 = cum2 .+ bt_mw

    p = plot(
        duration_h,
        cum1,
        label = "HeatPump -> Demand",
        xlabel = "Dauer [h]",
        ylabel = "Leistung [MW]",
        linewidth = 1.2,
        color = :blue,
        fillrange = 0.0,
        fillalpha = 0.85,
    )
    plot!(
        p,
        duration_h,
        cum2,
        label = "ElectrodeBoiler -> Demand",
        linewidth = 1.2,
        color = :orange,
        fillrange = cum1,
        fillalpha = 0.85,
    )
    plot!(
        p,
        duration_h,
        cum3,
        label = "BufferTank -> Demand",
        linewidth = 1.2,
        color = :green,
        fillrange = cum2,
        fillalpha = 0.85,
    )
    savefig(p, output_new_png)

    println("Jahresdauerlinien geschrieben:")
    println("  Alt (getrennt) CSV: $output_old_csv")
    println("  Alt (getrennt) PNG: $output_old_png")
    println("  Neu (gestapelt) CSV: $output_new_csv")
    println("  Neu (gestapelt) PNG: $output_new_png")
end

main()
