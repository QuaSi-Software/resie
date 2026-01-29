using CSV
using DataFrames
using Dates

println("== Start Script ==")

# ------------------------------------------------------------
# Einstellungen
# ------------------------------------------------------------
data_dir = raw"C:/Users/jenter/Documents/Backup Masterarbeit/CBMP"  # <-- anpassen
out_file = joinpath(data_dir, "aFRR_prices_2024_year_no_leapday.csv")

col_isp  = Symbol("ISP (CET/CEST)")
col_area = :Area
col_up   = Symbol("Price Up (EUR/MWh)")
col_down = Symbol("Price Down (EUR/MWh)")

isp_dt_format = dateformat"dd/mm/yyyy HH:MM:SS"

# ------------------------------------------------------------
# Hilfsfunktionen
# ------------------------------------------------------------
function day_from_filename(fname::AbstractString)::Date
    m = match(r"-(\d{8})\d{4}\.csv$", fname)
    m === nothing && error("Konnte Datum nicht aus Dateiname extrahieren: $fname")
    return Date(m.captures[1], dateformat"yyyymmdd")
end

function parse_isp_start(isp::AbstractString)::DateTime
    # robust: nimmt nur den ersten Timestamp am Anfang
    m = match(r"^\s*(\d{2}/\d{2}/\d{4} \d{2}:\d{2}:\d{2})", isp)
    m === nothing && error("Konnte Startzeit nicht aus ISP extrahieren: $(repr(isp))")
    return DateTime(m.captures[1], isp_dt_format)
end

# ------------------------------------------------------------
# Dateien sammeln
# ------------------------------------------------------------
println("Data dir: ", data_dir)

files = filter(f -> endswith(lowercase(f), ".csv"), readdir(data_dir))
println("Gefundene CSVs: ", length(files))
isempty(files) && error("Keine CSV-Dateien im Ordner gefunden: $data_dir")

sort!(files; by=day_from_filename)
println("Erste Datei: ", first(files))
println("Letzte Datei: ", last(files))

# ------------------------------------------------------------
# Einlesen & zusammenfügen (nur 2024-Zeilen)
# ------------------------------------------------------------
dfs = DataFrame[]

for (i, f) in enumerate(files)
    fp = joinpath(data_dir, f)
    println("Lese Datei $i / $(length(files)): $f")

    df = CSV.read(fp, DataFrame;
        delim = ',',
        quotechar = '"',
        ignorerepeated = true,
        silencewarnings = true,
        pool = false,
        missingstring = ["", "NA", "NaN", "N/A"],
        types = Dict(
            col_isp  => Union{Missing,String},
            col_area => Union{Missing,String},
            col_up   => Union{Missing,Float64},
            col_down => Union{Missing,Float64},
        ),
    )

    dropmissing!(df, col_isp)

    # ISP -> Startzeitpunkt
    df[!, :_dt] = parse_isp_start.(df[!, col_isp])

    # Nur Jahr 2024 behalten
    df = df[year.(df[!, :_dt]) .== 2024, :]

    if nrow(df) == 0
        continue
    end

    # Area ignorieren
    df2 = select(df, col_isp, col_up, col_down, :_dt)
    push!(dfs, df2)
end

isempty(dfs) && error("Nach Filterung auf Jahr 2024 sind keine Daten übrig geblieben.")

println("Concatenate…")
df_year = vcat(dfs...)
sort!(df_year, :_dt)

# ------------------------------------------------------------
# Schalttag entfernen + timestamp so berechnen, als wäre 2024 kein Schaltjahr
# ------------------------------------------------------------

# 1) Schalttag (29.02.2024) entfernen (falls doch vorhanden)
df_year = df_year[Date.(df_year[!, :_dt]) .!= Date(2024, 2, 29), :]

# 2) timestamp: Sekunden ab 01.01.2024 00:00:00,
#    aber für alle Zeitpunkte ab 01.03.2024 einen Tag (86400 s) abziehen
t0    = DateTime(2024, 1, 1, 0, 0, 0)
t_cut = DateTime(2024, 3, 1, 0, 0, 0)

df_year.timestamp = Int.(Dates.value.(df_year[!, :_dt] .- t0) .÷ 1000)

mask = df_year[!, :_dt] .>= t_cut
df_year.timestamp[mask] .-= 86400

# Spalten final
df_year = select(df_year, col_isp, :timestamp, col_up, col_down)

# Plausibilität: 4s Raster ohne Sprünge
Δt = diff(df_year.timestamp)
if any(Δt .!= 4)
    bad = unique(Δt[Δt .!= 4])
    @warn "Nicht alle Zeitschritte sind 4 Sekunden. Abweichungen (unique): $(bad[1:min(end, 10)])"
end

println("Schreibe Output: ", out_file)
CSV.write(out_file, df_year; delim = ',', quotechar = '"')

println("== Fertig ==")
println("Zeilen: ", nrow(df_year))
println("Start dt: ", first(df_year[!, col_isp]))
println("Ende dt: ", last(df_year[!, col_isp]))
println("timestamp Ende (s): ", last(df_year.timestamp))
