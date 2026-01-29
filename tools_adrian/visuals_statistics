############################################################
# Parameterstudie-Auswertung (CSV) – 4 Einzelparameter-Plots
#
# Ziel:
#  - CSV einlesen
#  - Runs mit Balance-Warnings (balance_power > 1 oder balance_heat > 1) rauswerfen
#  - Für den Basisfall (annuity_no …) je einen Plot erstellen:
#      y = "annuity_no A_total / €"
#      x = jeweils einer der 4 Anlagenparameter
#    Dabei müssen die anderen 3 Parameter auf den Fixwerten liegen.
#
# Fixwerte:
#  Hp_Power = 6.1e6 W
#  Boiler_Power = 3.5e6 W
#  BufferTank_Capacity = 6.5e7 Wh
#  Battery_Capacity = 200e3 Wh
############################################################

using CSV
using DataFrames
using Statistics
using Printf
using Plots

# ----------------------------
# User-Einstellungen
# ----------------------------
const CSV_PATH = "C:/Users/jenter/Documents/resie/output/parameterstudy/results_240runs_260127_113848.csv"   # <- anpassen
const OUTDIR   = "C:/Users/jenter/Documents/resie/output/parameterstudy/plots"                               # Ordner für PNGs

# Fixpunkt (Baseline), an dem "die anderen Parameter" festgehalten werden
const FIX = Dict(
    "Hp_Power / W"             => 6.1e6,
    "Boiler_Power / W"         => 3.5e6,
    "BufferTank_Capacity / Wh" => 6.5e7,
    "Battery_Capacity / Wh"    => 200e3,
)

# Zielgröße: Gesamtannuität im Basisfall
const YCOL = "annuity_no A_total / €"

# Balance-Warning-Regel: Werte > 1 sollen ausgeschlossen werden
const BAL_COLS = ["balance_power", "balance_heat"]
const BAL_THRESHOLD = 1.0

# Toleranz für Float-Vergleiche (wichtig, falls Werte aus Berechnungen kommen)
# Wenn deine CSV wirklich exakt die Fixwerte schreibt, kannst du eps kleiner machen.
const REL_TOL = 1e-9
const ABS_TOL = 1e-6

# ----------------------------
# Hilfsfunktionen
# ----------------------------

# "Nahe Gleichheit" für Floats (robust gegen Rundungsartefakte)
isapprox_num(a::Real, b::Real; rtol=REL_TOL, atol=ABS_TOL) = isapprox(a, b; rtol=rtol, atol=atol)

# Prüft, ob eine Zeile bei allen Fixparametern (außer dem variierenden) passt
function row_matches_fix(row, varying_col::String)
    for (col, val) in FIX
        col == varying_col && continue
        # Missing abfangen
        if ismissing(row[col]) || ismissing(val)
            return false
        end
        # Robust vergleichen
        if !(row[col] isa Real) || !(val isa Real)
            return false
        end
        if !isapprox_num(row[col], val)
            return false
        end
    end
    return true
end

# Filtert DataFrame auf:
#  - keine Balance-Warnings
#  - Fixpunkt für alle Parameter außer varying_col
function subset_for_single_param(df::DataFrame, varying_col::String)
    # 1) Balance-Warnings raus
    df2 = df
    for bcol in BAL_COLS
        @assert hasproperty(df2, Symbol(bcol)) || (bcol in names(df2)) "Spalte '$bcol' fehlt in der CSV."
        df2 = df2[.!ismissing.(df2[!, bcol]) .&& (df2[!, bcol] .<= BAL_THRESHOLD), :]
    end

    # 2) Fixpunkt (andere Parameter festhalten)
    keep = Vector{Bool}(undef, nrow(df2))
    for i in 1:nrow(df2)
        keep[i] = row_matches_fix(df2[i, :], varying_col)
    end
    df3 = df2[keep, :]

    return df3
end

# Macht einen Plot (x = varying_col, y = YCOL), sortiert nach x, speichert als PNG
function plot_single_param(df::DataFrame, varying_col::String; outdir::String=OUTDIR)
    @assert varying_col in names(df) "Spalte '$varying_col' nicht gefunden."
    @assert YCOL in names(df) "Spalte '$YCOL' nicht gefunden."

    # Nur Zeilen mit gültigem y und x
    dfp = df[.!ismissing.(df[!, varying_col]) .&& .!ismissing.(df[!, YCOL]), :]

    # Sortieren nach x (damit Linienplot sauber ist)
    sort!(dfp, varying_col)

    # Kurze Statistik (optional hilfreich)
    n = nrow(dfp)
    if n == 0
        @warn "Keine Datenpunkte für '$varying_col' (nach Filtern)."
        return nothing
    end

    ymin, iy = findmin(dfp[!, YCOL])
    xmin_at_ymin = dfp[iy, varying_col]

    @printf("\n[%s] Punkte: %d | min(annuity_no) = %.3e € bei %s = %.3e\n",
            varying_col, n, ymin, varying_col, xmin_at_ymin)

    # Plot: Punkte + Linie
    p = plot(
        dfp[!, varying_col], dfp[!, YCOL],
        seriestype = :line,
        linewidth = 2,
        marker = :circle,
        markersize = 4,
        xlabel = varying_col,
        ylabel = YCOL,
        title  = "Gesamtannuität (no) vs. $(varying_col)\n(Fix: übrige Parameter am Basispunkt, Balance ≤ $(BAL_THRESHOLD))",
        legend = false,
    )

    # Min markieren
    scatter!(
        p, [xmin_at_ymin], [ymin],
        marker = :star5,
        markersize = 9,
        label = false
    )

    # Speichern
    isdir(outdir) || mkpath(outdir)
    fname = replace(varying_col, r"[^\w]+" => "_")  # Dateiname säubern
    savefig(p, joinpath(outdir, "annuity_no_vs_$(fname).png"))

    return p
end

# ----------------------------
# Main
# ----------------------------

# CSV einlesen (wichtig: Strings als Spaltennamen beibehalten)
df = CSV.read(CSV_PATH, DataFrame)

# Sicherstellen, dass die wichtigsten Spalten da sind
needed = vcat(collect(keys(FIX)), [YCOL], BAL_COLS)
missing_cols = setdiff(needed, names(df))
if !isempty(missing_cols)
    error("CSV fehlt Spalten: $(missing_cols)")
end

# Vier Plots erzeugen (jeweils genau ein Parameter variiert)
PARAMS = [
    "Hp_Power / W",
    "Boiler_Power / W",
    "BufferTank_Capacity / Wh",
    "Battery_Capacity / Wh",
]

plots = Dict{String, Any}()

for varying in PARAMS
    sub = subset_for_single_param(df, varying)
    p = plot_single_param(sub, varying)
    plots[varying] = p
end

println("\nFertig. PNGs liegen in: $(OUTDIR)/")
println("Wenn bei einem Parameter 0 Punkte rauskommen: dann existiert in der CSV kein Run, der *genau* am Fixpunkt liegt (Float-Toleranz / Diskretisierung).")
