############################################################
# Parameterstudie-Auswertung (CSV) – NUR Szenario "no"
# ---------------------------------------------------
# Fixes gegenüber deiner letzten Version:
# 1) Balance-Warnings: ALLE Runs mit |balance| > 1 werden gefiltert (Betrag!)
# 2) Spearman: Statistics.cor unterstützt kein keyword "method"
#    → Spearman wird berechnet als Pearson-Korrelation der Ränge (inkl. Tie-Handling)
# 3) World-Age-Warnung (Julia 1.12 + Revise): main wird über invokelatest gestartet
############################################################

using CSV                                                # CSV-Dateien einlesen
using DataFrames                                         # Tabellenverarbeitung
using Statistics                                         # mean/std/cor
using Printf                                             # formatierte Ausgaben
using Plots                                              # Plots

############################################################
# Konfiguration (hier anpassen)
############################################################
const CSV_PATH = "c:/Users/jenter/Documents/resie/output/parameterstudy/results_1050runs_260130_090710.csv"             # <- Pfad zur CSV anpassen
const OUTDIR   = "c:/Users/jenter/Documents/resie/output/parameterstudy/plots/"                                 # Plot-Ausgabeordner

############################################################
# Spaltennamen (wie in deiner CSV)
############################################################
const XCOLS = [                                          # Einflussgrößen (Anlagenparameter)
    "Hp_Power / W",                                       # Wärmepumpenleistung
    "Boiler_Power / W",                                   # Boiler/Heizstab-Leistung
    "BufferTank_Capacity / Wh",                           # Wärmespeichergröße
    "Battery_Capacity / Wh",                              # Batteriegröße
]

const YCOL_NO  = "annuity_no A_total / €"                # Zielgröße (Basisfall "no")
const BAL_COLS = ["balance_power", "balance_heat"]       # Balance-Warning-Spalten
const BAL_THRESHOLD = 1.0                                 # alle |balance| > 1 werden verworfen

############################################################
# Fixpunkt (für 1D-Slices)
############################################################
const FIX = Dict(                                        # Fixwerte für 1D-Plots
    "Hp_Power / W"             => 6.5e6,                  # W
    "Boiler_Power / W"         => 3.6e6,                  # W
    "BufferTank_Capacity / Wh" => 6.5e7,                  # Wh
    "Battery_Capacity / Wh"    => 350e3,                  # Wh
)

############################################################
# Toleranzen für Fixpunktvergleich (isapprox)
############################################################
const REL_TOL = 1e-9                                      # relative Toleranz
const ABS_TOL = 1e-6                                      # absolute Toleranz

############################################################
# Hilfsfunktion: Zellwert → Float64 (oder missing)
############################################################
function to_float(x)                                      # beliebigen Zellwert konvertieren
    if x === missing                                      # Missing bleibt missing
        return missing                                    # Rückgabe: missing
    elseif x isa Real                                     # bereits numerisch?
        return Float64(x)                                 # in Float64 umwandeln
    elseif x isa AbstractString                            # String?
        s = strip(x)                                      # Whitespaces entfernen
        isempty(s) && return missing                       # leere Strings → missing
        s2 = replace(s, "," => ".")                        # Dezimal-Komma → Punkt
        v = tryparse(Float64, s2)                          # sicher parsen
        return v === nothing ? missing : v                 # parsefehler → missing
    else
        return missing                                    # unbekannter Typ → missing
    end
end

############################################################
# Hilfsfunktion: mehrere Spalten auf Float64 zwingen
############################################################
function coerce_columns!(df::DataFrame, cols::Vector{String})  # DataFrame + Spaltenliste
    for c in cols                                          # über gewünschte Spalten iterieren
        if c in names(df)                                  # nur vorhandene Spalten bearbeiten
            df[!, c] = to_float.(df[!, c])                  # elementweise Umwandlung
        end
    end
    return df                                              # DataFrame zurückgeben
end

############################################################
# isapprox-Wrapper für Zahlen (Fixpunktvergleich)
############################################################
function isapprox_num(a::Real, b::Real)                    # zwei Real-Zahlen vergleichen
    return isapprox(a, b; rtol=REL_TOL, atol=ABS_TOL)       # mit Toleranzen
end

############################################################
# Prüft, ob eine Zeile dem Fixpunkt entspricht (außer varying_col)
############################################################
function row_matches_fix(row, fix::Dict{String,Float64}, varying_col::String)  # Zeile + Fixpunkt + variierende Spalte
    for (col, val) in fix                                  # über alle Fixparameter iterieren
        col == varying_col && continue                     # variierenden Parameter überspringen
        r = row[col]                                       # Wert in der Zeile holen
        if r === missing                                   # fehlt?
            return false                                   # dann passt Fixpunkt nicht
        end
        if !(r isa Real)                                   # nicht numerisch?
            return false                                   # dann passt Fixpunkt nicht
        end
        if !isapprox_num(r, val)                            # nicht nahe genug?
            return false                                   # dann passt Fixpunkt nicht
        end
    end
    return true                                            # alle Bedingungen erfüllt
end

############################################################
# Filter: Balance-Warnings entfernen (Betrag!)
# Regel: alles mit |balance| > threshold wird verworfen
############################################################
function filter_balance(df::DataFrame, bal_cols::Vector{String}, threshold::Float64)  # DataFrame + Balance-Spalten + Schwelle
    dfv = df                                               # Arbeitskopie
    for bcol in bal_cols                                   # über Balance-Spalten iterieren
        @assert bcol in names(dfv) "Spalte '$bcol' fehlt in der CSV."  # Existenz prüfen
        dfv = dfv[
            .!ismissing.(dfv[!, bcol]) .&&                 # missing raus
            (abs.(dfv[!, bcol]) .<= threshold),            # Betrag <= Schwelle behalten
            :
        ]
    end
    return dfv                                             # gefilterten DataFrame zurückgeben
end

############################################################
# Subset für 1D-Plot: Balance-Filter + Fixpunktfilter + Diagnose
############################################################
function subset_for_single_param(df::DataFrame, varying_col::String;
                                 fix::Dict{String,Float64},
                                 bal_cols::Vector{String},
                                 bal_threshold::Float64)

    n_start = nrow(df)                                     # Startzeilen
    df2 = filter_balance(df, bal_cols, bal_threshold)       # Balance-Warnings filtern
    n_after_balance = nrow(df2)                             # Zeilen nach Balance-Filter

    keep = Vector{Bool}(undef, nrow(df2))                   # Bool-Vektor für Fixpunktfilter
    for i in 1:nrow(df2)                                    # über alle Zeilen iterieren
        keep[i] = row_matches_fix(df2[i, :], fix, varying_col)  # Fixpunkt prüfen
    end
    df3 = df2[keep, :]                                      # Fixpunktsubset bilden
    n_after_fix = nrow(df3)                                 # Zeilen nach Fixpunktfilter

    println("--------------------------------------------------")  # Trennlinie
    println("Parameter: ", varying_col)                       # aktueller Parameter
    @printf("Start:                 %5d Zeilen\n", n_start)    # Startausgabe
    @printf("nach Balance-Filter:   %5d (-%d)\n",              # Balanceausgabe
            n_after_balance, n_start - n_after_balance)       # Differenz
    @printf("nach Fixpunkt-Filter:  %5d (-%d)\n",              # Fixpunktausgabe
            n_after_fix, n_after_balance - n_after_fix)       # Differenz

    return df3                                               # Subset zurückgeben
end

############################################################
# Plot: y (Annuität) vs x (ein Parameter)
############################################################
function plot_single_param(df::DataFrame, varying_col::String, ycol::String; outdir::String=OUTDIR)

    @assert varying_col in names(df) "Spalte '$varying_col' nicht gefunden."  # x muss existieren
    @assert ycol in names(df) "Spalte '$ycol' nicht gefunden."                # y muss existieren

    dfp = df[
        .!ismissing.(df[!, varying_col]) .&&                # x muss gültig sein
        .!ismissing.(df[!, ycol]),                          # y muss gültig sein
        :
    ]

    sort!(dfp, varying_col)                                 # nach x sortieren
    n = nrow(dfp)                                            # Anzahl Punkte

    if n == 0                                                # keine Punkte?
        @warn "Keine gültigen Punkte für $varying_col"        # Warnung ausgeben
        return nothing                                       # nichts plotten
    end

    ymin, iy = findmin(dfp[!, ycol])                         # Minimum von y finden
    xmin = dfp[iy, varying_col]                              # zugehöriges x

    @printf("Minimum: %.3e € bei %s = %.3e\n",                # Minimum ausgeben
            ymin, varying_col, xmin)                         # Format

    p = plot(                                                # Plotobjekt erstellen
        dfp[!, varying_col],                                 # x-Daten
        dfp[!, ycol],                                        # y-Daten
        seriestype = :line,                                  # Linienplot
        marker = :circle,                                    # Marker
        linewidth = 2,                                       # Linienbreite
        markersize = 4,                                      # Markergröße
        xlabel = varying_col,                                # x-Label
        ylabel = ycol,                                       # y-Label
        title  = "Gesamtannuität (no) vs. $(varying_col)\n(Fixpunkt, |Balance| ≤ $(BAL_THRESHOLD))",  # Titel
        legend = false,                                      # keine Legende
    )

    scatter!(p, [xmin], [ymin], marker=:star5, markersize=9, label=false)  # Minimum markieren

    isdir(outdir) || mkpath(outdir)                          # Ordner anlegen falls nötig
    fname = replace(varying_col, r"[^\w]+" => "_")           # Dateiname säubern
    savefig(p, joinpath(outdir, "annuity_no_vs_$(fname).png"))# Plot speichern

    return p                                                 # Plot zurückgeben
end

############################################################
# Globales Optimum: Minimum der Annuität über alle Runs
############################################################
function global_optimum(df::DataFrame, ycol::String)
    dfv = df                                                 # Arbeitskopie
    dfv = dfv[.!ismissing.(dfv[!, ycol]), :]                  # gültige y-Werte behalten
    ymin, idx = findmin(dfv[!, ycol])                         # Minimum + Index
    best = dfv[idx, :]                                        # Optimumzeile

    println("\n================ GLOBALOPTIMUM ================") # Überschrift
    @printf("Minimalwert %s: %.3e €\n\n", ycol, ymin)          # Minimalwert ausgeben

    for c in names(best)                                      # alle Spalten ausgeben
        println(rpad(c, 35), " = ", best[c])                  # Name + Wert
    end

    return best                                               # Optimum zurückgeben
end

############################################################
# Sensitivität 1: standardisierte lineare Sensitivität (β)
############################################################
function linear_sensitivity(df::DataFrame, xcols::Vector{String}, ycol::String)

    dfv = df[.!ismissing.(df[!, ycol]), :]                    # nur gültige y-Werte

    X = reduce(hcat, [dfv[!, c] for c in xcols])              # Matrix X (n x 4)
    y = dfv[!, ycol]                                          # Zielvariable y (n)

    Xs = (X .- mean(X, dims=1)) ./ std(X, dims=1)             # X standardisieren
    ys = (y .- mean(y)) ./ std(y)                             # y standardisieren

    β = Xs \ ys                                               # Least-Squares β

    println("\n=========== LINEARE SENSITIVITÄT (β) ===========")# Überschrift
    for (c, b) in zip(xcols, β)                               # Parameter + β
        @printf("%-30s β = %+6.3f\n", c, b)                   # Ausgabe
    end

    return β                                                  # β zurückgeben
end

############################################################
# Ranking-Funktion (für Spearman), mit Tie-Handling
# - gleiche Werte erhalten den Mittelwert ihrer Rangpositionen
# - Ränge starten bei 1
############################################################
function rank_average_ties(v::Vector{Float64})                # Input: Float64-Vektor
    n = length(v)                                             # Länge
    order = sortperm(v)                                       # Index-Reihenfolge sortierter Werte
    r = similar(v)                                            # Rang-Vektor anlegen
    i = 1                                                     # Startindex im sortierten Array
    while i <= n                                              # über alle sortierten Werte laufen
        j = i                                                 # j startet bei i
        while j < n && v[order[j]] == v[order[j + 1]]         # solange Gleichheit (Tie)
            j += 1                                            # Tie-Block erweitern
        end
        avg_rank = (i + j) / 2                                # mittlerer Rang für Tie-Block
        for k in i:j                                          # jedem Element im Tie-Block Rang zuweisen
            r[order[k]] = avg_rank                            # Rang zurück auf Originalindex
        end
        i = j + 1                                             # zum nächsten Block springen
    end
    return r                                                  # Ränge zurückgeben
end

############################################################
# Spearman-Korrelation = Pearson-Korrelation der Ränge
############################################################
function spearman_corr(x::Vector{Float64}, y::Vector{Float64})
    rx = rank_average_ties(x)                                 # Ränge von x
    ry = rank_average_ties(y)                                 # Ränge von y
    return cor(rx, ry)                                        # Pearson-Korrelation der Ränge
end

############################################################
# Sensitivität 2: Spearman-Rangkorrelation (ρ)
############################################################
function spearman_sensitivity(df::DataFrame, xcols::Vector{String}, ycol::String)

    println("\n=========== SPEARMAN-SENSITIVITÄT (ρ) ===========")# Überschrift

    for c in xcols                                            # über Parameter iterieren
        sub = df[
            .!ismissing.(df[!, c]) .&&                        # x gültig
            .!ismissing.(df[!, ycol]),                        # y gültig
            :
        ]
        x = Vector{Float64}(sub[!, c])                        # x als Float64-Vektor
        y = Vector{Float64}(sub[!, ycol])                     # y als Float64-Vektor
        ρ = spearman_corr(x, y)                               # Spearman berechnen
        @printf("%-30s ρ = %+6.3f\n", c, ρ)                   # Ausgabe
    end
end

############################################################
# Optional: Werte-Check (zeigt, ob Fixwerte im Grid vorkommen)
############################################################
function nearest_value_report(df::DataFrame, col::String, target::Float64; max_show::Int=12)
    vals = unique(skipmissing(df[!, col]))                    # eindeutige Werte
    vals = sort(collect(vals))                                # sortieren
    if isempty(vals)                                          # falls leer
        println("[$col] keine Werte (alles missing?)")        # Hinweis
        return                                                # abbrechen
    end
    idx = findmin(abs.(vals .- target))[2]                    # Index nächster Wert
    nearest = vals[idx]                                       # nächster Wert
    println("[$col] target = $(target), nearest = $(nearest), diff = $(nearest - target)") # Report
    println("    Beispiel-Werte: ",                           # Beispielwerte
            join(vals[1:min(end, max_show)], ", "),           # join
            length(vals) > max_show ? ", ..." : "")           # ggf. mehr anzeigen
end

############################################################
# MAIN
############################################################
function main()

    df = CSV.read(CSV_PATH, DataFrame)                        # CSV einlesen
    println("Gesamtanzahl Zeilen CSV: ", nrow(df))            # Zeilenanzahl ausgeben

    numeric_cols = vcat(XCOLS, [YCOL_NO], BAL_COLS)           # benötigte numerische Spalten
    coerce_columns!(df, numeric_cols)                         # Konvertierung anwenden

    println("\n--- Werte-Check: Fixpunkt vs. vorhandene Grid-Werte ---")  # Überschrift
    for (col, target) in FIX                                  # über Fixwerte iterieren
        nearest_value_report(df, col, target)                 # report ausgeben
    end

    df_valid = filter_balance(df, BAL_COLS, BAL_THRESHOLD)    # technisch gültige Runs (|Balance| ≤ 1)
    println("\nTechnisch gültige Runs (|Balance| ≤ 1): ", nrow(df_valid)) # Count ausgeben

    ########################################################
    # 1) 1D-Plots am Fixpunkt (Basisfall, no)
    ########################################################
    for varying in XCOLS                                      # jeden Parameter separat
        sub = subset_for_single_param(                        # Subset inkl. Diagnose
            df,                                                # Originaldaten
            varying;                                           # variierende Spalte
            fix=FIX,                                           # Fixpunkt
            bal_cols=BAL_COLS,                                 # Balance-Spalten
            bal_threshold=BAL_THRESHOLD                        # Schwelle (Betrag)
        )
        plot_single_param(sub, varying, YCOL_NO; outdir=OUTDIR)# Plot speichern
    end

    ########################################################
    # 2) Globales Optimum (Basisfall no)
    ########################################################
    best_no = global_optimum(df_valid, YCOL_NO)                # Optimum aus gültigen Runs

    ########################################################
    # 3) Sensitivitäten (Basisfall no)
    ########################################################
    β_no = linear_sensitivity(df_valid, XCOLS, YCOL_NO)        # lineare Sensitivität
    spearman_sensitivity(df_valid, XCOLS, YCOL_NO)             # Spearman-Sensitivität

    println("\nFertig. Plots liegen in: $(OUTDIR)/")           # Abschlussmeldung
end

Base.invokelatest(main)                                       # Start (robust gegen World-Age/Revise)


