using CSV                                                # CSV-Dateien einlesen
using DataFrames                                         # Tabellenverarbeitung
using Statistics                                         # mean/std/cor
using Printf                                             # formatierte Ausgaben
using Plots                                              # Plots
import PlotlyJS

############################################################
# Konfiguration (hier anpassen)
############################################################
const CSV_PATH = "C:/Users/jenter/Documents/resie/output/parameterstudy/results_2640runs_260210_181232.csv"             # <- Pfad zur CSV anpassen
const OUTDIR   = "c:/Users/jenter/Documents/resie/output/parameterstudy/plots/"                                 # Plot-Ausgabeordner

############################################################
# Spaltennamen (wie in deiner CSV)
############################################################
const XCOLS = [                                          # Einflussgrößen (Anlagenparameter)
    "HeatPump_Power/W",                                       # Wärmepumpenleistung
    "ElectrodeBoiler_Power/W",                                   # ElectrodeBoiler/Heizstab-Leistung
    "BufferTank_Capacity/Wh",                           # Wärmespeichergröße
    "Battery_Capacity/Wh",                              # Batteriegröße
]

const YCOL_NO  = "annuity_no A_total/€"                # Zielgröße (Basisfall "no")
const BAL_COLS = ["balance_power", "balance_heat"]       # Balance-Warning-Spalten
const BAL_THRESHOLD = 1.0                                 # alle |balance| > 1 werden verworfen

############################################################
# Fixpunkt (für 1D-Slices)
############################################################
const FIX = Dict(                                        # Fixwerte für 1D-Plots
    "HeatPump_Power/W"             => 6.0e6,                  # W
    "ElectrodeBoiler_Power/W"         => 2.0e6,                  # W
    "BufferTank_Capacity/Wh" => 3.0e7,                  # Wh
    "Battery_Capacity/Wh"    => 0.0,                  # Wh
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
# Hilfsfunktion: Parameter-Label konvertieren (W/Wh → MW/MWh)
############################################################
function convert_param_label(col::String)                  # Spaltennamen konvertieren
    if contains(col, "/ W")                                # Watt-Einheit?
        return replace(col, "/ W" => "/ MW")               # in MW umwandeln
    elseif contains(col, "/ Wh")                           # Wh-Einheit?
        return replace(col, "/ Wh" => "/ MWh")             # in MWh umwandeln
    else
        return col                                         # sonst unverändert
    end
end

############################################################
# Hilfsfunktion: Annuitäts-Label konvertieren (€ → Mio. €)
############################################################
function convert_annuity_label(col::String)                # Spaltennamen für Annuitäten
    if contains(col, "/€")                               # Euro-Einheit?
        return replace(col, "/€" => "/Mio. €")         # in Mio. € umwandeln
    else
        return col                                         # sonst unverändert
    end
end

############################################################
# Hilfsfunktion: Wert-Konvertierung (W/Wh → MW/MWh, € → Mio. €)
############################################################
function convert_param_value(col::String, value::Real)     # Spaltennamen + Wert
    if contains(col, "/W") || contains(col, "/Wh")       # Watt oder Wh?
        return value / 1e6                                 # durch 1e6 teilen
    elseif contains(col, "/€")                           # Euro?
        return value / 1e6                                 # durch 1e6 teilen
    else
        return value                                       # sonst unverändert
    end
end

############################################################
# Hilfsfunktion: Wert formatieren (1 Nachkommastelle)
############################################################
function format_param_value(value::Real; decimals::Int=1)  # Wert mit Dezimalstellen
    rounded = round(value; digits=decimals)                # runden
    return string(rounded)                                 # als String zurückgeben
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
# Fixpunktfilter (OHNE Balance), damit du "grundsätzlich" siehst,
# wie viele Runs es am Fixpunkt gibt.
############################################################
function filter_fix_only(df::DataFrame, varying_col::String, fix::Dict{String,Float64})
    keep = Vector{Bool}(undef, nrow(df))                    # Bool-Vektor anlegen
    for i in 1:nrow(df)                                     # über alle Zeilen iterieren
        keep[i] = row_matches_fix(df[i, :], fix, varying_col) # Fixpunkt prüfen
    end
    return df[keep, :]                                      # gefilterten DataFrame zurückgeben
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
# Subset für 1D-Plot:
# - erst Fixpunktfilter (ohne Balance)
# - dann Balancefilter
# + Diagnoseausgabe für beide Schritte
############################################################
function subset_for_single_param(df::DataFrame, varying_col::String;
                                 fix::Dict{String,Float64},
                                 bal_cols::Vector{String},
                                 bal_threshold::Float64)

    n_start = nrow(df)                                     # Startzeilen

    df_fix = filter_fix_only(df, varying_col, fix)          # Fixpunktfilter zuerst
    n_after_fix = nrow(df_fix)                              # Anzahl am Fixpunkt (ohne Balance)

    df_fix_bal = filter_balance(df_fix, bal_cols, bal_threshold) # dann Balancefilter
    n_after_balance = nrow(df_fix_bal)                      # Anzahl am Fixpunkt (mit Balance)

    varying_col_label = convert_param_label(varying_col)    # Label konvertieren

    println("--------------------------------------------------")  # Trennlinie
    println("Parameter: ", varying_col_label)                # aktueller Parameter (konvertiert)
    @printf("Start:                      %5d Zeilen\n", n_start)           # Startausgabe
    @printf("nach Fixpunkt-Filter:       %5d (-%d)\n", n_after_fix, n_start - n_after_fix)  # Fix-Ausgabe
    @printf("nach Balance-Filter (auf Fix): %5d (-%d)\n", n_after_balance, n_after_fix - n_after_balance) # Balance-Ausgabe

    return df_fix_bal                                        # Subset für Plot zurückgeben
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

    xmin_conv = convert_param_value(varying_col, xmin)       # x-Wert konvertieren
    xmin_str = format_param_value(xmin_conv)                 # formatieren
    xlabel_conv = convert_param_label(varying_col)           # Label konvertieren
    ylabel_conv = convert_annuity_label(ycol)                # Annuitäts-Label konvertieren
    ymin_conv = convert_param_value(ycol, ymin)              # y-Minimum konvertieren
    ymin_str = format_param_value(ymin_conv)                 # formatieren

    @printf("Minimum: %s Mio. € bei %s = %s\n",              # Minimum ausgeben
            ymin_str, xlabel_conv, xmin_str)                 # Format

    # Konvertiere x-Daten (alle x-Werte)
    x_data = Vector{Float64}(dfp[!, varying_col])
    if contains(varying_col, "/ W") || contains(varying_col, "/ Wh")
        x_data = x_data ./ 1e6                               # In MW/MWh umwandeln
    end

    # Konvertiere y-Daten (Annuität in Mio. €)
    y_data = Vector{Float64}(dfp[!, ycol]) ./ 1e6            # In Mio. € umwandeln

    p = Plots.plot(                                          # Plotobjekt erstellen
        x_data,                                              # x-Daten konvertiert
        y_data,                                              # y-Daten konvertiert
        seriestype = :line,                                  # Linienplot
        marker = :circle,                                    # Marker
        linewidth = 2,                                       # Linienbreite
        markersize = 4,                                      # Markergröße
        xlabel = xlabel_conv,                                # x-Label konvertiert
        ylabel = ylabel_conv,                                # y-Label konvertiert
        title  = "Gesamtannuität (no) vs. $(xlabel_conv)\n(Fixpunkt zuerst, dann |Balance| ≤ $(BAL_THRESHOLD))", # Titel
        legend = false,                                      # keine Legende
    )

    Plots.scatter!(p, [xmin_conv], [ymin_conv], marker=:star5, markersize=9, label=false)  # Minimum markieren (konvertiert!)

    isdir(outdir) || mkpath(outdir)                          # Ordner anlegen falls nötig
    fname = replace(varying_col, r"[^\w]+" => "_")           # Dateiname säubern
    Plots.savefig(p, joinpath(outdir, "annuity_no_vs_$(fname).png"))# Plot speichern

    return p                                                 # Plot zurückgeben
end


############################################################
# Top-N Globaloptima (nur Szenario "no")
# - Filtert technisch gültige Runs (|balance| ≤ threshold)
# - Sortiert nach annuity_no A_total
# - Gibt Parameter + ALLE no-Annuitäten aus
############################################################
function topN_optima(df::DataFrame,
                     xcols::Vector{String},
                     ycol::String;
                     bal_cols::Vector{String}=BAL_COLS,
                     bal_threshold::Float64=BAL_THRESHOLD,
                     N::Int=10)

    # -------- technisch gültige Runs --------
    dfv = filter_balance(df, bal_cols, bal_threshold)

    # -------- nur gültige Annuitäten --------
    dfv = dfv[.!ismissing.(dfv[!, ycol]), :]

    # -------- nach Gesamtannuität sortieren --------
    sort!(dfv, ycol)

    # -------- Top-N auswählen --------
    nshow = min(N, nrow(dfv))
    bestN = dfv[1:nshow, :]

    # -------- Spaltenauswahl --------
    annuity_cols_no = [
        "annuity_no A_cap/€",
        "annuity_no A_misc/€",
        "annuity_no A_op/€",
        "annuity_no A_energy/€",
        "annuity_no A_rev_control/€",
        "annuity_no A_rev_feed/€",
        "annuity_no A_total/€",
    ]

    # -------- Konvertiere Leistungs-/Kapazitätswerte in bestN (für plot_top10_annuity_breakdown) --------
    for col in xcols
        if contains(col, "/ W") || contains(col, "/ Wh")
            bestN[!, col] = bestN[!, col] ./ 1e6  # In MW/MWh umwandeln
        end
    end

    # -------- Konvertiere Annuitätswerte in bestN (für plot_top10_annuity_breakdown) --------
    for col in annuity_cols_no
        bestN[!, col] = bestN[!, col] ./ 1e6     # In Mio. € umwandeln
    end

    # -------- Spaltenauswahl für Anzeige --------
    keep_cols = vcat(xcols, annuity_cols_no)
    bestN_small = bestN[:, keep_cols]

    # -------- Konvertiere Spaltenüberschriften für Anzeige --------
    rename_pairs = Pair{String, String}[]
    for col in names(bestN_small)
        if contains(col, "/W") || contains(col, "/Wh")
            push!(rename_pairs, col => convert_param_label(col))
        elseif contains(col, "/€")
            push!(rename_pairs, col => convert_annuity_label(col))
        end
    end
    if !isempty(rename_pairs)
        rename!(bestN_small, rename_pairs)
    end

    # -------- Ausgabe --------
    println("\n================ TOP-$(nshow) GLOBALOPTIMA (no) ================")
    println("Hinweis: gefiltert auf |Balance| ≤ $(bal_threshold), Leistungen in MW, Kapazitäten in MWh, Annuitäten in Mio. €\n")

    show(stdout, MIME"text/plain"(), bestN_small; allrows=true, allcols=true)
    println()

    return bestN  # Rückgabe der kompletten DataFrame mit konvertierten Werten aber Original-Labels
end


############################################################
# Plot: Top-10 Annuitäts-Zerlegung (gestapelte Balken) – GR-sicher
# - Kosten positiv (oben)
# - Erlöse negativ (unten)
# - A_total als Punkt
############################################################
function plot_top10_annuity_breakdown(best10::DataFrame; outdir::String=OUTDIR)

    # ----------------------------
    # Spalten (No-Szenario)
    # ----------------------------
    col_cap    = "annuity_no A_cap/€"            # CAPEX-Anteil
    col_misc   = "annuity_no A_misc/€"           # Sonstige Kosten
    col_op     = "annuity_no A_op/€"             # Betriebskosten
    col_energy = "annuity_no A_energy/€"         # Energiekosten
    col_rev_c  = "annuity_no A_rev_control/€"    # Erlöse Regelenergie (wird abgezogen)
    col_rev_f  = "annuity_no A_rev_feed/€"       # Erlöse Einspeisung (wird abgezogen)
    col_total  = "annuity_no A_total/€"          # Gesamtannuität

    # ----------------------------
    # x-Achse
    # ----------------------------
    n = nrow(best10)                               # Anzahl Runs
    x = collect(1:n)                               # 1..n

    # ----------------------------
    # Daten (robust nach Float)
    # Werte sind bereits in Mio. € konvertiert von topN_optima()
    # ----------------------------
    cap    = to_float.(best10[!, col_cap])         # A_cap (bereits in Mio. €)
    misc   = to_float.(best10[!, col_misc])        # A_misc (bereits in Mio. €)
    op     = to_float.(best10[!, col_op])          # A_op (bereits in Mio. €)
    energy = to_float.(best10[!, col_energy])      # A_energy (bereits in Mio. €)
    rev_c  = to_float.(best10[!, col_rev_c])       # A_rev_control (bereits in Mio. €)
    rev_f  = to_float.(best10[!, col_rev_f])       # A_rev_feed (bereits in Mio. €)
    total  = to_float.(best10[!, col_total])       # A_total (bereits in Mio. €)

    # ----------------------------
    # Missing-Check (sonst Plot-Fehler)
    # ----------------------------
    for (name, v) in [("cap",cap), ("misc",misc), ("op",op), ("energy",energy), ("rev_c",rev_c), ("rev_f",rev_f), ("total",total)]
        @assert all(.!ismissing.(v)) "Missing in $(name) – prüfe CSV/Parsing"
    end

    # ----------------------------
    # Positive Stacks: Kosten nach oben
    # Wir bauen den "Bottom" kumulativ selbst
    # ----------------------------
    bottom_pos = zeros(Float64, n)                 # Start bei 0
    p = Plots.plot(                                # leeren Plot initialisieren
        xlabel = "Top-Run (1 = bestes A_total)",
        ylabel = "Annuitäten/Mio. €",
        title  = "Top-$n: Zerlegung der Annuitäten (no)\nKosten oben, Erlöse unten",
        legend = :topright
    )

    # Hilfsfunktion: eine Stack-Schicht als Balken zeichnen
    function stack_layer!(p, x, bottom, height; label::String)
        top = bottom .+ height                     # obere Kante
        Plots.bar!(p, x, top;                      # Balken bis "top"
            fillrange = bottom,                    # Start bei "bottom" (macht stacking!)
            label = label)
        return bottom .+ height                    # neuen Bottom zurückgeben
    end

    # Kosten stapeln
    bottom_pos = stack_layer!(p, x, bottom_pos, cap;    label="A_cap")
    bottom_pos = stack_layer!(p, x, bottom_pos, misc;   label="A_misc")
    bottom_pos = stack_layer!(p, x, bottom_pos, op;     label="A_op")
    bottom_pos = stack_layer!(p, x, bottom_pos, energy; label="A_energy")

    # ----------------------------
    # Negative Stacks: Erlöse nach unten
    # (wir stapeln negative Höhen)
    # ----------------------------
    bottom_neg = zeros(Float64, n)                 # Start bei 0
    bottom_neg = stack_layer!(p, x, bottom_neg, -rev_c; label="-A_rev_control")
    bottom_neg = stack_layer!(p, x, bottom_neg, -rev_f; label="-A_rev_feed")

    # ----------------------------
    # A_total als Marker
    # ----------------------------
    Plots.scatter!(p, x, total; label="A_total", markersize=5)

    # ----------------------------
    # Achsenticks
    # ----------------------------
    Plots.xticks!(p, x, string.(x))

    # ----------------------------
    # Speichern
    # ----------------------------
    isdir(outdir) || mkpath(outdir)
    Plots.savefig(p, joinpath(outdir, "top10_annuity_breakdown_no.png"))

    return p
end



############################################################
# 3D-Plot: HeatPump_Power vs. ElectrodeBoiler_Power vs. BufferTank_Capacity
# mit farblicher Codierung der Annuität (grün=niedrig, rot=hoch)
############################################################
function plot_3d_parameter_space(df::DataFrame; outdir::String=OUTDIR)

    # ----------------------------
    # Spaltennamen
    # ----------------------------
    col_HeatPump     = "HeatPump_Power/W"
    col_ElectrodeBoiler = "ElectrodeBoiler_Power/W"
    col_buffer = "BufferTank_Capacity/Wh"
    col_annuity = "annuity_no A_total/€"

    # ----------------------------
    # Balance-Filter anwenden
    # ----------------------------
    dfv = filter_balance(df, BAL_COLS, BAL_THRESHOLD)
    println("\n---------- 3D-Parameter-Space ----------")
    println("Gesamtanzahl Runs: $(nrow(df))")
    println("Nach Balance-Filter: $(nrow(dfv))")

    # ----------------------------
    # Nur gültige Zeilen selektieren
    # ----------------------------
    dfv = dfv[
        .!ismissing.(dfv[!, col_HeatPump]) .&&
        .!ismissing.(dfv[!, col_ElectrodeBoiler]) .&&
        .!ismissing.(dfv[!, col_buffer]) .&&
        .!ismissing.(dfv[!, col_annuity]),
        :
    ]

    n = nrow(dfv)
    if n == 0
        @warn "Keine gültigen Datenpunkte für 3D-Plot"
        return nothing
    end

    println("Mit allen 4 Parametern gültig: $n")

    # ----------------------------
    # Daten extrahieren und konvertieren
    # ----------------------------
    HeatPump      = Vector{Float64}(dfv[!, col_HeatPump]) ./ 1e6        # in MW konvertieren
    ElectrodeBoiler  = Vector{Float64}(dfv[!, col_ElectrodeBoiler]) ./ 1e6    # in MW konvertieren
    buffer  = Vector{Float64}(dfv[!, col_buffer]) ./ 1e6    # in MWh konvertieren
    annuity = Vector{Float64}(dfv[!, col_annuity]) ./ 1e6    # in Mio. € konvertieren

    # ----------------------------
    # Labelnamen konvertieren
    # ----------------------------
    col_HeatPump_label     = convert_param_label(col_HeatPump)
    col_ElectrodeBoiler_label = convert_param_label(col_ElectrodeBoiler)
    col_buffer_label = convert_param_label(col_buffer)
    col_annuity_label = convert_annuity_label(col_annuity)

    # ----------------------------
    # Annuität-Statistik
    # ----------------------------
    annuity_min = minimum(annuity)
    annuity_max = maximum(annuity)
    annuity_mean = mean(annuity)
    @printf("Annuität min: %.2f Mio. €\n", annuity_min)
    @printf("Annuität max: %.2f Mio. €\n", annuity_max)
    @printf("Annuität mean: %.2f Mio. €\n", annuity_mean)

     # ----------------------------
    # Interaktiver 3D-Scatterplot (PlotlyJS)
    # ----------------------------
    trace = PlotlyJS.scatter3d(
        x = HeatPump,
        y = ElectrodeBoiler,
        z = buffer,
        mode = "markers",
        marker = PlotlyJS.attr(
            size = 4,
            color = annuity,                        # Farbcodierung nach Annuität
            colorscale = "Viridis",
            opacity = 0.9,
            colorbar = PlotlyJS.attr(title = "Gesamtannuität [Mio. €]")
        ),
        name = col_annuity_label
    )

    layout = PlotlyJS.Layout(
        title = "3D-Parameterraum: Parameter vs. Gesamtannuität",
        scene = PlotlyJS.attr(
            xaxis = PlotlyJS.attr(title = col_HeatPump_label),
            yaxis = PlotlyJS.attr(title = col_ElectrodeBoiler_label),
            zaxis = PlotlyJS.attr(title = col_buffer_label)
        ),
        margin = PlotlyJS.attr(l = 0, r = 0, b = 0, t = 60)
    )

    p = PlotlyJS.Plot(trace, layout)

    # ----------------------------
    # Speichern als interaktive HTML-Datei
    # ----------------------------
    isdir(outdir) || mkpath(outdir)
    html_path = joinpath(outdir, "3d_parameter_space_interactive.html")
    PlotlyJS.savefig(p, html_path)
    println("Interaktiver Plot gespeichert: $html_path")

    return p
end

############################################################
# Factor-Screening Sensitivität (SRC = standardisierte Regression)
# - z-standardisiert x und y
# - lineares Modell: y = b0 + sum(bi * xi)
# - |bi| zeigt relativen Einfluss
############################################################
function factor_screening_sensitivity(df::DataFrame, xcols::Vector{String}, ycol::String)

    println("\n=========== FACTOR SCREENING (SRC) ===========")

    # Nur vollständige Zeilen
    mask = trues(nrow(df))
    for c in xcols
        mask .&= .!ismissing.(df[!, c])
    end
    mask .&= .!ismissing.(df[!, ycol])
    sub = df[mask, :]

    n = nrow(sub)
    p = length(xcols)
    if n == 0
        @warn "Keine gültigen Zeilen für Factor Screening"
        return nothing
    end

    # Designmatrix und Ziel
    X = Matrix{Float64}(undef, n, p)
    for (j, c) in enumerate(xcols)
        X[:, j] = Vector{Float64}(sub[!, c])
    end
    y = Vector{Float64}(sub[!, ycol])

    # Standardisieren (z-Score)
    Xz = similar(X)
    valid_cols = trues(p)
    for j in 1:p
        xj = X[:, j]
        sj = std(xj)
        if sj == 0.0
            valid_cols[j] = false
            Xz[:, j] .= 0.0
        else
            Xz[:, j] = (xj .- mean(xj)) ./ sj
        end
    end
    sy = std(y)
    if sy == 0.0
        @warn "Zielvariable hat keine Variation"
        return nothing
    end
    yz = (y .- mean(y)) ./ sy

    # Lineare Regression mit Achsenabschnitt
    Xreg = hcat(ones(n), Xz)
    beta = Xreg \ yz
    coeffs = beta[2:end]

    # Ausgabe
    for (j, c) in enumerate(xcols)
        c_label = convert_param_label(c)
        if !valid_cols[j]
            @printf("%-30s SRC =   NaN  (keine Variation)\n", c_label)
        else
            @printf("%-30s SRC = %+6.3f\n", c_label, coeffs[j])
        end
    end
    return coeffs
end


############################################################
# MAIN
############################################################
function main()

    df = CSV.read(CSV_PATH, DataFrame)                        # CSV einlesen
    println("Gesamtanzahl Zeilen CSV: ", nrow(df))            # Zeilenanzahl ausgeben

    annuity_cols_no = [
    "annuity_no A_cap/€",
    "annuity_no A_cap_incentive/€",
    "annuity_no A_misc/€",
    "annuity_no A_op/€",
    "annuity_no A_energy/€",
    "annuity_no A_rev_control/€",
    "annuity_no A_rev_feed/€",
    "annuity_no A_total/€",
    "annuity_no A_total_incentive/€",
    ]
    numeric_cols = vcat(XCOLS, annuity_cols_no, BAL_COLS)      # benötigte numerische Spalten
    coerce_columns!(df, numeric_cols)                         # Konvertierung anwenden

    ########################################################
    # 1) 1D-Plots: Diagnose Fixpunkt zuerst, dann Balance
    ########################################################
    for varying in XCOLS                                      # jeden Parameter separat
        sub = Base.invokelatest(                              # Subset inkl. Diagnose
            getfield(@__MODULE__, :subset_for_single_param),
            df,                                                # Originaldaten
            varying;                                           # variierende Spalte
            fix=FIX,                                           # Fixpunkt
            bal_cols=BAL_COLS,                                 # Balance-Spalten
            bal_threshold=BAL_THRESHOLD                        # Schwelle (Betrag)
        )
        Base.invokelatest(
            getfield(@__MODULE__, :plot_single_param),
            sub, varying, YCOL_NO; outdir=OUTDIR
        )                                                      # Plot speichern
    end

    ########################################################
    # 2) 3D-Plot: HeatPump_Power vs. ElectrodeBoiler_Power vs. BufferTank
    #    mit Farbcodierung der Annuität
    ########################################################
    Base.invokelatest(getfield(@__MODULE__, :plot_3d_parameter_space), df; outdir=OUTDIR)

    ########################################################
    # 3) Top-10 Globaloptima (nur Parameter + annuity_no total)
    ########################################################
    best10 = Base.invokelatest(getfield(@__MODULE__, :topN_optima), df, XCOLS, YCOL_NO; N=10)
    Base.invokelatest(getfield(@__MODULE__, :plot_top10_annuity_breakdown), best10; outdir=OUTDIR)

    # 4) Sensitivität (Factor Screening) auf technisch gültigen Runs
    ########################################################
    df_valid = Base.invokelatest(getfield(@__MODULE__, :filter_balance), df, BAL_COLS, BAL_THRESHOLD)
    Base.invokelatest(getfield(@__MODULE__, :factor_screening_sensitivity), df_valid, XCOLS, YCOL_NO)

    println("\nFertig. Plots liegen in: $(OUTDIR)/")           # Abschlussmeldung
end

Base.invokelatest(getfield(@__MODULE__, :main))              # Start (robust gegen World-Age/Revise)
