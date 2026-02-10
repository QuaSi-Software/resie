using CSV                                                # CSV-Dateien einlesen
using DataFrames                                         # Tabellenverarbeitung
using Statistics                                         # mean/std/cor
using Printf                                             # formatierte Ausgaben
using Plots                                              # Plots

############################################################
# Konfiguration (hier anpassen)
############################################################
const CSV_PATH = "C:/Users/jenter/Documents/resie/output/parameterstudy/results_2975runs_260209_184340.csv"             # <- Pfad zur CSV anpassen
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
    "Hp_Power / W"             => 7.0e6,                  # W
    "Boiler_Power / W"         => 1.0e6,                  # W
    "BufferTank_Capacity / Wh" => 6.0e7,                  # Wh
    "Battery_Capacity / Wh"    => 0.0,                  # Wh
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

    println("--------------------------------------------------")  # Trennlinie
    println("Parameter: ", varying_col)                       # aktueller Parameter
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
        title  = "Gesamtannuität (no) vs. $(varying_col)\n(Fixpunkt zuerst, dann |Balance| ≤ $(BAL_THRESHOLD))", # Titel
        legend = false,                                      # keine Legende
    )

    scatter!(p, [xmin], [ymin], marker=:star5, markersize=9, label=false)  # Minimum markieren

    isdir(outdir) || mkpath(outdir)                          # Ordner anlegen falls nötig
    fname = replace(varying_col, r"[^\w]+" => "_")           # Dateiname säubern
    savefig(p, joinpath(outdir, "annuity_no_vs_$(fname).png"))# Plot speichern

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
        "annuity_no A_cap / €",
        "annuity_no A_misc / €",
        "annuity_no A_op / €",
        "annuity_no A_energy / €",
        "annuity_no A_rev_control / €",
        "annuity_no A_rev_feed / €",
        "annuity_no A_total / €",
    ]

    keep_cols = vcat(xcols, annuity_cols_no)
    bestN_small = bestN[:, keep_cols]

    # -------- Ausgabe --------
    println("\n================ TOP-$(nshow) GLOBALOPTIMA (no) ================")
    println("Hinweis: gefiltert auf |Balance| ≤ $(bal_threshold)\n")

    show(stdout, MIME"text/plain"(), bestN_small; allrows=true, allcols=true)
    println()

    return bestN_small
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
    col_cap    = "annuity_no A_cap / €"            # CAPEX-Anteil
    col_misc   = "annuity_no A_misc / €"           # Sonstige Kosten
    col_op     = "annuity_no A_op / €"             # Betriebskosten
    col_energy = "annuity_no A_energy / €"         # Energiekosten
    col_rev_c  = "annuity_no A_rev_control / €"    # Erlöse Regelenergie (wird abgezogen)
    col_rev_f  = "annuity_no A_rev_feed / €"       # Erlöse Einspeisung (wird abgezogen)
    col_total  = "annuity_no A_total / €"          # Gesamtannuität

    # ----------------------------
    # x-Achse
    # ----------------------------
    n = nrow(best10)                               # Anzahl Runs
    x = collect(1:n)                               # 1..n

    # ----------------------------
    # Daten (robust nach Float)
    # ----------------------------
    cap    = to_float.(best10[!, col_cap])         # A_cap
    misc   = to_float.(best10[!, col_misc])        # A_misc
    op     = to_float.(best10[!, col_op])          # A_op
    energy = to_float.(best10[!, col_energy])      # A_energy
    rev_c  = to_float.(best10[!, col_rev_c])       # A_rev_control (positiv in CSV)
    rev_f  = to_float.(best10[!, col_rev_f])       # A_rev_feed (positiv in CSV)
    total  = to_float.(best10[!, col_total])       # A_total

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
    p = plot(                                      # leeren Plot initialisieren
        xlabel = "Top-Run (1 = bestes A_total)",
        ylabel = "Annuitäten / €",
        title  = "Top-$n: Zerlegung der Annuitäten (no)\nKosten oben, Erlöse unten",
        legend = :topright
    )

    # Hilfsfunktion: eine Stack-Schicht als Balken zeichnen
    function stack_layer!(p, x, bottom, height; label::String)
        top = bottom .+ height                     # obere Kante
        bar!(p, x, top;                            # Balken bis "top"
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
    scatter!(p, x, total; label="A_total", markersize=5)

    # ----------------------------
    # Achsenticks
    # ----------------------------
    xticks!(p, x, string.(x))

    # ----------------------------
    # Speichern
    # ----------------------------
    isdir(outdir) || mkpath(outdir)
    savefig(p, joinpath(outdir, "top10_annuity_breakdown_no.png"))

    return p
end



############################################################
# Ranking-Funktion (für Spearman), mit Tie-Handling
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
# Sensitivität: Spearman (ρ)
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
        if length(unique(x)) < 2                               # keine Variation?
            @printf("%-30s ρ =   NaN  (keine Variation)\n", c) # Hinweis ausgeben
            continue                                           # weiter
        end
        ρ = spearman_corr(x, y)                               # Spearman berechnen
        @printf("%-30s ρ = %+6.3f\n", c, ρ)                   # Ausgabe
    end
end

############################################################
# MAIN
############################################################
function main()

    df = CSV.read(CSV_PATH, DataFrame)                        # CSV einlesen
    println("Gesamtanzahl Zeilen CSV: ", nrow(df))            # Zeilenanzahl ausgeben

    annuity_cols_no = [
    "annuity_no A_cap / €",
    "annuity_no A_cap_incentive / €",
    "annuity_no A_misc / €",
    "annuity_no A_op / €",
    "annuity_no A_energy / €",
    "annuity_no A_rev_control / €",
    "annuity_no A_rev_feed / €",
    "annuity_no A_total / €",
    "annuity_no A_total_incentive / €",
    ]
    numeric_cols = vcat(XCOLS, annuity_cols_no, BAL_COLS)      # benötigte numerische Spalten
    coerce_columns!(df, numeric_cols)                         # Konvertierung anwenden

    ########################################################
    # 1) 1D-Plots: Diagnose Fixpunkt zuerst, dann Balance
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
    # 2) Top-10 Globaloptima (nur Parameter + annuity_no total)
    ########################################################
    best10 = Base.invokelatest(topN_optima, df, XCOLS, YCOL_NO; N=10)
    Base.invokelatest(plot_top10_annuity_breakdown, best10; outdir=OUTDIR)

    ########################################################
    # 3) Sensitivität (Spearman) auf technisch gültigen Runs
    ########################################################
    df_valid = filter_balance(df, BAL_COLS, BAL_THRESHOLD)     # technisch gültige Runs
    spearman_sensitivity(df_valid, XCOLS, YCOL_NO)             # Spearman auf gültigen Runs

    println("\nFertig. Plots liegen in: $(OUTDIR)/")           # Abschlussmeldung
end

Base.invokelatest(main)                                       # Start (robust gegen World-Age/Revise)
