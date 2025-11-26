# Parameterstudie-Driver für ReSiE
#
# Zweck:
#  - Automatisiert vollständiges Durchlaufen aller Kombinationen von drei
#    Anlagenparametern (HeatPump Leistung, HeatingRod Leistung, BufferTank Kapazität).
#  - Erzeugt für jeden Run eine eigene input-JSON (auf Basis von resie/input_file.json),
#    setzt eindeutige absolute Output-Pfade in der JSON (csv + auxiliary info),
#    startet ReSiE als separaten Julia-Prozess und protokolliert Ergebnis in einer Logdatei.
#
# Annahmen / Einheiten:
#  - Alle Leistungswerte (Pth_HP, Pth_HR) werden in Watt (W) erwartet und vom Nutzer so eingegeben.
#  - BufferTank capacity wird in Wattstunden (Wh) erwartet und vom Nutzer so eingegeben.
#  - Das Skript setzt cfg["io_settings"]["csv_output_file"] auf einen absoluten Pfad,
#    sodass ReSiE direkt die erwartete CSV in das parameterstudy-Verzeichnis schreibt.
#
# Verwendung:
#  - In Projektverzeichnis ausführen:
#      julia --project=. parameterstudy.jl
#
# Hinweise:
#  - Bei vielen Kombinationen kann die Laufzeit erheblich sein. Testweise zuerst mit kleinen
#    Parameterlisten laufen lassen.
#  - Wenn ReSiE intern andere Feldnamen erwartet, passe die cfg["components"][...] Zuweisungen an.
#  - Dieses Skript startet für jeden Run einen separaten Julia-Prozess (sauberer Zustand).
#
using JSON
using Dates
using Resie
include("resie_logger.jl")
using .Resie_Logger
using OrderedCollections: OrderedDict
using UUIDs

# Basis-Input laden (Vorlage JSON)
# Die Datei resie/input_file.json dient als Template. Für jeden Run wird eine tiefe Kopie
# erzeugt und die gewünschten Parameter werden dort eingetragen.
# base_input = JSON.parsefile("./examples/AdJe_Bsp/input_file.json")
base_input = Resie.read_JSON("./examples/AdJe_Bsp/input_file.json")

# -----------------------
# Parameterdefinition
# -----------------------
# Die Benutzerwerte sind in W bzw. Wh angegeben (keine Umrechnung).
# Hier werden die Bereiche (Start, Stop, Step) festgelegt und in Vektoren gesammelt.
#
# Pth_HP = thermische Leistung der Wärmepumpe [W]
# Pth_HR = thermische Leistung des Heizstabs / Heizelements [W]
# Cap_Wh = Speicherkapazität des Puffers [Wh]
Pth_HP_lo = 10_000_000
Pth_HP_hi = 13_000_000
Pth_HP_step = 1_000_000
Pth_HP_vals = collect(range(Pth_HP_lo; step=Pth_HP_step, stop=Pth_HP_hi))

Pth_HR_lo = 500_000
Pth_HR_hi = 1_000_000
Pth_HR_step = 500_000
Pth_HR_vals = collect(range(Pth_HR_lo; step=Pth_HR_step, stop=Pth_HR_hi))

Cap_lo_Wh = 10_000
Cap_hi_Wh = 30_000
Cap_step_Wh = 10_000
Cap_vals_Wh = collect(range(Cap_lo_Wh; step=Cap_step_Wh, stop=Cap_hi_Wh))

# -----------------------
# Ausgabe-Ordner vorbereiten
# -----------------------
# outdir ist das zentrale Verzeichnis, in das alle erzeugten JSONs, CSVs und auxiliary-Files
# geschrieben werden. mkpath erstellt den Ordnerfalls erforderlich.
outdir = "./output/parameterstudy"
mkpath(outdir)

# -----------------------
# Hilfsfunktionen / Log
# -----------------------
# safe(x): Erzeugt eine dateinamensfreundliche Darstellung einer Zahl.
#            Ersetzt '.' durch 'p' (z.B. 1000000.0 -> "1000000p0").
#            Vermeidet Probleme mit Punkt im Dateinamen / Shell.
function safe(x)
    replace(string(x), "." => "p")
end

# Logdatei: enthält für jeden Run Zeitstempel, Parameter und Status (OK / ERROR + Info)
logfile = joinpath(outdir, "parameterstudy_log.csv")
if !isfile(logfile)
    open(logfile, "w") do f
        # Kopfzeile: zeit, index, gesamtanzahl, Parameterwerte, Status, zusätzliche Info (z.B. Fehlermeldung)
        println(f, "timestamp,runidx,total_runs,Pth_HP_W,Pth_HR_W,Cap_Wh,status,info")
    end
end

# -----------------------
# Kernfunktion: run_resie_variant
# -----------------------
# Zweck (detaillierter):
#  1) Erzeugt aus base_input eine lauf-spezifische Konfiguration (tiefe Kopie).
#  2) Trägt die drei Parameter in cfg["components"] ein (HeatPump, HeatingRod, BufferTank).
#  3) Legt eindeutige absolute Output-Pfade (csv + auxiliary_info) in cfg["io_settings"] fest,
#     damit ReSiE die Dateien direkt mit eindeutigen Namen erzeugt.
#  4) Schreibt die Lauf-JSON in outdir (zur Nachvollziehbarkeit gespeichert).
#  5) Startet ReSiE als externen Julia-Prozess: `julia --project=. src/run_resie.jl <input-json>`
#  6) Protokolliert Erfolg/Fehler in logfile.
#
# Eingaben:
#  - outdir: Verzeichnis für Outputs
#  - base_input: Template-JSON als Dict
#  - Pth_HP, Pth_HR: Leistungswerte in W
#  - Cap_Wh: Kapazität in Wh
#  - runidx, total_runs: numerische Indizes für Fortschrittsanzeige / Log
#
# Rückgabe:
#  - Tuple (status, info) — status: "OK" oder "ERROR", info: Fehlerbeschreibung falls vorhanden
function run_resie_variant(outdir::AbstractString, base_input::Union{Dict, OrderedDict},
                           Pth_HP::Number, Pth_HR::Number, Cap_Wh::Number,
                           runidx::Int, total_runs::Int, write_output::Bool=false)
    # --- Timestamp für Dateinamen und Log ---
    timestamp = Dates.format(now(), "yyyy-mm-dd_HH:MM:SS")

    # --- Tiefe Kopie der Vorlage ---
    # deepcopy stellt sicher, dass base_input unverändert bleibt und jeder Lauf mit sauberer Vorlage beginnt.
    cfg = deepcopy(base_input)

    # --- Parameter in die Konfiguration eintragen ---
    # ACHTUNG: Die hier verwendeten Schlüssel müssen mit der Struktur deiner input_file.json übereinstimmen.
    #          Falls die Component-Namen oder Feldnamen anders lauten, passe sie hier an.
    cfg["components"]["HeatPump"]["power_th"] = Pth_HP
    cfg["components"]["HeatingRod"]["power_th"] = Pth_HR
    cfg["components"]["BufferTank"]["capacity"] = Cap_Wh

    # --- Sicherstellen, dass io_settings vorhanden ist ---
    # Wenn die Vorlage kein io_settings-Objekt hat, legen wir eines an, damit wir die Outputpfade
    # gezielt setzen können.
    if !haskey(cfg, "io_settings") || cfg["io_settings"] === nothing
        cfg["io_settings"] = Dict{String,Any}()
    end

    # Schreibe output nur als Dateien, falls die Daten benötigt werden, da dies den Prozess 
    # deutlich verlangsamt
    if write_output
        # --- Eindeutige Output-Dateinamen direkt in cfg setzen ---
        # Wir verwenden absolute Pfade, damit ReSiE unabhängig vom aktuellen Arbeitsverzeichnis 
        # die Datei an den gewünschten Ort schreibt.
        csv_name = "out_HP$(safe(Pth_HP))W_HR$(safe(Pth_HR))W_Cap$(safe(Cap_Wh))Wh.csv"
        aux_name = "aux_HP$(safe(Pth_HP))W_HR$(safe(Pth_HR))W_Cap$(safe(Cap_Wh))Wh.md"
        csv_abs = abspath(joinpath(outdir, csv_name))
        aux_abs = abspath(joinpath(outdir, aux_name))
        cfg["io_settings"]["csv_output_file"] = csv_abs
        cfg["io_settings"]["auxiliary_info_file"] = aux_abs

    else
        cfg["io_settings"]["csv_output_keys"] = "nothing"
        cfg["io_settings"]["auxiliary_info"] = false
        cfg["io_settings"]["sankey_plot"] = "nothing"
    end

    # --- Lauf-spezifische Input-JSON speichern ---
    # Die Input-JSON bleibt im outdir als Beleg erhalten (nützlich für Reproduzierbarkeit).
    fname = joinpath(outdir, "input_HP$(safe(Pth_HP))W_HR$(safe(Pth_HR))W_Cap$(safe(Cap_Wh))Wh.json")
    open(fname, "w") do io
        JSON.print(io, cfg)
    end

    # --- ReSiE-Prozess starten ---
    # Wir starten ReSiE extern, damit jeder Lauf in einem frischen Julia-Prozess läuft.
    # Bei Fehlern wird eine Exception ausgelöst und abgefangen.
    println("[$runidx/$total_runs] Starte Lauf: HP=$(Pth_HP) W, HR=$(Pth_HR) W, Cap=$(Cap_Wh) Wh -> JSON: $(basename(fname))")
    status = "OK"
    info = ""
    sim_output = Dict()

    log_to_console = true
    log_to_file = false
    general_logfile_path = nothing
    balanceWarn_logfile_path = nothing
    min_log_level = Resie_Logger.Logging.Warn
    _, _ = Resie_Logger.start_logger(log_to_console,
                                        log_to_file,
                                        general_logfile_path,
                                        balanceWarn_logfile_path,
                                        min_log_level,
                                        fname)

    try
        _, sim_output = Resie.load_and_run(fname)
    catch e
        # Fehler auffangen, kurz protokollieren und im Log ablegen
        status = "ERROR"
        info = sprint(showerror, e)
        @warn "Run fehlgeschlagen: $info"
    end

    # try
    #     run(`julia --project=. src/resie-cli.jl run $fname --exit_after_run`)
    # catch e
    #     # Fehler auffangen, kurz protokollieren und im Log ablegen
    #     status = "ERROR"
    #     info = sprint(showerror, e)
    #     @warn "Run fehlgeschlagen: $info"
    # end

    # --- Logeintrag schreiben ---
    # Jede Zeile enthält Zeitstempel, Run-Index, Parametrierung und Status --> nützlich für Auswertung
    open(logfile, "a") do f
        # Info-Feld wird in Anführungszeichen gesetzt, damit Kommas in Fehlermeldungen nicht die CSV-Spalte trennen.
        println(f, "$(timestamp),$(runidx),$(total_runs),$(Pth_HP),$(Pth_HR),$(Cap_Wh),$(status),\"$(info)\"")
    end

    return sim_output
end

# -----------------------
# Hauptschleife: vollständiges kartesisches Produkt
# -----------------------
# Diese Schleife erzeugt und startet für jede mögliche Kombination (Pth_HP, Pth_HR, Cap_Wh) einen Run.
total_runs = length(Pth_HP_vals) * length(Pth_HR_vals) * length(Cap_vals_Wh)
println("Starte Parameterstudie: $total_runs Läufe (vollständiges kartesisches Produkt)")

runidx = 0
for (Pth_HP, Pth_HR, Cap_Wh) in Iterators.product(Pth_HP_vals, Pth_HR_vals, Cap_vals_Wh)
    global runidx += 1
    # Delegiere die eigentliche Erstellung der JSON, das Starten von ReSiE und das Logging an die Funktion
    sim_output = run_resie_variant(outdir, base_input, Pth_HP, Pth_HR, Cap_Wh, runidx, total_runs)
    
end