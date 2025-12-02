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
#  - BufferTank capacity wird in Wattstunden (Wh)  erwartet und vom Nutzer so eingegeben.
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
include("vdi2067.jl")
using .VDI2067

# Basis-Input laden (Vorlage JSON)
# Default: kann per Kommandozeilen-Argument überschrieben werden, z.B.:
#    julia --project=. parameterstudy.jl ./examples/AdJe_Bsp/input_file_V1a.json
base_input_path = length(ARGS) > 0 ? ARGS[1] : "./inputfiles/inputfile_V3a.json"
base_input = Resie.read_JSON(base_input_path)

# Variante aus Dateiname (oder JSON-Feld) erkennen
function detect_variant_from_path(path::AbstractString, cfg)  # accept any Dict-like (OrderedDict etc.)
    # try to find _Vn or Vn pattern in filename
    fname_uc = uppercase(basename(path))
    m = match(r"_?V\d+[A-Z]?", fname_uc)
    if m !== nothing
        v = replace(m.match, "_" => "")
        return v
    end
    # fallback: try JSON field "variant"
    if haskey(cfg, "variant") && cfg["variant"] !== nothing
        return uppercase(string(cfg["variant"]))
    end
    return "UNKNOWN"
end

variant = detect_variant_from_path(base_input_path, base_input)
println("Gefundene Variante aus Dateiname / JSON: $variant")

# Wenn V0 dann Abbruch (nicht parametrisiert)
if startswith(variant, "V0")
    println("Variante V0 (abbruch): wird für diese Parameterstudie übersprungen bzw. nicht parametrisiert.")
    exit(0)
end

# -----------------------
# Komponenten- / Variantenlogik
# -----------------------
# Der Buffer wird in allen Varianten parametrisiert, außer in V0 und allen V5-Varianten (V5*, v5a..v5d).
# Die Batterie wird nur parametrisiert, wenn Variantengruppe V5 (startswith("V5")) aktiv ist.
# HP / HR / Preise werden nur parametrisiert, wenn die Komponente in der Vorlage existiert.
has_hp = haskey(base_input["components"], "HeatPump")
has_hr = haskey(base_input["components"], "HeatingRod")
has_buffer = haskey(base_input["components"], "BufferTank")
has_battery = haskey(base_input["components"], "Battery")

# Variante-spezifische Parametrisierung:
buffer_param = has_buffer && !startswith(variant, "V5")   # Buffer nur wenn vorhanden und nicht V5*
battery_param = has_battery && startswith(variant, "V5") # Battery nur param. bei V5*
hp_param = has_hp  # Wenn im Template vorhanden -> parametrierbar (weitere Variantenausnahmen können ergänzt werden)
hr_param = has_hr

# Preis-Parameter: nur parametrisieren, wenn das control_module ein limit_price bzw. reserve_price / reserve-Feld enthält
function detect_price_params(cfg)
    use_stock = false
    use_reserve = false
    bp = get(cfg, "components", nothing)
    if bp !== nothing
        busp = get(bp, "BUS_Power", nothing)
        if busp !== nothing
            cms = get(busp, "control_modules", nothing)
            if cms !== nothing
                for cm in cms
                    if isa(cm, Dict)
                        if haskey(cm, "price_profile_path") || haskey(cm, "limit_price")
                            use_stock = true
                        end
                        if haskey(cm, "reserve_price") || occursin("reserve", lowercase(string(get(cm, "name", ""))))
                            use_reserve = true
                        end
                    end
                end
            end
        end
    end
    return use_stock, use_reserve
end

use_p_stock_param, use_p_reserve_param = detect_price_params(base_input)

println("Eingesetzte Komponenten/Parametrierung in Variant '$variant':")
println("  HeatPump vorhanden: $has_hp ; parametrierbar: $hp_param")
println("  HeatingRod vorhanden: $has_hr ; parametrierbar: $hr_param")
println("  Buffer vorhanden: $has_buffer ; parametrierbar: $buffer_param")
println("  Battery vorhanden: $has_battery ; parametrierbar: $battery_param")
println("  Marktpreis parametrierbar: $use_p_stock_param")
println("  Regelreserve parametrierbar: $use_p_reserve_param")

# -----------------------
# Hilfsfunktionen: defaults aus Vorlage
# -----------------------
function get_existing_component_value(cfg::AbstractDict, comp::AbstractString, key::AbstractString)
    if haskey(cfg["components"], comp) && isa(cfg["components"][comp], Dict)
        return get(cfg["components"][comp], key, nothing)
    end
    return nothing
end
hp_default = get_existing_component_value(base_input, "HeatPump", "power_th")
hr_default = get_existing_component_value(base_input, "HeatingRod", "power_th")
buf_default = get_existing_component_value(base_input, "BufferTank", "capacity")
batt_default = get_existing_component_value(base_input, "Battery", "capacity")
# Preis-Parameter (Grenzpreise für wirtschaftliche Steuerung)
p_stock_lo = 80.0
p_stock_hi = 80.0
p_stock_step = 10.0
p_stock_marginal_vals = collect(range(p_stock_lo; step=p_stock_step, stop=p_stock_hi))

p_reserve_lo = 10.0
p_reserve_hi = 10.0
p_reserve_step = 1.0
p_reserve_marginal_vals = collect(range(p_reserve_lo; step=p_reserve_step, stop=p_reserve_hi))

# helper functions for defaults to avoid soft-scope warnings
function get_default_limit_price(cfg)
    local p = 80.0
    if haskey(cfg["components"], "BUS_Power") && haskey(cfg["components"]["BUS_Power"], "control_modules")
        for cm in cfg["components"]["BUS_Power"]["control_modules"]
            if isa(cm, Dict) && haskey(cm, "limit_price")
                return cm["limit_price"]
            end
        end
    end
    return p
end
function get_default_reserve_price(cfg)
    local r = 0.0
    if haskey(cfg["components"], "BUS_Power") && haskey(cfg["components"]["BUS_Power"], "control_modules")
        for cm in cfg["components"]["BUS_Power"]["control_modules"]
            if isa(cm, Dict) && haskey(cm, "reserve_price")
                return cm["reserve_price"]
            end
        end
    end
    return r
end

stock_default = get_default_limit_price(base_input)
reserve_default = get_default_reserve_price(base_input)

# -----------------------
# Parameterdefinition
# -----------------------
# Die Benutzerwerte bleiben intern in W bzw. Wh (keine Änderung der internen Einheiten).
# Für die Dateinamen und Ausgaben werden die Werte verkürzt (W -> MW, Wh -> MWh).
#
# Pth_HP = thermische Leistung der Wärmepumpe [W]
# Pth_HR = thermische Leistung des Heizstabs / Heizelements [W]
# Cap_Wh = Speicherkapazität des Puffers [Wh]
# BattCap_Wh = Speicherkapazität der Batterie [Wh]
# p_stock_marginal = Grenzpreis Börsenstrom [EUR/MWh]
# p_reserve_marginal = Grenzpreis Regelreservevergütung [EUR/MW]
Pth_HP_lo = 10_000_000
Pth_HP_hi = 10_000_000
Pth_HP_step = 1_000_000
Pth_HP_vals = collect(range(Pth_HP_lo; step=Pth_HP_step, stop=Pth_HP_hi))

Pth_HR_lo = 500_000
Pth_HR_hi = 500_000
Pth_HR_step = 500_000
Pth_HR_vals = collect(range(Pth_HR_lo; step=Pth_HR_step, stop=Pth_HR_hi))

Cap_lo_Wh = 10_000
Cap_hi_Wh = 10_000
Cap_step_Wh = 10_000
Cap_vals_Wh = collect(range(Cap_lo_Wh; step=Cap_step_Wh, stop=Cap_hi_Wh))

BattCap_lo_Wh = 200e3
BattCap_hi_Wh = 200e3
BattCap_step_Wh = 50e3
BattCap_vals_Wh = collect(range(BattCap_lo_Wh; step=BattCap_step_Wh, stop=BattCap_hi_Wh))

# Preis-Parameter (Grenzpreise für wirtschaftliche Steuerung)
p_stock_lo = 80.0
p_stock_hi = 80.0
p_stock_step = 10.0
p_stock_marginal_vals = collect(range(p_stock_lo; step=p_stock_step, stop=p_stock_hi))

p_reserve_lo = 10.0
p_reserve_hi = 10.0
p_reserve_step = 1.0
p_reserve_marginal_vals = collect(range(p_reserve_lo; step=p_reserve_step, stop=p_reserve_hi))

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

# Anzeige-Hilfsfunktionen (nur zur Beschriftung / Dateinamenkonvention)
# Intern werden die Werte unverändert in W / Wh weiterverwendet; diese Funktionen
# erzeugen die kompakte Darstellung (MW / MWh) für Dateinamen etc.
function display_MW(W::Number)
    round(W / 1e6; digits=3)
end
function display_MWh(Wh::Number)
    round(Wh / 1e6; digits=3)
end

# Logdatei: enthält für jeden Run Zeitstempel, Parameter und Status (OK / ERROR + Info)
logfile = joinpath(outdir, "parameterstudy_log.csv")
if !isfile(logfile)
    open(logfile, "w") do f
        # Kopfzeile erweitert: variant
        println(f, "timestamp,runidx,total_runs,variant,Pth_HP_W,Pth_HR_W,Cap_Wh,BattCap_Wh,p_stock_marginal,p_reserve_marginal,status,info")
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
#  - variant: Variantenbezeichnung (z.B. V1a)
#  - buffer_param: Bool, ob Buffer-Parameter gesetzt werden (abweichend von Vorlage)
#  - battery_param: Bool, ob Batterie-Parameter gesetzt werden (abweichend von Vorlage)
#
# Rückgabe:
#  - Tuple (status, info) — status: "OK" oder "ERROR", info: Fehlerbeschreibung falls vorhanden
function run_resie_variant(outdir::AbstractString, base_input::Union{Dict, OrderedDict},
                           Pth_HP::Number, Pth_HR::Number, Cap_Wh::Number, BattCap_Wh::Number,
                           p_stock_marginal::Number, p_reserve_marginal::Number,
                           runidx::Int, total_runs::Int, variant::String,
                           buffer_param::Bool=true, battery_param::Bool=false; write_output::Bool=false)
    # --- Timestamp für Dateinamen und Log ---
    timestamp = Dates.format(now(), "yyyy-mm-dd_HH:MM:SS")

    cfg = deepcopy(base_input)

    # nur vorhandene Komponenten anpassen, und nur falls sie für die Variantengruppe parametrisiert werden
    if haskey(cfg["components"], "HeatPump") && haskey(cfg["components"]["HeatPump"], "power_th")
        cfg["components"]["HeatPump"]["power_th"] = Pth_HP
    end

    if haskey(cfg["components"], "HeatingRod") && haskey(cfg["components"]["HeatingRod"], "power_th")
        cfg["components"]["HeatingRod"]["power_th"] = Pth_HR
    end

    # Buffer nur setzen, wenn buffer_param True (Buffer wird in V5 nicht parametriert)
    if buffer_param && haskey(cfg["components"], "BufferTank") && haskey(cfg["components"]["BufferTank"], "capacity")
        cfg["components"]["BufferTank"]["capacity"] = Cap_Wh
    end

    # Battery nur für Variante V5 parametrisieren (battery_param True)
    if battery_param
        if haskey(cfg["components"], "Battery")
            cfg["components"]["Battery"]["capacity"] = BattCap_Wh
        else
            cfg["components"]["Battery"] = Dict("type" => "Battery",
                                               "output_refs" => ["BUS_Power"],
                                               "medium" => "m_power",
                                               "capacity" => BattCap_Wh,
                                               "load" => 175e3)
        end
    end

    # Steuerungs- / Preis-Parameter (nur falls control_modules existieren)
    bp = get(cfg, "components", nothing)
    if bp !== nothing
        busp = get(bp, "BUS_Power", nothing)
        if busp !== nothing
            cms = get(busp, "control_modules", nothing)
            if cms !== nothing
                # candidate targets for control module linkage
                targets = default_control_targets(cfg)
                for cm in cms
                    if isa(cm, Dict)
                        # Marktpreis (limit)
                        if haskey(cm, "limit_price") && !isnan(p_stock_marginal)
                            cm["limit_price"] = p_stock_marginal
                            # economic_control requires new_connections_below_limit/above_limit when limit_price is used:
                            # - if the template already defines meaningful lists, keep them
                            # - otherwise set defaults (non-empty), so the control can be applied
                            if !haskey(cm, "new_connections_below_limit") || get(cm, "new_connections_below_limit", nothing) === nothing || isempty(cm["new_connections_below_limit"])
                                cm["new_connections_below_limit"] = copy(targets)
                            end
                            if !haskey(cm, "new_connections_above_limit") || get(cm, "new_connections_above_limit", nothing) === nothing
                                cm["new_connections_above_limit"] = []
                            end
                        end
                        # Reserve-Preis (z.B. aFRR)
                        if haskey(cm, "reserve_price") && !isnan(p_reserve_marginal)
                            cm["reserve_price"] = p_reserve_marginal
                            if !haskey(cm, "new_connections_when_reserve_bought") || get(cm, "new_connections_when_reserve_bought", nothing) === nothing || isempty(cm["new_connections_when_reserve_bought"])
                                cm["new_connections_when_reserve_bought"] = copy(targets)
                            end
                            if !haskey(cm, "new_connections_when_reserve_sold") || get(cm, "new_connections_when_reserve_sold", nothing) === nothing
                                cm["new_connections_when_reserve_sold"] = []
                            end
                        end
                    end
                end
            end
        end
    end

    # --- Sicherstellen, dass io_settings vorhanden ist ---
    # Wenn die Vorlage kein io_settings-Objekt hat, legen wir eines an, damit wir die Outputpfade
    # gezielt setzen können.
    if !haskey(cfg, "io_settings") || cfg["io_settings"] === nothing
        cfg["io_settings"] = Dict{String,Any}()
    end

    # Schreibe output nur als Dateien, falls die Daten benötigt werden, da dies den Prozess 
    # deutlich verlangsamt
    if write_output
        # --- Eindeutige Output-Dateinamen (Anzeige in MW / MWh) ---
        csv_name = "out_HP$(safe(display_MW(Pth_HP)))MW_HR$(safe(display_MW(Pth_HR)))MW_Buf$(safe(display_MWh(Cap_Wh)))MWh_Batt$(safe(display_MWh(BattCap_Wh)))MWh_Pmkt$(safe(p_stock_marginal))EUR_Pres$(safe(p_reserve_marginal))EUR.csv"
        aux_name = "aux_HP$(safe(display_MW(Pth_HP)))MW_HR$(safe(display_MW(Pth_HR)))MW_Buf$(safe(display_MWh(Cap_Wh)))MWh_Batt$(safe(display_MWh(BattCap_Wh)))MWh_Pmkt$(safe(p_stock_marginal))EUR_Pres$(safe(p_reserve_marginal))EUR.md"
        csv_abs = abspath(joinpath(outdir, csv_name))
        aux_abs = abspath(joinpath(outdir, aux_name))
        cfg["io_settings"]["csv_output_file"] = csv_abs
        cfg["io_settings"]["auxiliary_info_file"] = aux_abs

    else
        cfg["io_settings"]["csv_output_keys"] = "nothing"
        cfg["io_settings"]["auxiliary_info"] = false
        cfg["io_settings"]["sankey_plot"] = "nothing"
        cfg["io_settings"]["output_plot"] = "nothing"
    end

    # --- Lauf-spezifische Input-JSON speichern (Dateiname ebenfalls in kompakten Einheiten) ---
    fname = joinpath(outdir, "input_HP$(safe(display_MW(Pth_HP)))MW_HR$(safe(display_MW(Pth_HR)))MW_Buf$(safe(display_MWh(Cap_Wh)))MWh_Batt$(safe(display_MWh(BattCap_Wh)))MWh_Pmkt$(safe(p_stock_marginal))EUR_Pres$(safe(p_reserve_marginal))EUR.json")
    open(fname, "w") do io
        JSON.print(io, cfg)
    end

    # --- ReSiE-Prozess starten ---
    # Wir starten ReSiE extern, damit jeder Lauf in einem frischen Julia-Prozess läuft.
    # Bei Fehlern wird eine Exception ausgelöst und abgefangen.
    println("[$runidx/$total_runs] Starte Lauf: HP=$(display_MW(Pth_HP)) MW, HR=$(display_MW(Pth_HR)) MW, Buf=$(display_MWh(Cap_Wh)) MWh, Batt=$(display_MWh(BattCap_Wh)) MWh, p_mkt=$(p_stock_marginal) EUR, p_res=$(p_reserve_marginal) EUR -> JSON: $(basename(fname))")
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
    _, sim_output = Resie.load_and_run(fname)
    # try
    #     _, sim_output = Resie.load_and_run(fname)
    # catch e
    #     # Fehler auffangen, kurz protokollieren und im Log ablegen
    #     status = "ERROR"
    #     info = sprint(showerror, e)
    #     @warn "Run fehlgeschlagen: $info"
    # end

    # try
    #     run(`julia --project=. src/resie-cli.jl run $fname --exit-after-run`)
    # catch e
    #     # Fehler auffangen, kurz protokollieren und im Log ablegen
    #     status = "ERROR"
    #     info = sprint(showerror, e)
    #     @warn "Run fehlgeschlagen: $info"
    # end

    # --- Logeintrag schreiben ---
    # Jede Zeile enthält Zeitstempel, Run-Index, Parametrierung und Status --> nützlich für Auswertung
    open(logfile, "a") do f
        println(f, "$(timestamp),$(runidx),$(total_runs),$(variant),$(Pth_HP),$(Pth_HR),$(Cap_Wh),$(BattCap_Wh),$(p_stock_marginal),$(p_reserve_marginal),$(status),\"$(info)\"")
    end

    return sim_output
end

# Helper: default targets for economic/price control modules
function default_control_targets(cfg::AbstractDict)
    names = String[]
    comps = get(cfg, "components", nothing)
    if comps === nothing
        return names
    end

    # Prefer commonly used keys
    for name in ("HeatPump", "HeatingRod", "Battery", "P2H", "ElectricHeater", "DirectElectricHeater", "BufferTank")
        if haskey(comps, name)
            push!(names, name)
        end
    end

    # If nothing yet, pick components that look like consumers/producers (heuristic)
    if isempty(names)
        for (k, v) in comps
            if isa(v, Dict) && (haskey(v, "output_refs") || haskey(v, "load") || haskey(v, "power_th") || haskey(v, "capacity"))
                push!(names, String(k))
            end
        end
    end

    # Last resort: use all top-level component keys
    if isempty(names)
        for k in keys(comps)
            push!(names, String(k))
        end
    end

    return unique(names)
end

# Hauptschleife: vollständiges kartesisches Produkt 
# -----------------------
# Diese Schleife erzeugt und startet für jede mögliche Kombination (Pth_HP, Pth_HR, Cap_Wh) einen Run.
# total_runs = length(Pth_HP_vals) * length(Pth_HR_vals) * length(Cap_vals_Wh) * length(BattCap_vals_Wh) * length(p_stock_marginal_vals) * length(p_reserve_marginal_vals)
# println("Starte Parameterstudie: $total_runs Läufe (vollständiges kartesisches Produkt)")
# 
# runidx = 0
# sim_output = OrderedDict()
# for (Pth_HP, Pth_HR, Cap_Wh, BattCap_Wh, p_stock_marginal, p_reserve_marginal) in Iterators.product(Pth_HP_vals, Pth_HR_vals, Cap_vals_Wh, BattCap_vals_Wh, p_stock_marginal_vals, p_reserve_marginal_vals)
#     global runidx += 1
#     # Delegiere die eigentliche Erstellung der JSON, das Starten von ReSiE und das Logging an die Funktion
#     print(@elapsed global sim_output = run_resie_variant(outdir, base_input, Pth_HP, Pth_HR, Cap_Wh, BattCap_Wh, p_stock_marginal, p_reserve_marginal, runidx, total_runs))
# end

# Die frühere, vollständige Product-Loop wurde entfernt. main() baut die Listen variantenabhängig auf
# und führt nur das angepasste kartesische Produkt aus.

function main()
    # Parameterlisten vorbereiten: nur für die registrierten Parameter werden Variationen erzeugt.
    Pth_HP_loop_vals = hp_param ? Pth_HP_vals : [hp_default === nothing ? 0.0 : hp_default]
    Pth_HR_loop_vals = hr_param ? Pth_HR_vals : [hr_default === nothing ? 0.0 : hr_default]
    Cap_loop_vals = buffer_param ? Cap_vals_Wh : [buf_default === nothing ? 0.0 : buf_default]
    Batt_loop_vals = battery_param ? BattCap_vals_Wh : [batt_default === nothing ? 0.0 : batt_default]
    p_stock_loop_vals = use_p_stock_param ? p_stock_marginal_vals : [stock_default]
    p_reserve_loop_vals = use_p_reserve_param ? p_reserve_marginal_vals : [reserve_default]

    total_runs = length(Pth_HP_loop_vals) * length(Pth_HR_loop_vals) * length(Cap_loop_vals) * length(Batt_loop_vals) * length(p_stock_loop_vals) * length(p_reserve_loop_vals)
    println("Starte Parameterstudie Version $variant: $total_runs Läufe (angepasstes kartesisches Produkt)")

    runidx = 0
    sim_output = OrderedDict()
    # TODO Komponentenkosten definieren; VDIParams evtl. anpassen
    components = [Component(), Component()]
    params = VDIParams()
    for (Pth_HP, Pth_HR, Cap_Wh, BattCap_Wh, p_stock_marginal, p_reserve_marginal) in Iterators.product(Pth_HP_loop_vals, Pth_HR_loop_vals, Cap_loop_vals, Batt_loop_vals, p_stock_loop_vals, p_reserve_loop_vals)
        runidx += 1
        uuid = UUIDs.uuid4()
        print(@elapsed sim_output[uuid] = run_resie_variant(outdir, base_input,
                                                                    Pth_HP, Pth_HR,
                                                                    Cap_Wh, BattCap_Wh,
                                                                    p_stock_marginal, p_reserve_marginal,
                                                                    runidx, total_runs,
                                                                    variant,
                                                                    buffer_param, battery_param;
                                                                    write_output=false))
        # TODO Output korrekt zuweisen
        vdi_input = Dict()
        vdi_input["Grid_IN"] = sim_output[uuid]["GridConnection:IN"]
        vdi2067_annuity(vdi_input, components, params)
    end
end

# main() nur am Ende aufrufen (verhindert Zugriffe auf Bindings vor Definition)
main()
