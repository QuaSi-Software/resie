############################################################
# parameterstudy.jl — neue Version ohne Varianten
# -----------------------------------------------------------
# - Alle Komponenten sind immer vorhanden.
# - Eine Komponente wird deaktiviert, indem ihr Parameter = 0 ist.
# - Simulationsausgabe + VDI2067-Wirtschaftlichkeit pro UUID.
############################################################

using JSON
using Dates
using Resie
include("resie_logger.jl")
using .Resie_Logger
using OrderedCollections: OrderedDict
using UUIDs

include("vdi2067.jl")
using .VDI2067


############################################################
# Basis-Input laden
############################################################

base_input_path = length(ARGS) > 0 ? ARGS[1] : "inputfiles/input_file_base.json"


############################################################
# Parameterdefinition: untere Grenze, obere Grenze, Schrittweite
############################################################

# Wärmepumpe (HeatPump)
Pth_HP_lo   = 0e6
Pth_HP_hi   = 20e6
Pth_HP_step = 5e6
Pth_HP_vals = collect(Pth_HP_lo:Pth_HP_step:Pth_HP_hi)

# Boiler (Elektrodenkessel)
Pth_Boiler_lo   = 0e6
Pth_Boiler_hi   = 10e6
Pth_Boiler_step = 5e6
Pth_Boiler_vals = collect(Pth_Boiler_lo:Pth_Boiler_step:Pth_Boiler_hi)

# Pufferspeicher
Cap_lo_Wh   = 0
Cap_hi_Wh   = 200e6
Cap_step_Wh = 50e6
Cap_vals_Wh = collect(Cap_lo_Wh:Cap_step_Wh:Cap_hi_Wh)

# Batterie
Batt_lo_Wh   = 0
Batt_hi_Wh   = 200e3
Batt_step_Wh = 100e3
BattCap_vals_Wh = collect(Batt_lo_Wh:Batt_step_Wh:Batt_hi_Wh)

# Marktpreisgrenzen
p_stock_lo   = 80.0
p_stock_hi   = 80.0
p_stock_step = 10.0
p_stock_vals = collect(p_stock_lo:p_stock_step:p_stock_hi)

# Regelenergiepreise
p_res_lo   = 10.0
p_res_hi   = 10.0
p_res_step = 1.0
p_reserve_vals = collect(p_res_lo:p_res_step:p_res_hi)


############################################################
# Hilfsfunktionen
############################################################

function safe(x)
    replace(string(x), "." => "p")
end

function display_MW(W)
    round(W / 1e6; digits=3)
end

function display_MWh(Wh)
    round(Wh / 1e6; digits=3)
end


############################################################
# Simulation eines einzelnen Laufs
############################################################

function run_resie_variant(
        outdir::String,
        base_input::OrderedDict,
        Pth_HP::Number,
        Pth_Boiler::Number,
        Cap_Wh::Number,
        BattCap_Wh::Number,
        p_stock::Number,
        p_reserve::Number,
        runidx::Int,
        total_runs::Int;
        write_output::Bool=true
    )

    timestamp = Dates.format(now(), "yyyy-mm-dd_HH:MM:SS")
    cfg = deepcopy(base_input)

    ########################################################
    # Komponenten setzen (0 = deaktivieren)
    ########################################################
    cfg["components"]["HeatPump"]["power_th"] = Pth_HP
    cfg["components"]["Boiler"]["power_th"]   = Pth_Boiler
    cfg["components"]["BufferTank"]["capacity"] = Cap_Wh
    cfg["components"]["Battery"]["capacity"]    = BattCap_Wh

    ########################################################
    # Economic control Preise setzen (falls verfügbar)
    ########################################################
    if haskey(cfg["components"], "BUS_Power") &&
       haskey(cfg["components"]["BUS_Power"], "control_modules")

        for cm in cfg["components"]["BUS_Power"]["control_modules"]
            if haskey(cm, "limit_price")
                cm["limit_price"] = p_stock
            end
            if haskey(cm, "reserve_price")
                cm["reserve_price"] = p_reserve
            end
        end
    end

    ########################################################
    # Output-Dateien
    ########################################################

    if write_output
        csv_name = "out_HP$(safe(display_MW(Pth_HP)))MW_" *
                   "BOI$(safe(display_MW(Pth_Boiler)))MW_" *
                   "BUF$(safe(display_MWh(Cap_Wh)))MWh_" *
                   "BAT$(safe(display_MWh(BattCap_Wh)))MWh_" *
                   "P$(safe(p_stock))EUR_R$(safe(p_reserve))EUR.csv"

        aux_name = replace(csv_name, ".csv" => ".md")

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

    ########################################################
    # Lauf-spezifische JSON speichern
    ########################################################

    fname = joinpath(
        outdir,
        "input_HP$(safe(display_MW(Pth_HP)))MW_" *
        "BOI$(safe(display_MW(Pth_Boiler)))MW_" *
        "BUF$(safe(display_MWh(Cap_Wh)))MWh_" *
        "BAT$(safe(display_MWh(BattCap_Wh)))MWh.json"
    )

    open(fname, "w") do io
        JSON.print(io, cfg)
    end

    ########################################################
    # Simulation durchführen
    ########################################################

    println("[$runidx/$total_runs] → Starte Simulation: $(basename(fname))")

    logger = Resie_Logger.start_logger(true, false, nothing, nothing,
                                       Resie_Logger.Logging.Warn, fname)

    success, sim_output = Resie.load_and_run(fname, UUID(runidx))

    Resie_Logger.close_logger(logger)

    return sim_output
end



############################################################
# Hauptschleife
############################################################

function main(base_input_path)
    base_input = Resie.read_JSON(base_input_path)
    outdir = "./output/parameterstudy"
    mkpath(outdir)

    total_runs =
        length(Pth_HP_vals) *
        length(Pth_Boiler_vals) *
        length(Cap_vals_Wh) *
        length(BattCap_vals_Wh) *
        length(p_stock_vals) *
        length(p_reserve_vals)

    println("Parameterstudy startet: $total_runs Läufe")

    sim_output = OrderedDict()
    runidx = 0


    for (Pth_HP, Pth_Boiler, Cap_Wh, BattCap_Wh, p_stock, p_reserve) in
        Iterators.product(Pth_HP_vals, Pth_Boiler_vals, Cap_vals_Wh,
                          BattCap_vals_Wh, p_stock_vals, p_reserve_vals)

        runidx += 1

        ####################################################
        # Simulation
        ####################################################
        raw_sim = run_resie_variant(
            outdir, base_input,
            Pth_HP, Pth_Boiler, Cap_Wh, BattCap_Wh,
            p_stock, p_reserve,
            runidx, total_runs;
            write_output=true
        )

        # SimOutput wieder im alten Format
        sim_output[runidx] = Dict("sim" => raw_sim)


        ####################################################
        # Wirtschaftlichkeit VDI 2067
        ####################################################
        A0_HP     = 300 * (Pth_HP / 1e3)
        A0_Boiler = 80  * (Pth_Boiler / 1e3)
        A0_Buffer = 25  * (Cap_Wh / 1e3)
        A0_Batt   = 150 * (BattCap_Wh / 1e3)

        components = VDIComponent[
            heatpump_component(A0_HP,     20),
            boiler_component(A0_Boiler,   20),
            buffertank_component(A0_Buffer, 30),
            battery_component(A0_Batt,    15)
        ]

        sim_output[runidx]["VDI_NO"] =
            vdi2067_annuity(raw_sim, components, VDI_SCENARIO_NO)

        sim_output[runidx]["VDI_MOD"] =
            vdi2067_annuity(raw_sim, components, VDI_SCENARIO_MOD)

        sim_output[runidx]["VDI_PRO"] =
            vdi2067_annuity(raw_sim, components, VDI_SCENARIO_PRO)
    end

    println("✔ Parameterstudy abgeschlossen.")

    return sim_output
end


main(base_input_path)

