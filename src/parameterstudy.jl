############################################################
# parameterstudy.jl - a parameterstudy tool for resie
# -----------------------------------------------------------
# - All components are defined here
# - A component is deactivated by setting power parameters = 0
# - Parameterstudy calls for vdi2067.jl and outputs annuity and  heat price per MWh
############################################################

using JSON
using Dates
using Resie
include("resie_logger.jl")
using .Resie_Logger
using OrderedCollections: OrderedDict
using UUIDs
using CSV

include("vdi2067.jl")
using .VDI2067

include("profiles/base.jl")
using .Profiles
using Infiltrator

############################################################
# Load base input
############################################################

base_input_path = length(ARGS) > 0 ? ARGS[1] : "inputfiles/inputfile_base.json"


############################################################
# Definition of variable component parameters: lower limit, upper limit, step size
############################################################

# HeatPump power (W)
Pth_HP_lo   = 0.0e6         # lower limit
Pth_HP_hi   = 9.0e6         # upper limit
Pth_HP_step = 3.0e6         # step size
Pth_HP_vals = collect(Pth_HP_lo:Pth_HP_step:Pth_HP_hi)  # array of values

# Boiler power (W)
Pth_Boiler_lo   = 0.0e6     # lower limit
Pth_Boiler_hi   = 7.0e6     # upper limit
Pth_Boiler_step = 3.5e6     # step size
Pth_Boiler_vals = collect(Pth_Boiler_lo:Pth_Boiler_step:Pth_Boiler_hi)  # creates an array of values

# BufferTank capacity (Wh)
Cap_lo_Wh   = 0.0e6         # lower limit
Cap_hi_Wh   = 90.0e6        # upper limit
Cap_step_Wh = 30.0e6        # step size
Cap_vals_Wh = collect(Cap_lo_Wh:Cap_step_Wh:Cap_hi_Wh)  # creates an array of values

# Battery capacity (Wh)
Batt_lo_Wh   = 0            # lower limit
Batt_hi_Wh   = 200e3        # upper limit
Batt_step_Wh = 100e3        # step size
BattCap_vals_Wh = collect(Batt_lo_Wh:Batt_step_Wh:Batt_hi_Wh)   # creates an array of values

# limit for energy stock prices for economic_control.jl (EUR/MWh)
# if no grid price limit is to be considered, limit is set "towards infinity"
# TODO adjust limits
p_stock_lo   = 20.0         # lower limit
p_stock_hi   = 100.0        # upper limit
p_stock_step = 20.0         # step size
p_stock_vals = collect(p_stock_lo:p_stock_step:p_stock_hi)  # creates an array of values

# benchmark (smallest accepted value) for control reserve revenue per 4 hour time slot for economic_control.jl (EUR/4h-slot)
# if control energy is not to be considerd, benchmark is set "towards infinity"
# TODO adjust limits
p_res_lo   = 10.0           # lower limit
p_res_hi   = 10.0           # upper limit
p_res_step = 1.0            # step size
p_reserve_vals = collect(p_res_lo:p_res_step:p_res_hi)  # creates an array of values


############################################################
# helper functions
############################################################

function safe(x)            
    replace(string(x), "." => "p")  # creates a filename-safe string
end

function display_MW(W)      
    round(W / 1e6; digits=3)        # displays MW
end

function display_MWh(Wh)    
    round(Wh / 1e6; digits=3)       # displays MWh
end


############################################################
# Simulation of a single run
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
    # TODO variable timestamp has beed assigned but not used?
    timestamp = Dates.format(now(), "yyyy-mm-dd_HH:MM:SS")
    cfg = deepcopy(base_input)

    ########################################################
    # Set components
    ########################################################
    cfg["components"]["HeatPump"]["power_th"] = Pth_HP
    cfg["components"]["Boiler"]["power_th"]   = Pth_Boiler
    cfg["components"]["BufferTank"]["capacity"] = Cap_Wh
    cfg["components"]["Battery"]["capacity"]    = BattCap_Wh
    
    # set price profiles paths
    price_profile_path_stock = "./profiles/MA/boersenpreis_EUR_kWh.prf"
    price_profile_path_reserve_power = "./profiles/MA/aFRR_neg_cap_EUR_MW_h.prf"
    price_profile_path_reserve_energy = "./profiles/MA/aFRR_neg_energy_EUR_MWh.prf"
    price_profile_path_market_value_pv = "./profiles/MA/market_value_PV_EUR_kWh.prf"
    price_profile_path_market_value_wind = "./profiles/MA/market_value_Wind_EUR_kWh.prf"

    ########################################################
    # Set limits / benchmark for economic_control.jl
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
    # Output file
    # TODO write all simulation outputs into one single file to reduce amount of files
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
    # safe run specific inputfile_"xxx".json
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
    # run simulation
    ########################################################

    println("[$runidx/$total_runs] → start simulation: $(basename(fname))")

    logger = Resie_Logger.start_logger(true, false, nothing, nothing,
                                       Resie_Logger.Logging.Warn, fname)
    Resie_Logger.Logging.global_logger(logger)

    run_ID = UUID(runidx)
    success, sim_output = Resie.load_and_run(fname, run_ID)
    Resie.close_run(run_ID)

    Resie_Logger.close_logger(logger)

    # create correct profiles from price_profile_paths and add the to sim_output
    sim_params, _, _ = Resie.prepare_inputs(cfg, run_ID)
    price_profile_stock = Profiles.Profile(price_profile_path_stock, sim_params)
    price_profile_reserve_power = Profiles.Profile(price_profile_path_reserve_power, sim_params)
    price_profile_reserve_energy = Profiles.Profile(price_profile_path_reserve_energy, sim_params)
    price_profile_market_value_pv = Profiles.Profile(price_profile_path_market_value_pv, sim_params)
    price_profile_market_value_wind = Profiles.Profile(price_profile_path_market_value_wind, sim_params)

    stock_values = []
    reserve_power_values = []
    reserve_energy_values = []
    market_value_pv_values = []
    market_value_wind_values = []
    for dt in keys(price_profile_stock.data)
        push!(stock_values, price_profile_stock.data[dt])
        push!(reserve_power_values, price_profile_reserve_power.data[dt])
        push!(reserve_energy_values, price_profile_reserve_energy.data[dt])
        push!(market_value_pv_values, price_profile_market_value_pv.data[dt])
        push!(market_value_wind_values, price_profile_market_value_wind.data[dt])
    end 
    sim_output["Stock_Price"] = stock_values
    sim_output["Reserve_Power_Price"] = reserve_power_values
    sim_output["Reserve_Energy_Price"] = reserve_energy_values
    sim_output["Market_Price_PV"] = market_value_pv_values
    sim_output["Market_Price_Wind"] = market_value_wind_values

    return sim_output
end



############################################################
# Main Loop
############################################################

function main(base_input_path, write_output)
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

    println("Parameterstudy starts: $total_runs runs")

    out_file_path = outdir * "/results_$(total_runs)runs_" * Dates.format(now(), "yymmdd_HHMMSS") * ".csv"
    touch(out_file_path)
    header = join(["Hp_Power / W", "Boiler_Power / W", "BufferTank_Capacity / Wh", 
                   "Battery_Capacity / Wh", "Stock_Price / €/MWh", "Reserve_Price / €/4h", 
                   "annuity_no / €", "annuity_mod / €", "annuity_pro / €", 
                   "balance_power", "balance_heat"], ';') * "\n"
    open(out_file_path, "a") do file_handle
        write(file_handle, header)
    end

    println("Results file at $out_file_path created.")

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
            write_output=write_output
        )

        # SimOutput in old format again
        sim_output[runidx] = Dict{String, Any}("sim" => raw_sim)
        sim_output[runidx]["Balance_heat"] = sum(raw_sim["BUS_Heat Balance"])
        sim_output[runidx]["Balance_power"] = sum(raw_sim["BUS_Power Balance"])


        ####################################################
        # Econcomy based on VDI 2067 principles
        ####################################################
        # function for component investment costs based on installed capacity (EUR/kW)
        # TODO Adjust factors in front
        A0_HP     = 300 * (Pth_HP / 1e3)        
        A0_Boiler = 80  * (Pth_Boiler / 1e3)
        A0_Buffer = 25  * (Cap_Wh / 1e3)
        A0_Batt   = 150 * (BattCap_Wh / 1e3)
        # component lifetimes # TODO Adjust values in front and double check years
        TN_HP = 20
        TN_Boiler = 20
        TN_Buffer = 20
        TN_Batt = 12

        components = VDI2067.VDIComponent[
            VDI2067.heatpump_component(A0_HP, TN_HP),
            VDI2067.boiler_component(A0_Boiler, TN_Boiler),
            VDI2067.buffertank_component(A0_Buffer, TN_Buffer),
            VDI2067.battery_component(A0_Batt, TN_Batt)
        ]

        sim_output[runidx]["VDI_NO"] = 1
            # VDI2067.vdi2067_annuity(raw_sim, components, VDI2067.VDI_SCENARIO_NO)

        sim_output[runidx]["VDI_MOD"] = 2
            # VDI2067.vdi2067_annuity(raw_sim, components, VDI2067.VDI_SCENARIO_MOD)

        sim_output[runidx]["VDI_PRO"] = 3
            # VDI2067.vdi2067_annuity(raw_sim, components, VDI2067.VDI_SCENARIO_PRO)

        # write important results to seperate file
        parameters = [Pth_HP, Pth_Boiler, Cap_Wh, BattCap_Wh, p_stock, p_reserve]
        annuities = collect(getindex.(Ref(sim_output[runidx]), ("VDI_NO", "VDI_MOD", "VDI_PRO")))
        balances = collect(getindex.(Ref(sim_output[runidx]), ("Balance_heat", "Balance_power")))
        row = join(vcat(parameters, annuities, balances), ';') * "\n"
        open(out_file_path, "a") do file_handle
            write(file_handle, row)
        end
    end

    println("✔ Parameterstudy completed.")

    return sim_output
end


main(base_input_path, false)

