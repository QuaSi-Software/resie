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
import Resie.Resie_Logger: start_logger, Logging, close_logger
using Resie.Profiles 
using OrderedCollections: OrderedDict
using UUIDs
using Statistics: mean

include("vdi2067.jl")
import .VDI2067
using Infiltrator
############################################################
# Load base input
############################################################

base_input_path = length(ARGS) > 0 ? ARGS[1] : "inputfiles/inputfile_base_no_ems.json"


############################################################
# Definition of variable component parameters: lower limit, upper limit, step size
############################################################

# HeatPump power (W)
Pth_HP_lo   = 5.5e6         # lower limit
Pth_HP_hi   = 6.5e6         # upper limit
Pth_HP_step = 0.2e6         # step size
Pth_HP_vals = collect(Pth_HP_lo:Pth_HP_step:Pth_HP_hi)  # array of values

# Boiler power (W)
Pth_Boiler_lo   = 3.1e6     # lower limit
Pth_Boiler_hi   = 3.9e6     # upper limit
Pth_Boiler_step = 0.2e6     # step size
Pth_Boiler_vals = collect(Pth_Boiler_lo:Pth_Boiler_step:Pth_Boiler_hi)  # creates an array of values

# BufferTank capacity (Wh)
Cap_lo_Wh   = 61.0e6        # lower limit
Cap_hi_Wh   = 69.0e6        # upper limit
Cap_step_Wh = 2.0e6         # step size
Cap_vals_Wh = collect(Cap_lo_Wh:Cap_step_Wh:Cap_hi_Wh)  # creates an array of values

# Battery capacity (Wh)
Batt_lo_Wh   = 0e3          # lower limit
Batt_hi_Wh   = 450e3        # upper limit
Batt_step_Wh = 150e3        # step size
BattCap_vals_Wh = collect(Batt_lo_Wh:Batt_step_Wh:Batt_hi_Wh)   # creates an array of values

# define adjustments to the different price profiles in the order of
# [stock_price, reserve_power_neg_price, reserve_energy_neg_price, reserve_power_pos_price, 
#  reserve_energy_pos_price, market_value_pv, market_value_wind, co2_value_grid]
# TODO adjust values
# Stock Price Addon consists for Grid Fees of 40 €/MWh and Taxes of 65 €/MWh
profile_addons = [105.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
profile_multipliers = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]

############################################################
# define states in control module
############################################################

# if price_profile_stock > p_stock
state_save_power = Dict{String,Any}(
             "input_order"=> [
                 "NegControlReserve",
                 "Photovoltaic",
                 "WindFarm",
                 "Grid_IN",
                 "Battery"
             ],
             "output_order"=> [
                 "PosControlReserve",
                 "Demand_Power",
                 "HeatPump",
                 "Boiler",
                 "Battery",
                 "Grid_OUT"
             ],
             "energy_flow"=> [
                [0, 0, 0, 0, 0, 0],
                [1, 1, 1, 1, 1, 1],
                [1, 1, 1, 1, 1, 1],
                [1, 1, 1, 0, 0, 1],
                [1, 0, 0, 0, 0, 1]
             ])

# if price_profile_stock <= p_stock && price_profile_stock <= 0 (Eigenerzeugung)
state_earn_power =Dict{String,Any}(
            "input_order"=> [
                "NegControlReserve",
                "Photovoltaic",
                "WindFarm",
                "Grid_IN",
                "Battery"
            ],
            "output_order"=> [
                "PosControlReserve",
                "Demand_Power",
                "HeatPump",
                "Boiler",
                "Battery",
                "Grid_OUT"
            ],
            "energy_flow"=> [
                [0, 0, 0, 0, 0, 0],
                [1, 1, 1, 1, 1, 1],
                [1, 1, 1, 1, 1, 1],
                [1, 1, 1, 0, 0, 1],
                [1, 0, 0, 0, 0, 1]
            ])

# if price_profile_stock > p_stock
state_save_heat = Dict{String,Any}(
             "input_order"=> [
                "HeatPump",
                "Boiler",
                "BufferTank",
                "HeatPump#secondary",
                "Boiler#secondary"
            ],
            "output_order"=> [
                "Demand_Heat",
                "BufferTank"
            ],
            "energy_flow"=> [
                [1, 1],
                [1, 1],
                [1, 0],
                [1, 0],
                [1, 0]
            ])

# if price_profile_stock <= p_stock && price_profile_stock <= 0 (Eigenerzeugung)
state_earn_heat =Dict{String,Any}(
            "input_order"=> [
                "HeatPump",
                "Boiler",
                "HeatPump#secondary",
                "Boiler#secondary",
                "BufferTank"
            ],
            "output_order"=> [
                "Demand_Heat",
                "BufferTank"
            ],
            "energy_flow"=> [
                [1, 1],
                [1, 1],
                [1, 1],
                [1, 1],
                [1, 0]
            ])


states_heat_bus = [state_save_heat, state_earn_heat]
states_power_bus = [state_save_power, state_earn_power]

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

function display_kWh(Wh)    
    round(Wh / 1e3; digits=3)       # displays kWh
end

function save_to_prf(timestamps::Array{Int,1}, values::Array{Float64,1}, filepath::String)
    header_variables = ["# data_type:", "# time_definition:", "# profile_start_date:", 
                        "# profile_start_date_format:", "# timestamp_format:", 
                        "# interpolation_type:"]
    header_values = ["intensive", "startdate_timestamp", "01.01.2024 00:00", 
                     "dd.mm.yyyy HH:MM", "seconds", "stepwise"]
    open(filepath, "w") do file_handle
        for (var, val) in zip(header_variables, header_values)
            write(file_handle, var * "\t" * val * "\n")
        end
        for (ts, val) in zip(timestamps, values)
            write(file_handle, string(ts) * ";" * string(val) * "\n")
        end      
    end

    println("Profile file at $filepath created.")
end

function save_to_prf(dates::Array{DateTime,1}, values::Array{Float64,1}, filepath::String)
    header_variables = ["# data_type:", "# time_definition:", "# timestamp_format:", 
                        "# time_zone:", "# interpolation_type:"]
    header_values = ["intensive", "datestamp", "dd.mm.yyyy HH:MM", 
                     "Europe/Berlin", "stepwise"]
    open(filepath, "w") do file_handle
        for (var, val) in zip(header_variables, header_values)
            write(file_handle, var * "\t" * val * "\n")
        end
        for (ts, val) in zip(dates, values)
            write(file_handle, string(ts) * ";" * string(val) * "\n")
        end      
    end

    println("Profile file at $filepath created.")
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
        states_power_bus::Array{Dict{String,Any}},
        states_heat_bus::Array{Dict{String,Any}},
        profile_addons,
        profile_multipliers,
        runidx::Int,
        total_runs::Int;
        write_output::Bool=true
        )
    cfg = deepcopy(base_input)

    ########################################################
    # Set components
    ########################################################
    cfg["components"]["HeatPump"]["power_th"] = Pth_HP
    cfg["components"]["Boiler"]["power_th"] = Pth_Boiler
    cfg["components"]["BufferTank"]["capacity"] = Cap_Wh
    cfg["components"]["Battery"]["capacity"] = BattCap_Wh
    
    ########################################################
    # Output file
    ########################################################

    if write_output
        csv_name = "out_HP$(safe(display_MW(Pth_HP)))MW_" *
                   "BOI$(safe(display_MW(Pth_Boiler)))MW_" *
                   "BUF$(safe(display_MWh(Cap_Wh)))MWh_" *
                   "BAT$(safe(display_kWh(BattCap_Wh)))kWh.csv"

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

    # set price profiles paths
    price_profile_path_stock = "./profiles/MA/boersenpreis_EUR_MWh.prf"
    price_profile_path_reserve_power_neg = "./profiles/MA/aFRR_neg_cap_EUR_MW_h.prf"
    price_profile_path_reserve_energy_neg = "./profiles/MA/cbmp_down_mean_15min.prf"
    price_profile_path_reserve_power_pos = "./profiles/MA/aFRR_pos_cap_EUR_MW_h.prf"
    price_profile_path_reserve_energy_pos = "./profiles/MA/cbmp_up_mean_15min.prf"
    price_profile_path_market_value_pv = "./profiles/MA/MW_Solar.prf"
    price_profile_path_market_value_wind = "./profiles/MA/MW_Wind.prf"
    co2_profile_path_grid = "./profiles/MA/CO2_life_g_kWh.prf"

    profile_paths = [price_profile_path_stock, price_profile_path_reserve_power_neg, price_profile_path_reserve_energy_neg,
                     price_profile_path_reserve_power_pos, price_profile_path_reserve_energy_pos,
                     price_profile_path_market_value_pv, price_profile_path_market_value_wind,
                     co2_profile_path_grid]

    # create correct profiles from price_profile_paths and add the to sim_output.
    # profiles will be overwritten for every run to make sure the profile_multipliers and 
    # profile_addons calculated correctly.
    # If multiple threads are used each thread gets their own profile.
    run_ID = UUID(runidx)
    sim_params, _, _ = Resie.prepare_inputs(cfg, run_ID)
    
    profile_dir = "./profiles/MA/profiles_parameterstudy"
    mkpath(profile_dir)

    profile_id = ifelse(Threads.nthreads() > 1, Threads.threadid(), 0)
    date_range = remove_leap_days(collect(sim_params["start_date"]:Second(sim_params["time_step_seconds"]):sim_params["end_date"]))
    new_paths = Array{String}(undef, length(profile_paths)) 
    profile_values = Array{Float64}(undef, length(profile_paths), length(date_range)) 

    for (p_idx, path) in enumerate(profile_paths)
        if profile_multipliers[p_idx] != 1 && profile_addons[p_idx] != 0 && profile_id == 0
            new_paths[p_idx] = path
        else
            profile = Profile(path, sim_params)
            values = [profile.data[dt] .* profile_multipliers[p_idx] .+ profile_addons[p_idx] for dt in date_range]
            new_path = profile_dir * "/" * split(path[1:end-4], '/')[end] * "_$profile_id.prf" 
            save_to_prf(collect(date_range), values, new_path)
            new_paths[p_idx] = path
            profile_values[p_idx, :] = values
        end
    end

    ########################################################
    # Set limits / benchmark for economic_control.jl
    ########################################################
    if haskey(cfg["components"], "BUS_Power") &&
       haskey(cfg["components"]["BUS_Power"], "control_modules")
        cm = cfg["components"]["BUS_Power"]["control_modules"][1]
        cm["price_profile_paths"] = new_paths[1:5]
        cm["bus_uacs"] = ["BUS_Power", "BUS_Heat"]
        cm["new_connections_below_limits"] = Dict("BUS_Power" => states_power_bus, 
                                                  "BUS_Heat" => states_heat_bus)
    end

    ########################################################
    # safe run specific inputfile_"xxx".json
    ########################################################

    fname = joinpath(
        outdir,
        "input_HP$(safe(display_MW(Pth_HP)))MW_" *
        "BOI$(safe(display_MW(Pth_Boiler)))MW_" *
        "BUF$(safe(display_MWh(Cap_Wh)))MWh_" *
        "BAT$(safe(display_kWh(BattCap_Wh)))kWh.json"
    )

    open(fname, "w") do io
        JSON.print(io, cfg)
    end

    ########################################################
    # run simulation
    ########################################################

    println("[$runidx/$total_runs] → start simulation: $(basename(fname))")

    success, sim_output = Resie.load_and_run(fname, run_ID)
    Resie.close_run(run_ID)

    # remove addon for used grid electricity to went into PosControlReserve
    energy_with_addon = max.(sim_output["Grid_IN m_power OUT"] .- sim_output["PosControlReserve m_power IN"], 0)
    energy_without_addon = min.(sim_output["Grid_IN m_power OUT"], sim_output["PosControlReserve m_power IN"])
    profile_values[1, :] .-= profile_addons[1] .* energy_without_addon ./ (energy_with_addon .+ energy_without_addon)


    # add the profiles to the sim_output to be used in vdi2067 calculation
    sim_output["Stock_Price"] = profile_values[1, :]    #TODO is this with or without grid add ons? -> with grid_addons for all profiles
    sim_output["Reserve_Power_Price_Neg"] = profile_values[2, :]
    sim_output["Reserve_Energy_Price_Neg"] = profile_values[3, :]
    sim_output["Reserve_Power_Price_Pos"] = profile_values[4, :]
    sim_output["Reserve_Energy_Price_Pos"] = profile_values[5, :]
    sim_output["Market_Price_PV"] = profile_values[6, :]
    sim_output["Market_Price_Wind"] = profile_values[7, :]
    sim_output["CO2_Grid"] = profile_values[8, :]

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
        length(BattCap_vals_Wh)

    println("Parameterstudy starts: $total_runs runs on $(Threads.nthreads()) parallel Threads")

    out_file_path = outdir * "/results_$(total_runs)runs_" * Dates.format(now(), "yymmdd_HHMMSS") * ".csv"
    touch(out_file_path)
    header = join([
                #component parameters
                "Hp_Power / W", "Boiler_Power / W", "BufferTank_Capacity / Wh", "Battery_Capacity / Wh", #compontet parameters

                # no cost escalation
                "annuity_no A_cap / €", "annuity_no A_cap_incentive / €", "annuity_no A_misc / €", "annuity_no A_op / €", "annuity_no A_energy / €", 
                "annuity_no A_rev_control / €", "annuity_no A_rev_feed / €", "annuity_no A_total / €", "annuity_no A_total_incentive / €",

                # moderate cost escalation
                "annuity_mod A_cap / €", "annuity_mod A_cap_incentive / €", "annuity_mod A_misc / €", "annuity_mod A_op / €", "annuity_mod A_energy / €", 
                "annuity_mod A_rev_control / €", "annuity_mod A_rev_feed / €", "annuity_mod A_total / €", "annuity_mod A_total_incentive / €",

                # progressive cost escalation
                "annuity_pro A_cap / €", "annuity_pro A_cap_incentive / €", "annuity_pro A_misc / €", "annuity_pro A_op / €", "annuity_pro A_energy / €", 
                "annuity_pro A_rev_control / €", "annuity_pro A_rev_feed / €", "annuity_pro A_total / €", "annuity_pro A_total_incentive / €",

                # balance warnings
                "balance_power", "balance_heat", 

                # yearly CO2-emissions
                # "CO2_yearly"
                ], ';') * "\n"

    open(out_file_path, "a") do file_handle
        write(file_handle, header)
    end

    println("Results file at $out_file_path created.")

    output_lock = ReentrantLock()

    runidx_global = Threads.Atomic{Int}(1)

    logger = start_logger(true, false, nothing, nothing,
                          Logging.Warn, nothing)
    Logging.global_logger(logger)

    Threads.@threads for (Pth_HP, Pth_Boiler, Cap_Wh, BattCap_Wh) in
        collect(Iterators.product(Pth_HP_vals, Pth_Boiler_vals, Cap_vals_Wh,
                                  BattCap_vals_Wh))
        
        runidx = Threads.atomic_add!(runidx_global, 1)
        sim_output = OrderedDict()
        start_time = now()
        ####################################################
        # Simulation
        ####################################################
        raw_sim = run_resie_variant(
            outdir, base_input,
            Pth_HP, Pth_Boiler, Cap_Wh, BattCap_Wh,
            states_power_bus, states_heat_bus,
            profile_addons, profile_multipliers,
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
        # TODO Adjust factors in front (EUR/kW or EUR/kWh)
        A0_HP     = 700 * (Pth_HP / 1e3)        #   
        A0_Boiler = 285  * (Pth_Boiler / 1e3)   # circa 1 Mio. € per 3.5 MW #TODO Christian fragen
        A0_Buffer = 39  * (Cap_Wh / 1e3)        # from FACT document
        A0_Batt   = 375 * (BattCap_Wh / 1e3)    # TODO Jule fragen
       
        components = VDI2067.VDIComponent[
            VDI2067.heatpump_component(A0_HP),
            VDI2067.boiler_component(A0_Boiler),
            VDI2067.buffertank_component(A0_Buffer),
            VDI2067.battery_component(A0_Batt)
        ]

        sim_output[runidx]["VDI_NO"] = 
            VDI2067.vdi2067_annuity(raw_sim, components, VDI2067.VDI_SCENARIO_NONE)

        sim_output[runidx]["VDI_MOD"] = 
            VDI2067.vdi2067_annuity(raw_sim, components, VDI2067.VDI_SCENARIO_MOD)

        sim_output[runidx]["VDI_PRO"] = 
            VDI2067.vdi2067_annuity(raw_sim, components, VDI2067.VDI_SCENARIO_PRO)

        # write important results to seperate file
        parameters = [Pth_HP, Pth_Boiler, Cap_Wh, BattCap_Wh]
        annuities_no = collect(values(sim_output[runidx]["VDI_NO"]))
        annuities_mod = collect(values(sim_output[runidx]["VDI_MOD"]))
        annuities_pro = collect(values(sim_output[runidx]["VDI_PRO"]))
        balances = collect(getindex.(Ref(sim_output[runidx]), ("Balance_heat", "Balance_power")))
        # co2 = collect(values(sim_output[runidx])) #TODO
        row = join(vcat(parameters, annuities_no, annuities_mod, annuities_pro, balances), ';') * "\n"
        row = replace(row, '.' => ',')
        # Lock the file writing
        Threads.lock(output_lock)
        try
            open(out_file_path, "a") do file_handle
                write(file_handle, row)
            end
        finally
            Threads.unlock(output_lock)
        end

        runtime = round(Int, Dates.seconds(now() - start_time))
        eta = round(Int, (total_runs - runidx) * runtime / Threads.nthreads() / 60)
        println("[$runidx/$total_runs] → completed in $runtime s. ETA: $eta min")
    end

    close_logger(logger)
    println("✔ Parameterstudy completed.")
end

main(base_input_path, false)
