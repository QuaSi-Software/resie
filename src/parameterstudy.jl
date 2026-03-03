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

base_input_path = length(ARGS) > 0 ? ARGS[1] : "inputfiles/inputfile_base_ems.json"


############################################################
# Definition of variable component parameters: lower limit, upper limit, step size
############################################################

# HeatPump power (W)
Pth_HeatPump_lo   = 4.0e6         # lower limit
Pth_HeatPump_hi   = 10.0e6         # upper limit
Pth_HeatPump_step = 1.0e6         # step size
Pth_HeatPump_vals = collect(Pth_HeatPump_lo:Pth_HeatPump_step:Pth_HeatPump_hi)  # array of values

# ElectrodeBoiler power (W)
Pth_ElectrodeBoiler_lo   = 0.0e6     # lower limit
Pth_ElectrodeBoiler_hi   = 5.0e6     # upper limit
Pth_ElectrodeBoiler_step = 1.0e6     # step size
Pth_ElectrodeBoiler_vals = collect(Pth_ElectrodeBoiler_lo:Pth_ElectrodeBoiler_step:Pth_ElectrodeBoiler_hi)  # creates an array of values

# BufferTank capacity (Wh)
Cap_lo_Wh   = 20.0e6        # lower limit
Cap_hi_Wh   = 100.0e6        # upper limit
Cap_step_Wh = 10.0e6         # step size
Cap_vals_Wh = collect(Cap_lo_Wh:Cap_step_Wh:Cap_hi_Wh)  # creates an array of values

# Battery capacity (Wh)
Batt_lo_Wh   = 1.0e3          # lower limit
Batt_hi_Wh   = 1.0e3        # upper limit
Batt_step_Wh = 0.5e3        # step size
BattCap_vals_Wh = collect(Batt_lo_Wh:Batt_step_Wh:Batt_hi_Wh)   # creates an array of values

# define adjustments to the different price profiles in the order of
# [grid_price, reserve_power_neg_price, reserve_energy_neg_price, reserve_power_pos_price, 
#  reserve_energy_pos_price, market_value_pv, market_value_wind, co2_value_grid]
# TODO adjust values
# Grid Price Addon consists for Grid Fees of 40 €/MWh and Taxes of 65 €/MWh
profile_addons = [105.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]   #TODO CHANGED VALUES FROM 105.0
profile_multipliers = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]

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


function save_to_prf(timestamps::Array{Int,1},
                     values::Array{Float64,1},
                     filepath::String;
                     profile_start_date::DateTime=DateTime(2024, 1, 1, 0, 0))
    header_variables = ["# data_type:", "# time_definition:", "# profile_start_date:", 
                        "# profile_start_date_format:", "# timestamp_format:", 
                        "# interpolation_type:"]
    header_values = ["intensive", "startdate_timestamp", Dates.format(profile_start_date, "dd.mm.yyyy HH:MM"), 
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
            write(file_handle, Dates.format(ts, "dd.mm.yyyy HH:MM") * ";" * string(val) * "\n")
        end      
    end

    println("Profile file at $filepath created.")
end


############################################################
# Simulation of a single run
############################################################

function create_variant(
        outdir::String,
        base_input::OrderedDict,
        Pth_HeatPump::Number,
        Pth_ElectrodeBoiler::Number,
        Cap_Wh::Number,
        BattCap_Wh::Number,
        profile_addons,
        profile_multipliers,
        run_ID::UUID;
        write_output::Bool=true
        )
    cfg = deepcopy(base_input)

    ########################################################
    # Set components
    ########################################################
    cfg["components"]["HeatPump"]["power_th"] = Pth_HeatPump
    cfg["components"]["ElectrodeBoiler"]["power_th"] = Pth_ElectrodeBoiler
    cfg["components"]["BufferTank"]["capacity"] = Cap_Wh
    cfg["components"]["Battery"]["capacity"] = BattCap_Wh
    
    ########################################################
    # Output file
    ########################################################

    if write_output
        csv_name = "out_HP$(safe(display_MW(Pth_HeatPump)))MW_" *
                   "EB$(safe(display_MW(Pth_ElectrodeBoiler)))MW_" *
                   "BUF$(safe(display_MWh(Cap_Wh)))MWh_" *
                   "BAT$(safe(display_MWh(BattCap_Wh)))MWh.csv"

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
    price_profile_path_grid = "./profiles/MA/boersenpreis_EUR_MWh.prf"
    price_profile_path_reserve_power_neg = "./profiles/MA/aFRR_neg_cap_EUR_MW_h.prf"
    price_profile_path_reserve_energy_neg = "./profiles/MA/cbmp_down_mean_15min.prf"
    profile_path_reserve_call_neg = "./profiles/MA/reserve_call_neg.prf" #TODO Abrufprofil negativ hinterlegen
    price_profile_path_reserve_power_pos = "./profiles/MA/aFRR_pos_cap_EUR_MW_h.prf"
    price_profile_path_reserve_energy_pos = "./profiles/MA/cbmp_up_mean_15min.prf"
    profile_path_reserve_call_pos = "./profiles/MA/reserve_call_pos.prf" #TODO Abrufprofil positiv hinterlegen
    price_profile_path_market_value_pv = "./profiles/MA/MW_Solar.prf"
    price_profile_path_market_value_wind = "./profiles/MA/MW_Wind.prf"
    co2_profile_path_grid = "./profiles/MA/CO2_life_g_kWh.prf"

    profile_paths = [price_profile_path_grid, 
                     price_profile_path_reserve_power_neg, 
                     price_profile_path_reserve_energy_neg, profile_path_reserve_call_neg,
                     price_profile_path_reserve_power_pos, 
                     price_profile_path_reserve_energy_pos, profile_path_reserve_call_pos,
                     price_profile_path_market_value_pv, 
                     price_profile_path_market_value_wind,
                     co2_profile_path_grid]

    # create correct profiles from price_profile_paths and add the to sim_output.
    # profiles will be overwritten for every run to make sure the profile_multipliers and 
    # profile_addons calculated correctly.
    # If multiple threads are used each thread gets their own profile.
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
            timestamps = collect(0:Int(sim_params["time_step_seconds"]):
                                 Int(sim_params["time_step_seconds"]) * (length(values) - 1))
            new_path = profile_dir * "/" * split(path[1:end-4], '/')[end] * "_$profile_id.prf" 
            save_to_prf(timestamps,
                        values,
                        new_path;
                        profile_start_date=sim_params["start_date"])
            new_paths[p_idx] = new_path
            profile_values[p_idx, :] = values
        end
    end

    ########################################################
    # Set limits / benchmark for ems.jl
    ########################################################
    if haskey(cfg["components"], "BUS_Power") &&
       haskey(cfg["components"]["BUS_Power"], "control_modules")
        cm = cfg["components"]["BUS_Power"]["control_modules"][1]
        cm["grid_price_path"] = new_paths[1]
        cm["reserve_power_price_neg_path"] = new_paths[2]
        cm["reserve_energy_price_neg_path"] = new_paths[3]
        cm["reserve_call_neg_path"] = new_paths[4]
        cm["reserve_power_price_pos_path"] = new_paths[5]
        cm["reserve_energy_price_pos_path"] = new_paths[6]
        cm["reserve_call_pos_path"] = new_paths[7]
        # cm["price_profile_paths"] = new_paths[1:5]
        # cm["bus_uacs"] = ["BUS_Power", "BUS_Heat"]
        # cm["new_connections_below_limits"] = Dict("BUS_Power" => states_power_bus, 
        #                                           "BUS_Heat" => states_heat_bus)
    end

    

    ########################################################
    # safe run specific inputfile_"xxx".json
    ########################################################

    fname = joinpath(
        outdir,
        "inputfiles",
        "input_HP$(safe(display_MW(Pth_HeatPump)))MW_" *
        "EB$(safe(display_MW(Pth_ElectrodeBoiler)))MW_" *
        "BUF$(safe(display_MWh(Cap_Wh)))MWh_" *
        "BAT$(safe(display_MWh(BattCap_Wh)))MWh.json"
    )

    open(fname, "w") do io
        JSON.print(io, cfg)
    end

    return fname, profile_values, Int(sim_params["time_step_seconds"])
end

function run_resie_variant(input_file::String, run_ID::UUID, profile_values, profile_addons)

    ########################################################
    # run simulation
    ########################################################
    success, sim_output = Resie.load_and_run(input_file, run_ID)

    # add the profiles to the sim_output to be used in vdi2067 calculation
    sim_output["Grid_Price"] = profile_values[1, :]    #TODO is this with or without grid add ons? -> with grid_addons for all profiles
    sim_output["Reserve_Power_Price_Neg"] = profile_values[2, :]
    sim_output["Reserve_Energy_Price_Neg"] = profile_values[3, :]
    sim_output["Reserve_Call_Neg"] = profile_values[4, :]
    sim_output["Reserve_Power_Price_Pos"] = profile_values[5, :]
    sim_output["Reserve_Energy_Price_Pos"] = profile_values[6, :]
    sim_output["Reserve_Call_Pos"] = profile_values[7, :]
    sim_output["Market_Price_PV"] = profile_values[8, :]
    sim_output["Market_Price_Wind"] = profile_values[9, :]
    sim_output["CO2_Grid"] = profile_values[10, :]
    # TODO activate for parameterstudy with ems
    #cfg = Resie.read_JSON(input_file)
    #cm = cfg["components"]["BUS_Power"]["control_modules"][1]
    #sim_output["Offer_only_freebids"] = get(cm, "offer_only_freebids", false)

    return sim_output
end

############################################################
# Main Loop
############################################################

function main(base_input_path::String, write_output::Bool=false, save_input_files::Bool=true)
    base_input = Resie.read_JSON(base_input_path)
    outdir = "./output/parameterstudy"
    mkpath(outdir)
    mkpath(joinpath(outdir, "inputfiles"))

    all_combinations = collect(Iterators.product(
        Pth_HeatPump_vals,
        Pth_ElectrodeBoiler_vals,
        Cap_vals_Wh,
        BattCap_vals_Wh
    ))
    valid_combinations = filter(all_combinations) do (Pth_HeatPump, Pth_ElectrodeBoiler, _, _)
        (Pth_HeatPump + Pth_ElectrodeBoiler) > 5.0e6
    end
    skipped_runs = length(all_combinations) - length(valid_combinations)
    total_runs = length(valid_combinations)
    

    println("Parameterstudy starts: $total_runs runs on $(Threads.nthreads()) parallel Threads")
    println("Skipped $skipped_runs invalid runs (HeatPump + ElectrodeBoiler <= 5 MW).")
    println("Input file: $(abspath(base_input_path))")

    out_file_path = outdir * "/results_$(total_runs)runs_" * Dates.format(now(), "yymmdd_HHMMSS") * ".csv"
    touch(out_file_path)
    #component parameters
    header_parameters = ["HeatPump_Power / W", "ElectrodeBoiler_Power / W", "BufferTank_Capacity / Wh", 
                         "Battery_Capacity / Wh"]
    # no cost escalation
    header_annuity_no = ["annuity_no A_cap / €", "annuity_no A_cap_incentive / €", 
                         "annuity_no A_misc / €", "annuity_no A_op / €", 
                         "annuity_no A_energy / €", "annuity_no A_rev_control / €", 
                         "annuity_no A_rev_feed / €", "annuity_no A_total / €", 
                         "annuity_no A_total_incentive / €"]
    # moderate cost escalation
    header_annuity_mod = ["annuity_mod A_cap / €", "annuity_mod A_cap_incentive / €", 
                          "annuity_mod A_misc / €", "annuity_mod A_op / €", 
                          "annuity_mod A_energy / €", "annuity_mod A_rev_control / €", 
                          "annuity_mod A_rev_feed / €", "annuity_mod A_total / €", 
                          "annuity_mod A_total_incentive / €"]
    # progressive cost escalation
    header_annuity_pro = ["annuity_pro A_cap / €", "annuity_pro A_cap_incentive / €", 
                          "annuity_pro A_misc / €", "annuity_pro A_op / €", 
                          "annuity_pro A_energy / €", "annuity_pro A_rev_control / €", 
                          "annuity_pro A_rev_feed / €", "annuity_pro A_total / €", 
                          "annuity_pro A_total_incentive / €"]
    vdi_annuity_keys = [
        "A_cap", "A_cap_incentive", "A_misc", "A_op", "A_energy",
        "A_rev_control", "A_rev_feed", "A_total", "A_total_incentive"
    ]
    # balance warnings
    header_balances = ["balance_power", "balance_heat", "Errors"]
    # yearly_CO2-emissions
    header_co2 = ["CO2_yearly / t/a"]
    header = join(vcat(header_parameters, header_annuity_no, header_annuity_mod, 
                       header_annuity_pro, header_co2, header_balances), ';') * "\n"

    open(out_file_path, "a") do file_handle
        write(file_handle, header)
    end

    println("Results file at $out_file_path created.")

    output_lock = ReentrantLock()

    runidx_global = Threads.Atomic{Int}(1)
    erridx_global = Threads.Atomic{Int}(0)

    logger = start_logger(true, false, nothing, nothing,
                          Logging.Warn, nothing)
    Logging.global_logger(logger)

    Threads.@threads for (Pth_HeatPump, Pth_ElectrodeBoiler, Cap_Wh, BattCap_Wh) in valid_combinations
        
        runidx = Threads.atomic_add!(runidx_global, 1)
        run_ID = uuid4()
        sim_output = OrderedDict()
        start_time = now()
        ####################################################
        # Simulation
        ####################################################
        input_file, 
        profile_values, 
        timestep = create_variant(outdir, base_input, 
                                  Pth_HeatPump, Pth_ElectrodeBoiler, Cap_Wh, BattCap_Wh,
                                  profile_addons, profile_multipliers,
                                  run_ID; 
                                  write_output=write_output)

        println("[$runidx/$total_runs] → start simulation: $(basename(input_file))")

        raw_sim = nothing
        try
            raw_sim = run_resie_variant(input_file, run_ID, profile_values, profile_addons)

            sim_output[runidx] = Dict{String, Any}("sim" => raw_sim)
            sim_output[runidx]["Balance_heat"] = sum(raw_sim["BUS_Heat Balance"])
            sim_output[runidx]["Balance_power"] = sum(raw_sim["BUS_Power Balance"])
            sim_output[runidx]["time_step_seconds"] = timestep
            sim_output[runidx]["Errors"] = ""
        catch e
            sim_output[runidx] = Dict{String, Any}("sim" => missing)
            sim_output[runidx]["Balance_heat"] = missing
            sim_output[runidx]["Balance_power"] = missing
            sim_output[runidx]["time_step_seconds"] = missing

            # save excact error message to output file
            error_message = sprint(showerror, e)
            full_error_message = error_message * "\n" * sprint(Base.show_backtrace, catch_backtrace())
            println(full_error_message)
            sim_output[runidx]["Errors"] =  "\"" * replace(full_error_message, "\"" => "\"\"") * "\"\n"
            Threads.atomic_add!(erridx_global, 1)

            # save input file that threw error seperately
            error_path = joinpath(outdir, "error_inputfiles")
            mkpath(error_path)
            cp(input_file, joinpath(error_path, splitpath(input_file)[end]), force=true)
        end
        if !save_input_files 
            rm(input_file) 
        end

        ####################################################
        # Econcomy based on VDI 2067 principles
        ####################################################
        # function for component investment costs based on installed capacity (EUR/kW)
        # TODO Adjust factors in front (EUR/kW or EUR/kWh)
        A0_HeatPump     = 700 * (Pth_HeatPump / 1e3)        #   alter Wert 700
        A0_ElectrodeBoiler = 285  * (Pth_ElectrodeBoiler / 1e3)   # circa 1 Mio. € per 3.5 MW #TODO Christian fragen  # alter Wert: 285
        A0_Buffer = 39  * (Cap_Wh / 1e3)        # from FACT document
        A0_Batt   = 375 * (BattCap_Wh / 1e3)    # TODO Jule fragen
       
        components = VDI2067.VDIComponent[
            VDI2067.heatpump_component(A0_HeatPump),
            VDI2067.ElectrodeBoiler_component(A0_ElectrodeBoiler),
            VDI2067.buffertank_component(A0_Buffer),
            VDI2067.battery_component(A0_Batt)
        ]

        if !isnothing(raw_sim)
            sim_output[runidx]["VDI_NO"] =
                Base.invokelatest(VDI2067.vdi2067_annuity, raw_sim, components, VDI2067.VDI_SCENARIO_NONE)

            sim_output[runidx]["VDI_MOD"] =
                Base.invokelatest(VDI2067.vdi2067_annuity, raw_sim, components, VDI2067.VDI_SCENARIO_MOD)

            sim_output[runidx]["VDI_PRO"] =
                Base.invokelatest(VDI2067.vdi2067_annuity, raw_sim, components, VDI2067.VDI_SCENARIO_PRO)

            # get values for writing into output file
            annuities_no = [sim_output[runidx]["VDI_NO"][k] for k in vdi_annuity_keys]
            annuities_mod = [sim_output[runidx]["VDI_MOD"][k] for k in vdi_annuity_keys]
            annuities_pro = [sim_output[runidx]["VDI_PRO"][k] for k in vdi_annuity_keys]
            co2_yearly = sim_output[runidx]["VDI_NO"]["CO2_yearly"]
        else
            annuities_no = fill(missing, length(header_annuity_no))
            annuities_mod = fill(missing, length(header_annuity_mod))
            annuities_pro = fill(missing, length(header_annuity_pro))
            co2_yearly = missing
        end
        # write important results to seperate file
        parameters = [Pth_HeatPump, Pth_ElectrodeBoiler, Cap_Wh, BattCap_Wh]
        balances = collect(getindex.(Ref(sim_output[runidx]), ("Balance_power", "Balance_heat", "Errors")))
        row = join(vcat(parameters, annuities_no, annuities_mod, annuities_pro, co2_yearly, balances), ';') * "\n"
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
    println("✔ Parameterstudy completed with $(erridx_global[]) Errors.")
end

Base.invokelatest(main, base_input_path, false, true)
