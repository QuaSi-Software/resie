module Resie

using Printf
using Dates: now, seconds
using UUIDs

"""
Contains the parameters, instantiated components and the order of operations for a simulation run.

In short, it bundles all the data required to perform a simulation from start to finish.
Through its fields it contains a lot of complexity in terms of hierarchical data structures.
However the struct should not be used as argument for complex calculations. Instead it is
intended to be used for the run registry only.

Because the types used to represent the field are defined in modules that depend upon the
definition of the struct in the first place, the fields are defined with the generic 'Any'
type.
"""
mutable struct SimulationRun
    parameters::Dict{String,Any}
    components::Dict{String,Any}
    order_of_operations::Vector{Any}
end

# this registry should be the only global state in the package and contains the state for
# ongoing or paused simulation runs
current_runs::Dict{UUID,SimulationRun} = Dict{UUID,SimulationRun}()

"""
    get_run(id)

Get the simulation run container for the given ID.

# Args
- `id::UUID`: The ID of the run
# Returns
- `SimulationRun`: The simulation run container
"""
function get_run(id::UUID)::SimulationRun
    return current_runs[id]
end

"""
    close_run(id)

Closes the given run, removing it from the run registry.

# Args
- `id:UUID`: The ID of the run
"""
function close_run(id::UUID)
    delete!(current_runs, id)
end

"""
Custom exception `InputError` used to signify that an input was not correctly set up,
outside the allowed range, etc.
Call with `throw(InputError("msg"))` or `throw(InputError())`.
"""
struct InputError <: Exception
    msg::Union{AbstractString,Nothing}
end
InputError() = InputError(nothing)

# note: includes that contain their own module, which have to be submodules of the Resie
# module, are included first, then can be accessed with the "using" keyword. files that
# contain code that is intended to be used in-place of their include statement (as part
# of the Resie module), are included after the "using" statements have been declared.
# this is done so the latter files can access the symbols of the submodules the same as
# if the code was inside this file.

include("profiles/solar_irradiance.jl")
using .SolarIrradiance
include("profiles/base.jl")
using .Profiles
include("profiles/weatherdata.jl")
using .Weatherdata

include("energy_systems/base.jl")
using .EnergySystems

include("project_loading.jl")
include("file_output.jl")

include("resie_logger.jl")
using .Resie_Logger

using PlotlyJS
using ColorSchemes
using Colors
using Interpolations
using JSON
using Dates

const HOURS_PER_SECOND::Float64 = 1.0 / 3600.0
const SECONDS_PER_HOUR::Float64 = 3600.0

"""
    get_simulation_params(project_config)

Constructs the dictionary of simulation parameters.

# Arguments
-`project_config::Dict{AbstractString,Any}`: The project config
# Returns
-`Dict{String,Any}`: The simulation parameter dictionary
"""
function get_simulation_params(project_config::AbstractDict{AbstractString,Any})::Dict{String,Any}
    time_step, start_date, start_date_output, end_date, nr_of_steps, nr_of_steps_output = get_timesteps(project_config["simulation_parameters"])

    sim_params = Dict{String,Any}(
        "time" => 0,
        "time_since_output" => 0,
        "current_date" => start_date,
        "time_step_seconds" => time_step,
        "number_of_time_steps" => nr_of_steps,
        "number_of_time_steps_output" => nr_of_steps_output,
        "start_date" => start_date,
        "start_date_output" => start_date_output,
        "end_date" => end_date,
        "epsilon" => default(project_config["simulation_parameters"], "epsilon", 1e-9),
        "latitude" => default(project_config["simulation_parameters"], "latitude", nothing),
        "longitude" => default(project_config["simulation_parameters"], "longitude", nothing),
        "timezone" => default(project_config["simulation_parameters"], "time_zone", nothing),
        "step_info_interval" => default(project_config["io_settings"],
                                        "step_info_interval",
                                        Integer(floor(nr_of_steps / 20))),
        "force_profiles_to_repeat" => default(project_config["simulation_parameters"], "force_profiles_to_repeat",
                                              false),
    )

    # add helper functions to convert power to work and vice-versa. this uses the time step
    # of the simulation as the duration required for the conversion.
    sim_params["watt_to_wh"] = function (watts::Float64)
        return watts * time_step * HOURS_PER_SECOND
    end
    sim_params["wh_to_watts"] = function (wh::Float64)
        return wh * SECONDS_PER_HOUR / time_step
    end

    # add helper function for using paths, absolute or relative to the run base path
    if haskey(project_config["io_settings"], "base_path")
        run_base_path = abspath(project_config["io_settings"]["base_path"])
    else
        run_base_path = abspath(joinpath(dirname(@__FILE__), ".."))
    end
    sim_params["run_path"] = function (path)
        return isabspath(path) ? path : abspath(joinpath(run_base_path, path))
    end

    # load weather profiles accesible for all components
    weather_file_path = default(project_config["simulation_parameters"],
                                "weather_file_path",
                                nothing)

    if weather_file_path !== nothing
        weather_interpolation_type_solar = default(project_config["simulation_parameters"],
                                                   "weather_interpolation_type_solar", "linear_solar_radiation")
        weather_interpolation_type_general = default(project_config["simulation_parameters"],
                                                     "weather_interpolation_type_general", "linear_classic")
        # WeatherData() writes the lat and long to sim_params if they are not given in the input file
        sim_params["weather_data"] = WeatherData(weather_file_path,
                                                 sim_params,
                                                 weather_interpolation_type_solar,
                                                 weather_interpolation_type_general)
    end

    return sim_params
end

"""
    prepare_inputs(project_config)

Construct and prepare parameters, energy system components and the order of operation.

# Arguments
-`project_config::Dict{AbstractString,Any}`: The project config
# Returns
-`Dict{String,Any}`: Simulation parameters
-`Grouping`: The constructed energy system components
-`OrderOfOperations`: Order of operations
"""
function prepare_inputs(project_config::AbstractDict{AbstractString,Any}, run_ID::UUID)
    sim_params = get_simulation_params(project_config)
    sim_params["run_ID"] = run_ID

    components = load_components(project_config["components"], sim_params)

    if haskey(project_config, "order_of_operation") && length(project_config["order_of_operation"]) > 0
        operations = load_order_of_operations(project_config["order_of_operation"], components)
        @info "The order of operations was successfully imported from the input file.\n" *
              "Note that the order of operations has a major impact on the simulation " *
              "result and should only be changed by experienced users!"
    else
        operations = calculate_order_of_operations(components)
    end

    return sim_params, components, operations
end

"""
    run_simulation_loop()

Performs the simulation as loop over time steps and records outputs.

# Arguments
-`project_config::Dict{AbstractString,Any}`: The project config
-`sim_params::Dict{String,Any}`: Simulation parameters
-`components::Grouping`: The energy system components
-`operations::OrderOfOperations`:: Order of operations
"""
function run_simulation_loop(project_config::AbstractDict{AbstractString,Any},
                             sim_params::Dict{String,Any},
                             components::Grouping,
                             operations::OrderOfOperations)
    # get list of requested output keys for lineplot and csv export
    output_keys_lineplot, output_keys_to_CSV = get_output_keys(project_config["io_settings"], components)
    weather_data_keys = get_weather_data_keys(sim_params)
    do_create_plot_data = output_keys_lineplot !== nothing
    do_create_plot_weather = weather_data_keys !== nothing &&
                             default(project_config["io_settings"], "plot_weather_data", false)
    do_write_CSV_weather = weather_data_keys !== nothing &&
                           default(project_config["io_settings"], "csv_output_weather", false)
    weather_CSV_keys = do_write_CSV_weather ? weather_data_keys : nothing
    do_write_CSV = output_keys_to_CSV !== nothing || do_write_CSV_weather
    csv_output_file_path = default(project_config["io_settings"],
                                   "csv_output_file",
                                   "./output/out.csv")
    csv_time_unit = default(project_config["io_settings"], "csv_time_unit", "seconds")
    if !(csv_time_unit in ["seconds", "minutes", "hours", "date"])
        @error "The `csv_time_unit` has to be one of: seconds, minutes, hours, date!"
        throw(InputError())
    end

    # Initialize the array for output plots
    output_data_lineplot = do_create_plot_data ?
                           zeros(Float64, sim_params["number_of_time_steps_output"], 1 + length(output_keys_lineplot)) :
                           nothing
    output_weather_lineplot = do_create_plot_weather ?
                              zeros(Float64, sim_params["number_of_time_steps_output"], 1 + length(weather_data_keys)) :
                              nothing

    # reset CSV file
    if do_write_CSV
        reset_file(sim_params["run_path"](csv_output_file_path),
                   output_keys_to_CSV,
                   weather_CSV_keys,
                   csv_time_unit)
    end

    # check if sankey should be plotted
    do_create_sankey = haskey(project_config["io_settings"], "sankey_plot") &&
                       project_config["io_settings"]["sankey_plot"] !== "nothing"
    if do_create_sankey
        # get information about all interfaces for Sankey
        nr_of_interfaces,
        medium_of_interfaces,
        output_sourcenames_sankey,
        output_targetnames_sankey = get_interface_information(components)
        # preallocate for speed: Matrix with data of interfaces in every timestep
        output_interface_values = zeros(Float64, sim_params["number_of_time_steps_output"], nr_of_interfaces)
    end

    # export order of operation and other additional info like optional plots
    dump_auxiliary_outputs(project_config, components, operations, sim_params)

    @info "-- Start time step loop"
    if sim_params["start_date_output"] == sim_params["start_date"]
        @info "-- No preheating activated. Starting output from beginning."
    else
        @info "-- Starting with preheating, no output will be written until completed."
    end
    start = now()
    for steps in 1:sim_params["number_of_time_steps"]
        # check if data should be output
        do_output = sim_params["current_date"] >= sim_params["start_date_output"]
        output_steps = Int(max(1,
                               Int(steps) - (Int(sim_params["number_of_time_steps"]) -
                                             Int(sim_params["number_of_time_steps_output"]))))
        if output_steps == 1 && do_output && sim_params["start_date_output"] != sim_params["start_date"]
            @info "-- Preheating completed. Starting output now."
        end

        operations_adjusted = reorder_operations_in_time_step(components, operations, sim_params)
        perform_operations(components, operations_adjusted, sim_params)

        if do_output
            # check if any component was not balanced
            warnings = check_balances(components, sim_params["epsilon"])
            if length(warnings) > 0
                for (key, balance) in warnings
                    @balanceWarn "Balance for component $key was not zero in timestep " *
                                 "$(sim_params["time_since_output"]): $balance"
                end
            end

            # write requested output data of the components to CSV-file
            # This is currently done in every time step to keep data even if 
            # an error occurs.
            if do_write_CSV
                write_to_file(sim_params["run_path"](csv_output_file_path),
                              output_keys_to_CSV,
                              weather_CSV_keys,
                              sim_params,
                              csv_time_unit)
            end

            # get the energy transported through each interface in every timestep for Sankey
            if do_create_sankey
                output_interface_values[output_steps, :] = collect_interface_energies(components, nr_of_interfaces)
            end

            # gather output data of each component for line plot
            if do_create_plot_data
                output_data_lineplot[output_steps, :] = gather_output_data(output_keys_lineplot,
                                                                           sim_params["time_since_output"])
            end
            if do_create_plot_weather
                output_weather_lineplot[output_steps, :] = gather_weather_data(weather_data_keys, sim_params)
            end

            # simulation update
            sim_params["time_since_output"] += Int(sim_params["time_step_seconds"])
        end

        # simulation update
        sim_params["time"] += Int(sim_params["time_step_seconds"])
        sim_params["current_date"] = add_ignoring_leap_days(sim_params["current_date"],
                                                            Second(sim_params["time_step_seconds"]))
        # progress report
        if sim_params["step_info_interval"] > 0 && steps % sim_params["step_info_interval"] == 0
            eta = ((sim_params["number_of_time_steps"] - steps)
                   *
                   seconds(now() - start)
                   /
                   sim_params["step_info_interval"])
            @info "Progress: $(steps)/$(sim_params["number_of_time_steps"]) in" *
                  " $(seconds(now() - start)) s. ETA: $(@sprintf("%.4f", eta)) s"
            start = now()
        end
    end
    @info "-- Finished time step loop"

    # create profile line plot
    if do_create_plot_data || do_create_plot_weather
        create_profile_line_plots(output_data_lineplot,
                                  output_keys_lineplot,
                                  output_weather_lineplot,
                                  weather_data_keys,
                                  project_config,
                                  sim_params)
        filepath = default(project_config["io_settings"],
                           "output_plot_file",
                           "./output/output_plot.html")
        @info "Line plot created and saved to $(sim_params["run_path"](filepath))"
    end

    # create Sankey diagram
    if do_create_sankey
        create_sankey(output_sourcenames_sankey,
                      output_targetnames_sankey,
                      output_interface_values,
                      medium_of_interfaces,
                      nr_of_interfaces,
                      project_config["io_settings"],
                      sim_params)
        filepath = default(project_config["io_settings"],
                           "sankey_plot_file",
                           "./output/output_sankey.html")
        @info "Sankey created and saved to $(sim_params["run_path"](filepath))"
    end

    if do_write_CSV
        @info "CSV-file with outputs written to $(sim_params["run_path"](csv_output_file_path))"
    end

    # plot additional figures potentially available from components after simulation
    if default(project_config["io_settings"], "auxiliary_plots", false)
        component_list = []
        output_path = default(project_config["io_settings"], "auxiliary_plots_path", "./output/")
        for component in components
            if plot_optional_figures_end(component[2], sim_params, output_path)
                push!(component_list, component[2].uac)
            end
        end
        if length(component_list) > 0
            @info "(Further) auxiliary plots are saved to folder " *
                  "$(sim_params["run_path"](output_path)) for the following components: " *
                  "$(join(component_list, ", "))"
        end
    end
end

"""
    load_and_run(filepath)

Load a project from the given file and run the simulation with it.

# Arguments
- `filepath::String`: Filepath to the project config file.
- `UUID`: The run ID used in the run registry
# Returns
- `Bool`: `true` if the simulation was successful, `false` otherwise.
"""
function load_and_run(filepath::String, run_ID::UUID)::Bool
    start = now()
    @info "---- Simulation setup ----"
    @info "-- Starting simulation at $(start)"
    @info "-- Now reading project config"

    project_config = nothing

    try
        # we can't use run_path() here because that is only defined after loading the config
        # in the first place. so we hope we've been given a valid path and forbid upwards
        # path traversal for security reasons
        if occursin("..", filepath)
            @error "Project config filepath must not contain .. path traversal"
            return false
        end
        project_config = read_JSON(abspath(filepath))
    catch exc
        if isa(exc, MethodError)
            @error "Could not parse project config file at $(abspath(filepath))"
            return false
        end
    end

    if project_config === nothing
        @error "Could not find or parse project config file at $(abspath(filepath))"
        return false
    end

    @info "-- Now preparing inputs"
    if !isnothing(project_config["optimizer"])
        create_output_file()
        if project_config["optimizer"] == "parameterstudy"
            opt_params = project_config["optimizer"]["optim_params"]

            output_lock = ReentrantLock()
            runidx_global = Threads.Atomic{Int}(1)
            erridx_global = Threads.Atomic{Int}(0)

            Threads.@threads for sample_params in collect(Iterators.product(opt_params))
                run_sample(outdir, project_config, sample_params)
            end
        elseif project_config["optimizer"] == "monte_carlo_anealing"
            objective_function = # ???
            optimize(run_single_sim)
        end

    else
        run_single_sim(project_config_run, run_ID)
    end
    return true
end

function run_single_sim(project_config, run_ID)
    sim_params, components, operations = prepare_inputs(project_config, run_ID)
    current_runs[run_ID] = SimulationRun(sim_params, components, operations)
    @info "-- Simulation setup complete in $(seconds(now() - start)) s"

    start = now()
    @info "---- Simulation loop ----"
    sim_output = run_simulation_loop(project_config, sim_params, components, operations)
    
    if !isnothing(project_config["economy"]) || !isnothing(project_config["emissions"])
        years_sim = sim_params["number_of_time_steps_output"] * sim_params["time_step_seconds"] / 3600*8760
        years_economy_data = length(price_profile) / sim_params["time_step_seconds"]
        if years < 1
            # repeat whole sim_output until 1 year full
            sim_output
        elseif  years > 1
            # repeat last year until project_config["economy"]["economy_parms"]["eval_time_a"]
            sim_output
        end
        if !isnothing(project_config["economy"])
            annuities, total_costs = calc_economy(sim_output) # annuity for each value for each component
        end
        if !isnothing(project_config["emissions"])
            # analog zu economy
        end
    end
    @info "-- Simulation loop complete in $(seconds(now() - start)) s"
    return sim_output, annuities, total_costs, emissions
end

function monte_carlo_anealing()
    # combined monte carlo and simulated annealing
    # the temperature determines if a completely random or existing sample is used as starting
    # point, determines the size of the neighborhood and the number of results (sorted by)
    # global measure, from which a new sample is drawn
    for idx in range(1, NR_TRIES)
        # temperature schedule is simple inverse logistic curve
        temperature = 1.0 - 1.0 / (1.0 + exp(-8.0 * (idx / NR_TRIES - 0.5)))

        if length(all_results) == 0 || rand() < temperature
            # set parameters to equally distributed random values across whole parameter space
            hp_power = rand(BOUNDS["TST_HP_01"][2]:BOUNDS["TST_HP_01"][3])
            pv_power = rand(BOUNDS["TST_PV_01"][2]:BOUNDS["TST_PV_01"][3])
            bat_cap = rand(BOUNDS["TST_BAT_01"][2]:BOUNDS["TST_BAT_01"][3])
            bt_cap = rand(BOUNDS["TST_BFT_01"][2]:BOUNDS["TST_BFT_01"][3])
            sample_params =

            print(". ")
        else
            # set parameters to neighborhood of existing result, drawn from the top results
            # by global measure, where temperature determines the results pool and size of
            # neighborhood
            sample_idx = rand(1:max(1, Int(round(length(all_results) * temperature))))
            print("$sample_idx ")
            sample = sample_idx >= 1 && sample_idx <= length(all_results) ? all_results[sample_idx] : all_results[1]

            range = NBH_SCALE * temperature * (BOUNDS["TST_HP_01"][3] - BOUNDS["TST_HP_01"][2])
            hp_power = sample["hp_power"] + rand((-0.5 * range):(0.5 * range))
            hp_power = max(BOUNDS["TST_HP_01"][2], min(BOUNDS["TST_HP_01"][3], hp_power))

            range = NBH_SCALE * temperature * (BOUNDS["TST_PV_01"][3] - BOUNDS["TST_PV_01"][2])
            pv_power = sample["pv_power"] + rand((-0.5 * range):(0.5 * range))
            pv_power = max(BOUNDS["TST_PV_01"][2], min(BOUNDS["TST_PV_01"][3], pv_power))

            range = NBH_SCALE * temperature * (BOUNDS["TST_BAT_01"][3] - BOUNDS["TST_BAT_01"][2])
            bat_cap = sample["bat_cap"] + rand((-0.5 * range):(0.5 * range))
            bat_cap = max(BOUNDS["TST_BAT_01"][2], min(BOUNDS["TST_BAT_01"][3], bat_cap))

            range = NBH_SCALE * temperature * (BOUNDS["TST_BFT_01"][3] - BOUNDS["TST_BFT_01"][2])
            bt_cap = sample["bt_cap"] + rand((-0.5 * range):(0.5 * range))
            bt_cap = max(BOUNDS["TST_BFT_01"][2], min(BOUNDS["TST_BFT_01"][3], bt_cap))

            sample_params = 
        end

        # run sim and record results
        run_results = run_sample(outdir, project_config, sample_params)
        obj_results = objective_function(run_results)
        result = Dict(
            "gm" => 0.0,
            "cost" => obj_results.lcc,
            "emissions" => obj_results.emissions,
            "hp_power" => hp_power,
            "hp_capex_per_kw" => inputs["components"]["TST_HP_01"]["capex_per_kW"],
            "pv_power" => pv_power,
            "pv_capex_per_kw" => inputs["components"]["TST_PV_01"]["capex_per_kW"],
            "bat_cap" => bat_cap,
            "bat_capex_per_kwh" => inputs["components"]["TST_BAT_01"]["capex_per_kWh"],
            "bt_cap" => bt_cap,
            "bt_capex_per_kwh" => inputs["components"]["TST_BFT_01"]["capex_per_kWh"],
            "base_cop" => inputs["components"]["TST_HP_01"]["base_cop"],
            "base_grid_price" => inputs["components"]["TST_BAT_01"]["base_grid_price"],
        )
        push!(all_results, result)
        min_lcc = obj_results.lcc < min_lcc ? obj_results.lcc : min_lcc
        min_emissions = obj_results.emissions < min_emissions ? obj_results.emissions : min_emissions

        # calculate global measure and sort by it
        for res in all_results
            res["gm"] = sqrt(((res["cost"]) / min_lcc - 1.0)^2 + (res["emissions"] / min_emissions - 1.0)^2)
        end
        sort!(all_results; by=x -> x["gm"])

        if idx % 5 == 0
            print("| ")
        end
    end
end

function run_sample(outdir, project_config, sample_params)
    runidx = Threads.atomic_add!(runidx_global, 1)
    run_ID = uuid4()
    combined_output = OrderedDict()
    start_time = now()
    ####################################################
    # Simulation
    ####################################################
    # TODO input_file wird nicht neu geschrieben, sondern nur profile
    project_config_run = create_variant(outdir, project_config, sample_params,
                                        run_ID; 
                                        write_output=write_output)


    @info("[$runidx/$total_runs] → start simulation: $(basename(input_file))")

    raw_sim = nothing
    try
        raw_sim, annuities, total_costs, emissions = run_single_sim(project_config_run, run_ID)

    catch e
        # save excact error message to output file
        error_message = sprint(showerror, e)
        full_error_message = error_message * "\n" * sprint(Base.show_backtrace, catch_backtrace())
        @info full_error_message
        combined_output[runidx]["Errors"] =  "\"" * replace(full_error_message, "\"" => "\"\"") * "\"\n"
        Threads.atomic_add!(erridx_global, 1)

        # save input file that threw error seperately
        error_path = joinpath(outdir, "error_inputfiles")
        mkpath(error_path)
        cp(input_file, joinpath(error_path, splitpath(input_file)[end]), force=true)
    end
    if !save_input_files 
        rm(input_file) 
    end
    return raw_sim, annuities, total_costs, emissions
end

end # module
