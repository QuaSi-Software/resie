# this file contains functionality pertaining to loading a project's metadata and the
# energy system components from the project config file, as well as constructing certain
# helpful information data structures from the inputs in the config
using JSON: JSON
using OrderedCollections: OrderedDict

const HOURS_PER_SECOND::Float64 = 1.0 / 3600.0
const SECONDS_PER_HOUR::Float64 = 3600.0

#! format: off
const IO_SETTINGS_DEF = Dict{String,Any}(
    "base_path" => (
        default=nothing,
        description="If given, this path will be used as the base path for all relative " *
                    "paths used in the config file. If not given it defaults to the " *
                    "current working directory for the Julia process running ReSiE, which " *
                    "in almost all cases is the directory from which ReSiE is started.",
        display_name="Base path",
        required=false,
        type=String,
        json_type="string",
        unit="-"
    ),
    "plot_weather_data" => (
        default=false,
        description="Toggle if the weather data read in from the given weather file " *
                    "should be included in the line plot",
        display_name="Plot weather data?",
        required=false,
        type=Bool,
        json_type="boolean",
        unit="-"
    ),
    "csv_output_weather" => (
        default=false,
        description="Toggle if the weather data read in from the given weather file " *
                    "should be included in the CSV output",
        display_name="Weather data in CSV?",
        required=false,
        type=Bool,
        json_type="boolean",
        unit="-"
    ),
    "write_csv_continuously" => (
        default=false,
        description="Toggle if CSV output will be written continuously, meaning in every " *
                    "time step. Activating this functionality will ensure partial output " *
                    "if the simulation fails during execution, however it also incurs a " *
                    "substantial performance penalty due to frequent file access.",
        display_name="Write CSV continously?",
        required=false,
        type=Bool,
        json_type="boolean",
        unit="-"
    ),
    "write_summary_CSV" => (
        default=true,
        description="Toggle if a CSV summary output with sum/mean values should be created " *
                    "additionally to the timestep-wise CSV output.",
        display_name="Write summary CSV?",
        required=false,
        type=Bool,
        json_type="boolean",
        unit="-"
    ),
    "csv_output_file" => (
        default="./output/out.csv",
        description="File path to where the CSV output will be written",
        display_name="CSV output file",
        required=false,
        type=String,
        json_type="string",
        unit="-"
    ),
    "csv_time_unit" => (
        default="date",
        description="Time unit for the timestamp of the CSV output",
        display_name="CSV time unit",
        required=false,
        options=["seconds", "minutes", "hours", "date"],
        type=String,
        json_type="string",
        unit="-"
    ),
    "output_plot_file" => (
        default="./output/output_plot.html",
        description="File path to where the output line plot will be written",
        display_name="Output plot file",
        required=false,
        type=String,
        json_type="string",
        unit="-"
    ),
    "output_plot_time_unit" => (
        default="date",
        description="Unit for x-axis of the output plot. Note that the plotted energies " *
                    "always refer to the simulation time step and not to the unit " *
                    "specified here!",
        display_name="Output plot time unit",
        required=false,
        options=["seconds", "minutes", "hours", "date"],
        type=String,
        json_type="string",
        unit="-"
    ),
    "sankey_plot_file" => (
        default="./output/output_sankey.html",
        description="File path to where the Sankey plot will be written",
        display_name="Sankey plot file",
        required=false,
        type=String,
        json_type="string",
        unit="-"
    ),
    "show_detailed_errors" => (
        default=false,
        description="Toggle if errors should show a more detailed message. Only affects " *
                    "some errors.",
        display_name="Show detailed errors?",
        required=false,
        type=Bool,
        json_type="boolean",
        unit="-"
    ),
    "auxiliary_info" => (
        default=false,
        description="Toggle if auxiliary info about the current run should be written to " *
                    "markdown file",
        display_name="Write auxiliary info?",
        required=false,
        type=Bool,
        json_type="boolean",
        unit="-"
    ),
    "auxiliary_info_file" => (
        default="./output/auxiliary_info.md",
        description="File path to where the auxiliary information will be written",
        display_name="Auxiliary info file",
        required=false,
        type=String,
        json_type="string",
        unit="-"
    ),
    "auxiliary_plots" => (
        default=false,
        description="Toggle if additional plots of components, if they are available, " *
                    "are created",
        display_name="Create auxiliary plots?",
        required=false,
        type=Bool,
        json_type="boolean",
        unit="-"
    ),
    "auxiliary_plots_path" => (
        default="./output/",
        description="Directory path to where the additional plots will be saved",
        display_name="Auxiliary plots path",
        required=false,
        type=String,
        json_type="string",
        unit="-"
    ),
    "auxiliary_plots_formats" => (
        default=["png"],
        description="Multiple selection of which file formats are used to create the " *
                    "auxiliary plots. Allowed formats are: html, pdf, png, ps, svg",
        display_name="Auxiliary plots formats",
        required=false,
        type=Vector{String},
        json_type="list",
        unit="-"
    ),
    "plot_economic_cashflows" => (
        default=true,
        description="Toggle if the economic results (as cashflows) should be plotted",
        display_name="Plot economic cashflows?",
        required=false,
        type=Bool,
        json_type="boolean",
        unit="-"
    ),
    "economic_plot_cashflows_file_path" => (
        default="./output/economic_results_cashflows.html",
        description="File path to where the economic cashflow plots are written",
        display_name="Economic cashflows plot file path",
        required=false,
        type=String,
        json_type="string",
        unit="-"
    ),
    "plot_economic_present_values" => (
        default=true,
        description="Toggle if the economic results (as present values) should be plotted",
        display_name="Plot economic present values?",
        required=false,
        type=Bool,
        json_type="boolean",
        unit="-"
    ),
    "economic_plot_present_values_file_path" => (
        default="./output/economic_results_present_values.html",
        description="File path to where the economic present value plots are written",
        display_name="Economic present value plot file path",
        required=false,
        type=String,
        json_type="string",
        unit="-"
    ),
    "output_economic_CSV" => (
        default=true,
        description="Toggle if a CSV with the economic results should be created",
        display_name="Output economic CSV?",
        required=false,
        type=Bool,
        json_type="boolean",
        unit="-"
    ),
    "economic_CSV_file_path" => (
        default="./output/economic_results.csv",
        description="File path to where the economic results are written to CSV",
        display_name="Economic results CSV file path",
        required=false,
        type=String,
        json_type="string",
        unit="-"
    ),
    "plot_emission_results" => (
        default=true,
        description="Toggle if the emission results should be plotted",
        display_name="Plot emission results?",
        required=false,
        type=Bool,
        json_type="boolean",
        unit="-"
    ),
    "emissions_plot_file_path" => (
        default="./output/emissions_plot.html",
        description="File path to where the emissions plot file is written",
        display_name="Emissions plot file path",
        required=false,
        type=String,
        json_type="string",
        unit="-"
    ),
    "output_emissions_CSV" => (
        default=true,
        description="Toggle if a CSV with the emission results should be created",
        display_name="Output emissions CSV?",
        required=false,
        type=Bool,
        json_type="boolean",
        unit="-"
    ),
    "emissions_CSV_file_path" => (
        default="./output/emissions_results.csv",
        description="File path to where the emissions are written to CSV",
        display_name="Emissions CSV file path",
        required=false,
        type=String,
        json_type="string",
        unit="-"
    ),
    "plot_price_and_emission_profiles" => (
        default=false,
        description="Toggle if a plot with the utilized price and emission profiles should be created",
        display_name="Plot price and emission profiles?",
        required=false,
        type=Bool,
        json_type="boolean",
        unit="-"
    ),
    "price_and_emission_profile_file_path" => (
        default="./output/price_and_emissions_profiles.html",
        description="File path to where the plot with price and emissions profiles should be written",
        display_name="Price and Emissions profiles plot file path",
        required=false,
        type=String,
        json_type="string",
        unit="-"
    ),
    "step_info_interval" => (
        default=nothing,
        description="Defines how often a progress report on the loop over the timesteps " *
                    "of the simulation is logged to the info channel. This is useful to " *
                    "get an estimation of how much longer the simulation requires " *
                    "(albeit that such estimation is always inaccurate). If no value is " *
                    "given, automatically sets a value such that 20 reports are printed " *
                    "over the course of the simulation. To deactivate these reports, set " *
                    "this to 0.",
        display_name="Step info interval",
        required=false,
        validations=[("self", "value_gte_num", 0.0)],
        type=Integer,
        json_type="int",
        unit="-"
    ),
    "sankey_plot" => (
        default="default",
        description="Sets the mode of the sankey plot output, switching between default " *
                    "and custom behaviour as well an option of not creating a sankey " *
                    "plot file.",
        display_name="Sankey plot mode",
        options=["default", "custom", "nothing"],
        required=false,
        type=String,
        json_type="string",
        unit="-"
    ),
    "sankey_plot_spec" => (
        default=nothing,
        description="Specification for the sankey plot in custom mode. See documentation " *
                    "for how this needs to be structured.",
        display_name="Sankey plot specification",
        required=false,
        conditionals=[("sankey_plot", "is", "custom")],
        type=Dict{String,Any},
        json_type="object",
        unit="-"
    ),
    "output_plot" => (
        default="all_incl_flows",
        description="Sets the mode of the output plot, switching between several default " *
                    "and custom behaviour modes as well an option of not creating a plot " *
                    "file at all.",
        display_name="Output plot mode",
        options=["custom", "all_excl_flows", "all_incl_flows", "nothing"],
        required=false,
        type=String,
        json_type="string",
        unit="-"
    ),
    "output_plot_spec" => (
        default=nothing,
        description="Specification for the output plot in custom mode. See documentation " *
                    "for how this needs to be structured.",
        display_name="Output plot specification",
        required=false,
        conditionals=[("output_plot", "is", "custom")],
        type=Dict{String,Any},
        json_type="object",
        unit="-"
    ),
    "csv_output" => (
        default="nothing",
        description="Sets the mode of the CSV output, switching between several default " *
                    "and custom behaviour modes as well an option of not creating a " *
                    "CSV file at all.",
        display_name="CSV output mode",
        options=["custom", "all_excl_flows", "all_incl_flows", "nothing"],
        required=false,
        type=String,
        json_type="string",
        unit="-"
    ),
    "csv_output_keys" => (
        default=nothing,
        description="Specification for the CSV output in custom mode. See documentation " *
                    "for how this needs to be structured.",
        display_name="CSV output specification",
        required=false,
        conditionals=[("csv_output", "is", "custom")],
        type=Dict{String,Any},
        json_type="object",
        unit="-"
    ),
    "fixed_output_precision" => (
        default=0,
        description="If given a non-zero value, uses this many digits as the fixed " *
                    "precision for float outputs in CSV and plot files. It is not " *
                    "recommended to use this setting. It is intended for making the " *
                    "output perfectly repeatable, which is useful for testing but not in " *
                    "normal simuation.",
        display_name="Fixed output precision",
        required=false,
        type=Integer,
        json_type="number",
        unit="-"
    ),
)

const SIMULATION_PARAMETERS_DEF = Dict{String,Any}(
    "start" => (
        description="Start time of the simulation as datetime format",
        display_name="Start datetime",
        required=true,
        type=String,
        json_type="string",
        unit="-"
    ),
    "start_output" => (
        default=nothing,
        description="The start time as datetime format at which the simulation begins " *
                    "to output the simulation results. Has to be equal or later than " *
                    "`start`. Can be used to perform heat-up simulation ahead of the " *
                    "actual simulation. Note that during heat-up, no warnings are " *
                    "printed. The energies in the various output files are starting at " *
                    "the time specified in `start_output`.",
        display_name="Start output",
        required=false,
        type=String,
        json_type="string",
        unit="-"
    ),
    "end" => (
        description="End time (inclusive) of the simulation as datetime format. Will be " *
                    "rounded down to the nearest multiple of the time step.",
        display_name="End datetime",
        required=true,
        type=String,
        json_type="string",
        unit="-"
    ),
    "start_end_unit" => (
        description="Datetime format specifier for parameters `start`, `start_output` " *
                    "and `end`",
        display_name="Start/end format",
        required=true,
        type=String,
        json_type="string",
        unit="-"
    ),
    "time_step" => (
        default=900,
        description="Time step for the simulation. The parameter `time_step_unit` " *
                    "determines what the value in `time_step` means.",
        display_name="Time step",
        required=false,
        validations=[("self", "value_gte_num", 1)],
        type=Integer,
        json_type="number",
        unit="-"
    ),
    "time_step_unit" => (
        default="seconds",
        description="Unit for the value given in parameter `time_step`.",
        display_name="Time step unit",
        options=["seconds", "minutes", "hours"],
        required=false,
        type=String,
        json_type="string",
        unit="-"
    ),
    "weather_file_path" => (
        default=nothing,
        description="File path to the project-wide weather file. Can either be an " *
                    "EnergyPlus Weather File (EPW, time step has to be one hour, without " *
                    "leap day or DST) or a .dat file from the DWD (see " *
                    "https://kunden.dwd.de/obt/, free registration is required). See the " *
                    "component parameters on how to link weather file data to a component.",
        display_name="Weather file path",
        required=false,
        type=String,
        json_type="string",
        unit="-"
    ),
    "weather_interpolation_type_general" => (
        default="linear_classic",
        description="Interpolation type for weather data from weather file, except for " *
                    "solar radiation data. See the documentation for more details.",
        display_name="Weather interpol. general",
        options=["stepwise", "linear_classic", "linear_time_preserving", "linear_solar_radiation"],
        required=false,
        conditionals=[("weather_file_path", "is_not_nothing")],
        type=String,
        json_type="string",
        unit="-"
    ),
    "weather_interpolation_type_solar" => (
        default="linear_solar_radiation",
        description="Interpolation method for solar radiation data from weather file. " *
                    "See the documentation for more details.",
        display_name="Weather interpol. solar",
        options=["stepwise", "linear_classic", "linear_time_preserving", "linear_solar_radiation"],
        required=false,
        conditionals=[("weather_file_path", "is_not_nothing")],
        type=String,
        json_type="string",
        unit="-"
    ),
    "latitude" => (
        default=nothing,
        description="The latitude of the location in WGS84. If given, it overwrites the " *
                    "coordinates read out of the weather file!",
        display_name="Latitude",
        required=false,
        validations=[
            ("self", "value_gte_num_or_nothing", -90.0),
            ("self", "value_lte_num_or_nothing", 90.0)
        ],
        type=Float64,
        json_type="number",
        unit="°"
    ),
    "longitude" => (
        default=nothing,
        description="The longitude of the location in WGS84. If given, it overwrites the " *
                    "coordinates read out of the weather file!",
        display_name="Latitude",
        required=false,
        validations=[
            ("self", "value_gte_num_or_nothing", -180.0),
            ("self", "value_lte_num_or_nothing", 180.0)
        ],
        type=Float64,
        json_type="number",
        unit="°"
    ),
    "time_zone" => (
        default=nothing,
        description="The time zone offset used in the current simulation in relation to " *
                    "UTC. If given, it overwrites the coordinates read out of the " *
                    "weather file! DWD-dat files are assumed to be in GMT+1.",
        display_name="Time zone offset",
        required=false,
        type=Float64,
        json_type="number",
        unit="h"
    ),
    "epsilon" => (
        default=1e-9,
        description="The absolute tolerance for many floating-point comparisons in the " *
                    "simulation. Two values whose difference falls below this threshold " *
                    "are treated as equal.",
        display_name="Epsilon",
        required=false,
        validations=[
            ("self", "value_gt_num", 0.0),
            ("self", "value_lte_num", 1e-4)
        ],
        type=Float64,
        json_type="number",
        unit=""
    ),
    "force_profiles_to_repeat" => (
        default=false,
        description="If set to true, all utilized profiles are allowed to be repeated, " *
                    "even if denied or not specified in the profile header! Attention: " *
                    "This parameter disables the profile parameter in the profile header!",
        display_name="Force profiles to repeat?",
        required=false,
        type=Bool,
        json_type="boolean",
        unit="-"
    ),
)

ECONOMIC_PARAMETERS_DEF = Dict{String,Any}(
    "calculate_economy" => (
        default=false,
        description="If set to true, performs the calculation of economic results, " *
                    "requiring the input parameters for all components.",
        display_name="Calculate economy?",
        required=false,
        type=Bool,
        json_type="boolean",
        unit="-"
    ),
    "observation_period_in_years" => (
        default=20,
        description="The period in consideration for the calculation of economic results.",
        display_name="Observation period",
        required=false,
        type=Int,
        json_type="number",
        unit="a"
    ),
    "interest_rate" => (
        default=0.02,
        description="Interest rate for the calculation of annuities.",
        display_name="Interest rate",
        required=false,
        type=Float64,
        json_type="number",
        unit="-"
    ),
    "labour_costs_per_hour" => (
        default=100,
        description="Cost of labour for the operation of components.",
        display_name="Labour costs",
        required=false,
        type=Float64,
        json_type="number",
        unit="€/h"
    ),
    "labour_costs_price_change_rate_per_year" => (
        default=0.035,
        description="Rate of change of labour costs per year.",
        display_name="Labour costs change rate",
        required=false,
        type=Float64,
        json_type="number",
        unit="-"
    ),
    "repeat_method" => (
        default="last_year",
        description="Defines which period of the result data is repeated to fill up the " *
                    "remainder of the observation period. This can be equal or less than " *
                    "the simulation period, for example simulating three years, but only " *
                    "repeating the last year for the entire observation period.",
        display_name="Repeat method",
        required=false,
        options=["all", "last_year", "last_month", "last_week"],
        type=String,
        json_type="string",
        unit="-"
    ),
)

EMISSIONS_PARAMATERS_DEF = Dict{String,Any}(
    "calculate_emissions" => (
        default=false,
        description="If set to true, performs the calculation of GHG emissions results, " *
                    "requiring the input parameters for all components.",
        display_name="Calculate emissions?",
        required=false,
        type=Bool,
        json_type="boolean",
        unit="-"
    ),
    "observation_period_in_years" => (
        default=20,
        description="The period in consideration for the calculation of GHG emissions.",
        display_name="Observation period",
        required=false,
        type=Int,
        json_type="number",
        unit="a"
    ),
    "include_embodied_emissions" => (
        default=true,
        description="If set to true, includes embodied emissions in the calculation.",
        display_name="Include embodied emissions?",
        required=false,
        type=Bool,
        json_type="boolean",
        unit="-"
    ),
    "repeat_method" => (
        default="last_year",
        description="Defines which period of the result data is repeated to fill up the " *
                    "remainder of the observation period. This can be equal or less than " *
                    "the simulation period, for example simulating three years, but only " *
                    "repeating the last year for the entire observation period.",
        display_name="Repeat method",
        required=false,
        options=["all", "last_year", "last_month", "last_week"],
        type=String,
        json_type="string",
        unit="-"
    ),
)

OUTPUT_SPECIFICATION_SETTINGS = [
    ("output_plot", "output_plot_spec"),
    ("sankey_plot", "sankey_plot_spec"),
    ("csv_output", "csv_output_keys"),
]
#! format: on

"""
    read_JSON(filepath)

Read and parse the JSON-encoded Dict in the given file.
"""
function read_JSON(filepath::String)::OrderedDict{AbstractString,Any}
    open(filepath, "r") do file_handle
        content = read(file_handle, String)
        return JSON.parse(content; dicttype=OrderedDict)
    end
end

"""
    get_io_settings(project_config::AbstractDict{String,Any})

Constructs the dictionary of IO settings from the given config, considering default values.

# Arguments
-`project_config::AbstractDict{String,Any}`: The project config
# Returns
-`Dict{String,Any}`: The IO settings dictionary
"""
function get_io_settings(project_config::AbstractDict{String,Any})::Dict{String,Any}
    # extract most parameter values using the extract function on the base type component
    # (even though we are not checking components...)
    io_settings = Dict{String,Any}()
    for (name, param_def) in pairs(IO_SETTINGS_DEF)
        if name in ["sankey_plot_spec", "output_plot_spec", "csv_output_keys"]
            continue
        end
        io_settings[name] = EnergySystems.extract_parameter(EnergySystems.Component, project_config["io_settings"],
                                                            name, param_def, Dict{String,Any}(), "IO settings")
    end

    # special parameters for output specifications
    for (mode, spec) in OUTPUT_SPECIFICATION_SETTINGS
        if io_settings[mode] == "custom"
            if haskey(project_config["io_settings"], spec)
                io_settings[spec] = Dict{String,Any}(pairs(project_config["io_settings"][spec]))
            else
                throw(InputError("IO setting $mode was set to custom, but no custom specification was provided"))
            end
        elseif mode == "sankey_plot" && io_settings[mode] == "default"
            io_settings[spec] = Dict{String,Any}(
                "m_h_w_lt1" => "steelblue1",
                "m_h_w_lt2" => "steelblue1",
                "m_h_w_ht1" => "darkred",
                "m_e_ac_230v" => "darkgoldenrod1",
                "m_c_g_natgas" => "purple3",
                "m_c_g_h2" => "green3",
                "m_c_g_o2" => "firebrick1",
                "Losses" => "grey40",
                "Gains" => "grey40",
            )
        end
    end

    # special handling for some other parameters
    if haskey(io_settings, "step_info_interval") && isnothing(io_settings["step_info_interval"])
        delete!(io_settings, "step_info_interval")
    end

    if haskey(io_settings, "base_path") && !isnothing(io_settings["base_path"])
        io_settings["base_path"] = abspath(io_settings["base_path"])
    else
        io_settings["base_path"] = abspath(joinpath(dirname(@__FILE__), ".."))
    end

    # again, we use the component validation function as it avoids duplicate code
    EnergySystems.validate_config(EnergySystems.Component, io_settings, "IO settings",
                                  Dict{String,Any}(), IO_SETTINGS_DEF)

    return io_settings
end

"""
    get_simulation_params(project_config, io_settings)

Constructs the dictionary of simulation parameters.

# Arguments
-`project_config::Dict{String,Any}`: The project config
-`io_settings::Dict{String,Any}`: IO settings, already extracted from the project config
# Returns
-`Dict{String,Any}`: The simulation parameter dictionary
"""
function get_simulation_params(project_config::AbstractDict{String,Any},
                               io_settings::Dict{String,Any})::Dict{String,Any}
    # load time and step info directly, bypassing extraction and validation
    time_step,
    start_date,
    start_date_output,
    end_date,
    nr_of_steps,
    nr_of_steps_output = get_timesteps(project_config["simulation_parameters"])

    # extract most parameter values using the extract function on the base type component
    # (even though we are not checking components...)
    sim_params = Dict{String,Any}()
    for (name, param_def) in pairs(SIMULATION_PARAMETERS_DEF)
        if name in ["start", "start_output", "end", "start_end_unit", "time_step", "time_step_unit"]
            continue
        end
        sim_params[name] = EnergySystems.extract_parameter(EnergySystems.Component,
                                                           project_config["simulation_parameters"],
                                                           name, param_def, Dict{String,Any}(), "Sim params")
    end

    # again, we use the component validation function as it avoids duplicate code
    EnergySystems.validate_config(EnergySystems.Component, sim_params, "Sim params",
                                  Dict{String,Any}(), SIMULATION_PARAMETERS_DEF)

    sim_params = merge(sim_params,
                       Dict{String,Any}(
                           "time" => 0,
                           "time_since_output" => 0,
                           "current_date" => start_date,
                           "time_step_seconds" => time_step,
                           "number_of_time_steps" => nr_of_steps,
                           "number_of_time_steps_output" => nr_of_steps_output,
                           "start_date" => start_date,
                           "start_date_output" => start_date_output,
                           "end_date" => end_date,
                           "step_info_interval" => default(io_settings, "step_info_interval",
                                                           Integer(floor(nr_of_steps / 20))),
                           "show_detailed_errors" => io_settings["show_detailed_errors"],
                       ))

    sim_params["economic_parameters"] = get_economic_parameters(project_config, sim_params)
    sim_params["emissions_parameters"] = get_emissions_parameters(project_config, sim_params)

    # add helper functions to convert power to work and vice-versa. this uses the time step
    # of the simulation as the duration required for the conversion.
    sim_params["watt_to_wh"] = function (watts::Float64)
        return watts * time_step * HOURS_PER_SECOND
    end
    sim_params["wh_to_watts"] = function (wh::Float64)
        return wh * SECONDS_PER_HOUR / time_step
    end

    # add helper function for using paths, absolute or relative to the run base path
    sim_params["run_path"] = function (path)
        return isabspath(path) ? path : abspath(joinpath(io_settings["base_path"], path))
    end

    # load weather profiles accessible for all components
    weather_file_path = sim_params["weather_file_path"]
    if weather_file_path !== nothing
        # WeatherData() writes the latitude and longitude to sim_params if either of them is
        # nothing at this point
        sim_params["weather_data"] = WeatherData(sim_params["run_path"](weather_file_path),
                                                 sim_params,
                                                 guess_file_format(sim_params["run_path"](weather_file_path)),
                                                 sim_params["weather_interpolation_type_solar"],
                                                 sim_params["weather_interpolation_type_general"])
    end

    return sim_params
end

"""
    prepare_inputs(project_config, run_ID)

Construct and prepare parameters, energy system components and the order of operation.

# Arguments
-`project_config::AbstractDict{String,Any}`: The project config
-`run_ID::UUID`: The run ID used in the run registry
# Returns
-`Dict{String,Any}`: Simulation parameters
-`Dict{String,Any}`: IO settings
-`Grouping`: The constructed energy system components
-`OrderOfOperations`: Order of operations
"""
function prepare_inputs(project_config::AbstractDict{String,Any}, run_ID::UUID)
    io_settings = get_io_settings(project_config)
    sim_params = get_simulation_params(project_config, io_settings)
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

    return sim_params, io_settings, components, operations
end

"""
    load_control_module_class_mapping()

Loads the control modules' classes index by their name as used in the input file.

Returns:
-`Dict{String, Any}`: The mapping from name (String) to the module class, which is probably
    of type `Symbol`, however `getproperty` does not specify the return type. In any case
    the entry can be used for calling the constructor of the control module as a function.
"""
function load_control_module_class_mapping()::Dict{String,Any}
    mapping = Dict{String,Any}()

    for name in names(EnergySystems; all=true)
        if startswith(String(name), "CM_")
            symbol = Symbol(String(name))
            unit_class = getproperty(EnergySystems, symbol)

            if unit_class <: EnergySystems.ControlModule
                module_name = nothing

                # type-level accessor function implemented per control module as following:
                # control_module_name(::Type{CM_ModuleTypeName})::String = "module_type_name"
                try
                    module_name = EnergySystems.control_module_name(unit_class)
                    mapping[module_name] = unit_class
                catch
                    @error("Control module type $name does not have a method defined for " *
                           "function control_module_name.")
                end
            end
        end
    end

    return mapping
end

"""
load_components(config, sim_params)

Construct instances of components from the given config.

The config must have the structure:
```
{
"UAC key": {
    "type": "PVPlant",
    ...
},
...
}
```

The required config to construct a component from one entry in the config must match what is
required for the particular component. The `type` parameter must be present and must match
the symbol of the component class exactly. The structure is described in more detail in the
accompanying documentation on the project file.
"""
function load_components(config_ordered::AbstractDict{String,Any}, sim_params::Dict{String,Any})::Grouping
    # convert OrderedDict to normal Dict to have a normal dict in all components as they do not
    # require any sorting
    to_dict(x) = x
    to_dict(x::OrderedDict) = Dict{String,Any}(k => to_dict(v) for (k, v) in x)
    to_dict(x::AbstractVector) = map(to_dict, x)
    config = to_dict(config_ordered)

    components = Grouping()

    # create instances
    for (unit_key, entry) in pairs(config)
        default_dict = Dict{String,Any}()
        unit_config = Base.merge(default_dict, entry)

        symbol = Symbol(String(unit_config["type"]))
        unit_class = getproperty(EnergySystems, symbol)
        if unit_class <: EnergySystems.Component
            instance = unit_class(unit_key, unit_config, sim_params)
            components[unit_key] = instance
        end
    end

    # link inputs/outputs
    for (unit_key, entry) in pairs(config)
        if String(entry["type"]) != "Bus" && haskey(entry, "output_refs") && length(entry["output_refs"]) > 0
            if isa(entry["output_refs"], AbstractDict)
                # components with multiple outputs should enter the output_refs as Dict to achieve uniqueness 
                media_keys = collect(keys(entry["output_refs"]))
                target_components = [Grouping(uac => components[uac])
                                     for uac in [entry["output_refs"][key] for key in media_keys]]
                media_sym = Symbol[]
                for medium in media_keys
                    if hasproperty(components[unit_key], Symbol(medium))
                        push!(media_sym, getproperty(components[unit_key], Symbol(medium)))
                    else
                        @error "For component $unit_key, the given key `$medium` in the `output_refs` is not a valid key!"
                        throw(InputError())
                    end
                end
                link_output_with(components[unit_key], target_components; given_media=media_sym)
            else
                if length(entry["output_refs"]) > 1
                    @warn "The component $unit_key has more than one output interface, but the `output_refs` are not " *
                          "specified explicitly! This can work, but it can also cause wrong interconnection between " *
                          "components! Consider using a mapping of the media to the target components "
                end
                target_components = Grouping(uac => components[uac] for uac in entry["output_refs"])
                link_output_with(components[unit_key], target_components)
            end
        elseif String(entry["type"]) == "Bus" && haskey(entry, "connections") && length(entry["connections"]) > 0
            target_components = Grouping(uac => components[uac] for uac in entry["connections"]["output_order"])
            link_output_with(components[unit_key], target_components)
        end
    end

    # add control modules to components
    mapping = load_control_module_class_mapping()
    for (unit_key, entry) in pairs(config)
        unit = components[unit_key]

        for module_config in default(entry, "control_modules", [])
            if !haskey(mapping, module_config["name"])
                @warn("Unknown control module type $(module_config["name"]) while loading " *
                      "unit $(unit.uac)")
                continue
            end
            module_class = mapping[module_config["name"]]
            push!(unit.controller.modules, module_class(module_config, components, sim_params, unit.uac))
        end
    end

    # the input/output interfaces of busses are constructed in the order of appearance in
    # the config, so after all components are loaded they need to be reordered to match
    # the input/output priorities
    components = reorder_interfaces_of_busses!(components)

    # other type-specific initialisation
    EnergySystems.initialise_components(components, sim_params)

    # create proxy busses from bus chains
    chains = find_chains(values(components), EnergySystems.sf_bus)
    EnergySystems.merge_bus_chains(chains, components, sim_params)

    return components
end

"""
    reorder_interfaces_of_busses!(components)

Calls reorder_interfaces_of_bus!() for all busses in the given grouping of components.

Args:
-`components::Grouping`: The components
Return:
-`Grouping`: The components with busses having their interfaces reordered
"""
function reorder_interfaces_of_busses!(components::Grouping)::Grouping
    for unit in each(components)
        if unit.sys_function == EnergySystems.sf_bus
            reorder_interfaces_of_bus!(unit)
        end
    end
    return components
end

"""
    reorder_interfaces_of_bus!(bus)

Reorder the input and output interfaces of busses according to their input and output
priorities given in the connectivity matrix.

Args:
-`bus::EnergySystems.Bus`: The bus for which to reorder interfaces
"""
function reorder_interfaces_of_bus!(bus::EnergySystems.Bus)
    # get correct order according to connectivity matrix
    output_order = bus.connectivity.output_order
    input_order = bus.connectivity.input_order

    # check for misconfigured bus (it should have at least one input and at least
    # one output)
    if length(input_order) == 0 || length(output_order) == 0
        return
    end

    # Create a dictionary to map 'uac' to its correct position
    output_order_dict = Dict(uac => idx for (idx, uac) in enumerate(output_order))
    input_order_dict = Dict(uac => idx for (idx, uac) in enumerate(input_order))

    # Get the permutation indices that would sort the 'source'/'target' field by
    # 'uac' order
    output_perm_indices = sortperm([output_order_dict[bus.output_interfaces[i].target.uac]
                                    for i in 1:length(bus.output_interfaces)])

    # Input side: include is_secondary_interface in the key
    component_key(uac::AbstractString, is_secondary_interface::Bool) = is_secondary_interface ?
                                                                       string(uac, "#secondary") : uac
    input_perm_indices = sortperm([input_order_dict[component_key(bus.input_interfaces[i].source.uac,
                                                                  bus.input_interfaces[i].is_secondary_interface)]
                                   for i in eachindex(bus.input_interfaces)])

    # Reorder the input and output interfaces using the permutation indices
    bus.output_interfaces = bus.output_interfaces[output_perm_indices]
    bus.input_interfaces = bus.input_interfaces[input_perm_indices]
end

"""
get_timesteps(simulation_parameters)

Function to read in the time step information from the input file.
If no information is given in the input file, the following defaults 
will be set:
time_step = 900 s
"""
function get_timesteps(simulation_parameters::AbstractDict{String,Any})
    start_date = DateTime(0)
    start_date_output = DateTime(0)
    end_date = DateTime(0)
    try
        start_date = Dates.DateTime(simulation_parameters["start"], simulation_parameters["start_end_unit"])
        end_date = Dates.DateTime(simulation_parameters["end"], simulation_parameters["start_end_unit"])
        if haskey(simulation_parameters, "start_output")
            start_date_output = Dates.DateTime(simulation_parameters["start_output"],
                                               simulation_parameters["start_end_unit"])
        else
            start_date_output = start_date
        end
    catch e
        @error("Time given 'start_end_unit' of the simulation parameters does not fit to the data.\n" *
               "'start_end_unit' has to be a daytime format, e.g. 'dd-mm-yyyy HH:MM:SS'.\n" *
               "'start_end_unit' is `$(simulation_parameters["start_end_unit"])` which does not fit to the start" *
               "and end time given: `$(simulation_parameters["start"])` and `$(simulation_parameters["end"])`.\n" *
               "The following error occurred: $e")
        throw(InputError())
    end
    if start_date_output < start_date
        @error "The start date of the output can not be prior to the start date of the simulation!"
        throw(InputError())
    end
    if simulation_parameters["time_step_unit"] == "seconds"
        time_step = simulation_parameters["time_step"]
    elseif simulation_parameters["time_step_unit"] == "minutes"
        time_step = simulation_parameters["time_step"] * 60
    elseif simulation_parameters["time_step_unit"] == "hours"
        time_step = simulation_parameters["time_step"] * 60 * 60
    else
        time_step = 900
        @info("The simulation time step is set to 900 s as default, as it could not be found in the input file" *
              "(`time_step` and `time_step_unit` have to be given!).")
    end

    nr_of_steps = UInt(max(0, floor(Dates.value(Second(sub_ignoring_leap_days(end_date, start_date))) / time_step)) + 1)
    nr_of_steps_output = UInt(max(0,
                                  floor(Dates.value(Second(sub_ignoring_leap_days(end_date, start_date_output))) /
                                        time_step)) + 1)

    # set end_date to be integer dividable by the timestep
    end_date = add_ignoring_leap_days(start_date, (nr_of_steps - 1) * Second(time_step))

    if (month(start_date) == 2 && day(start_date) == 29) ||
       (month(end_date) == 2 && day(end_date) == 29) ||
       (month(start_date_output) == 2 && day(start_date_output) == 29)
        @error "The simulation start and end date and the start date of the output can not be at a leap day!"
        throw(InputError())
    end
    return UInt(time_step), start_date, start_date_output, end_date, nr_of_steps, nr_of_steps_output
end

# calculation of the order of operations has its own include files due to its complexity
include("order_of_operations.jl")


function parse_optimizer_function(eff_def::String)::Function
    splitted = split(eff_def, ":")

    if length(splitted) > 1
        method = lowercase(splitted[1])
        data = splitted[2]

        #TODO get more function definitions analog to different optimization packages
        if method == "poly-2"
            params = parse.(Float64, split(data, ","))
            return function (x, y)
                return params[1] +
                       params[2] * x +
                       params[3] * y +
                       params[4] * x * x +
                       params[5] * x * y +
                       params[6] * y * y +
                       params[7] * x * x * x +
                       params[8] * x * x * y +
                       params[9] * x * y * y +
                       params[10] * y * y * y
            end
        end
    end

    @error "Cannot parse 2-dimensional function from: $eff_def"
    return (x, y) -> 0.0
end