# this file contains functionality pertaining to loading a project's metadata and the
# energy system components from the project config file, as well as constructing certain
# helpful information data structures from the inputs in the config
using JSON: JSON
using OrderedCollections: OrderedDict
using Logging

const HOURS_PER_SECOND::Float64 = 1.0 / 3600.0
const SECONDS_PER_HOUR::Float64 = 3600.0

"""
Shared preparation cache used during parameter studies.

The cache stores expensive, effectively read-only setup results that may be reused between
simulation samples. Components themselves are intentionally not cached because component
instances are mutable during a simulation run.
"""
mutable struct PreparationCache
    lock::ReentrantLock
    profiles::Dict{Any,Any}
    weather_data::Dict{Any,Any}
    operations::Dict{Any,Any}
end

PreparationCache() = PreparationCache(ReentrantLock(),
                                      Dict{Any,Any}(),
                                      Dict{Any,Any}(),
                                      Dict{Any,Any}())

function operation_cache_allowed(sim_params::Dict{String,Any})::Bool
    parameter_study = get(sim_params, "parameter_study", Dict{String,Any}())

    runtime = get(parameter_study, "runtime", Dict{String,Any}())

    if !haskey(runtime, "parameter_keys")
        return true
    end

    # Conservative blocklist. If one of these parameters is varied, the component graph
    # or the operation order may change, so the operation cache is disabled.
    structural_params = Set(["type",
                             "medium",
                             "m_el_in",
                             "m_heat_in",
                             "m_heat_out",
                             "m_heat_out_secondary",
                             "has_secondary_interface",
                             "primary_el_sources",
                             "secondary_el_sources",
                             "input_refs",
                             "output_refs",
                             "input_order",
                             "output_order",
                             "connections",
                             "energy_flow"])

    for key in runtime["parameter_keys"]
        parts = split(key, " ")
        if parts[end] in structural_params
            return false
        end
    end

    return true
end

function operation_cache_key(project_config::AbstractDict{String,Any},
                             sim_params::Dict{String,Any})::String
    parameter_study = get(sim_params, "parameter_study", Dict{String,Any}())
    runtime = get(parameter_study, "runtime", Dict{String,Any}())

    component_cfg = deepcopy(project_config["components"])

    # Non-structural varied component values should not invalidate the
    # operation-order cache.
    if haskey(runtime, "parameter_keys")
        for key in runtime["parameter_keys"]
            uac, param_key = split(key, " ")
            if haskey(component_cfg, uac) && haskey(component_cfg[uac], param_key)
                component_cfg[uac][param_key] = "__PARAMETER_STUDY_NONSTRUCTURAL_VALUE__"
            end
        end
    end

    order_cfg = get(project_config, "order_of_operation", Any[])

    return JSON.json(Dict(
                         "components" => component_cfg,
                         "order_of_operation" => order_cfg,
                     ))
end

function build_operations(project_config::AbstractDict{String,Any},
                          components::Grouping)::OrderOfOperations
    if haskey(project_config, "order_of_operation") && length(project_config["order_of_operation"]) > 0
        operations = load_order_of_operations(project_config["order_of_operation"], components)
        @info "The order of operations was successfully imported from the input file.\n" *
              "Note that the order of operations has a major impact on the simulation " *
              "result and should only be changed by experienced users!"
        return operations
    else
        return calculate_order_of_operations(components)
    end
end

function get_operations(project_config::AbstractDict{String,Any},
                        components::Grouping,
                        sim_params::Dict{String,Any},
                        preparation_cache::Union{Nothing,PreparationCache})::OrderOfOperations
    if preparation_cache === nothing || !operation_cache_allowed(sim_params)
        return build_operations(project_config, components)
    end

    key = operation_cache_key(project_config, sim_params)

    cached = lock(preparation_cache.lock) do
        get(preparation_cache.operations, key, nothing)
    end

    if cached !== nothing
        # Return a fresh vector so reorderings of the run-local operation vector cannot
        # mutate the cached template.
        return deepcopy(cached)
    end

    operations = build_operations(project_config, components)

    lock(preparation_cache.lock) do
        if !haskey(preparation_cache.operations, key)
            preparation_cache.operations[key] = deepcopy(operations)
        end
    end

    return operations
end


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
        display_name="Write CSV continuously?",
        required=false,
        type=Bool,
        json_type="boolean",
        unit="-"
    ),
    "write_summary_csv" => (
        default=true,
        description="Toggle if a CSV summary output with sum/mean values should be created " *
                    "additionally to the timestep-wise CSV output.",
        display_name="Write summary CSV?",
        required=false,
        type=Bool,
        json_type="boolean",
        unit="-"
    ),
    "csv_output_file_path" => (
        default="./output/out.csv",
        description="File path to where the CSV output will be written",
        display_name="CSV output file path",
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
    "output_plot_file_path" => (
        default="./output/output_plot.html",
        description="File path to where the output line plot will be written",
        display_name="Output plot file path",
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
    "sankey_plot_file_path" => (
        default="./output/output_sankey.html",
        description="File path to where the Sankey plot will be written",
        display_name="Sankey plot file path",
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
    "output_economic_csv" => (
        default=true,
        description="Toggle if a CSV with the economic results should be created",
        display_name="Output economic CSV?",
        required=false,
        type=Bool,
        json_type="boolean",
        unit="-"
    ),
    "economic_csv_file_path" => (
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
    "output_emissions_csv" => (
        default=true,
        description="Toggle if a CSV with the emission results should be created",
        display_name="Output emissions CSV?",
        required=false,
        type=Bool,
        json_type="boolean",
        unit="-"
    ),
    "emissions_csv_file_path" => (
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
        description="If given a non-zero value, uses this many significant digits as the " *
                    "fixed precision for float outputs in CSV and plot files. It is not " *
                    "recommended to use this setting in normal simulation. It's intended " *
                    "for making the output perfectly repeatable, which is useful for " *
                    "testing but changes the results.",
        display_name="Fixed output precision",
        required=false,
        type=Integer,
        json_type="number",
        unit="-"
    ),
    "output_parameter_study_csv" => (
        default=true,
        description="Toggle if CSV file(s) with the parameter-study results should be created",
        display_name="Output parameter-study CSV?",
        required=false,
        type=Bool,
        json_type="boolean",
        unit="-"
    ),
    "parameter_study_csv_path" => (
        default="./output/parameter_study",
        description="Directory path to where the parameter-study CSV result files are written to",
        display_name="Parameter-study CSV directory path",
        required=false,
        type=String,
        json_type="string",
        unit="-"
    ),
    "write_parameter_study_csv_continuously" => (
        default=false,
        description="Toggle if parameter-study CSV output will be written continuously, " * 
                    "meaning after every run. Activating this functionality will ensure " * 
                    "partial output if the parameter study is stopped during execution. " *
                    "It incurs a slight performance penalty depending on the run time of " *
                    "one simulation and the number of parallel runs, since the threads " *
                    "might have to wait for write access.", 
        display_name="Write parameter-study CSV continuously?",
        required=false,
        type=Bool,
        json_type="boolean",
        unit="-"
    ),   
    "output_parameter_study_plots" => (
        default=true,
        description="Toggle if plots with the parameter-study results should be created",
        display_name="Output parameter-study plots?",
        required=false,
        type=Bool,
        json_type="boolean",
        unit="-"
    ),
    "parameter_study_plots_path" => (
        default="./output/parameter_study_plots",
        description="Directory path to where the parameter-study result plots will be written",
        display_name="Parameter-study plots directory path",
        required=false,
        type=String,
        json_type="string",
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
    ("csv_output", "csv_output_keys")
]

const OPTIMISER_LIMITS_DEF = Dict{String,Any}(
    "max_runs" => (
        default=nothing,
        description="Maximum number of objective evaluations executed by the optimiser.",
        display_name="Maximum runs",
        required=false,
        validations=[("self", "value_gte_num_or_nothing", 1)],
        type=Int64,
        json_type="number",
        unit="-"
    ),
    "max_time" => (
        default=nothing,
        description="Maximum optimisation runtime in seconds.",
        display_name="Maximum time",
        required=false,
        validations=[("self", "value_gte_num_or_nothing", 0.0)],
        type=Float64,
        json_type="number",
        unit="s"
    ),
    "x_tol_abs" => (
        default=nothing,
        description="Tolerance for normalised optimisation parameters in the range " *
                    "[0, 1]. It is translated to the closest backend-specific stopping " *
                    "criterion, which may be absolute or relative; support depends on " *
                    "the selected algorithm.",
        display_name="Parameter tolerance",
        required=false,
        conditionals=[("type", "is_one_of", ("Optim", "NLopt", "NOMAD"))],
        validations=[("self", "value_gt_num_or_nothing", 0.0)],
        type=Float64,
        json_type="number",
        unit="-"
    ),
    "f_tol_abs" => (
        default=nothing,
        description="Tolerance for changes or spread in the objective function. It is " *
                    "translated to the closest backend-specific stopping criterion; " *
                    "support depends on the selected algorithm.",
        display_name="Objective tolerance",
        required=false,
        conditionals=[("type", "is_one_of", ("Optim", "BlackBoxOptim", "NLopt"))],
        validations=[("self", "value_gt_num_or_nothing", 0.0)],
        type=Float64,
        json_type="number",
        unit="-"
    ),
)

const REFINEMENT_OPTIMISATION_DEF = Dict{String,Any}(
    "type" => (
        default=nothing,
        description="Optimisation backend used for the refinement stage.",
        display_name="Refinement optimisation type",
        required=true,
        options=["Optim", "BlackBoxOptim", "Metaheuristics", "NLopt", "NOMAD"],
        type=String,
        json_type="string",
        unit="-"
    ),
    "algorithm" => (
        default=nothing,
        description="Algorithm used by the selected refinement optimisation backend.",
        display_name="Refinement algorithm",
        required=true,
        type=String,
        json_type="string",
        unit="-"
    ),
    "max_runs" => OPTIMISER_LIMITS_DEF["max_runs"],
    "max_time" => OPTIMISER_LIMITS_DEF["max_time"],
    "x_tol_abs" => OPTIMISER_LIMITS_DEF["x_tol_abs"],
    "f_tol_abs" => OPTIMISER_LIMITS_DEF["f_tol_abs"],
    "optim_kwargs" => (
        default=nothing,
        description="Additional keyword arguments passed to the selected refinement backend.",
        display_name="Refinement keyword arguments",
        required=false,
        type=Dict{String,Any},
        json_type="object",
        unit="-"
    ),
)

const PARAMETER_STUDY_PARAMETER_DEF = Dict{String,Any}(
    "bounds_min" => (
        default=nothing,
        description="Lower bound used by optimisation and global sensitivity.",
        display_name="Minimum bound",
        required=false,
        type=Float64,
        json_type="number",
        unit="-"
    ),
    "bounds_max" => (
        default=nothing,
        description="Upper bound used by optimisation and global sensitivity.",
        display_name="Maximum bound",
        required=false,
        type=Float64,
        json_type="number",
        unit="-"
    ),
    "start" => (
        default=nothing,
        description="Initial optimisation value or standalone local-sensitivity reference.",
        display_name="Start value",
        required=false,
        type=Float64,
        json_type="number",
        unit="-"
    ),
    "sensitivity_lower" => (
        default=nothing,
        description="Explicit lower evaluation value for local sensitivity.",
        display_name="Lower sensitivity value",
        required=false,
        type=Float64,
        json_type="number",
        unit="-"
    ),
    "sensitivity_upper" => (
        default=nothing,
        description="Explicit upper evaluation value for local sensitivity.",
        display_name="Upper sensitivity value",
        required=false,
        type=Float64,
        json_type="number",
        unit="-"
    ),
    "values" => (
        default=nothing,
        description="Explicit discrete values used for parameter variation.",
        display_name="Parameter variation values",
        required=false,
        type=Vector,
        json_type="list",
        unit="-"
    ),
    "range_start" => (
        default=nothing,
        description="First value of a generated parameter-variation range.",
        display_name="Range start",
        required=false,
        type=Float64,
        json_type="number",
        unit="-"
    ),
    "range_stop" => (
        default=nothing,
        description="Last value of a generated parameter-variation range.",
        display_name="Range stop",
        required=false,
        type=Float64,
        json_type="number",
        unit="-"
    ),
    "range_step" => (
        default=nothing,
        description="Step size of a generated parameter-variation range.",
        display_name="Range step",
        required=false,
        type=Float64,
        json_type="number",
        unit="-"
    ),
    "range_length" => (
        default=nothing,
        description="Number of values in a generated parameter-variation range.",
        display_name="Range length",
        required=false,
        validations=[("self", "value_gte_num_or_nothing", 2)],
        type=Integer,
        json_type="number",
        unit="-"
    ),
)

const PARAMETER_STUDY_DEF = Dict{String,Any}(
    "parameters" => (
        default=nothing,
        description="Defines the project parameters used by the parameter study. Each " *
                    "selected parameter can define bounds, an initial value, explicit " *
                    "or generated variation values and local-sensitivity values.",
        display_name="Parameter study parameters",
        required=true,
        type=Dict{String,Any},
        json_type="object",
        unit="-"
    ),
    "objective_params" => (
        default=nothing,
        description="Defines which parameters are set as the objective. Definition " *
                    "follows the definition of components. See the documentation for " *
                    "more details.",
        display_name="Objective parameters",
        required=true,
        type=Dict{String,Any},
        json_type="object",
        unit="-"
    ),
    "objective_function" => (
        default="sum",
        description="Defines how objective_params are combined. Supported values are " *
                    "`sum`, `linear` and `multi-objective`. For `linear`, coefficients " *
                    "must be supplied by objective_factors using the flattened objective " *
                    "parameter names as keys.",
        display_name="Objective function",
        required=false,
        options=["sum", "linear", "multi-objective"],
        type=String,
        json_type="string",
        unit="-"
    ),
    "objective_factors" => (
        default=nothing,
        description="Named coefficients for objective_function=`linear`. Every flattened " *
                    "objective parameter key must occur exactly once. Example: " *
                    "{\"economic total_annuity\": 1.0, " *
                    "\"sum GridOut m_e_ac_230v IN\": -2.0}.",
        display_name="Linear objective factors",
        required=false,
        conditionals=[("objective_function", "is", "linear")],
        type=Dict{String,Any},
        json_type="object",
        unit="-"
    ),
    "objective_senses" => (
        default=nothing,
        description="Defines whether each objective is minimised or maximised " *
                    "for objective_function=`multi-objective`. Keys must exactly " *
                    "match the flattened objective parameter names. Allowed " *
                    "values are `min` and `max`. If omitted, all objectives are " *
                    "minimised.",
        display_name="Objective directions",
        required=false,
        conditionals=[("objective_function", "is", "multi-objective")],
        type=Dict{String,Any},
        json_type="object",
        unit="-",
    ),
    "disable_all_simulation_outputs" => (
        default=true,
        description="Disables simulation outputs written to the hard drive during " *
                    "parameter studies.",
        display_name="Disable all simulation outputs",
        required=false,
        type=Bool,
        json_type="boolean",
        unit="-"
    ),
    "parameter_variation" => (
        default=nothing,
        description="Configuration of a parameter variation.",
        display_name="Parameter variation",
        required=false,
        type=Dict{String,Any},
        json_type="object",
        unit="-"
    ),
    "optimisation" => (
        default=nothing,
        description="Configuration of an optimisation.",
        display_name="Optimisation",
        required=false,
        type=Dict{String,Any},
        json_type="object",
        unit="-"
    ),
    "sensitivity_analysis" => (
        default=nothing,
        description="Configuration of global and local sensitivity analyses.",
        display_name="Sensitivity analysis",
        required=false,
        type=Dict{String,Any},
        json_type="object",
        unit="-"
    ),
)

const PARAMETER_VARIATION_DEF = Dict{String,Any}(
    "run_parameter_variation" => (
        default=false,
        description="If set to true, executes the configured parameter variation.",
        display_name="Run parameter variation?",
        required=false,
        type=Bool,
        json_type="boolean",
        unit="-"
    ),
    "algorithm" => (
        default="product",
        description="Combination method for parameter values. Supported values are " *
                    "product, zip and random_*, where * is the number of samples.",
        display_name="Parameter variation algorithm",
        required=false,
        type=String,
        json_type="string",
        unit="-"
    ),
)

const PARAMETER_STUDY_OPTIMISATION_DEF = Dict{String,Any}(
    "run_optimisation" => (
        default=false,
        description="If set to true, runs the configured optimisation.",
        display_name="Run optimisation?",
        required=false,
        type=Bool,
        json_type="boolean",
        unit="-"
    ),
    "type" => (
        default="NLopt",
        description="Selects the optimisation backend.",
        display_name="Optimisation type",
        required=true,
        conditionals=[("run_optimisation", "is", true)],
        options=["Optim", "BlackBoxOptim", "Metaheuristics", "NLopt", "NOMAD"],
        type=String,
        json_type="string",
        unit="-"
    ),
    "algorithm" => (
        default="LN_SBPLX",
        description="Algorithm used by the selected optimisation backend.",
        display_name="Optimisation algorithm",
        required=false,
        conditionals=[("run_optimisation", "is", true)],
        type=String,
        json_type="string",
        unit="-"
    ),
    "max_runs" => OPTIMISER_LIMITS_DEF["max_runs"],
    "max_time" => OPTIMISER_LIMITS_DEF["max_time"],
    "x_tol_abs" => OPTIMISER_LIMITS_DEF["x_tol_abs"],
    "f_tol_abs" => OPTIMISER_LIMITS_DEF["f_tol_abs"],
    "optim_kwargs" => (
        default=nothing,
        description="Additional keyword arguments passed to the selected optimisation backend.",
        display_name="Optimisation keyword arguments",
        required=false,
        type=Dict{String,Any},
        json_type="object",
        unit="-"
    ),
    "refinement" => (
        default=nothing,
        description="Optional second optimisation stage starting from the best result " *
                    "of the primary optimisation. Its nested fields are exposed as " *
                    "refinement_optimisation by all_general_parameters().",
        display_name="Refinement optimiser",
        required=false,
        type=Dict{String,Any},
        json_type="object",
        unit="-"
    ),
)

const SENSITIVITY_ANALYSIS_DEF = Dict{String,Any}(
    "run_global_sensitivity" => (
        default=false,
        description="If set to true, calculates global sensitivity indices. Existing " *
                    "results from a parameter variation or optimisation are reused.",
        display_name="Run global sensitivity analysis?",
        required=false,
        type=Bool,
        json_type="boolean",
        unit="-"
    ),
    "run_local_sensitivity" => (
        default=false,
        description="If set to true, calculates local sensitivity around configured " *
                    "start values or the best optimisation result.",
        display_name="Run local sensitivity analysis?",
        required=false,
        type=Bool,
        json_type="boolean",
        unit="-"
    ),
    "local_reference" => (
        default="start",
        description="Reference point for local sensitivity. Use `start` for configured " *
                    "parameter start values or `best_result` for the best result so far " *
                    "from optimisation or global sensitivity analysis.",
        display_name="Local sensitivity reference",
        required=false,
        options=["start", "best_result"],
        type=String,
        json_type="string",
        unit="-"
    ),
    "local_variation" => (
        default=0.1,
        description="Relative variation applied in positive and negative direction when " *
                    "no explicit sensitivity_lower and sensitivity_upper values are set.",
        display_name="Local relative variation",
        required=false,
        validations=[("self", "value_gt_num", 0.0)],
        type=Float64,
        json_type="number",
        unit="-"
    ),
    "include_local_sensitivity_results_in_figures" => (
        default=true,
        description="If set to true, simulation results generated exclusively for local " *
                    "sensitivity analysis are additionally included in the general " *
                    "parameter-study result figures.",
        display_name="Include local sensitivity results in figures?",
        required=false,
        type=Bool,
        json_type="boolean",
        unit="-"
    ),
    "max_runs" => (
        default=nothing,
        description="Maximum number of (additional) simulation performed for global " *
                    "sensitivity. If omitted, twice the optimisation max_runs value are used.",
        display_name="Maximum sensitivity runs",
        required=false,
        validations=[("self", "value_gte_num_or_nothing", 0)],
        type=Int64,
        json_type="number",
        unit="-"
    ),
)
#! format: on

"""
    read_JSON(filepath)

Read and parse the JSON-encoded Dict in the given file.
"""
function read_JSON(filepath::String)::OrderedDict{String,Any}
    open(filepath, "r") do file_handle
        content = read(file_handle, String)
        return JSON.parse(content; dicttype=OrderedDict)
    end
end

"""
    get_min_log_level(project_config, logger)

Determine log level:
- for parameter-study runs, allow only logs >= GlobalInfo
- for single simulations, set min log level to Info level

Therefore, we need early access to the parameter-study run flags.

# Arguments
- `project_config::OrderedDict{String,Any}`: Base project configuration used to generate
  simulation variants.
- `logger::Union{Nothing,Resie_Logger.CustomLogger}`: Logger used for ReSiE
# Returns
- `min_log_level::LogLevel`: minimum log level used in ReSiE for logging
"""
function get_min_log_level(project_config, logger)
    parameter_study = get(project_config, "parameter_study", nothing)
    parameter_variation = parameter_study isa AbstractDict ?
                          something(get(parameter_study, "parameter_variation", nothing), Dict{String,Any}()) :
                          Dict{String,Any}()
    optimisation = parameter_study isa AbstractDict ?
                   something(get(parameter_study, "optimisation", nothing), Dict{String,Any}()) :
                   Dict{String,Any}()
    sensitivity_analysis = parameter_study isa AbstractDict ?
                           something(get(parameter_study, "sensitivity_analysis", nothing), Dict{String,Any}()) :
                           Dict{String,Any}()
    run_parameter_study = get(parameter_variation, "run_parameter_variation", false) == true ||
                          get(optimisation, "run_optimisation", false) == true ||
                          get(sensitivity_analysis, "run_global_sensitivity", false) == true ||
                          get(sensitivity_analysis, "run_local_sensitivity", false) == true

    if logger !== nothing
        if run_parameter_study
            # for parameter studies, allow only logs >= GlobalInfo
            min_log_level = Logging.LogLevel(600)
        else
            # for single simulations, set min log level to Info level
            min_log_level = Logging.Info
        end
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
        if name in last.(OUTPUT_SPECIFICATION_SETTINGS)
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
                               io_settings::Dict{String,Any};
                               preparation_cache::Union{Nothing,PreparationCache}=nothing)::Dict{String,Any}
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

    if preparation_cache !== nothing
        sim_params["preparation_cache"] = preparation_cache
    end

    sim_params["economic_parameters"] = get_economic_parameters(project_config, sim_params)
    sim_params["emissions_parameters"] = get_emissions_parameters(project_config, sim_params)
    sim_params["parameter_study"] = get_parameter_study_parameters(project_config, sim_params)

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
        weather_path_abs = sim_params["run_path"](weather_file_path)

        if preparation_cache === nothing
            # WeatherData() writes the latitude and longitude to sim_params if either of them is
            # nothing at this point
            @globalInfo "Loading weather data..."
            sim_params["weather_data"] = WeatherData(weather_path_abs,
                                                     sim_params,
                                                     guess_file_format(weather_path_abs),
                                                     sim_params["weather_interpolation_type_solar"],
                                                     sim_params["weather_interpolation_type_general"])
        else
            # WeatherData may infer latitude, longitude, and time_zone from the file and write
            # them into sim_params. Cache those inferred values together with the weather data.
            key = (weather_path_abs,
                   sim_params["start_date"],
                   sim_params["end_date"],
                   sim_params["time_step_seconds"],
                   sim_params["weather_interpolation_type_solar"],
                   sim_params["weather_interpolation_type_general"],
                   sim_params["latitude"],
                   sim_params["longitude"],
                   sim_params["time_zone"])

            entry = lock(preparation_cache.lock) do
                get(preparation_cache.weather_data, key, nothing)
            end

            if entry === nothing
                @globalInfo "Loading weather data..."
                weather_data = WeatherData(weather_path_abs,
                                           sim_params,
                                           guess_file_format(weather_path_abs),
                                           sim_params["weather_interpolation_type_solar"],
                                           sim_params["weather_interpolation_type_general"])

                entry = (weather_data=weather_data,
                         latitude=sim_params["latitude"],
                         longitude=sim_params["longitude"],
                         time_zone=sim_params["time_zone"])

                lock(preparation_cache.lock) do
                    if !haskey(preparation_cache.weather_data, key)
                        preparation_cache.weather_data[key] = entry
                    else
                        entry = preparation_cache.weather_data[key]
                    end
                end
            end

            sim_params["weather_data"] = entry.weather_data
            sim_params["latitude"] = entry.latitude
            sim_params["longitude"] = entry.longitude
            sim_params["time_zone"] = entry.time_zone
        end
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
function prepare_inputs(project_config::AbstractDict{String,Any},
                        run_ID::UUID;
                        preparation_cache::Union{Nothing,PreparationCache}=nothing)
    io_settings = get_io_settings(project_config)
    sim_params = get_simulation_params(project_config, io_settings;
                                       preparation_cache=preparation_cache)
    sim_params["run_ID"] = run_ID

    components = load_components(project_config["components"], sim_params)
    operations = get_operations(project_config, components, sim_params, preparation_cache)

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

        if !isdefined(EnergySystems, symbol)
            @error "The component type `$(string(unit_config["type"]))` of component `$(unit_key)` is " *
                   "not a supported component type by ReSiE."
            throw(InputError())
        end
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
    function component_key(uac::AbstractString, is_secondary_interface::Bool)
        is_secondary_interface ?
        string(uac, "#secondary") : uac
    end
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

"""
    get_economic_parameters(project_config, sim_params)

Extract economic parameters form input file.

Args:
-`project_config::Dict{String,Any}`: The project config data
-`sim_params::Dict{String,Any}`: Simulation parameters, that have been extracted so far
Return:
-`Dict{String,Any}`: The economic parameters from the input file. If none are given,
    calculate_economy will be set to false.
"""
function get_economic_parameters(project_config::AbstractDict{String,Any},
                                 sim_params::Dict{String,Any})::Dict{String,Any}
    if !haskey(project_config, "economic_parameters")
        return Dict{String,Any}(
            "calculate_economy" => false,
        )
    end

    economic_parameters = Dict{String,Any}()
    for (name, param_def) in pairs(ECONOMIC_PARAMETERS_DEF)
        economic_parameters[name] = EnergySystems.extract_parameter(EnergySystems.Component,
                                                                    project_config["economic_parameters"],
                                                                    name,
                                                                    param_def,
                                                                    sim_params,
                                                                    "Economic parameters")
    end

    EnergySystems.validate_config(EnergySystems.Component, economic_parameters, "Economic parameters",
                                  sim_params, ECONOMIC_PARAMETERS_DEF)

    economic_parameters["repeat_period"], economic_parameters["repeat_method"] = get_repeat_period(economic_parameters["repeat_method"],
                                                                                                   sim_params,
                                                                                                   "economic")

    return economic_parameters
end

"""
    get_repeat_period(repeat_method::String, sim_params::Dict{String,Any}, type::String)

Get the repeat period as Dates module object (Day or similar) from the given repeat method.

# Arguments
-`repeat_method::String`: The repeat method from the settings
-`sim_params::Dict{String,Any}`: Simulation parameters
-`type::String`: The type (economic or emissions), used for error messages.
# Returns
-`Any`: An object of the Dates module, typically Day
"""
function get_repeat_period(repeat_method::String, sim_params::Dict{String,Any}, type::String)
    if repeat_method == "all"
        # With repeat method "all", the whole  profile is taken as it is and it is repeated until the 
        # observation_period_in_years is reached.
        repeat_period = sub_ignoring_leap_days(sim_params["end_date"], sim_params["start_date_output"]) +
                        Millisecond(Second(sim_params["time_step_seconds"]))
    elseif repeat_method == "last_year"
        repeat_period = Day(365)
    elseif repeat_method == "last_month"
        repeat_period = Day(30)
    elseif repeat_method == "last_week"
        repeat_period = Day(7)
    end
    if sub_ignoring_leap_days(sim_params["end_date"], sim_params["start_date_output"]) +
       Second(sim_params["time_step_seconds"]) < repeat_period
        @warn "In $type calculation, the repeat_method that defines the method for energy and price profile repetitions " *
              "was set to 'all' as the given repeat_method is longer than the simulation period."
        repeat_period = sub_ignoring_leap_days(sim_params["end_date"], sim_params["start_date_output"]) +
                        Millisecond(Second(sim_params["time_step_seconds"]))
        repeat_method = "all"
    end
    return repeat_period, repeat_method
end

"""
    get_emissions_parameters(project_config)

Extract emissions parameters form input file.

Args:
-`project_config::Dict{String,Any}`: The project config data
-`sim_params::Dict{String,Any}`: Simulation parameters, that have been extracted so far
Return:
-`Dict{String,Any}`: The emissions parameters from the input file. If none are given,
    calculate_emissions will be set to false.
"""
function get_emissions_parameters(project_config::AbstractDict{String,Any},
                                  sim_params::Dict{String,Any})::Dict{String,Any}
    if !haskey(project_config, "emissions_parameters")
        return Dict{String,Any}(
            "calculate_emissions" => false,
        )
    end

    emissions_parameters = Dict{String,Any}()
    for (name, param_def) in pairs(EMISSIONS_PARAMATERS_DEF)
        emissions_parameters[name] = EnergySystems.extract_parameter(EnergySystems.Component,
                                                                     project_config["emissions_parameters"],
                                                                     name,
                                                                     param_def,
                                                                     sim_params,
                                                                     "Emissions parameters")
    end

    EnergySystems.validate_config(EnergySystems.Component, emissions_parameters, "Emissions parameters",
                                  sim_params, EMISSIONS_PARAMATERS_DEF)

    emissions_parameters["repeat_period"], emissions_parameters["repeat_method"] = get_repeat_period(emissions_parameters["repeat_method"],
                                                                                                     sim_params,
                                                                                                     "emissions")

    return emissions_parameters
end

function extract_parameter_section(config::AbstractDict{String,Any},
                                   definitions::Dict{String,Any},
                                   sim_params::Dict{String,Any},
                                   display_name::String)::Dict{String,Any}
    extracted = Dict{String,Any}()
    for (name, param_def) in pairs(definitions)
        extracted[name] = EnergySystems.extract_parameter(EnergySystems.Component,
                                                          config,
                                                          name,
                                                          param_def,
                                                          sim_params,
                                                          display_name)
    end

    EnergySystems.validate_config(EnergySystems.Component, extracted, display_name,
                                  sim_params, definitions)

    return extracted
end

function get_parameter_study_parameters(project_config::AbstractDict{String,Any},
                                        sim_params::Dict{String,Any})::Dict{String,Any}
    if haskey(project_config, "parameter_study")
        return load_parameter_study_parameters(project_config["parameter_study"],
                                               sim_params,
                                               project_config)
    end

    return Dict{String,Any}(
        "runtime" => Dict{String,Any}("enabled" => false),
    )
end

"""
    all_general_parameters()::Dict{String,Any}

Lists all general, non-component parameters for the simulation engine.

# Returns
- `Dict{String,Any}`: The parameter definitions.
"""
function all_general_parameters()::Dict{String,Any}
    return Dict{String,Any}(
        "simulation" => SIMULATION_PARAMETERS_DEF,
        "io_settings" => IO_SETTINGS_DEF,
        "economic" => ECONOMIC_PARAMETERS_DEF,
        "emissions" => EMISSIONS_PARAMATERS_DEF,
        "parameter_study" => PARAMETER_STUDY_DEF,
        "parameter_study_parameter" => PARAMETER_STUDY_PARAMETER_DEF,
        "parameter_variation" => PARAMETER_VARIATION_DEF,
        "parameter_study_optimisation" => PARAMETER_STUDY_OPTIMISATION_DEF,
        "refinement_optimisation" => REFINEMENT_OPTIMISATION_DEF,
        "sensitivity_analysis" => SENSITIVITY_ANALYSIS_DEF,
    )
end

# calculation of the order of operations has its own include files due to its complexity
include("order_of_operations.jl")

function load_parameter_study_parameters(parameter_study_config::AbstractDict{String,Any},
                                         sim_params::Dict{String,Any},
                                         project_config::AbstractDict{String,Any})::Dict{String,Any}
    parameter_study = extract_parameter_section(parameter_study_config,
                                                PARAMETER_STUDY_DEF,
                                                sim_params,
                                                "Parameter study")

    parameter_study["parameter_variation"] = extract_parameter_section(something(parameter_study["parameter_variation"],
                                                                                 Dict{String,Any}()),
                                                                       PARAMETER_VARIATION_DEF,
                                                                       sim_params,
                                                                       "Parameter variation")
    parameter_study["optimisation"] = extract_parameter_section(something(parameter_study["optimisation"],
                                                                          Dict{String,Any}()),
                                                                PARAMETER_STUDY_OPTIMISATION_DEF,
                                                                sim_params,
                                                                "Optimisation")
    parameter_study["sensitivity_analysis"] = extract_parameter_section(something(parameter_study["sensitivity_analysis"],
                                                                                  Dict{String,Any}()),
                                                                        SENSITIVITY_ANALYSIS_DEF,
                                                                        sim_params,
                                                                        "Sensitivity analysis")

    parameter_variation = parameter_study["parameter_variation"]
    optimisation = parameter_study["optimisation"]
    sensitivity_analysis = parameter_study["sensitivity_analysis"]

    optimisation["refinement"] = get_refinement_optimiser_config(optimisation["refinement"],
                                                                 sim_params)
    if isnothing(optimisation["optim_kwargs"])
        optimisation["optim_kwargs"] = Dict{String,Any}()
    end

    run_parameter_variation = parameter_variation["run_parameter_variation"]
    run_optimisation = optimisation["run_optimisation"]
    run_global_sensitivity = sensitivity_analysis["run_global_sensitivity"]
    run_local_sensitivity = sensitivity_analysis["run_local_sensitivity"]
    enabled = run_parameter_variation || run_optimisation || run_global_sensitivity || run_local_sensitivity

    parameter_study["runtime"] = Dict{String,Any}(
        "enabled" => enabled,
        "run_primary_study" => run_parameter_variation || run_optimisation,
    )

    if !enabled
        return parameter_study
    end

    if run_parameter_variation && run_optimisation
        @error "Parameter variation and optimisation cannot be run in the same parameter study."
        throw(InputError())
    end

    if run_local_sensitivity && sensitivity_analysis["local_reference"] == "best_result" &&
       (!run_optimisation && !run_global_sensitivity && !run_parameter_variation)
        @error "Local sensitivity with local_reference='best_result' requires an optimisation, " *
               "parameter variation or a global sensitivity analysis activated."
        throw(InputError())
    end

    if (run_global_sensitivity || run_local_sensitivity) && parameter_study["objective_function"] == "multi-objective"
        @error "Sensitivity analysis currently supports only single-objective calculations."
        throw(InputError())
    end

    return prepare_parameter_study!(parameter_study, sim_params, project_config)
end

function prepare_parameter_study!(parameter_study::Dict{String,Any},
                                  sim_params::Dict{String,Any},
                                  project_config::AbstractDict{String,Any})::Dict{String,Any}
    parameter_variation = parameter_study["parameter_variation"]
    optimisation = parameter_study["optimisation"]
    sensitivity_analysis = parameter_study["sensitivity_analysis"]
    runtime = parameter_study["runtime"]

    runtime["iterator"] = [1]
    runtime["parameter_keys"] = String[]
    runtime["parameter_values"] = Vector{Vector{Float64}}()
    runtime["objective_output_spec"] = Dict{String,Any}()
    runtime["objective_params_keys"] = String[]

    sensitivity_lower_values = Float64[]
    sensitivity_upper_values = Float64[]
    parameter_bounds = Array{Float64}(undef, 0, 3)

    for (uac, params) in pairs(sort(parameter_study["parameters"]; by=lowercase))
        for (key_param, raw_def) in pairs(sort(params; by=lowercase))
            key = uac * " " * key_param
            parameter_exists = if haskey(project_config["components"], uac)
                haskey(project_config["components"][uac], key_param)
            elseif haskey(project_config, uac)
                haskey(project_config[uac], key_param)
            else
                false
            end

            if !parameter_exists
                @error "$key of parameters is not a parameter in the input file. " *
                       "Check for spelling errors or the description in the documentation."
                throw(InputError())
            end

            if !(raw_def isa AbstractDict)
                @error "The definition of parameter $key has to be an object."
                throw(InputError())
            end

            def = Dict{String,Any}(String(k) => v for (k, v) in pairs(raw_def))
            allowed_keys = Set(keys(PARAMETER_STUDY_PARAMETER_DEF))
            unsupported_keys = sort(collect(setdiff(Set(keys(def)), allowed_keys)); by=lowercase)
            if !isempty(unsupported_keys)
                @error "Unsupported parameter definition keys for $key: " *
                       join(unsupported_keys, ", ")
                throw(InputError())
            end

            for (name, param_def) in pairs(PARAMETER_STUDY_PARAMETER_DEF)
                if haskey(def, name)
                    def[name] = EnergySystems.extract_parameter(EnergySystems.Component,
                                                                def,
                                                                name,
                                                                param_def,
                                                                sim_params,
                                                                "Parameter study parameter $key")
                end
            end
            parameter_study["parameters"][uac][key_param] = def
            push!(runtime["parameter_keys"], key)

            lower_bound = NaN
            upper_bound = NaN
            start_value = haskey(def, "start") ? Float64(def["start"]) : NaN
            values = Float64[]

            if haskey(def, "start") && !isfinite(start_value)
                @error "The start value of parameter $key must be finite."
                throw(InputError())
            end

            has_bounds_min = haskey(def, "bounds_min")
            has_bounds_max = haskey(def, "bounds_max")
            if has_bounds_min != has_bounds_max
                @error "Parameter $key must define both bounds_min and bounds_max."
                throw(InputError())
            end
            if has_bounds_min
                lower_bound = Float64(def["bounds_min"])
                upper_bound = Float64(def["bounds_max"])
                if !(isfinite(lower_bound) && isfinite(upper_bound))
                    @error "The bounds_min and bounds_max of parameter $key must be finite."
                    throw(InputError())
                end
                if lower_bound >= upper_bound
                    @error "The bounds_max of parameter $key must be greater than bounds_min."
                    throw(InputError())
                end
            end

            has_values = haskey(def, "values")
            range_keys = ["range_start", "range_stop", "range_step", "range_length"]
            has_range = any(haskey(def, range_key) for range_key in range_keys)
            if has_values && has_range
                @error "Parameter $key cannot define values and a range definition at the same time."
                throw(InputError())
            end

            if has_values
                raw_values = collect(def["values"])
                if isempty(raw_values)
                    @error "The values list of parameter $key must not be empty."
                    throw(InputError())
                end
                if !all(value -> value isa Real, raw_values)
                    @error "The values list of parameter $key must contain only numbers."
                    throw(InputError())
                end
                values = Float64.(raw_values)
                if any(.!isfinite.(values))
                    @error "The values list of parameter $key must contain only finite numbers."
                    throw(InputError())
                end
            elseif has_range
                has_range_start = haskey(def, "range_start")
                has_range_stop = haskey(def, "range_stop")
                has_range_step = haskey(def, "range_step")
                has_range_length = haskey(def, "range_length")
                if !(has_range_start && has_range_stop)
                    @error "Parameter $key must define both range_start and range_stop."
                    throw(InputError())
                end
                if has_range_step == has_range_length
                    @error "Parameter $key must define exactly one of range_step and range_length."
                    throw(InputError())
                end

                range_start = Float64(def["range_start"])
                range_stop = Float64(def["range_stop"])
                if !(isfinite(range_start) && isfinite(range_stop))
                    @error "The range_start and range_stop of parameter $key must be finite."
                    throw(InputError())
                end
                if range_start >= range_stop
                    @error "The range_stop of parameter $key must be greater than range_start."
                    throw(InputError())
                end

                if has_range_step
                    range_step = Float64(def["range_step"])
                    if !isfinite(range_step)
                        @error "The range_step of parameter $key must be finite."
                        throw(InputError())
                    end
                    if range_step <= 0.0
                        @error "The range_step of parameter $key must be greater than zero."
                        throw(InputError())
                    end
                    values = Float64.(collect(range(; start=range_start,
                                                    stop=range_stop,
                                                    step=range_step)))
                else
                    range_length = Int(def["range_length"])
                    if range_length < 2
                        @error "The range_length of parameter $key must be at least 2."
                        throw(InputError())
                    end
                    values = Float64.(collect(range(; start=range_start,
                                                    stop=range_stop,
                                                    length=range_length)))
                end
            end

            has_lower = haskey(def, "sensitivity_lower")
            has_upper = haskey(def, "sensitivity_upper")
            if has_lower != has_upper
                @error "Parameter $key must define both sensitivity_lower and sensitivity_upper."
                throw(InputError())
            end

            local_lower = has_lower ? Float64(def["sensitivity_lower"]) : NaN
            local_upper = has_upper ? Float64(def["sensitivity_upper"]) : NaN
            if has_lower && !(isfinite(local_lower) && isfinite(local_upper))
                @error "The sensitivity_lower and sensitivity_upper of parameter $key must be finite."
                throw(InputError())
            end
            if has_lower && local_lower >= local_upper
                @error "The sensitivity_lower of parameter $key must be smaller than sensitivity_upper."
                throw(InputError())
            end

            push!(runtime["parameter_values"], values)
            parameter_bounds = vcat(parameter_bounds, [lower_bound upper_bound start_value])
            push!(sensitivity_lower_values, local_lower)
            push!(sensitivity_upper_values, local_upper)
        end
    end

    isempty(runtime["parameter_keys"]) && begin
        @error "At least one parameter must be defined for a parameter study."
        throw(InputError())
    end

    needs_bounds = optimisation["run_optimisation"] ||
                   sensitivity_analysis["run_global_sensitivity"]
    if needs_bounds &&
       (any(.!isfinite.(parameter_bounds[:, 1])) ||
        any(.!isfinite.(parameter_bounds[:, 2])))
        @error "Optimisation and global sensitivity require bounds_min and bounds_max " *
               "for every parameter."
        throw(InputError())
    end

    if needs_bounds
        ranges = parameter_bounds[:, 2] .- parameter_bounds[:, 1]
        if any(ranges .<= 0.0)
            @error "The bounds_max of every parameter must be greater than bounds_min."
            throw(InputError())
        end
    end

    if optimisation["run_optimisation"]
        if any(.!isfinite.(parameter_bounds[:, 3])) ||
           any(parameter_bounds[:, 3] .< parameter_bounds[:, 1]) ||
           any(parameter_bounds[:, 3] .> parameter_bounds[:, 2])
            @error "The start value of every optimisation parameter must lie within its bounds."
            throw(InputError())
        end
    end

    if parameter_variation["run_parameter_variation"] &&
       any(isempty.(runtime["parameter_values"]))
        @error "Parameter variation requires values or range_start, range_stop and " *
               "either range_step or range_length for every parameter."
        throw(InputError())
    end

    if sensitivity_analysis["run_local_sensitivity"] &&
       sensitivity_analysis["local_reference"] == "start"
        if any(.!isfinite.(parameter_bounds[:, 3]))
            @error "Local sensitivity with local_reference='start' requires a start value for every parameter."
            throw(InputError())
        end

        for (i, key) in enumerate(runtime["parameter_keys"])
            lower_value = sensitivity_lower_values[i]
            upper_value = sensitivity_upper_values[i]
            if isfinite(lower_value) &&
               !(lower_value < parameter_bounds[i, 3] < upper_value)
                @error "Parameter $key must satisfy sensitivity_lower < start < sensitivity_upper."
                throw(InputError())
            end
        end
    end

    runtime["parameter_bounds"] = parameter_bounds
    runtime["start_values"] = copy(parameter_bounds[:, 3])
    runtime["sensitivity_lower_values"] = sensitivity_lower_values
    runtime["sensitivity_upper_values"] = sensitivity_upper_values

    sensitivity_max_runs = sensitivity_analysis["max_runs"]
    runtime["sensitivity_max_runs"] = if !isnothing(sensitivity_max_runs)
        sensitivity_max_runs
    elseif optimisation["run_optimisation"] && !isnothing(optimisation["max_runs"])
        2 * optimisation["max_runs"]
    else
        200
    end

    normalised_bounds = nothing
    if optimisation["run_optimisation"]
        ranges = parameter_bounds[:, 2] .- parameter_bounds[:, 1]
        normalised_bounds = hcat(zeros(size(parameter_bounds, 1)),
                                 ones(size(parameter_bounds, 1)),
                                 (parameter_bounds[:, 3] .- parameter_bounds[:, 1]) ./ ranges)
    end

    for (func, raw_spec) in pairs(parameter_study["objective_params"])
        spec = raw_spec
        if func == "sum" || func == "mean"
            for (output_group, output_entries) in pairs(raw_spec)
                group_key = String(output_group)
                existing_entries = get!(runtime["objective_output_spec"],
                                        group_key,
                                        Any[])
                append!(existing_entries, output_entries)
            end
            spec = parse_outkeys(raw_spec)
        elseif func != "economic" && func != "emissions"
            @error "Objective parameter {$func: $raw_spec} could not be read. $func has " *
                   "to be one of 'sum', 'mean', 'economic', 'emissions'."
            throw(InputError())
        end

        for key in spec
            push!(runtime["objective_params_keys"], func * " " * key)
        end

        if func == "economic" && !sim_params["economic_parameters"]["calculate_economy"]
            @error "For the parameter study, the objective parameter `economic` is chosen, " *
                   "but the flag `calculate_economy` is set to false. Activate it in " *
                   "the `economic_parameters` section."
            throw(InputError())
        elseif func == "emissions" && !sim_params["emissions_parameters"]["calculate_emissions"]
            @error "For the parameter study, the objective parameter `emissions` is chosen, " *
                   "but the flag `calculate_emissions` is set to false. Activate it " *
                   "in the `emissions_parameters` section."
            throw(InputError())
        end
    end

    for output_entries in values(runtime["objective_output_spec"])
        unique!(output_entries)
    end
    sort!(runtime["objective_params_keys"]; by=lowercase)

    runtime["objective_function"],
    runtime["objective_function_name"],
    runtime["objective_factors"] = parse_objective_function(parameter_study["objective_function"],
                                                            runtime["objective_params_keys"],
                                                            parameter_study["objective_factors"])

    if runtime["objective_function_name"] == "multi-objective"
        runtime["N_obj"] = length(runtime["objective_params_keys"])
        runtime["objective_senses"],
        runtime["objective_signs"] = parse_objective_senses(runtime["objective_params_keys"],
                                                            parameter_study["objective_senses"])
    else
        runtime["N_obj"] = 1
        runtime["objective_senses"] = Dict{String,Symbol}()
        runtime["objective_signs"] = Float64[]
    end

    if parameter_variation["run_parameter_variation"]
        algorithm = parameter_variation["algorithm"]
        if algorithm == "product"
            runtime["iterator"] = Iterators.product(runtime["parameter_values"]...)
        elseif algorithm == "zip"
            runtime["iterator"] = zip(runtime["parameter_values"]...)
        elseif startswith(algorithm, "random_")
            parts = split(algorithm, "_")
            if length(parts) != 2 || isnothing(tryparse(Int, parts[2]))
                @error "Algorithm $algorithm is not supported for parameter variation. " *
                       "Use product, zip or random_*, where * is an integer."
                throw(InputError())
            end
            iterator = Iterators.product(runtime["parameter_values"]...)
            n_samples = min(parse(Int, parts[2]), length(iterator))
            runtime["iterator"] = rand(collect(iterator), n_samples)
        else
            @error "Algorithm $algorithm is not supported for parameter variation. " *
                   "Use product, zip or random_*, where * is an integer."
            throw(InputError())
        end
    elseif optimisation["run_optimisation"]
        configure_optimiser_backend!(runtime, optimisation, normalised_bounds)
    end

    return parameter_study
end

function get_refinement_optimiser_config(config::Union{Nothing,AbstractDict},
                                         sim_params::Dict{String,Any})
    isnothing(config) && return nothing

    refinement = extract_parameter_section(config,
                                           REFINEMENT_OPTIMISATION_DEF,
                                           sim_params,
                                           "Refinement optimisation")
    if isnothing(refinement["optim_kwargs"])
        refinement["optim_kwargs"] = Dict{String,Any}()
    end

    return refinement
end

function configure_optimiser_backend!(optimiser, optimiser_config, normalised_bounds)
    optim_kwargs = something(get(optimiser_config, "optim_kwargs", nothing), Dict{String,Any}())

    # A configured zero disables or prevents convergence in the supported backends.
    function reject_duplicate_tolerance(common_name::String,
                                        backend_name::String,
                                        backend::String)::Nothing
        if !isnothing(get(optimiser_config, common_name, nothing)) && haskey(optim_kwargs, backend_name)
            @error "$backend option '$backend_name' cannot be set in optim_kwargs together " *
                   "with the corresponding common setting '$common_name'."
            throw(InputError())
        end
        return nothing
    end

    if optimiser_config["type"] == "Optim"
        if optimiser["objective_function_name"] == "multi-objective"
            @error "Objective function multi-objective not supported for algorithms from " *
                   "package 'Optim'"
            throw(InputError())
        end
        optimiser["args"] = []

        #TODO most Optim algorithms ignore bounds which can be supposedly added with wrapper 
        # Optim.Fminbox() but it doesn't work
        if optimiser_config["algorithm"] == "NelderMead"
            if !isnothing(optimiser_config["x_tol_abs"])
                @error "x_tol_abs is not supported by Optim.NelderMead. " *
                       "Use max_runs or max_time instead."
                throw(InputError())
            end
            reject_duplicate_tolerance("f_tol_abs", "g_abstol", "Optim")
            alg = Optim.NelderMead()

        elseif optimiser_config["algorithm"] == "SAMIN"
            samin_kwargs = Dict{Symbol,Any}()
            if !isnothing(optimiser_config["x_tol_abs"])
                samin_kwargs[:x_tol] = optimiser_config["x_tol_abs"]
            end
            if !isnothing(optimiser_config["f_tol_abs"])
                samin_kwargs[:f_tol] = optimiser_config["f_tol_abs"]
            end
            alg = Optim.SAMIN(; samin_kwargs...)
            push!(optimiser["args"], normalised_bounds[:, 1])
            push!(optimiser["args"], normalised_bounds[:, 2])

        elseif optimiser_config["algorithm"] == "ParticleSwarm"
            if !isnothing(optimiser_config["x_tol_abs"]) ||
               !isnothing(optimiser_config["f_tol_abs"])
                @error "x_tol_abs and f_tol_abs are not supported by " *
                       "Optim.ParticleSwarm because this algorithm does not assess " *
                       "convergence. Use max_runs or max_time instead."
                throw(InputError())
            end
            alg = Optim.ParticleSwarm(; upper=normalised_bounds[:, 2], lower=normalised_bounds[:, 1])
        else
            @error "For optimisation type 'Optim' the algorithm has to be one of " *
                   "'NelderMead', 'SAMIN', 'ParticleSwarm'."
            throw(InputError())
        end

        push!(optimiser["args"], normalised_bounds[:, 3])
        push!(optimiser["args"], alg)

        optimiser["kwargs"] = Dict{Symbol,Any}()
        optimiser["kwargs"][:show_trace] = true
        for (keyword, val) in pairs(optim_kwargs)
            optimiser["kwargs"][Symbol(keyword)] = val
        end
        if !isnothing(optimiser_config["max_runs"])
            optimiser["kwargs"][:f_calls_limit] = optimiser_config["max_runs"]
        end
        if !isnothing(optimiser_config["max_time"])
            optimiser["kwargs"][:time_limit] = optimiser_config["max_time"]
        end
        if optimiser_config["algorithm"] == "NelderMead" && !isnothing(optimiser_config["f_tol_abs"])
            # Nelder-Mead uses the spread of objective values at the simplex vertices
            # as its convergence criterion, exposed through Optim.Options.g_abstol.
            optimiser["kwargs"][:g_abstol] = optimiser_config["f_tol_abs"]
        end

        push!(optimiser["args"], Optim.Options(; optimiser["kwargs"]...))

    elseif optimiser_config["type"] == "BlackBoxOptim"
        if !isnothing(optimiser_config["x_tol_abs"])
            @error "x_tol_abs is not supported by BlackBoxOptim. " *
                   "Use max_runs or max_time instead."
            throw(InputError())
        end
        reject_duplicate_tolerance("f_tol_abs", "MinDeltaFitnessTolerance", "BlackBoxOptim")

        alg = Symbol(optimiser_config["algorithm"])

        optimiser["args"] = [normalised_bounds[:, 3]]
        optimiser["kwargs"] = Dict{Symbol,Any}()

        if optimiser["objective_function_name"] == "multi-objective"
            if optimiser_config["algorithm"] != "borg_moea"
                @error "Optimisation algorithm '$(optimiser_config["algorithm"])' doesn't " *
                       "support multi-objective optimisation. Choose a different " *
                       "objective_function or algorithm 'borg_moea'."
                throw(InputError())
            end
            if !isnothing(optimiser_config["f_tol_abs"])
                @error "f_tol_abs is not supported for multi-objective BlackBoxOptim."
                throw(InputError())
            end
            optimiser["kwargs"][:FitnessScheme] = BlackBoxOptim.ParetoFitnessScheme{optimiser["N_obj"]}(;
                                                                                                        is_minimizing=true)
        end

        optimiser["kwargs"][:Method] = alg
        optimiser["kwargs"][:SearchRange] = Tuple.(eachrow(normalised_bounds[:, 1:2]))
        optimiser["kwargs"][:NumDimensions] = size(normalised_bounds, 1)
        n_evaluation_threads = Threads.nthreads(:default) - 1
        if n_evaluation_threads > 1
            optimiser["kwargs"][:NThreads] = n_evaluation_threads
            optimiser["kwargs"][:PopulationSize] = max(20, 4 * n_evaluation_threads)
        end
        if !isnothing(optimiser_config["max_runs"])
            optimiser["kwargs"][:MaxFuncEvals] = optimiser_config["max_runs"]
        end
        if !isnothing(optimiser_config["max_time"])
            optimiser["kwargs"][:MaxTime] = optimiser_config["max_time"]
        end
        for (keyword, val) in pairs(optim_kwargs)
            optimiser["kwargs"][Symbol(keyword)] = val
        end
        if !isnothing(optimiser_config["f_tol_abs"])
            # Stop when two consecutive best-fitness improvements differ by less
            # than the configured objective tolerance.
            optimiser["kwargs"][:MinDeltaFitnessTolerance] = optimiser_config["f_tol_abs"]
        end

    elseif optimiser_config["type"] == "Metaheuristics"
        if !isnothing(optimiser_config["x_tol_abs"]) || !isnothing(optimiser_config["f_tol_abs"])
            @error "x_tol_abs and f_tol_abs are not currently supported for " *
                   "Metaheuristics algorithms. Metaheuristics combines its parameter and " *
                   "objective tolerances with additional mandatory population-convergence " *
                   "conditions, so they cannot be translated to independent stopping " *
                   "criteria. Use max_runs or max_time instead."
            throw(InputError())
        end
        m_obj_algs = ["MOEAD_DE", "NSGA2", "NSGA3", "SMS_EMOA", "SPEA2", "CCMO"]
        if optimiser["objective_function_name"] == "multi-objective" &&
           !(optimiser_config["algorithm"] in m_obj_algs)
            @error "Optimisation algorithm '$(optimiser_config["algorithm"])' doesn't " *
                   "support multi-objective optimisation. Choose a different " *
                   "objective_function or algorithm, e.g. one of $(join(string.(m_obj_algs), ", "))."
            throw(InputError())
        end

        alg = getproperty(Metaheuristics, Symbol(optimiser_config["algorithm"]))

        optimiser["args"] = Any[[normalised_bounds[:, 1] normalised_bounds[:, 2]]']

        args_alg = Any[]
        kwargs_general = Dict{Symbol,Any}()
        kwargs_alg = Dict{Symbol,Any}()

        if optimiser_config["algorithm"] == "MOEAD_DE"
            if optimiser["N_obj"] == 1
                @error "Algorithm 'MOEAD_DE' can only be used with " *
                       "objective_function='multi-objective'"
                throw(InputError())
            end

            n_partitions = get(optim_kwargs, "n_partitions", 12)
            reference_directions = Metaheuristics.gen_ref_dirs(optimiser["N_obj"], n_partitions)

            push!(args_alg, reference_directions)
        end

        if !isnothing(optimiser_config["max_runs"]) && optimiser_config["max_runs"] < 100 && :N in fieldnames(alg)
            @warn "'max_runs' is smaller than the default population size of 100. This may lead to poor results."
            kwargs_alg[:N] = max(2, ceil(Int, optimiser_config["max_runs"] / 4))
        end

        if Threads.nthreads() > 1
            kwargs_general[:parallel_evaluation] = true
        end
        if !isnothing(optimiser_config["max_runs"])
            kwargs_general[:f_calls_limit] = optimiser_config["max_runs"]
        end
        if !isnothing(optimiser_config["max_time"])
            kwargs_general[:time_limit] = optimiser_config["max_time"]
        end

        for (keyword, val) in pairs(optim_kwargs)
            # This is consumed above when constructing MOEA/D weights.
            keyword == "n_partitions" && continue

            key = Symbol(keyword)

            if key in fieldnames(Metaheuristics.Options)
                kwargs_general[key] = val
            elseif key in fieldnames(alg)
                kwargs_alg[key] = val
            else
                @warn "Unknown Metaheuristics option '$keyword' for " *
                      "algorithm '$(optimiser_config["algorithm"])'."
            end
        end

        options = Metaheuristics.Options(; kwargs_general...)
        if optimiser_config["algorithm"] == "CCMO"
            push!(optimiser["args"], alg(Metaheuristics.NSGA2(args_alg...; kwargs_alg...); options=options))
        else
            push!(optimiser["args"], alg(args_alg...; kwargs_alg..., options=options))
        end

    elseif optimiser_config["type"] == "NLopt"
        configured_algorithm = optimiser_config["algorithm"]
        algorithm = startswith(configured_algorithm, "NLOPT_") ? configured_algorithm[7:end] : configured_algorithm

        if startswith(algorithm, "LD_") || startswith(algorithm, "GD_")
            @error "The chosen algorithm `$configured_algorithm` needs a gradient which " *
                   "is not supported in ReSiE"
            throw(InputError())
        end

        reject_duplicate_tolerance("x_tol_abs", "xtol_abs", "NLopt")
        reject_duplicate_tolerance("f_tol_abs", "ftol_abs", "NLopt")

        n_dim_opt = length(normalised_bounds[:, 1])
        alg = NLopt.Opt(Symbol(algorithm), n_dim_opt)
        optimiser["kwargs"] = Dict{Symbol,Any}()

        optimiser["kwargs"][:lower_bounds] = normalised_bounds[:, 1]
        optimiser["kwargs"][:upper_bounds] = normalised_bounds[:, 2]
        if !isnothing(optimiser_config["max_runs"])
            optimiser["kwargs"][:maxeval] = optimiser_config["max_runs"]
        end
        if !isnothing(optimiser_config["max_time"])
            optimiser["kwargs"][:maxtime] = optimiser_config["max_time"]
        end
        if !isnothing(optimiser_config["x_tol_abs"])
            optimiser["kwargs"][:xtol_abs] = optimiser_config["x_tol_abs"]
        end
        if !isnothing(optimiser_config["f_tol_abs"])
            optimiser["kwargs"][:ftol_abs] = optimiser_config["f_tol_abs"]
        end

        for (keyword, val) in pairs(optim_kwargs)
            if Symbol(keyword) in propertynames(alg)
                optimiser["kwargs"][Symbol(keyword)] = val
            elseif NLopt.nlopt_has_param(alg, keyword)
                NLopt.nlopt_set_param(alg, keyword, val)
            end
        end

        for (keyword, val) in pairs(optimiser["kwargs"])
            NLopt.setproperty!(alg, keyword, val)
        end
        optimiser["args"] = [alg, normalised_bounds[:, 3]]

    elseif optimiser_config["type"] == "NOMAD"
        if !isnothing(optimiser_config["f_tol_abs"])
            @error "f_tol_abs is not supported by NOMAD. " *
                   "Use max_runs or max_time instead."
            throw(InputError())
        end

        optimiser["kwargs"] = Dict{Symbol,Any}()
        optimiser["kwargs"][:lower_bound] = normalised_bounds[:, 1]
        optimiser["kwargs"][:upper_bound] = normalised_bounds[:, 2]

        if !isnothing(optimiser_config["x_tol_abs"])
            optimiser["kwargs"][:min_mesh_size] = fill(optimiser_config["x_tol_abs"], size(normalised_bounds, 1))
        end

        kwargs_general = Dict{Symbol,Any}()
        if !isnothing(optimiser_config["max_runs"])
            kwargs_general[:max_bb_eval] = optimiser_config["max_runs"]
        end
        if !isnothing(optimiser_config["max_time"])
            kwargs_general[:max_time] = optimiser_config["max_time"]
        end

        for (keyword, val) in pairs(optim_kwargs)
            if Symbol(keyword) in fieldnames(NOMAD.NomadOptions)
                kwargs_general[Symbol(keyword)] = val
            end
        end

        optimiser["kwargs"][:options] = NOMAD.NomadOptions(; kwargs_general...)
        optimiser["args"] = [size(normalised_bounds, 1), optimiser["N_obj"],
                             fill("OBJ", optimiser["N_obj"]), normalised_bounds[:, 3]]
    end
end

function load_refinement_optimiser(base_parameter_study::Dict{String,Any},
                                   stage_config::AbstractDict,
                                   start_values::AbstractVector{<:Real})
    runtime = base_parameter_study["runtime"]
    n_inputs = length(runtime["parameter_keys"])
    normalised_start = clamp.(Float64.(start_values), 0.0, 1.0)
    normalised_bounds = hcat(zeros(n_inputs), ones(n_inputs), normalised_start)

    stage = deepcopy(base_parameter_study)
    stage["optimisation"] = deepcopy(stage_config)
    stage["optimisation"]["run_optimisation"] = true
    stage["optimisation"]["refinement"] = nothing

    stage_runtime = stage["runtime"]
    pop!(stage_runtime, "args", nothing)
    pop!(stage_runtime, "kwargs", nothing)

    try
        configure_optimiser_backend!(stage_runtime,
                                     stage["optimisation"],
                                     normalised_bounds)
    catch e
        @globalInfo "Loading of refinement algorithm was not successful. The generation of the results will be " *
                    "continued without a refinement. Check the inputs.\n" *
                    "  The following error occurred: $e"
        return nothing
    end
    return stage
end

function parse_objective_function(eff_def::String,
                                  objective_keys::Vector{String},
                                  configured_factors)::Tuple{Function,String,Dict{String,Float64}}
    method = lowercase(strip(eff_def))
    has_configured_factors = configured_factors !== nothing && !isempty(configured_factors)

    parsed_factors = Dict{String,Float64}()

    if method == "sum"
        f = x -> sum(Float64.(x))
    elseif method == "linear"
        if !has_configured_factors
            @error("objective_function=\"linear\" requires objective_factors.")
            throw(InputError())
        end

        configured_keys = String.(collect(keys(configured_factors)))
        missing_keys = setdiff(objective_keys, configured_keys)
        unknown_keys = setdiff(configured_keys, objective_keys)

        if !isempty(missing_keys)
            @error("objective_factors is missing coefficients for: " *
                   join(missing_keys, ", "))
            throw(InputError())
        end

        if !isempty(unknown_keys)
            @error("objective_factors contains unknown objective keys: " *
                   join(unknown_keys, ", "))
            throw(InputError())
        end

        for key in objective_keys
            raw_factor = configured_factors[key]

            parsed_factors[key] = try
                raw_factor isa Number ?
                Float64(raw_factor) :
                parse(Float64, String(raw_factor))
            catch
                @error("The objective factor for \"$key\" must be numeric, got " *
                       "\"$raw_factor\".")
                throw(InputError())
            end
        end

        # The vector is derived from explicit names. Its order only follows the
        # canonical objective key order used by Resie for objective_values.
        ordered_factors = Float64[parsed_factors[key] for key in objective_keys]

        f = x -> begin
            values = Float64.(x)
            if length(values) != length(ordered_factors)
                throw(ArgumentError("Objective value count ($(length(values))) does not match " *
                                    "factor count ($(length(ordered_factors)))."))
            end
            sum(values .* ordered_factors)
        end

    elseif method == "multi-objective"
        f = x -> collect(Float64.(x))
    else
        @error("Cannot parse objective function from: $eff_def. Has to be one of " *
               "'sum', 'linear', 'multi-objective'.")
        throw(InputError())
    end

    return f, method, parsed_factors
end

function parse_objective_senses(objective_keys::Vector{String},
                                configured_senses)::Tuple{Dict{String,Symbol},Vector{Float64}}
    if configured_senses === nothing
        senses = Dict{String,Symbol}(
            key => :min
            for key in objective_keys
        )
    else
        normalized = Dict{String,Any}(
            String(key) => value
            for (key, value) in pairs(configured_senses)
        )

        configured_keys = collect(keys(normalized))
        missing_keys = setdiff(objective_keys, configured_keys)
        unknown_keys = setdiff(configured_keys, objective_keys)

        if !isempty(missing_keys)
            @error ("objective_senses is missing directions for: " * join(missing_keys, ", "))
            throw(InputError())
        end

        if !isempty(unknown_keys)
            @error ("objective_senses contains unknown objective keys: " * join(unknown_keys, ", "),)
            throw(InputError())
        end

        senses = Dict{String,Symbol}()

        for key in objective_keys
            sense = Symbol(lowercase(strip(String(normalized[key]))))

            if !(sense in (:min, :max))
                @error("The objective sense for \"$key\" must be " * "`min` or `max`, got \"$(normalized[key])\".",)
                throw(InputError())
            end

            senses[key] = sense
        end
    end

    signs = Float64[senses[key] == :min ? 1.0 : -1.0
                    for key in objective_keys]

    return senses, signs
end
