# this file contains functionality pertaining to loading a project's metadata and the
# energy system components from the project config file, as well as constructing certain
# helpful information data structures from the inputs in the config
using JSON: JSON
using OrderedCollections: OrderedDict
using Logging

const HOURS_PER_SECOND::Float64 = 1.0 / 3600.0
const SECONDS_PER_HOUR::Float64 = 3600.0

"""
Shared preparation cache used during optimisation.

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
    optimiser = get(sim_params, "optimisation", Dict{String,Any}())

    if !haskey(optimiser, "optim_params_keys")
        return true
    end

    # Conservative blocklist. If one of these parameters is optimised, the component graph
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

    for key in optimiser["optim_params_keys"]
        parts = split(key, " ")
        if parts[end] in structural_params
            return false
        end
    end

    return true
end

function operation_cache_key(project_config::AbstractDict{String,Any},
                             sim_params::Dict{String,Any})::String
    optimiser = get(sim_params, "optimisation", Dict{String,Any}())

    component_cfg = deepcopy(project_config["components"])

    # Non-structural optimised component values should not invalidate the
    # operation-order cache. 
    if haskey(optimiser, "optim_params_keys")
        for opt_key in optimiser["optim_params_keys"]
            uac, param_key = split(opt_key, " ")
            if haskey(component_cfg, uac) && haskey(component_cfg[uac], param_key)
                component_cfg[uac][param_key] = "__OPTIMISED_NONSTRUCTURAL_VALUE__"
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
    "output_optimisation_csv" => (
        default=true,
        description="Toggle if a csv with the optimisation results should be created",
        display_name="Output optimisation csv?",
        required=false,
        type=Bool,
        json_type="boolean",
        unit="-"
    ),
    "optimisation_csv_file_path" => (
        default="./output/optim_results.csv",
        description="File path to where the optimisation results are written to csv",
        display_name="optimisation csv file path",
        required=false,
        type=String,
        json_type="string",
        unit="-"
    ),
    "write_optimisation_csv_continuously" => (
        default=false,
        description="Toggle if csv output of optimisation will be written continuously, " * 
                    "meaning after every run. Activating this functionality will ensure " * 
                    "partial output if the optimisation is stopped during execution. " *
                    "It incurs a slight performance penalty depending on the run time of " *
                    "one simulation and the number of parallel runs, since the threads " *
                    "might have to wait for write access.", 
        display_name="Write optimisation csv continuously?",
        required=false,
        type=Bool,
        json_type="boolean",
        unit="-"
    ),   
    "matrix_plot" => (
        default="default",
        description="Sets the mode of the matrix plot, switching between several default " *
                    "and custom behaviour modes as well an option of not creating a plot " *
                    "file at all.",
        display_name="Matrix plot mode",
        options=["custom", "default", "nothing"],
        required=false,
        type=String,
        json_type="string",
        unit="-"
    ),
    "matrix_plot_spec" => (
        default=nothing,
        description="Specification of the objective in the matrix plot in custom mode. " *
                    "Has same structure as objective_params.",
        display_name="Matrix plot specification",
        required=false,
        conditionals=[("matrix_plot", "is", "custom")],
        type=Dict{String,Any},
        json_type="object",
        unit="-"
    ),
    "optim_plots_file_path" => (
        default="./output/optim_plots",
        description="File path to where the optimisation result plots will be written",
        display_name="Optim plots file path",
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
    ("csv_output", "csv_output_keys"),
    ("matrix_plot", "matrix_plot_spec")
]

OPTIMISATION_PARAMATERS_DEF = Dict{String,Any}(
    "run_optimisation" => (
        default=false,
        description="If set to true, executes multiple runs with the chosen optimisation" *
                    "algorithm.",
        display_name="Run optimisation?",
        required=false,
        type=Bool,
        json_type="boolean",
        unit="-"
    ),
    #TODO remove type and only use algorithm to choose package internally
    "type" => (
        default="nothing",
        description="Sets type of optimisation algorithm",
        display_name="optimisation type",
        required=false,
        conditionals=["run_optimisation", "is", true],
        options=["parametervariation", "Optim", "BlackBoxOptim", "Metaheuristics", "NLopt", "NOMAD"],
        type=String,
        json_type="string",
        unit="-"
    ),
    "optim_params" => (
        default=nothing,
        description="Defines which parameters of which components should be varied in " *
                    "the optimisation. Definition follows the definition of components. " *
                    "See the documentation for more details.",
        display_name="optimisation parameters",
        required=false,
        conditionals=[("type", "is_not_nothing")],
        type=Dict{String,Any},
        json_type="object",
        unit="-"
    ),
    "algorithm" => (
        default="NelderMead",
        description="Optimisation algorithm that is used from packages `Optim` and " *
                    "`BlackBoxOptim`",
        display_name="Optimisation algorithm",
        required=false,
        conditionals=[("type", "is_one_of", 
                       ("Optim", "BlackBoxOptim", "Metaheuristics", "NLopt", "parametervariation")
                       )],
        type=String,
        json_type="string",
        unit="-"
    ),
    #TODO discuss if this should get an option to choose results from a file from optimisation
    "run_sensitivity" => (
        default=false,
        description="If set to true, a sensitivity analysis will be run either partially " * 
                    "reusing results from optimisation if available or running seperately.",
        display_name="Run sensitivity analysis",
        required=false,
        conditionals=[("objective_function", "is_one_of", ("sum", "linear"))],
        type=Bool,
        json_type="boolean",
        unit="-"
    ),
    "objective_params" => (
        default=nothing,
        description="Defines which parameters are set as the objective. Definition " *
                    "follows the definition of components. See the documentation for " *
                    "more details.",
        display_name="Objective parameters",
        required=true,
        conditionals=[("run_optimisation", "is", true)],
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
        conditionals=[
            ("objective_function", "is", "multi-objective"),
        ],
        type=Dict{String,Any},
        json_type="object",
        unit="-",
    ),
    "disable_all_simulation_outputs" => (
        default=true,
        description="Disables all simulation outputs written to the hard drive during optimisation.",
        display_name="Disable all simulation outputs",
        required=false,
        type=Bool,
        json_type="boolean",
        unit="-"
    ),
    "max_runs" => (
        default=nothing,
        description="Set the maximum number of runs to be executed by the algorithm.",
        display_name="Max runs",
        required=false,
        type=Int64,
        json_type="number",
        unit="-"
    ),
    "max_time" => (
        default=nothing,
        description="Set the maximum time in seconds before the optimisation stops.",
        display_name="Max time",
        required=false,
        type=Int64,
        json_type="number",
        unit="-"
    ),
    "x_tol_abs" => (
        default=nothing,
        description="Absolute tolerance for the normalised optimisation parameters (`optim_params`) in the range [0,1]",
        display_name="Absolute tolerance of normalised optimisation parameters",
        required=false,
        conditionals=[("type", "is_one_of", ("Optim", "NLopt", "NOMAD"))],
        type=Float64,
        json_type="number",
        unit="-"
    ),
    "f_tol_abs" => (
        default=nothing,
        description="Absolute tolerance for the objective function",
        display_name="Absolute tolerance objective function",
        required=false,
        conditionals=[("type", "is_one_of", ("Optim", "NLopt"))],
        type=Float64,
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
- for optimisation runs, allow only logs >= GlobalInfo
- for single simulations, set min log level to Info level

Therefore, we need to have early-access to the flag "run_optimisation".

# Arguments
- `project_config::OrderedDict{String,Any}`: Base project configuration used to generate
  simulation variants.
- `logger::Union{Nothing,Resie_Logger.CustomLogger}`: Logger used for ReSiE
# Returns
- `min_log_level::LogLevel`: minimum log level used in ReSiE for logging
"""
function get_min_log_level(project_config, logger)
    optimisation_config = get(project_config, "optimisation_parameters", nothing)
    run_optimisation = optimisation_config isa AbstractDict &&
                       get(optimisation_config, "run_optimisation", false) == true
    if logger !== nothing
        if run_optimisation
            # for optimisation runs, allow only logs >= GlobalInfo
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
    sim_params["optimisation"] = get_optimisation_parameters(project_config, sim_params)

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

function get_optimisation_parameters(project_config::AbstractDict{String,Any},
                                     sim_params::Dict{String,Any})::Dict{String,Any}
    if !haskey(project_config, "optimisation_parameters") ||
       !get(project_config["optimisation_parameters"], "run_optimisation", false)
        # no optimisation
        return Dict{String,Any}("run_optimisation" => false)
    end
    optimiser_config = Dict{String,Any}()
    for (name, param_def) in pairs(OPTIMISATION_PARAMATERS_DEF)
        optimiser_config[name] = EnergySystems.extract_parameter(EnergySystems.Component,
                                                                 project_config["optimisation_parameters"],
                                                                 name,
                                                                 param_def,
                                                                 sim_params,
                                                                 "optimisation parameters")
    end

    EnergySystems.validate_config(EnergySystems.Component, optimiser_config, "Optimisation parameters",
                                  sim_params, OPTIMISATION_PARAMATERS_DEF)

    optimiser = load_optimiser(optimiser_config, sim_params, project_config)

    return optimiser
end

"""
    all_general_parameters()::Dict{String,Any}

Lists all general, non-component parameters for the simulation engine.

# Returns
-`Dict{String,Any}`: The parameter definitions
"""
function all_general_parameters()::Dict{String,Any}
    return Dict{String,Any}(
        "simulation" => SIMULATION_PARAMETERS_DEF,
        "io_settings" => IO_SETTINGS_DEF,
        "economic" => ECONOMIC_PARAMETERS_DEF,
        "emissions" => EMISSIONS_PARAMATERS_DEF,
        "optimisation" => OPTIMISATION_PARAMATERS_DEF,
    )
end

# calculation of the order of operations has its own include files due to its complexity
include("order_of_operations.jl")

function load_optimiser(optimiser_config::Dict{String,Any}, sim_params::Dict{String,Any},
                        project_config::AbstractDict{String,Any})::Dict{String,Any}
    optimiser = Dict{String,Any}()
    optimiser["type"] = optimiser_config["type"]
    optimiser["run_optimisation"] = true
    optimiser["iterator"] = [1]
    optimiser["run_sensitivity"] = optimiser_config["run_sensitivity"]
    optimiser["disable_all_simulation_outputs"] = optimiser_config["disable_all_simulation_outputs"]

    # read and parse optim_params in Arrays to preserve order
    optimiser["optim_params_keys"] = String[]
    optimiser["optim_params_values"] = []
    # Matrix with bounds with colums being lower_bound, upper_bound, start_value
    bounds = Array{Float64}(undef, 0, 3)
    for (uac, params) in pairs(sort(optimiser_config["optim_params"]; by=lowercase))
        for (key_param, def) in pairs(sort(params; by=lowercase))
            key = uac * " " * key_param
            if !(uac in keys(project_config["components"]) ||
                 uac in keys(project_config) ||
                 key_param in keys(project_config["components"][uac]) ||
                 key_param in keys(project_config[uac]))
                # end of expression
                @error "$key of optim_params is not a parameter in the input file. " *
                       "Check for spelling errors or the description in the documentation."
                throw(InputError())
            end
            push!(optimiser["optim_params_keys"], key)
            if haskey(def, "values")
                push!(optimiser["optim_params_values"], def["values"])
                bounds = vcat(bounds,
                              [minimum(def["values"]) maximum(def["values"]) (minimum(def["values"]) +
                                                                              maximum(def["values"])) / 2])
            elseif haskey(def, "min") && haskey(def, "max")
                values = range(; start=def["min"], stop=def["max"], length=100)
                push!(optimiser["optim_params_values"], values)
                start_val = ifelse(haskey(def, "start"), def["start"], (def["min"] + def["max"]) / 2)
                bounds = vcat(bounds, [def["min"] def["max"] start_val])
            else
                def = Dict(Symbol(k) => v for (k, v) in def)
                values = range(; def...)
                push!(optimiser["optim_params_values"], values)
                bounds = vcat(bounds, [minimum(values) maximum(values) (minimum(values) + maximum(values)) / 2])
            end
        end
    end

    # check start values
    if any(bounds[:, 3] .< bounds[:, 1]) || any(bounds[:, 3] .> bounds[:, 2])
        @error "The start value of every optimisation parameter must lie within its bounds."
        throw(InputError())
    end
    optimiser["bounds"] = bounds

    # normalise bounds and start values to [0,1]
    ranges = bounds[:, 2] .- bounds[:, 1]
    if any(ranges .<= 0.0)
        @error "The maximum of every optimisation parameter must be greater than its minimum."
        throw(InputError())
    end
    normalised_bounds = hcat(zeros(size(bounds, 1)), ones(size(bounds, 1)), (bounds[:, 3] .- bounds[:, 1]) ./ ranges)

    # read and parse objective_params
    optimiser["objective_keys_sum_mean"] = Dict{String,Any}()
    optimiser["objective_params_keys"] = String[]
    if !isnothing(optimiser_config["objective_params"])
        for (func, spec) in pairs(optimiser_config["objective_params"])
            if func == "sum" || func == "mean"
                spec = parse_outkeys(spec)
            elseif func != "economic" && func != "emissions"
                @error "Objective parameter {$func: $value} could not be read. $func has " *
                       "to be one of 'sum', 'mean', 'economic', 'emissions'."
                throw(InputError())
            end

            for key in spec
                push!(optimiser["objective_params_keys"], func * " " * key)
            end

            if func == "economic" && !sim_params["economic_parameters"]["calculate_economy"]
                @error "For optimisation, the objective parameter `economic` is chosen, " *
                       "but the flag `calculate_economy` is set to false. Activate it in " *
                       "the `economic_parameters` section."
                throw(InputError())
            elseif func == "emissions" && !sim_params["emissions_parameters"]["calculate_emissions"]
                @error "For optimisation, the objective parameter `emissions` is chosen, " *
                       "but the flag `calculate_emissions` is set to false. Activate it " *
                       "in the `emissions_parameters` section."
                throw(InputError())
            end
        end

        # Canonical internal order. Linear factor assignment does not depend on this
        # order because factors are matched by objective name.
        sort!(optimiser["objective_params_keys"]; by=lowercase)

        optimiser["objective_function"],
        optimiser["objective_function_name"],
        optimiser["objective_factors"] = parse_objective_function(optimiser_config["objective_function"],
                                                                  optimiser["objective_params_keys"],
                                                                  optimiser_config["objective_factors"])

        if optimiser["objective_function_name"] == "multi-objective"
            optimiser["N_obj"] = length(optimiser["objective_params_keys"])
            optimiser["objective_senses"],
            optimiser["objective_signs"] = parse_objective_senses(optimiser["objective_params_keys"],
                                                                  optimiser_config["objective_senses"])
        else
            optimiser["N_obj"] = 1
            optimiser["objective_senses"] = Dict{String,Symbol}()
            optimiser["objective_signs"] = Float64[]
        end
    end

    if !isnothing(optimiser_config["max_runs"])
        optimiser["max_runs"] = optimiser_config["max_runs"]
    end

    if optimiser_config["type"] == "parametervariation"
        if optimiser_config["algorithm"] == "product"
            optimiser["iterator"] = Iterators.product(optimiser["optim_params_values"]...)
        elseif optimiser_config["algorithm"] == "zip"
            optimiser["iterator"] = zip(optimiser["optim_params_values"]...)
        elseif split(optimiser_config["algorithm"], "_")[1] == "random"
            iter = Iterators.product(optimiser["optim_params_values"]...)
            n_samples = min(parse(Int, split(optimiser_config["algorithm"], "_")[2]),
                            length(iter))
            optimiser["iterator"] = rand(collect(iter), n_samples)
        else
            @error "Algorithm $(optimiser_config["algorithm"]) is not supported for type " *
                   "`parametervariation`. Has to be one of `product`, `zip` or " *
                   "`random_*`, where * is a integer]" *
                   throw(InputError())
        end

    elseif optimiser_config["type"] == "Optim"
        if optimiser["objective_function_name"] == "multi-objective"
            @error "Objective function multi-objective not supported for algorithms from " *
                   "package 'Optim'"
            throw(InputError())
        end
        optimiser["args"] = []

        #TODO most Optim algorithms ignore bounds which can be supposedly added with wrapper 
        # Optim.Fminbox() but it doesn't work
        if optimiser_config["algorithm"] == "NelderMead"
            alg = Optim.NelderMead()

        elseif optimiser_config["algorithm"] == "SAMIN"
            alg = Optim.SAMIN()
            push!(optimiser["args"], normalised_bounds[:, 1])
            push!(optimiser["args"], normalised_bounds[:, 2])

        elseif optimiser_config["algorithm"] == "ParticleSwarm"
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
        if haskey(optimiser_config, "optim_kwargs")
            for (keyword, val) in pairs(optimiser_config["optim_kwargs"])
                optimiser["kwargs"][Symbol(keyword)] = val
            end
        end
        if !isnothing(optimiser_config["max_runs"])
            optimiser["kwargs"][:f_calls_limit] = optimiser_config["max_runs"]
        end
        if !isnothing(optimiser_config["max_time"])
            optimiser["kwargs"][:time_limit] = optimiser_config["max_time"]
        end
        if !isnothing(optimiser_config["x_tol_abs"])
            optimiser["kwargs"][:x_abstol] = optimiser_config["x_tol_abs"]
        end
        if !isnothing(optimiser_config["f_tol_abs"])
            optimiser["kwargs"][:f_abstol] = optimiser_config["f_tol_abs"]
        end

        push!(optimiser["args"], Optim.Options(; optimiser["kwargs"]...))

    elseif optimiser_config["type"] == "BlackBoxOptim"
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
            optimiser["kwargs"][:FitnessScheme] = BlackBoxOptim.ParetoFitnessScheme{optimiser["N_obj"]}(;
                                                                                                        is_minimizing=true)
        end

        optimiser["kwargs"][:Method] = alg
        optimiser["kwargs"][:SearchRange] = Tuple.(eachrow(normalised_bounds[:, 1:2]))
        optimiser["kwargs"][:NumDimensions] = size(normalised_bounds, 1)
        optimiser["kwargs"][:NThreads] = Threads.nthreads() - 1
        if !isnothing(optimiser_config["max_runs"])
            optimiser["kwargs"][:MaxFuncEvals] = optimiser_config["max_runs"]
        end
        if !isnothing(optimiser_config["max_time"])
            optimiser["kwargs"][:MaxTime] = optimiser_config["max_time"]
        end
        if haskey(optimiser_config, "optim_kwargs")
            for (keyword, val) in pairs(optimiser_config["optim_kwargs"])
                optimiser["kwargs"][Symbol(keyword)] = val
            end
        end

    elseif optimiser_config["type"] == "Metaheuristics"
        m_obj_algs = ["MOEAD_DE", "NSGA2", "NSGA3", "SMS_EMOA", "SPEA2", "CCMO"]
        if optimiser["objective_function_name"] == "multi-objective" && !(optimiser_config["algorithm"] in m_obj_algs)
            @error "Optimisation algorithm '$(optimiser_config["algorithm"])' doesn't " *
                   "support multi-objective optimisation. Choose a different " *
                   "objective_function or algorithm."
            throw(InputError())
        end

        alg = getproperty(Metaheuristics, Symbol(optimiser_config["algorithm"]))

        optimiser["args"] = Any[[normalised_bounds[:, 1] normalised_bounds[:, 2]]']

        args_alg = []
        kwargs_general = Dict{Symbol,Any}()
        kwargs_alg = Dict{Symbol,Any}()

        if !isnothing(optimiser_config["max_runs"]) && optimiser_config["max_runs"] < 100
            @warn "'max_runs' of optimiser are smaller than algorithm default for " *
                  "one generation of 100. This may lead to poor results."
            kwargs_alg[:N] = ceil(optimiser_config["max_runs"]/4)
        end

        if optimiser_config["algorithm"] == "MOEAD_DE"
            if optimiser["N_obj"] == 1
                @error "Algorithm 'MOEAD_DE' can only be used with " *
                       "objective_function='multi-objective'"
                throw(InputError())
            end
            push!(args_alg, Metaheuristics.gen_ref_dirs(size(normalised_bounds, 1), population_size))
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

        if haskey(optimiser_config, "optim_kwargs")
            for (keyword, val) in pairs(optimiser_config["optim_kwargs"])
                if Symbol(keyword) in fieldnames(Metaheuristics.Options)
                    kwargs_general[Symbol(keyword)] = val
                elseif Symbol(keyword) in fieldnames(alg)
                    kwargs_alg[Symbol(keyword)] = val
                end
            end
        end
        options = Metaheuristics.Options(; kwargs_general...)
        if optimiser_config["algorithm"] == "CCMO"
            push!(optimiser["args"], alg(Metaheuristics.NSGA2(args_alg...; kwargs_alg...); options=options))
        else
            push!(optimiser["args"], alg(args_alg...; kwargs_alg..., options=options))
        end

    elseif optimiser_config["type"] == "NLopt"
        if occursin(r"LD_.*", optimiser_config["algorithm"])
            @error "The chosen algorithm `$(optimiser_config["algorithm"])` needs a " *
                   "gradient which is not supported in ReSiE"
            throw(InputError())
        end
        if startswith(optimiser_config["algorithm"], "NLOPT_")
            optimiser_config["algorithm"] = split(optimiser_config["algorithm"], "NLOPT_")[2]
        end

        n_dim_opt = length(normalised_bounds[:, 1])
        alg = NLopt.Opt(Symbol(optimiser_config["algorithm"]), n_dim_opt)
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

        if haskey(optimiser_config, "optim_kwargs")
            for (keyword, val) in pairs(optimiser_config["optim_kwargs"])
                if Symbol(keyword) in propertynames(alg)
                    optimiser["kwargs"][Symbol(keyword)] = val
                elseif NLopt.nlopt_has_param(alg, keyword)
                    NLopt.nlopt_set_param(alg, keyword, val)
                end
            end
        end

        for (keyword, val) in pairs(optimiser["kwargs"])
            NLopt.setproperty!(alg, keyword, val)
        end
        optimiser["args"] = [alg, normalised_bounds[:, 3]]

    elseif optimiser_config["type"] == "NOMAD"
        optimiser["kwargs"] = Dict{Symbol,Any}()
        optimiser["kwargs"][:lower_bound] = normalised_bounds[:, 1]
        optimiser["kwargs"][:upper_bound] = normalised_bounds[:, 2]

        if !isnothing(optimiser_config["x_tol_abs"])
            optimiser["kwargs"][:min_mesh_size] = fill(optimiser_config["x_tol_abs"], size(bounds, 1))
        end

        kwargs_general = Dict{Symbol,Any}()
        if !isnothing(optimiser_config["max_runs"])
            kwargs_general[:max_bb_eval] = optimiser_config["max_runs"]
        end
        if !isnothing(optimiser_config["max_time"])
            kwargs_general[:max_time] = optimiser_config["max_time"]
        end

        if haskey(optimiser_config, "optim_kwargs")
            for (keyword, val) in pairs(optimiser_config["optim_kwargs"])
                if Symbol(keyword) in fieldnames(NOMAD.NomadOptions)
                    kwargs_general[Symbol(keyword)] = val
                end
            end
        end

        optimiser["kwargs"][:options] = NOMAD.NomadOptions(; kwargs_general...)

        optimiser["args"] = [size(normalised_bounds, 1), optimiser["N_obj"],
                             fill("OBJ", optimiser["N_obj"]), normalised_bounds[:, 3]]

    else
        #TODO plot outputs like pareto front -> example see optimisation-cli.jl
    end

    return optimiser
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

        missing_keys = setdiff(objective_keys,
                               configured_keys)

        unknown_keys = setdiff(configured_keys,
                               objective_keys)

        if !isempty(missing_keys)
            @error("objective_senses is missing directions for: " *
                   join(missing_keys, ", "),)
            throw(InputError())
        end

        if !isempty(unknown_keys)
            @error("objective_senses contains unknown objective keys: " *
                   join(unknown_keys, ", "),)
            throw(InputError())
        end

        senses = Dict{String,Symbol}()

        for key in objective_keys
            sense = Symbol(lowercase(strip(String(normalized[key]))))

            if !(sense in (:min, :max))
                @error("The objective sense for \"$key\" must be " *
                       "`min` or `max`, got \"$(normalized[key])\".",)
                throw(InputError())
            end

            senses[key] = sense
        end
    end

    signs = Float64[senses[key] == :min ? 1.0 : -1.0
                    for key in objective_keys]

    return senses, signs
end
