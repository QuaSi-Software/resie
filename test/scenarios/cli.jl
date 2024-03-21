"""
Utility script to generate the output from predefined scenarios for testing purposes.

The script can be used like so:
    `julia --project=. ./test/scenarios/cli.jl generate_output multisector_district`

This generates the output for the scenario `multisector_district`. It is assumed that the
script is called from the directory two levels above this script, which is the main
directory of ReSiE.

Several commands exist for the script, which are:
    * `generate_output`: Generate the output for a scenario by running the simulation.
    * `set_reference`: Set the reference outputs for a scenario. Please note that this will
        not automatically commit the changes to the repository.
    * `compare_ooo`: Compare the order of operations between the reference and output.

Each of the commands can be given, as second argument, the name of a scenario.

If no scenario name is given, performs the command for all scenarios. Be aware that this
might take a while.
"""

include("../../src/resie_logger.jl")
using .Resie_Logger
using Resie

KNOWN_COMMANDS = Set([
    "generate_output",
    "set_reference",
    "compare_ooo"
])

"""
    setup_logger(subdir)

Setup the general and balance warning loggers for the given directory.

The loggers are setup to log into files `general.log` and `balanceWarn.log` within the
directory and to not output to the console.

# Arguments
- `subdir::String`: The directory in which to setup loggers
# Returns
-`IO`: The file stream for the general log
-`IO`: The file stream for the balance warning log
"""
function setup_logger(subdir)
    log_to_console = false
    log_to_file = true
    general_logfile_path = joinpath(subdir, "general.log")
    balanceWarn_logfile_path = joinpath(subdir, "balanceWarn.log")
    min_log_level = Resie_Logger.Logging.Debug
    log_file_general, log_file_balanceWarn = Resie_Logger.start_logger(
        log_to_console,
        log_to_file,
        general_logfile_path,
        balanceWarn_logfile_path,
        min_log_level
    )
    @info "Logging set up"
    return log_file_general, log_file_balanceWarn
end

"""
    compare_files(file_1, file_2)

List of line-by-line differences between the two given files.

If the number of lines differs between the files, adds a difference in line 0 noting the
number of lines for each file. The files are only compared line by line up until the length
of the shorter of the two files.

# Arguments
-`file_1::String`: Filepath to the first file
-`file_2::String`: Filepath to the first file
# Returns
-`List{Tuple{Integer,String}}`: A list of differences with the line number and a
    concatenation of the two lines, each surrounded by ' and seperated with |
"""
function compare_files(file_1, file_2)
    lines_1 = split(read(file_1, String), "\n")
    lines_2 = split(read(file_2, String), "\n")
    differences = []

    if length(lines_1) != length(lines_2)
        push!(
            differences,
            (0, "Different number of lines: $(length(lines_1))|$(length(lines_2))")
        )
    end

    for i in 1:min(length(lines_1), length(lines_2))
        if lines_1[i] != lines_2[i]
            push!(differences, (i, "'" * lines_1[i] * "' | '" * lines_2[i] * "'"))
        end
    end

    return differences
end

"""
    generate_output(name, subdir)

Generate the output for the given scenario.

# Arguments
-`name::String`: The name of the scenario
-`subdir::String`: Full path of the subdir for the scenario
"""
function generate_output(name::String, subdir::String)
    println("Generating scenario $name in dir $subdir")
    println("|Logr+|RdJsn|SmPrm|LdCmp|OrdOp|RnSim|Logr-|")

    log_file_general, log_file_balanceWarn = setup_logger(subdir)
    print("|  ✓  ")

    project_config = nothing
    try
        project_config = Resie.read_JSON(joinpath(subdir, "inputs.json"))
        print("|  ✓  ")
    catch error
        @error error
        print("|  X  ")
    end

    sim_params = nothing
    try
        if project_config !== nothing
            sim_params = Resie.get_simulation_params(project_config)
            print("|  ✓  ")
        else
            print("|     ")
        end
    catch error
        @error error
        print("|  X  ")
    end

    components = nothing
    try
        if sim_params !== nothing && project_config !== nothing
            components = Resie.load_components(
                project_config["components"],
                sim_params
            )
            print("|  ✓  ")
        else
            print("|     ")
        end
    catch error
        @error error
        print("|  X  ")
    end

    step_order = nothing
    try
        if project_config !== nothing && components !== nothing
            if (
                haskey(project_config, "order_of_operation")
                && length(project_config["order_of_operation"]) > 0
            )
                step_order = Resie.load_order_of_operations(
                    project_config["order_of_operation"],
                    components
                )
            else
                step_order = Resie.calculate_order_of_operations(components)
            end
            print("|  ✓  ")
        else
            print("|     ")
        end
    catch error
        @error error
        print("|  X  ")
    end

    if (
        project_config !== nothing
        && sim_params !== nothing
        && components !== nothing
        && step_order !== nothing
    )
        Resie.run_simulation_loop(project_config, sim_params, components, step_order)
        print("|  ✓  ")
    else
        print("|     ")
    end

    Resie_Logger.close_logger(log_file_general, log_file_balanceWarn)
    print("|  ✓  ")

    println("|")
end

"""
    set_reference(name, subdir)

Set the reference outputs by copying relevant output files.

# Arguments
-`name::String`: The name of the scenario
-`subdir::String`: Full path of the subdir for the scenario
"""
function set_reference(name, subdir)
    outputs_to_set = [
        "auxiliary_info.md",
        "balanceWarn.log",
        "out.csv",
        "output_plot.html",
        "sankey_plot.html"
    ]

    println("Setting reference outputs for scenario $name")

    for file_name in outputs_to_set
        output_file_path = joinpath(subdir, file_name)
        ref_file_path = joinpath(subdir, "ref_" * file_name)

        if isfile(output_file_path)
            cp(output_file_path, ref_file_path, force=true)
        end
    end
end

"""
    compare_ooo(name, subdir)

Compare the order of operations between the reference and output files.

# Arguments
-`name::String`: The name of the scenario
-`subdir::String`: Full path of the subdir for the scenario
"""
function compare_ooo(name, subdir)
    print("Comparing order of operations for scenario $name: ")

    differences = compare_files(
        joinpath(subdir, "auxiliary_info.md"),
        joinpath(subdir, "ref_auxiliary_info.md")
    )

    if length(differences) == 0
        print("✓\n")
    else
        print("X\n")
        for diff_line in differences
            println("Line $(diff_line[1]): $(diff_line[2])")
        end
    end
end

"""
    main()

Entry point for the script.

Parses arguments and performs the commands as described in the documentation above.
"""
function main()
    if length(ARGS) < 1
        println("No command given.")
        return
    end

    command = lowercase(ARGS[1])
    if !(command in KNOWN_COMMANDS)
        println("Unknown command.")
        return
    end

    if length(ARGS) >= 2 && ARGS[2] != ""
        scenario_name = ARGS[2]
        println("Given scenario: $scenario_name")
    else
        scenario_name = nothing
        println("No scenario given - performing command for all scenarios")
    end

    scenarios_dir = joinpath(".", "test", "scenarios")
    for name in readdir(scenarios_dir)
        if scenario_name !== nothing && scenario_name != name
            continue
        end

        subdir = joinpath(scenarios_dir, name)
        if !isdir(subdir)
            continue
        end

        if command == "generate_output"
            generate_output(name, subdir)
        elseif command == "set_reference"
            set_reference(name, subdir)
        elseif command == "compare_ooo"
            compare_ooo(name, subdir)
        end
    end
end

main()
