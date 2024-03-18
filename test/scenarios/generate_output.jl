"""
Utility script generate the output from predefined scenarios for testing purposes. The
script can be used like so:
    `julia --project=. ./test/scenarios/generate_output.jl "multisector_district"`

This generates the output for the scenario `multisector_district`. It is assumed that the
script is called from the directory two levels above this script, which is the main
directory of ReSiE.

If no argument is given, generates the output for all scenarios. Be aware that this might
take a while.
"""

include("../../src/resie_logger.jl")
using .Resie_Logger

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
    main()

Entry point for the script.

Parses arguments and performs the functions as described in the documentation above.
"""
function main()
    if length(ARGS) > 0 && ARGS[1] != ""
        scenario_name = strip(ARGS[1])
        println("Given scenario to generate: $scenario_name")
    else
        scenario_name = nothing
        println("No argument given - generating outputs for all scenarios")
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

        println("Generating scenario $name in dir $subdir")
        println("|Logr+|Logr-|")

        log_file_general, log_file_balanceWarn = setup_logger(subdir)
        print("|  x  ")

        Resie_Logger.close_logger(log_file_general, log_file_balanceWarn)
        print("|  x  ")

        println("|")
    end
end

main()
