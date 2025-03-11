using Resie

include("resie_logger.jl")
using .Resie_Logger

"""
    main()

Entry point of the CLI.

# Command line arguments
## Positional arguments
- `String`: Filepath to the project config file (see documentation on file format). Can be
a path relative to the CWD of the caller.
"""
function main()
    if length(ARGS) < 1
        @error "No project config file argument given."
        return
    end

    input_filepath = ARGS[1]
    if input_filepath === nothing || strip(input_filepath) == ""
        @error "Could not find or access project config file at $input_filepath"
        return
    end

    # set up logging
    log_to_console = true
    log_to_file = true
    general_logfile_path = "output/logfile_general.log"
    balanceWarn_logfile_path = "output/logfile_balanceWarn.log"
    min_log_level = Resie_Logger.Logging.Error
    log_file_general, log_file_balanceWarn = Resie_Logger.start_logger(log_to_console,
                                                                       log_to_file,
                                                                       general_logfile_path,
                                                                       balanceWarn_logfile_path,
                                                                       min_log_level,
                                                                       input_filepath)

    # run the simulation
    Resie.load_and_run(input_filepath)

    # close logging
    Resie_Logger.close_logger(log_file_general, log_file_balanceWarn)
    return
end

main()
