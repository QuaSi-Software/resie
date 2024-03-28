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
    # set up Resie_Logger
    log_to_console = false                                        # flag if logs should be printed to console
    log_to_file = false                                           # flag if logs should be written to a file
    general_logfile_path = "output/logfile_general.log"          # general log file path for debug, info, warn, error
    balanceWarn_logfile_path = "output/logfile_balanceWarn.log"  # balanceWarn log file path, only for balanceWarning
    min_log_level = Resie_Logger.Logging.Info                    # the minimal log level for the output. Can be one of Debug, Info, Warn, Error or
                                                                 # a custom LogLevel, like Logging.LogLevel(500) for BalanceWarning
    log_file_general, log_file_balanceWarn = Resie_Logger.start_logger(log_to_console, 
                                                                       log_to_file,
                                                                       general_logfile_path,
                                                                       balanceWarn_logfile_path,
                                                                       min_log_level)

    if length(ARGS) > 0
        filepath = ARGS[1]
        if filepath !== nothing && filepath != ""
            Resie.load_and_run(filepath)
            Resie_Logger.close_logger(log_file_general, log_file_balanceWarn)
            return
        end

        @error "Could not find or access project config file at $(ARGS[1])"
        Resie_Logger.close_logger(log_file_general, log_file_balanceWarn)
        return
    end

    @error "No project config file argument given."
    Resie_Logger.close_logger(log_file_general, log_file_balanceWarn)
end

main()