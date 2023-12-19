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
    log_to_console = true                        # flag if logs should be printed to console
    log_to_file = true                           # flag if logs should be written to a file
    logfile_path = "output/logfile.log"          # log file path
    min_log_level = Resie_Logger.Logging.Info    # the minimal log level for the output. Can be one of Debug, Info, Warn, Error or
                                                 # a custom LogLevel, like Logging.LogLevel(500) for BalanceWarning
    log_file = Resie_Logger.start_logger(log_to_console, log_to_file, logfile_path, min_log_level)

    if length(ARGS) > 0
        filepath = ARGS[1]
        if filepath !== nothing && filepath != ""
            Resie.load_and_run(filepath)
            Resie_Logger.close_logger(log_file)
            return
        end

        @error "Could not find or access project config file at $(ARGS[1])"
        Resie_Logger.close_logger(log_file)
        return
    end

    @error "No project config file argument given."
    Resie_Logger.close_logger(log_file)
end

main()