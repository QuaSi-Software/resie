using Infiltrator
using Revise
using Resie
Revise.track(Resie)
include("resie_logger.jl")
using .Resie_Logger

function run_repl(input_filepath::String)
    revise()
    log_to_console = true
    log_to_file = true
    general_logfile_path = "output/logfile_general.log"
    balanceWarn_logfile_path = "output/logfile_balanceWarn.log"
    min_log_level = Resie_Logger.Logging.Info
    log_file_general, log_file_balanceWarn = Resie_Logger.start_logger(log_to_console,
                                                                       log_to_file,
                                                                       general_logfile_path,
                                                                       balanceWarn_logfile_path,
                                                                       min_log_level,
                                                                       input_filepath)

    success = Resie.load_and_run(input_filepath)

    Resie_Logger.close_logger(log_file_general, log_file_balanceWarn)
end