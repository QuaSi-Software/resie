using Test

include("../src/resie_logger.jl")
using .Resie_Logger

# set up Resie_Logger
log_to_console = false                                       # flag if logs should be printed to console
log_to_file = false                                          # flag if logs should be written to a file
general_logfile_path = "output/logfile_general.log"          # general log file path for debug, info, warn, error
balanceWarn_logfile_path = "output/logfile_balanceWarn.log"  # balanceWarn log file path, only for balanceWarning
min_log_level = Resie_Logger.Logging.Debug                   # the minimal log level for the output. Can be one of Debug, Info, Warn, Error or
                                                             # a custom LogLevel, like Logging.LogLevel(500) for BalanceWarning
log_file_general, log_file_balanceWarn = Resie_Logger.start_logger(log_to_console,
                                                                   log_to_file,
                                                                   general_logfile_path,
                                                                   balanceWarn_logfile_path,
                                                                   min_log_level)
@testset "tests_are_working" begin
    @test true
end

# there is an order to the includes in the sense that later tests might use functionality
# that is being tested in earlier test sets, in particular the project loading methods
include("tests_project_loading.jl")
include("tests_control.jl")
include("tests_energy_systems.jl")

Resie_Logger.close_logger(log_file_general, log_file_balanceWarn)