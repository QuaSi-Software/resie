using Test

include("../src/resie_logger.jl")
using .Resie_Logger

# set up logging
log_to_console = false
log_to_file = false
general_logfile_path = "output/logfile_general.log"
balanceWarn_logfile_path = "output/logfile_balanceWarn.log"
min_log_level = Resie_Logger.Logging.Debug
log_file_general, log_file_balanceWarn = Resie_Logger.start_logger(
    log_to_console,
    log_to_file,
    general_logfile_path,
    balanceWarn_logfile_path,
    min_log_level
)

@testset "tests_are_working" begin
    @test true
end

# there is an order to the includes in the sense that later tests might use functionality
# that is being tested in earlier test sets, in particular the project loading methods
include("tests_project_loading.jl")
include("tests_control.jl")
include("tests_energy_systems.jl")

Resie_Logger.close_logger(log_file_general, log_file_balanceWarn)