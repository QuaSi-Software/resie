using Test
using Logging

include("../src/resie_logger.jl")
using .Resie_Logger

# set up logging
global log_to_console = false
global log_to_file = false
global general_logfile_path = nothing
global balanceWarn_logfile_path = nothing
global min_log_level = Resie_Logger.Logging.Debug
global logger = Resie_Logger.start_logger(log_to_console,
                                          log_to_file,
                                          general_logfile_path,
                                          balanceWarn_logfile_path,
                                          min_log_level)
global_logger(logger)

@testset "tests_are_working" begin
    @test true
end

# there is an order to the includes in the sense that later tests might use functionality
# that is being tested in earlier test sets, in particular the project loading methods
include("tests_project_loading.jl")
include("tests_control.jl")
include("tests_energy_systems.jl")

Resie_Logger.close_logger(logger)
