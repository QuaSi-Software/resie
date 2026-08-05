using Test
using Logging
using Resie
using Resie.Resie_Logger

# set up logging
log_to_console = false
log_to_file = false
general_logfile_path = nothing
balanceWarn_logfile_path = nothing
min_log_level = Logging.Debug
logger = Resie_Logger.start_logger(log_to_console,
                                   log_to_file,
                                   general_logfile_path,
                                   balanceWarn_logfile_path,
                                   min_log_level)

try
    with_logger(logger) do
        @testset "tests_are_working" begin
            @test true
        end

        # there is an order to the includes in the sense that later tests might use functionality
        # that is being tested in earlier test sets, in particular the project loading methods
        include("tests_project_loading.jl")
        include("tests_control.jl")
        include("tests_energy_systems.jl")
    end
finally
    Resie_Logger.close_logger(logger)
end
