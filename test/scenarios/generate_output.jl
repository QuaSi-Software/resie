include("../../src/resie_logger.jl")
using .Resie_Logger

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

function main()
    println("No argument given - generating outputs for all scenarios")

    scenarios_dir = joinpath(".", "test", "scenarios")
    for name in readdir(scenarios_dir)
        subdir = joinpath(scenarios_dir, name)
        if !isdir(subdir)
            continue
        end

        println("Generating scenario $name in dir $subdir")
        println("|Logr+|Logr-|")

        log_file_general, log_file_balanceWarn = setup_logger(subdir)
        print("|  x  ")

        Resie_Logger.close_logger(log_file_general, log_file_balanceWarn)
        print("|  x  |")
    end
end

main()
