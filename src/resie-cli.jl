using Resie

include("resie_logger.jl")
using .Resie_Logger

function run_cli_loop()
    """
    Main loop for the CLI. This function will run the CLI until the user exits.
    """
    while true
        println("Enter command or 'exit' to quit or 'help' for more info:")
        parts = split(strip(readline()), " ")
        command_input = lowercase(parts[1])

        if command_input == "exit"
            break
        end

        if command_input == "help"
            println("Commands:")
            println("  - 'exit': Exit the CLI.")
            println("  - 'help': Display this help message.")
            println("  - 'run': Run the simulation with arguments:")
            println("    - <file_path>: Project config file")
            println("")
            continue
        end

        if command_input == "run"
            if length(parts) < 2
                println("Insufficient arguments. Type 'help' for more info.")
                println("")
                continue
            end

            println("Running simulation...")
            try
                run(parts[2:end])
            catch
                println("An error occurred while running the simulation:")
                for (exception, backtrace) in current_exceptions()
                    showerror(stdout, exception, backtrace)
                    println(stdout)
                end
                println("")
                continue
            end
            println("Simulation complete.")
            println("")
            continue
        end

        println("Invalid command. Type 'help' for more info.")
        println("")
    end
end

"""
    run(arguments)

Runs the simulation with the given arguments.

The function has a single argument, which is an array of CLI-like arguments. The positional
arguments are as follows:
- `String`: Filepath to the project config file (see documentation on file format). Can be
a path relative to the CWD of the caller.

# Args
- `arguments::Array{String}`: Array of CLI-like arguments. See above for details.
"""
function run(arguments)
    if length(arguments) < 1
        @error "No project config file argument given."
        return
    end

    input_filepath = string(arguments[1])
    if input_filepath === nothing || strip(input_filepath) == ""
        @error "Given project config file argument is empty."
        return
    end

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

    Resie.load_and_run(input_filepath)

    Resie_Logger.close_logger(log_file_general, log_file_balanceWarn)
    return
end

run_cli_loop()
