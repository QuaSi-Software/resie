using UUIDs
using Resie

include("resie_logger.jl")
using .Resie_Logger

"""
    parse_arguments(input)

Parses a string into a vector of strings.

The input string is split by spaces, but strings enclosed in double quotes are treated as a
single argument. For example, the string `run "some file"` is parsed as
`["run", "some file"]`. Note that the double quotes are removed from the resulting strings.

# Args
- `input::String`: The string to parse.
# Returns
- `Vector{String}`: The parsed arguments.
"""
function parse_arguments(input::String)::Vector{String}
    arguments = []
    for match in eachmatch(r"\"[^\"]*\"|[^\s]+", input)
        push!(arguments, string(strip(match.match, '"')))
    end
    return arguments
end

"""
Main loop for the CLI. This function will run the CLI until the user exits.

If there are additional arguments passed to the script, they will be used as the first
command to run with the arguments following that as arguments to the command. If no
additional arguments are passed, the user will be prompted to enter a command.
"""
function run_cli_loop()
    # setting this to false causes Ctrl+c to throw InterruptException instead of quitting
    # Julia entirely
    Base.exit_on_sigint(false)

    is_first = true
    while true
        parts = []
        if is_first
            is_first = false
            if length(ARGS) > 0
                parts = ARGS
            end
        end

        if length(parts) == 0
            println("Enter command or 'exit' to quit or 'help' for more info:")
            input = string(strip(readline()))
            parts = parse_arguments(input)
            if length(parts) == 0
                continue
            end
        end

        command_input = lowercase(parts[1])
        deleteat!(parts, 1)

        if command_input == "exit"
            break
        end

        if command_input == "help"
            println("Commands:")
            println("  - 'exit': Exit the CLI.")
            println("  - 'help': Display this help message.")
            println("  - 'run': Run the simulation with arguments:")
            println("    - <file_path>: Project config file")
            println("    - --exit-after-run: (Optional) Exit the CLI after running the simulation")
            println("")
            continue
        end

        if command_input == "run"
            success = false
            exit_after_run = false

            try
                success, exit_after_run = run(map(string, parts))
            catch exc
                if exc isa InterruptException
                    println("Simulation was interrupted.")
                    success = true # not quite correct, but suppresses the next message
                else
                    println("An error occurred while running the simulation:")
                    for (exception, backtrace) in current_exceptions()
                        showerror(stdout, exception, backtrace)
                        println(stdout)
                    end
                end
            end

            if !success
                println("Simulation was not successful. Check logs for details.")
            end
            println("")

            if exit_after_run
                break
            end
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

Keyword arguments:
- `--exit-after-run`: If this argument is present, the CLI will exit after running the
simulation.

# Args
- `arguments::Array{String}`: Array of CLI-like arguments. See above for details.
# Returns
- `Bool`: `true` if the simulation was successful, `false` otherwise.
- `Bool`: `true` if the CLI should be quit after running the simulation, `false` otherwise.
"""
function run(arguments::Array{String})::Tuple{Bool,Bool}
    exit_after_run = false
    for (idx, arg) in enumerate(arguments)
        if lowercase(arg) == "--exit-after-run"
            exit_after_run = true
            deleteat!(arguments, idx)
        end
    end

    if length(arguments) < 1
        @error "No project config file argument given."
        return false, exit_after_run
    end

    input_filepath = string(arguments[1])
    if input_filepath === nothing || strip(input_filepath) == ""
        @error "Given project config file argument is empty."
        return false, exit_after_run
    end

    log_to_console = true
    log_to_file = true
    general_logfile_path = abspath(joinpath(dirname(@__FILE__), "..", "output", "logfile_general.log"))
    balanceWarn_logfile_path = abspath(joinpath(dirname(@__FILE__), "..", "output", "logfile_balanceWarn.log"))
    min_log_level = Resie_Logger.Logging.Info
    logger = Resie_Logger.start_logger(log_to_console,
                                       log_to_file,
                                       general_logfile_path,
                                       balanceWarn_logfile_path,
                                       min_log_level,
                                       input_filepath)
    global_logger(logger)

    success = false
    run_ID = uuid1()
    try
        success = Resie.load_and_run(input_filepath, run_ID)
    catch exc
        # exceptions should just be rethrown, but the finally block should also make sure
        # that the logger is closed and the run removed from the run registry.
        # it's not fully clear if this actually works since rethrow might quit out of the
        # execution entirely
        rethrow(exc)
    finally
        Resie.close_run(run_ID)
        Resie_Logger.close_logger(logger)
    end

    return success, exit_after_run
end

try
    run_cli_loop()
catch exc
    if exc isa InterruptException
        # nothing to do, just exit
    else
        # exceptions not caught inside the loop should be carried to the shell
        rethrow(exc)
    end
end
