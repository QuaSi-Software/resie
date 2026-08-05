# Suppress GR/GKS GUI windows for file-only plots.
# Windows uses "nul"; other systems use GR's no-output workstation type "100".
if Sys.iswindows()
    ENV["GKSwstype"] = "nul"
    ENV["GKS_WSTYPE"] = "nul"
else
    ENV["GKSwstype"] = "100"
    ENV["GKS_WSTYPE"] = "100"
end

using UUIDs
using Logging
using Resie
const Resie_Logger = Resie.Resie_Logger

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
            println("  - 'run': Run a local simulation with arguments:")
            println("    - <file_path>: Project config file")
            println("    - --exit-after-run: (Optional) Exit the CLI after running the simulation")
            println("  - 'run-confined': Run a confined simulation with arguments:")
            println("    - <file_path>: Project config file")
            println("    - <input_root>: Trusted input directory")
            println("    - <output_root>: Trusted output directory")
            println("    - --exit-after-run: (Optional) Exit the CLI after running the simulation")
            println("")
            continue
        end

        if command_input == "run" || command_input == "run-confined"
            success = false
            exit_after_run = false

            try
                if command_input == "run-confined"
                    success, exit_after_run = run_confined(map(string, parts))
                else
                    success, exit_after_run = run(map(string, parts))
                end
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
    parse_run_arguments(arguments, expected_positional)

Parse common CLI run arguments. Only `--exit-after-run` is accepted so the selected command
fully determines whether execution is local or confined.

# Args
- `arguments::Array{String}`: CLI-like arguments.
- `expected_positional::Int`: Required number of positional arguments.
# Returns
- `Vector{String}`: Positional arguments.
- `Bool`: Whether the CLI should exit after the run.
- `Bool`: Whether parsing was successful.
"""
function parse_run_arguments(arguments::Array{String}, expected_positional::Int)::Tuple{Vector{String},Bool,Bool}
    exit_after_run = false
    positional_arguments = String[]

    for argument in arguments
        if lowercase(argument) == "--exit-after-run"
            exit_after_run = true
        elseif startswith(argument, "--")
            @error "Unknown option: $argument"
            return positional_arguments, exit_after_run, false
        else
            push!(positional_arguments, argument)
        end
    end

    if length(positional_arguments) != expected_positional
        @error "Expected $expected_positional positional argument(s), received $(length(positional_arguments))."
        return positional_arguments, exit_after_run, false
    end

    if any(isempty(strip(argument)) for argument in positional_arguments)
        @error "Run arguments must not be empty."
        return positional_arguments, exit_after_run, false
    end

    return positional_arguments, exit_after_run, true
end

"""
    execute_run(input_filepath, exit_after_run; input_root=nothing, output_root=nothing)

Execute a parsed local or confined CLI run.
"""
function execute_run(input_filepath::String, exit_after_run::Bool;
                     input_root::Union{Nothing,String}=nothing,
                     output_root::Union{Nothing,String}=nothing)::Tuple{Bool,Bool}
    confined = input_root !== nothing || output_root !== nothing
    if confined && (input_root === nothing || output_root === nothing)
        @error "Confined execution requires both input and output roots."
        return false, exit_after_run
    end

    log_to_console = true
    log_to_file = true
    run_ID = uuid4()
    log_directory = confined ? abspath(output_root) : abspath(joinpath(dirname(@__FILE__), "..", "output"))
    mkpath(log_directory)
    general_logfile_path = joinpath(log_directory, "logfile_general_$(run_ID).log")
    balanceWarn_logfile_path = joinpath(log_directory, "logfile_balanceWarn_$(run_ID).log")
    min_log_level = Resie_Logger.Logging.Info
    logger = Resie_Logger.start_logger(log_to_console,
                                       log_to_file,
                                       general_logfile_path,
                                       balanceWarn_logfile_path,
                                       min_log_level,
                                       input_filepath)

    success = false
    with_logger(logger) do
        try
            if confined
                success = Resie.load_and_run_confined(input_filepath, run_ID;
                                                      logger=logger,
                                                      input_root=input_root,
                                                      output_root=output_root)
            else
                success = Resie.load_and_run(input_filepath, run_ID; logger=logger)
            end
        finally
            Resie.close_run(run_ID)
            Resie_Logger.close_logger(logger)
        end
    end

    return success, exit_after_run
end

"""
    run(arguments)

Run a local simulation. This command does not accept confinement options.

The positional argument is:
- `String`: Filepath to the project config file.

Options:
- `--exit-after-run`: Exit the CLI after running the simulation.
"""
function run(arguments::Array{String})::Tuple{Bool,Bool}
    positional_arguments, exit_after_run, success = parse_run_arguments(arguments, 1)
    success || return false, exit_after_run

    return execute_run(positional_arguments[1], exit_after_run)
end

"""
    run_confined(arguments)

Run a simulation through the confined-only server entry point.

The positional arguments are:
- `String`: Filepath to the project config file.
- `String`: Trusted input root.
- `String`: Trusted output root.

Options:
- `--exit-after-run`: Exit the CLI after running the simulation.

This command has no local-mode option. Its roots must be supplied by trusted caller code and
must not be copied from an uploaded project or public request parameter.
"""
function run_confined(arguments::Array{String})::Tuple{Bool,Bool}
    positional_arguments, exit_after_run, success = parse_run_arguments(arguments, 3)
    success || return false, exit_after_run

    return execute_run(positional_arguments[1], exit_after_run;
                       input_root=positional_arguments[2],
                       output_root=positional_arguments[3])
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
