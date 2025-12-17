module Resie_Logger

using Dates
using Logging

export @balanceWarn

"""
CustomLogger

sets up a custom logger to enable logging to file and/or to the console while using CustomLogger
logging formats.
    io_general          IO handler with opend general log file (or nothing if log_to_file==false)
    io_balanceWarnings  IO handler with opend balanceWarn log file (or nothing if log_to_file==false)   
    log_to_console      bool if loogig to console shoule be enabled
    log_to_file         bool if logging to file path should be enabled
    min_level           minimum logging level for output, required by logger


Avaiblable logging level:
 name           |  macro name   |  logging level
----------------|---------------|------------------------------------------------------
 Debug          | @debug        | Logging.LogLevel(-1000) 
 Info           | @info         | Logging.LogLevel(    0)
 BalanceWarning | @balanceWarn  | Logging.LogLevel(  500)  <-- this is a custom log level!
 Warn           | @warn         | Logging.LogLevel( 1000)
 Error          | @error        | Logging.LogLevel( 2000)
"""
struct CustomLogger <: Logging.AbstractLogger
    io_general::Union{IO,Nothing}
    io_balanceWarnings::Union{IO,Nothing}
    log_to_console::Bool
    log_to_file::Bool
    min_level::Logging.LogLevel
end

"""
CustomLevel

defines a custom warn level for balance errors
    level       defines the Logging.LogLevel (see documentation of Logging-module for details)
    name        name of the custom logger
"""
struct CustomLevel
    level::Int32
    name::String
end

"""
Adding custom level "BalanceWarning" of level 500.
This can be used with the @balanceWarn macro.

To add further custom levels, you need to define a new constant using the 
CustomLevel struct, the corresponding macro and export the macro. 
Also define the behaviour in Base.show() and Logging.handle_message() below.
"""
const BalanceWarning = CustomLevel(500, "BalanceWarning")
# const BalanceWarning = LogLevel(500)
macro balanceWarn(exprs...)
    quote
        @logmsg BalanceWarning $(map(x -> esc(x), exprs)...)
    end
end

# define basic functions to handle the CustomLevels
Base.isless(a::CustomLevel, b::LogLevel) = isless(a.level, b.level)
Base.isless(a::LogLevel, b::CustomLevel) = isless(a.level, b.level)
Base.convert(::Type{LogLevel}, level::CustomLevel) = LogLevel(level.level)
Base.show(io::IO, level::CustomLevel) =
    if level == BalanceWarning
        print(io, "BalanceWarning")
    else
        show(io, LogLevel(level))
    end

"""
handle CustomLogger and CustomLevel

functions to handle CustomLogger.
Current implementation:
- log is written to console with standard ConsoleLogger
- logs of BalanceWarning are written using a simplified custom logger
- all logs are written to a file without any formatting
"""
Logging.min_enabled_level(logger::CustomLogger) = logger.min_level
Logging.shouldlog(logger::CustomLogger, level, _module, group, id) = level >= logger.min_level
function Logging.handle_message(logger::CustomLogger, level, message, _module, group, id, file, line; kwargs...)
    # log to console
    if logger.log_to_console
        if level.level == BalanceWarning.level  # directly comparing the levels here as this is easier to implement
            handle_BalanceWarning_message(level, message)
        else
            default_logger = ConsoleLogger(stderr, logger.min_level)
            Logging.handle_message(default_logger, level, message, _module, group, id, file, line; kwargs...)
        end
    end

    # log to file
    if logger.log_to_file
        # create log message
        if level > Logging.LogLevel(BalanceWarning.level)
            log_message = string("[", level, "] ", message, " (", file, ":", line, ")")
        else # no julia file path for log level below BalanceWarning
            log_message = string("[", level, "] ", message)
        end

        # write log message to file
        if level.level == BalanceWarning.level
            println(logger.io_balanceWarnings, log_message)
        else
            println(logger.io_general, log_message)
        end

        # flush message to log files if an error occurs
        if level >= Logging.LogLevel(2000)  # error
            flush(logger.io_general)
            flush(logger.io_balanceWarnings)
        end
    end
end

"""
    handle_BalanceWarning_message()

defines the console-output of the BalanceWarn messages.
This is a very simple implementation of the general handle_message()
of the Logging module of Julia with only the basic functionalities.

Color can be one of  :normal, :default, :bold, :black, :blink, :blue, :cyan, 
:green, :hidden, :light_black, :light_blue, :light_cyan, :light_green, 
:light_magenta, :light_red, :light_white, :light_yellow, :magenta, :nothing, 
:red, :reverse, :underline, :white, or :yellow or an integer between 0 and 255 inclusive. 
Note that not all terminals support 256 colors..
"""
function handle_BalanceWarning_message(level, message)
    color = :light_magenta
    # prefix = "[ " * string(level.name, ':')
    prefix = "[ BalanceWarning:"

    stream = stderr
    buf = IOBuffer()
    iob = IOContext(buf, stream)
    printstyled(iob, prefix; bold=true, color=color)
    print(iob, " ", message, "\n")
    write(stream, take!(buf))
end

"""
    start_logger(
        log_to_console, log_to_file, general_logfile_path, balanceWarn_logfile_path,
        min_log_level, input_file
    )

Starts the general and balance warning loggers and opens the files if requested.

Note: This does not set the global logger to the created object. Use `global_logger` for this.

# Arguments
- `log_to_console::Bool`: If the log should be printed to the console
- `log_to_file::Bool`: If the log should be printed to a file
- `general_logfile_path::String`: Path to the general log file, if file logging is requested
- `balanceWarn_logfile_path::String`: Path to the balance warning log file, if file logging
    is requested
- `min_log_level::Logging.LogLevel`: The minimal log level for the output. Can be one of
    Debug, Info, Warn, Error or a custom LogLevel, like Logging.LogLevel(500) for
    BalanceWarning
# Returns
- `CustomLogger`: The created logger
"""
function start_logger(log_to_console::Bool,
                      log_to_file::Bool,
                      general_logfile_path::Union{String,Nothing},
                      balanceWarn_logfile_path::Union{String,Nothing},
                      min_log_level::Logging.LogLevel,
                      input_file::Union{String,Nothing}=nothing)::CustomLogger
    if log_to_file
        log_file_general = open(general_logfile_path, "w")
        log_file_balanceWarn = open(balanceWarn_logfile_path, "w")
    else
        log_file_general = nothing
        log_file_balanceWarn = nothing
    end
    logger = CustomLogger(log_file_general, log_file_balanceWarn, log_to_console, log_to_file, min_log_level)

    if log_to_file
        time_now = Dates.format(now(), "yyyy-mm-dd HH:MM:SS")
        println(log_file_general, "ReSiE general log file of simulation started at: $(time_now)")
        if input_file !== nothing
            println(log_file_general, "Input file: ", abspath(input_file))
        end
        println(log_file_general, "This log file contains general warnings, errors and information written by ReSiE")
        println(log_file_general, "---------------------------------")

        println(log_file_balanceWarn, "ReSiE balanceWarn log file of simulation started at: $(time_now)")
        if input_file !== nothing
            println(log_file_balanceWarn, "Input file: ", abspath(input_file))
        end
        println(log_file_balanceWarn, "This log file contains only balance warnings written by ReSiE.")
        println(log_file_balanceWarn, "---------------------------------")
    end
    return logger
end

"""
    close_logger(logger)

Closes the given logger and prints final statements if file logging is activated.

# Args:
- `logger::CustomLogger`: The logger to close
"""
function close_logger(logger)
    if logger.io_general !== nothing
        @info "general log saved to $(match(r"<file (.*?)>", logger.io_general.name).captures[1])"
        close(logger.io_general)
    end
    if logger.io_balanceWarnings !== nothing
        @info "balanceWarn log saved to $(match(r"<file (.*?)>", logger.io_balanceWarnings.name).captures[1])"
        close(logger.io_balanceWarnings)
    end
end

end
