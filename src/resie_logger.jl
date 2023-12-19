module Resie_Logger

using Dates
using Logging

export @balanceWarn

"""
CustomLogger

sets up a custom logger to enable logging to file and/or to the console while using CustomLogger
logging formats.
    io                  logging file path
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
    io::Union{IO, Nothing}
    log_to_console::Bool
    log_to_file::Bool
    min_level::Logging.LogLevel
end

"""
BalanceWarnLogLevel

defines a custom warn level for balance errors
    level       defines the Logging.LogLevel (see documentation of Logging-module for details)
    name        name of the custom logger
"""
struct BalanceWarnLogLevel
    level::Int32
    name::String
end
const BalanceWarning = BalanceWarnLogLevel(500, "BalanceWarning")

# define basic functions to handle the custom BalanceWarnLogLevel
Base.isless(a::BalanceWarnLogLevel, b::LogLevel) = isless(a.level, b.level)
Base.isless(a::LogLevel, b::BalanceWarnLogLevel) = isless(a.level, b.level)
Base.convert(::Type{LogLevel}, level::BalanceWarnLogLevel) = LogLevel(level.level)
Base.show(io::IO, level::BalanceWarnLogLevel) =
    if level == BalanceWarning
        print(io, "BalanceWarning")
    else
        show(io, LogLevel(level))
    end

"""
handle CustomLogger

functions to handle CustomLogger.
Current implementation:
- log is written to console with standard ConsoleLogger
- logs of BalanceWarnLogLevel are written using a simplified custom logger
- all logs are written to a file without any formatting
"""
Logging.min_enabled_level(logger::CustomLogger) = logger.min_level
Logging.shouldlog(logger::CustomLogger, level, _module, group, id) = level >= logger.min_level
function Logging.handle_message(logger::CustomLogger, level, message, _module, group, id, file, line; kwargs...)
    if logger.log_to_console
        if level.level == BalanceWarning.level  # directly comparing the levels here as this is easier to implement
            handle_BalanceWarning_message(level, message)
        else
            default_logger = ConsoleLogger(stderr, logger.min_level)
            Logging.handle_message(default_logger, level, message, _module, group, id, file, line; kwargs...)
        end
    end   
    if logger.log_to_file
        if level > Logging.LogLevel(BalanceWarning.level)
            log_message = string("[", level, "] ", message, " (", file, ":", line, ")")
        else # no message for log level below BalanceWarning
            log_message = string("[", level, "] ", message)
        end           
        println(logger.io, log_message)
        if level >= Logging.LogLevel(2000)  # error
            flush(logger.io)   # flush message to write log to file when an error occurs
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
    prefix = "[ " * string(level.name, ':')

    stream = stderr
    buf = IOBuffer()
    iob = IOContext(buf, stream)
    printstyled(iob, prefix, bold=true, color=color)
    print(iob, " ", message, "\n")
    write(stream, take!(buf))
end

# macro as alias for BalanceWarning: @balanceWarn = @logmsg BalanceWarning
macro balanceWarn(exprs...)
    quote
        @logmsg BalanceWarning $(map(x -> esc(x), exprs)...)
    end
end

"""
    start_logger()

starts the logger and opens the file if the log should be written to a file.
    log_to_console::Bool                bool if the log should be printed to the console
    log_to_file::Bool                   bool if the log should be printed to a file
    logfile_path::String                path to a the log file, if log_to_file==false, this can be "" or nothing.
    min_log_level::Logging.LogLevel     the minimal log level for the output. Can be one of Debug, Info, Warn, Error or
                                        a custom LogLevel, like Logging.LogLevel(500) for BalanceWarning
                                        
""" 
function start_logger(log_to_console::Bool, log_to_file::Bool, logfile_path::Union{String, Nothing}, min_log_level::Logging.LogLevel)
    if log_to_file
        log_file = open(logfile_path, "w")
    else
        log_file  = nothing
    end
    logger = CustomLogger(log_file, log_to_console, log_to_file, min_log_level)  
    global_logger(logger)

    if log_to_file
        println(log_file, "ReSiE log file of simulation started at: ", Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))
        println(log_file, "Input file: ", ARGS[1])
        println(log_file, "---------------------------------")
    end
    return log_file
end

"""
    close_logger()

Function to close the logger. 
Prints final statement and closes the logging file.
    log_file::Union{IO, Nothing}        IO hanlder for log file (if present), otherwise this should be nothing
"""
function close_logger(log_file::Union{IO, Nothing})
    if log_file !== nothing
        @info "log saved to ./$(match(r"<file (.*?)>", log_file.name).captures[1])"
        close(log_file)
    end
end

end