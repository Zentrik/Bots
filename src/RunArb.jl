using Pkg
Pkg.activate("RunBots")
using Revise, Bots # need to dev Bots, so don't activate Bots env. Also revise needs to be loaded first

using Logging, LoggingExtras, Dates

function not_Bots_message_filter(log)  
    if (@isdefined Bots)
        return log._module === Bots || parentmodule(log._module) == Bots || log._module === Main
    else
        return log._module === Main
    end
end

timestamp_logger(logger) = TransformerLogger(logger) do log
    merge(log, (; message = "$(Dates.format(now(), "yyyy-mm-dd HH:MM:SS.sss")) $(log.message)"))
end

const DIR = "Bots\\src\\logs"

global_logger(TeeLogger(
    EarlyFilteredLogger(not_Bots_message_filter, ConsoleLogger(stderr, Logging.Info)), 
    EarlyFilteredLogger(not_Bots_message_filter, MinLevelLogger(FileLogger("$DIR\\$(Dates.format(now(), "YYYY-mm-dd-HH-MM")).log"), Logging.Info)),
    timestamp_logger(EarlyFilteredLogger(not_Bots_message_filter, MinLevelLogger(FileLogger("$DIR\\Debug-$(Dates.format(now(), "YYYY-mm-dd-HH-MM")).log"), Logging.Debug))),
    timestamp_logger(MinLevelLogger(FileLogger("$DIR\\Verbose-$(Dates.format(now(), "YYYY-mm-dd-HH-MM")).log"), Logging.Debug))
))

ArbBot.retryProd()