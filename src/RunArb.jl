using Pkg
Pkg.activate("Bots")
using Bots

using Logging, LoggingExtras, Dates
timestamp_logger(logger) = TransformerLogger(logger) do log
    merge(log, (; message = "$(Dates.format(now(), "yyyy-mm-dd HH:MM:SS.sss")) $(log.message)"))
end
const DIR = "Bots\\src\\logs"
global_logger(TeeLogger(
    EarlyFilteredLogger(Bots.not_Bots_message_filter, ConsoleLogger(stderr, Logging.Info)), 
    EarlyFilteredLogger(Bots.not_Bots_message_filter, MinLevelLogger(FileLogger("$DIR\\$(Dates.format(now(), "YYYY-mm-dd")).log"), Logging.Info)),
    EarlyFilteredLogger(Bots.not_Bots_message_filter, MinLevelLogger(FileLogger("$DIR\\Debug-$(Dates.format(now(), "YYYY-mm-dd")).log"), Logging.Debug)),
    timestamp_logger(MinLevelLogger(FileLogger("$DIR\\Verbose-$(Dates.format(now(), "YYYY-mm-dd")).log"), Logging.Debug))
))

ArbBot.production(skip=true)