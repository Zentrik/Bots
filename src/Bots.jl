module Bots

export ArbBot, Extremities

not_Bots_message_filter(log) = log._module === Bots || parentmodule(log._module) == Bots || log._module === Main

module ArbBot
include("Arb.jl")
end

module Extremities
include("Extremities.jl")
end

end