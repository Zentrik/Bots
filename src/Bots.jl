module Bots

export ArbBot, Extremities

macro spawn_showerr(ex)
    esc(quote
        Threads.@spawn try
            eval($ex)
        catch err
            @error "Something went wrong" err
            # https://github.com/JuliaLang/julia/pull/48282#issuecomment-1426083522
            for line in ["[$ii] $frame" for (ii, frame) in enumerate(stacktrace(catch_backtrace()))]
                @error line
            end
            rethrow(err)
        end
    end)
end

macro async_showerr(ex)
    esc(quote
        @async try
            eval($ex)
        catch err
            @error "Something went wrong" err
            # https://github.com/JuliaLang/julia/pull/48282#issuecomment-1426083522
            for line in ["[$ii] $frame" for (ii, frame) in enumerate(stacktrace(catch_backtrace()))]
                @error line
            end
            rethrow()
        end
    end)
end

using SmartAsserts

macro smart_assert_showerr(ex, msg=nothing)
    esc(quote
        try
            @smart_assert eval($ex) eval($msg)
        catch err
            @error "Something went wrong" err
            # https://github.com/JuliaLang/julia/pull/48282#issuecomment-1426083522
            for line in ["[$ii] $frame" for (ii, frame) in enumerate(stacktrace(catch_backtrace()))]
                @error line
            end
            rethrow()
        end
    end)
end

function testLogging()
    @info("You won't see this")
    @warn("won't see this either")
    @error("You will only see this")
    @debug "test"
    @smart_assert_showerr 1 != 1
end

module ArbBot
include("Arb.jl")
end

module Extremities
include("Extremities.jl")
end

end