module Bots

export ArbBot, Extremities

# https://github.com/julia-vscode/julia-vscode/blob/ed2fb65b42d11b6beec51aaaa4703824cfa0fd1c/scripts/packages/VSCodeServer/src/repl.jl#L233
function display_error(err, bt)
    io = IOBuffer()
    st = stacktrace(crop_backtrace(bt))
    showerror(IOContext(io, :limit => true), err, st)
    @error String(take!(io))
end

# https://github.com/julia-vscode/julia-vscode/blob/ed2fb65b42d11b6beec51aaaa4703824cfa0fd1c/scripts/packages/VSCodeServer/src/misc.jl#L6
function crop_backtrace(bt)
    i = find_first_topelevel_scope(bt)
    return bt[1:(i === nothing ? end : i)]
end

function find_first_topelevel_scope(bt::Vector{<:Union{Base.InterpreterIP,Ptr{Cvoid}}})
    for (i, ip) in enumerate(bt)
        st = Base.StackTraces.lookup(ip)
        ind = findfirst(st) do frame
            linfo = frame.linfo
            if linfo isa Core.CodeInfo
                linetable = linfo.linetable
                if isa(linetable, Vector) && length(linetable) ≥ 1
                    lin = first(linetable)
                    if isa(lin, Core.LineInfoNode) && lin.method === Symbol("top-level scope")
                        return true
                    end
                end
            else
                return frame.func === Symbol("top-level scope")
            end
        end
        ind === nothing || return i
    end
    return
end

macro spawn_showerr(ex)
    esc(quote
        Threads.@spawn try
            eval($ex)
        catch err
            display_error(err, catch_backtrace())
            rethrow()
        end
    end)
end

macro async_showerr(ex)
    esc(quote
        @async try
            eval($ex)
        catch err
            display_error(err, catch_backtrace())
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
            display_error(err, catch_backtrace())
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