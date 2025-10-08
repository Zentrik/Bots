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
            rethrow(err)
        end
    end)
end

try
    @sync begin
        @async while true
            sleep(5)
            println("going")
            throw(MethodError)
        end
        throw(ArgumentError)
    end
catch
    println("Caught 1")
end

try
    @sync begin
        @async throw(ArgumentError)
    end
catch
    println("Caught 1")
end

try
    @sync begin
        try 
            @async throw(ArgumentError)
        catch
            rethrow()
        end
        @async while true
            sleep(5)
            # println("oops")
        end
    end
catch
    println("Caught 1")
end

try
    @sync begin
        @async_showerr throw(ArgumentError)
        @async while true
            sleep(5)
            # println("oops")
        end
    end
catch
    println("Caught 1")
    rethrow()
end

try 
    try
        @sync begin
            @async throw(ArgumentError)
            @async for i in 1:2
                sleep(1)
                println("oops")
            end
        end
    catch
        println("Caught 1")
        rethrow()
    end
catch err
    println("Caught")
    # rethrow()
finally
    println("done")
end

function wait_until(c, timeout::Real) # `c` is any object that is both wait-able and cancel-able (e.g. any IO or a Channel, etc.)
    timer = Timer(timeout) do t
        isready(c) || close(c)
    end
    try
        return wait(c)
    finally
        close(timer)
    end
end

function wait_until(c::Condition, timeout::Real) # `c` is any object that is both wait-able and cancel-able (e.g. any IO or a Channel, etc.)
    timer = Timer(timeout) do _
        # notify(c)
        notify(c, "Wait timed out", error=true)
    end
    try
        return wait(c)
    finally
        close(timer)
    end
end

function testCatch()
    try
        @sync begin
            # @async throw(ArgumentError)
            @async begin 
                Bots.wait_until(Condition(), 2)
            end
        end
    catch err
        # for e in err
        #     println(e)
        # end
        # println(err)
        return err
    end
end

function test()
    rerun = :A
    try
        @sync for i in 1:2
            x = Condition()
            # @async throw(ArgumentError)
            Bots.@async_showerr println("1 $i")
            Bots.@async_showerr begin 
                # @debug Bots.wait_until(x, 2)
                Bots.wait_until(x, 5)
            end
            Bots.@async_showerr println("2 $i")
        end
    catch err
        println("Caught 1")
        return :Fail
    end

    return rerun
end

function test()
    rerun = :A
    try
        @sync for i in 1:2
            x = Condition()
            # @async throw(ArgumentError)
            @async println("1 $i")
            @async begin 
                # @debug Bots.wait_until(x, 2)
                Bots.wait_until(x, 5)
            end
            @async println("2 $i")
        end
    catch err
        println("Caught 1")
        return :Fail
    end

    return rerun
end

try 
    try
        @sync begin
            x = Condition()
            # @async throw(ArgumentError)
            @async_showerr begin 
                wait_until(x, 2)
            end
        end
    catch
        println("Caught 1")
        rethrow()
    end
catch err
    println("Caught")
    rethrow()
finally
    println("done")
end

try 
    try
        Base.Experimental.@sync begin
            x = Condition()
            # @async throw(ArgumentError)
            @async begin 
                println(wait(x))
                println(x)
            end
            @async begin
                println("Sleeping")
                sleep(4)
                println("Awoken")
                notify(x)
            end
        end
    catch
        println("Caught 1")
        rethrow()
    end
catch err
    println("Caught")
    # rethrow()
finally
    println("done")
end

try 
    @sync let timer
        @async begin
            global timer = Timer(2) do _
                println("oops")
            end
        end

        try
            @sync begin
                @async throw(ArgumentError)
            end
        catch err
            println("Caught 1")
            close(timer)
            rethrow()
        end
    end
catch err
    println("Caught")
finally
    println("done")
end

try 
    @sync begin
        @async for i in 1:2
            sleep(2)
            println("oops")
        end
        throw(ArgumentError)
    end
catch err
    println("Caught")
finally
    println("done")
end


@sync begin 
    currentTask = @async sleep(2)

    for i in 1:2
        @async begin
            if !istaskdone(currentTask)
                println("1: $currentTask, $i")
                wait(currentTask)
                println("2: $currentTask, $i")
            end
            println("3: $currentTask, $i")
            currentTask = current_task()
            println("4: $currentTask, $i")
            # yield()
            println("5: $currentTask, $i")
            println(i)
            println("6: $currentTask, $i")
            sleep(1)
            println("7: $currentTask, $i")
            println(i)
        end
    end
end

@sync begin 
    currentTask = @async sleep(2)

    for i in 1:2
        @async begin
            while !istaskdone(currentTask)
                println("1: $currentTask, $i")
                wait(currentTask)
                println("2: $currentTask, $i")
            end
            println("3: $currentTask, $i")
            currentTask = current_task()
            println("4: $currentTask, $i")
            # yield()
            println("5: $currentTask, $i")
            println(i)
            println("6: $currentTask, $i")
            sleep(1)
            println("7: $currentTask, $i")
            println(i)
        end
    end
end

@sync try 
    # need to add check to see if we're actually receiving messages
    @sync for i in 1:2
        Bots.@async_showerr throw(ArgumentError(i))
    end
catch
    return :PostFailure
end

try
    throw(ArgumentError)
finally
    println(first(current_exceptions()).exception)
end

function testRethrow()
    try
        try
            throw(ArgumentError)
        finally
            println("Finally")
            println(err)
        end
    catch
        println("Caught")
    end
end

using HTTP.WebSockets
import ..Bots: display_error
function testSocket()
    error = false

    WebSockets.open("wss://demo.piesocket.com/v3/channel_123?api_key=VCXCEuvhGcBDP7XhiJJUDvR1e1D3eiVjgZ9VRiaV&notify_self") do socket
        try
            throw(ArgumentError)
        finally
            error = current_exceptions()
            println("Finally")
            if !WebSockets.isclosed(socket)
                println("Left Channel")
    
                close(socket)
            end
        end
    end
    
    return error
end

function testSocket()
    WebSockets.open("wss://demo.piesocket.com/v3/channel_123?api_key=VCXCEuvhGcBDP7XhiJJUDvR1e1D3eiVjgZ9VRiaV&notify_self", suppress_close_error=true) do socket
        try
            throw(MethodError)
        finally
            println("Finally")
            close(socket)
        end
    end
end

import ..Bots: display_error
function testSocketPrint()
    try
        testSocket()
    catch
        # println(typeof(err))
        # println(current_exceptions())
        # display_error(err, catch_backtrace())
        # return err
        for error in reverse(current_exceptions())
            display_error(error.exception, error.backtrace)
        end
    end
end

function shouldBreak(exception)
    if exception isa MethodError || exception isa InterruptException || exception isa DimensionMismatch
        return true
    elseif exception isa TaskFailedException
        return shouldBreak(exception.task.exception)
    elseif exception isa CompositeException
        return any(shouldBreak.(exception.exceptions))
    else
        return false
    end
end

function testBreak()
    try
        try
            # throw(MethodError("", ""))
            sleep("")
        finally
            println("Finally")
        end
    catch err
        println(err)
        for (exception, _) in current_exceptions()
            if shouldBreak(exception)
                println("Broke")
                return nothing
            end
        end
    end
end

function testAsyncBreak()
    try
        @sync try
            @async throw(MethodError("", ""))
        finally
            println("Finally")
        end
    catch err
        println(err)

        for (exception, _) in current_exceptions()
            if shouldBreak(exception)
                println("Broke")
                return nothing
            end
        end
    end
end

function testOrdering()
    task = @task println(1)
    task2 = @task begin
        sleep(5)
        println(2)
    end

    schedule(task2)
    sleep(0.001)
    schedule(task)

    wait(task)
    wait(task2)
end