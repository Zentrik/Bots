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

try 
    try
        @sync begin
            x = Condition()
            # @async throw(ArgumentError)
            @async begin 
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