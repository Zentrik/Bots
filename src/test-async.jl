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

@sync try
    try 
        @async throw(ArgumentError)
    catch
        rethrow()
    end
    @async while true
        sleep(5)
        # println("oops")
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
            # @async while true
            #     sleep(5)
            #     # println("oops")
            # end
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