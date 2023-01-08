using ManifoldMarkets, TOML, Optimization, OptimizationBBO, Dates, Combinatorics, ThreadsX
using Suppressor

struct Group
    name::String
    slugs::Vector{String}
    y_matrix::BitMatrix
    n_matrix::BitMatrix
    noMarkets::Int64

    function Group(name, GroupDict)
        name = name
        slugs = urlToSlug.(collect(keys(GroupDict)))
        actionVectors = reduce(hcat, values(GroupDict))
        y_matrix = actionVectors .== "YES"
        n_matrix = actionVectors .== "NO"
        noMarkets = length(slugs)
    
        return new(name, slugs, y_matrix, n_matrix, noMarkets)
    end
end

struct PlannedBet
    amount::Float64
    shares::Float64
    outcome::String
    market::Market
end

function Base.show(io::IO, plannedBet::PlannedBet)
    printstyled(io, "$(plannedBet.market.question)\n", color=:green)
    print(io, "Buy $(plannedBet.shares) $(plannedBet.outcome) shares for $(plannedBet.amount)")
end

function execute(bet, APIKEY)
    response = createBet(APIKEY, bet.market.id, bet.amount, bet.outcome)
    # need to check if returned info matches what we wanted to bet, i.e. if we got less shares than we wanted to. If we got more ig either moved or smth weird with limit orders.
    if response.shares ≉ bet.shares
        println(bet.market.question)
        println(bet.market.url)
        println(response)
        println(response.fills)
    end

    return response
end

function updateShares!(noSharesBySlug, yesSharesBySlug, slug, newBet)
    if newBet.outcome == "NO"
        noSharesBySlug[slug] += newBet.shares
    elseif newBet.outcome == "YES"
        yesSharesBySlug[slug] += newBet.shares
    end
end

function redeemShares!(noSharesBySlug, yesSharesBySlug)
    for slug in keys(noSharesBySlug)
        if yesSharesBySlug[slug] >= noSharesBySlug[slug]
            yesSharesBySlug[slug] -= noSharesBySlug[slug]
            noSharesBySlug[slug] = 0.
        elseif yesSharesBySlug[slug] < noSharesBySlug[slug]
            noSharesBySlug[slug] -= yesSharesBySlug[slug]
            yesSharesBySlug[slug] = 0.
        end
    end

    # @assert mapreduce(slug -> min(noSharesBySlug[slug], yesSharesBySlug[slug]), max, getSlugs(GROUPS)) ≈ 0
end

function f(betAmount, group, markets, limitOrdersBySlug, noSharesBySlug, yesSharesBySlug, bettableSlugsIndex)
    newProb = [markets[slug].probability for slug in group.slugs]
    noShares = zeros(group.noMarkets)
    yesShares = zeros(group.noMarkets)

    fees = 0.

    for (i, j) in enumerate(bettableSlugsIndex)
        slug = group.slugs[j]
        market = markets[slug]

        if abs(betAmount[i]) >= 1 
            shares = betToShares(market, limitOrdersBySlug[slug], betAmount[i]).shares
            newProb[j] = betToShares(market, limitOrdersBySlug[slug], betAmount[i]).probability

            fees += 0.1
            if betAmount[i] >= 1
                yesShares[j] += shares
            elseif betAmount[i] <= -1
                noShares[j] += shares
            end
        else
            betAmount[i] = 0
        end
    end

    profitsByEvent = group.y_matrix * (yesShares + getindex.(Ref(yesSharesBySlug), group.slugs)) + group.n_matrix * (noShares + getindex.(Ref(noSharesBySlug), group.slugs)) .- sum(abs.(betAmount)) .- fees # we shouldn't count resolved markets.

    return (profitsByEvent=profitsByEvent, noShares=noShares, yesShares=yesShares, newProbability=newProb)
end

function optimise(group, markets, limitOrdersBySlug, maxBetAmount, noSharesBySlug, yesSharesBySlug, bettableSlugsIndex)    
    profitF = OptimizationFunction((betAmount, _) -> -minimum( f(betAmount, group, markets, limitOrdersBySlug, noSharesBySlug, yesSharesBySlug, bettableSlugsIndex).profitsByEvent ))

    x0 = repeat([0.], length(bettableSlugsIndex))
    lb = repeat([-maxBetAmount], length(bettableSlugsIndex))
    ub = repeat([maxBetAmount], length(bettableSlugsIndex))

    problem = Optimization.OptimizationProblem(profitF, x0, lb=lb, ub=ub)

    sol = solve(problem, BBO_adaptive_de_rand_1_bin_radiuslimited(), maxtime=10.)

    bestSolution = repeat([0.], length(bettableSlugsIndex))
    maxRiskFreeProfit = f(bestSolution, group, markets, limitOrdersBySlug, noSharesBySlug, yesSharesBySlug, bettableSlugsIndex).profitsByEvent |> minimum

    nonZeroIndices = findall(!iszero, sol.u)

    for indices in powerset(nonZeroIndices)
        betAmount = copy(sol.u)
        betAmount[indices] .= 0

        riskFreeProfit = f(betAmount, group, markets, limitOrdersBySlug, noSharesBySlug, yesSharesBySlug, bettableSlugsIndex).profitsByEvent |> minimum

        if riskFreeProfit > maxRiskFreeProfit
            maxRiskFreeProfit = riskFreeProfit
            bestSolution = betAmount
        end
    end

    return bestSolution
end

function getMarkets(slugs)
    markets = Dict{String, Market}()
    @sync for slug in slugs
        @async try
            markets[slug] = getMarketBySlug(slug)
        catch err
            bt = catch_backtrace()
            println()
            showerror(stderr, err, bt)
        end
    end

    return markets
end

function getMarketsAndBets!(oldUserBalance, GROUPS, USERNAME)
    markets = Dict{String, Market}()
    betsBySlug = Dict{String, Vector{Bet}}()
    # betsByMe = Dict{String, Vector{Bet}}()
    botBalance = 0.

    @sync begin
        # println(Dates.format(now(), "HH:MM:SS.sss"))
        for group in values(GROUPS), url in keys(group)
            slug = urlToSlug(url)

            @async try
                markets[slug] = getMarketBySlug(slug)
            catch err
                bt = catch_backtrace()
                println()
                showerror(stderr, err, bt)
            end

            @async try
                betsBySlug[slug] = getBets(slug=slug, limit=200)
            catch err
                bt = catch_backtrace()
                println()
                showerror(stderr, err, bt)
            end
            # betsBySlug[slug] = []
        end

        # println(Dates.format(now(), "HH:MM:SS.sss"))
        for userId in keys(oldUserBalance)  # we need to drop users after not needed for a while
            @async try
                oldUserBalance[userId] = getUserById(userId).balance
            catch err
                bt = catch_backtrace()
                println()
                showerror(stderr, err, bt)
            end
        end 
        # println(Dates.format(now(), "HH:MM:SS.sss"))

        @async try
            botBalance = getUserByUsername(USERNAME).balance
        catch err
            bt = catch_backtrace()
            println()
            showerror(stderr, err, bt)
        end
    end
    # println(Dates.format(now(), "HH:MM:SS.sss"))

    return (markets=markets, betsBySlug=betsBySlug, botBalance=botBalance)
end

function getSlugs(GROUPS::Dict)
    return mapreduce(x -> urlToSlug.(x), vcat, keys.(values(GROUPS)))
    # for group in values(GROUPS), url in keys(group)
    #     slug = urlToSlug(url)
    # end
end

function getSlugs(groups::Vector{Group})
    return mapreduce(group -> group.slugs, vcat, groups)
end

isMarketClosingSoon(market) = market.isResolved || market.closeTime / 1000 < time() + 60 # if resolved or closing in 60 seconds

function processGroups!(GROUPS, markets)
    # Remove any group with no open markets

    for (name, group) in GROUPS
        allClosed = true

        for url in keys(group)
            slug = urlToSlug(url)

            market = markets[slug]
            if !isMarketClosingSoon(market)
                allClosed = false
            end
        end

        if allClosed
            printstyled("\nDeleted group $name\n", bold=true, blink=true, underline=true)
            pop!(GROUPS, name)
        end
    end
end

function getLimits!(userBalance, groups, betsBySlug)
    limitOrdersBySlug = Dict()
    limitOrdersByProb = Dict()
    splitBets = Dict()
    userIds = Set()

    for group in groups, slug in group.slugs
        tmp, userIdsSlug = sortLimitOrders!(userBalance, betsBySlug[slug])
        union!(userIds, userIdsSlug)
        limitOrdersByProb[slug] = tmp
    end

    @sync for userId in userIds
        if userId ∉ Set(keys(userBalance))
            @async userBalance[userId] = getUserById(userId).balance
        end
    end
    
    for group in groups, slug in group.slugs
        limitOrdersBySlug[slug], splitBets[slug] = getLimitOrders(limitOrdersByProb[slug], userBalance)
    end

    return (limitOrdersBySlug=limitOrdersBySlug, splitBets=splitBets)
end

function calculateMyShares(slugs, betsByMe)
    yesShares = Dict(slugs .=> 0.)
    noShares = Dict(slugs .=> 0.)

    for slug in slugs
        for bet in betsByMe[slug]
            if bet.outcome == "YES"
                yesShares[slug] += bet.shares
            elseif bet.outcome == "NO"
                noShares[slug] += bet.shares
            end
        end
    end

    redeemShares!(noShares, yesShares)
    return (noShares=noShares, yesShares=yesShares)
end


function arbitrage(GROUPS, APIKEY, USERNAME, noSharesBySlug, yesSharesBySlug, oldUserBalance, live=false, confirmBets=true, printDebug=true)
    fetchTime = time()

    markets, betsBySlug, botBalance = getMarketsAndBets!(oldUserBalance, GROUPS, USERNAME)
    
    processGroups!(GROUPS, markets)

    groups = Group.(keys(GROUPS), values(GROUPS))

    limitOrdersBySlug, splitBets = getLimits!(oldUserBalance, groups, betsBySlug)

    # totalNumberOfMarkets = mapreduce(group -> group.noMarkets, +, groups)
    maxNumberOfMarkets = mapreduce(group -> group.noMarkets, max, groups)
    maxBetAmount = botBalance / (3 + 1.5*maxNumberOfMarkets)

    # noSharesBySlug, yesSharesBySlug = calculateShares(groups, betsByMe)

    # sols = map(group -> optimise(group, markets, limitOrdersBySlug, maxBetAmount, noSharesBySlug, yesSharesBySlug, [i for (i, slug) in enumerate(group.slugs) if !isMarketClosingSoon(markets[slug])]), groups)
    # sols = ThreadsX.map(group -> optimise(group, markets, limitOrdersBySlug, maxBetAmount, noSharesBySlug, yesSharesBySlug, [i for (i, slug) in enumerate(group.slugs) if !isMarketClosingSoon(markets[slug])]), groups)

    executedBetPairs = []
    @sync for (groupNumber, group) in enumerate(groups)
        plannedBets = PlannedBet[]
        
        # for market in markets:
        #     skip = skip_market(market)
        #     if skip:
        #         print()
        #         print(skip)
        #         print("Skipping group.\n")
        #         return

        printedGroupName = false

        bettableSlugsIndex = [i for (i, slug) in enumerate(group.slugs) if !isMarketClosingSoon(markets[slug])]

        betAmounts = optimise(group, markets, limitOrdersBySlug, maxBetAmount, noSharesBySlug, yesSharesBySlug, bettableSlugsIndex)
        # betAmounts = sols[groupNumber]

        if time() - fetchTime > 20 # caching so needs to be more than 15, or sleep until then
            markets = getMarkets(getSlugs(groups))
            fetchTime = time()
        end

        newProfitsByEvent, noShares, yesShares, newProb = f(betAmounts, group, markets, limitOrdersBySlug, noSharesBySlug, yesSharesBySlug, bettableSlugsIndex)

        oldProfitsByEvent, _, _, oldProb = f(repeat([0.], length(bettableSlugsIndex)), group, markets, limitOrdersBySlug, noSharesBySlug, yesSharesBySlug, bettableSlugsIndex)

        profit = minimum(newProfitsByEvent) - minimum(oldProfitsByEvent)

        if printDebug
            printstyled("=== $(group.name) ===\n", color=:bold)
            printedGroupName = true

            println(bettableSlugsIndex)
            println()
            println(oldProb)
            println(newProb)
            println()
            println(oldProfitsByEvent)
            println(newProfitsByEvent)
            println(profit)
            println()
            println(betAmounts)
            println()
            println(group.y_matrix)
            println(group.n_matrix)
            println()
            println(getindex.(Ref(yesSharesBySlug), group.slugs))
            println(getindex.(Ref(noSharesBySlug), group.slugs))
            println(yesShares)
            println(noShares)
            println(group.y_matrix * yesShares)
            println(group.n_matrix * noShares)
            println(sum(abs.(betAmounts)))
            println()
        end

        skipMarket = false
        for (i, j) in enumerate(bettableSlugsIndex)
            amount = betAmounts[i]
            slug = group.slugs[j]

            if isapprox(amount, 0., atol=1e-6)
                continue
            elseif abs(amount) >= 1.
                bet = PlannedBet(abs(amount), yesShares[j] + noShares[j], amount > 0. ? "YES" : "NO", markets[slug])
                push!(plannedBets, bet)
            else
                if !printedGroupName
                    printstyled("=== $(group.name) ===\n", color=:bold)
                    printedGroupName = true
                end
                println("Bet amount: $(amount) is too small, $(slug), $(markets[slug].question)")
                skipMarket = true
                break
            end

            if !printedGroupName
                printstyled("=== $(group.name) ===\n", color=:bold)
                printedGroupName = true
            end

            printstyled("$slug\n", color=:red)
            println("Prior probs:     $(markets[slug].probability * 100)%")
            println("Posterior probs: $(newProb[j]*100)%")
        end

        if skipMarket
            continue
        end
                
        if profit <= .01 * length(plannedBets)
            # print()
            if !isempty(plannedBets)
                if !printedGroupName
                    printstyled("=== $(group.name) ===\n", color=:bold)
                    printedGroupName = true
                end

                println("Insufficient profit $(profit) for $(length(plannedBets)) bets")
            end
            continue
        end

        if !printedGroupName
            printstyled("=== $(group.name) ===\n", color=:bold)
            printedGroupName = true
        end
        printstyled("Profits:         $profit\n", color=:yellow)

        for bet in plannedBets
            println(bet)
        end

        if confirmBets
            println("Proceed? (y/n)") 
            if readline() !="y"
                continue
            end
        end

        if live
            for bet in plannedBets 
                @async try
                    push!(executedBetPairs, (bet, execute(bet, APIKEY)))
                catch err
                    bt = catch_backtrace()
                    println()
                    showerror(stderr, err, bt)
                end
            # print(f"Balance: {my_balance()}\n")
            end
        end
    end

    if live
        for (plannedBet, executedBet) in executedBetPairs
            updateShares!(noSharesBySlug, yesSharesBySlug, urlToSlug(plannedBet.market.url), executedBet)
        end
    end
end

function run(groupNames = nothing; live=false, confirmBets=true, printDebug=true)
    data = TOML.parsefile("ArbBot/Arb.toml")
    GROUPS = data["GROUPS"]
    APIKEY = data["APIKEY"]
    USERNAME = data["USERNAME"]

    if groupNames !== nothing
        GROUPS = Dict(name => data["GROUPS"][name] for name in groupNames)
    end
    
    # slugs = getSlugs(GROUPS)
    # markets = getMarkets(slugs)
    # processGroups!(GROUPS, markets)

    betsByMe = Dict{String, Vector{Bet}}()

    @sync for slug in getSlugs(GROUPS)
        @async try
            betsByMe[slug] = getAllBets(slug=slug, username=USERNAME) # If we don't have the exact position arbitrage gets fucked, so we need all bets
        catch err
            bt = catch_backtrace()
            println()
            showerror(stderr, err, bt)
        end
    end

    noSharesBySlug, yesSharesBySlug = calculateMyShares(getSlugs(GROUPS), betsByMe)

    userBalance = Dict{String, Float64}()

    arbitrage(GROUPS, APIKEY, USERNAME, noSharesBySlug, yesSharesBySlug, userBalance, live, confirmBets, printDebug)
end

# run()
# run(live=true)
# run(["Next Speaker of the House"])

function prod(groupNames = nothing; live=true, confirmBets=false, printDebug=false)
    data = TOML.parsefile("ArbBot/Arb.toml")
    GROUPS = data["GROUPS"]
    APIKEY = data["APIKEY"]
    USERNAME = data["USERNAME"]

    if groupNames !== nothing
        GROUPS = Dict(name => data["GROUPS"][name] for name in groupNames)
    end
    
    # slugs = getSlugs(GROUPS)
    # markets = getMarkets(slugs)
    # processGroups!(GROUPS, markets)

    betsByMe = Dict{String, Vector{Bet}}()

    @sync for slug in getSlugs(GROUPS)
        @async try
            betsByMe[slug] = getAllBets(slug=slug, username=USERNAME) # If we don't have the exact position arbitrage gets fucked, so we need all bets
        catch err
            bt = catch_backtrace()
            println()
            showerror(stderr, err, bt)
        end
    end

    noSharesBySlug, yesSharesBySlug = calculateMyShares(getSlugs(GROUPS), betsByMe)

    userBalance = Dict{String, Float64}()

    while true
        # oldTime = time()
        printstyled("Running at $(Dates.format(now(), "HH:MM:SS.sss"))\n"; color = :blue)

        arbitrage(GROUPS, APIKEY, USERNAME, noSharesBySlug, yesSharesBySlug, userBalance, live, confirmBets, printDebug)
        redeemShares!(noSharesBySlug, yesSharesBySlug) # just needs to be run periodically to prevent overflow

        printstyled("Sleeping at $(Dates.format(now(), "HH:MM:SS.sss"))\n"; color = :magenta)

        sleep(120 + 2*(rand()-.5) * 20) # - (time() - oldTime) # add some randomness so it can't be exploited based on predicability of betting time.
    end
end

# https://github.com/innerlee/PrintLog.jl
"""
    @printlog "file.log"
Enable `printlog`.
Outputs of `print` and `println` will be logged into file `file.log`.
It will create the log file if it does not exist.
All prints will be appended to the end of the file.
    @printlog "file.log" silent
Silent mode.
Do not output `info`.
"""
macro printlog(file, silent=false)
    file = "$(@__DIR__)/$file"
    @eval begin
        import Base.println, Base.print, Base.printstyled
        @suppress Base.println(xs...) =
            open(f -> (println(f, xs...); println(stdout, xs...)), $file, "a")
        @suppress Base.print(xs...) =
            open(f -> (print(f, xs...); print(stdout, xs...)), $file, "a")
        @suppress Base.printstyled(xs...; kwargs...) =
            open(f -> (printstyled(f, xs...; kwargs...); printstyled(stdout, xs..., ;  kwargs...)), $file, "a")
    end
    silent != :silent &&
        @info("`print`, `println` and `printstyled` will be logged into file `$file`")
    nothing
end

"""
    @noprintlog
Disable `printlog`.
    @noprintlog silent
Silent mode.
"""
macro noprintlog(silent=false)
    @eval begin
        import Base.println, Base.print, Base.printstyled
        @suppress Base.println(xs...) = println(stdout, xs...)
        @suppress Base.print(xs...) = print(stdout, xs...)
        @suppress Base.printstyled(xs...; kwargs...) = printstyled(stdout, xs...; kwargs...)
    end
    silent != :silent &&
        @info("`print`, `println` and `printstyled` are resumed.")
    nothing
end

@printlog "log.txt"
# prod()