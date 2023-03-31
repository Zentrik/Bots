using ManifoldMarkets, TOML, Optimization, OptimizationBBO, Dates, Combinatorics, LinearAlgebra, Suppressor

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

mutable struct MarketData 
    Shares::Dict{Symbol, F}
    limitOrders::Dict{String, Vector{Float64}}
end

struct BotData @deftype String
    APIKEY
    USERNAME
    USERID
end

struct GroupData
    groups::Vector{Group}
    contractIdSet::Set{String}
    contractIdToGroupIndex
    contractIdToSlug
end

@with_kw struct Arguments @deftype Bool
    live=false
    confirmBets=true
    printDebug=true
end

function Base.show(io::IO, plannedBet::PlannedBet)
    printstyled(io, "\e]8;;$(plannedBet.market.url)\e\\$(plannedBet.market.question)\e]8;;\e\\\n", color=:green) # hyperlink
    print(io, "Buy $(plannedBet.shares) $(plannedBet.outcome) shares for $(plannedBet.amount)")
end

function execute(bet, APIKEY)
    ohno = false
    response = createBet(APIKEY, bet.market.id, bet.amount, bet.outcome)
    # need to check if returned info matches what we wanted to bet, i.e. if we got less shares than we wanted to. If we got more ig either moved or smth weird with limit orders.
    if response.shares ≉ bet.shares
        printstyled("\e]8;;$(bet.market.url)\e\\$(bet.market.question)\e]8;;\e\\\n", color=:green) # hyperlink
        println(bet.market.url)
        println(response)
        println(response.fills)

        ohno = true
    end

    return response, ohno
end

function updateShares!(noShares, yesShares, newBet)
    if newBet.outcome == "NO"
        noShares += newBet.shares
    elseif newBet.outcome == "YES"
        yesShares += newBet.shares
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

function f(betAmount, group, markets, limitOrdersBySlug, sortedLimitProbs, currentNoShares, currentYesShares, bettableSlugsIndex)
    A = zeros(size(group.y_matrix)[1])
    B = zeros(size(group.n_matrix)[1])
    profitsByEvent = zeros(size(group.y_matrix)[1])

    # newProb = zeros(group.noMarkets) # Makes it obvious which markets we don't bet on. We can print this manually, but this hides errors in fetching market probabilities
    newProb = [markets[slug].probability::Float64 for slug in group.slugs] # So we return the correct results for markets we don't bet on  and closing soon markets

    return f!(betAmount, group, markets, limitOrdersBySlug, sortedLimitProbs, currentNoShares, currentYesShares, bettableSlugsIndex, copy(newProb), A, B, profitsByEvent)
end

function f!(betAmount, group, markets, limitOrdersBySlug, sortedLimitProbs, currentNoShares, currentYesShares, bettableSlugsIndex, newProb, A, B, profitsByEvent)
    noShares = copy(currentNoShares)
    yesShares = copy(currentYesShares)

    fees = 0.

    for (i, j) in enumerate(bettableSlugsIndex)
        slug = group.slugs[j]
        market = markets[slug]

        if abs(betAmount[i]) >= 1.
            pool = market.pool
            shares = 0.
            shares, newProb[j] = betToShares(market.p::Float64, pool, market.probability::Float64, limitOrdersBySlug[slug], sortedLimitProbs[slug], betAmount[i])

            # fees += 0.1
            if betAmount[i] >= 1.
                yesShares[j] += shares
            elseif betAmount[i] <= -1.
                noShares[j] += shares
            end
        else
            betAmount[i] = 0.
        end
    end

    profitsByEvent .= mul!(A, group.y_matrix, yesShares) .+ mul!(B, group.n_matrix, noShares) .- sum(abs.(betAmount)) .- fees # we shouldn't count resolved markets.

    return (profitsByEvent=profitsByEvent, noShares=noShares, yesShares=yesShares, newProbability=newProb)
end

function optimise(group, markets, limitOrdersBySlug, sortedLimitProbs, maxBetAmount, noSharesBySlug, yesSharesBySlug, bettableSlugsIndex)
    newProb = zeros(group.noMarkets)
    noShares = [noSharesBySlug[slug] for slug in group.slugs]
    yesShares = [yesSharesBySlug[slug] for slug in group.slugs]

    A = zeros(size(group.y_matrix)[1])
    B = zeros(size(group.n_matrix)[1])
    profitsByEvent = zeros(size(group.y_matrix)[1])

    profitF = OptimizationFunction((betAmount, _) -> -minimum( f!(betAmount, group, markets, limitOrdersBySlug, sortedLimitProbs, noShares, yesShares, bettableSlugsIndex, newProb, A, B, profitsByEvent).profitsByEvent ))

    x0 = repeat([0.], length(bettableSlugsIndex))
    lb = repeat([-maxBetAmount], length(bettableSlugsIndex))
    ub = repeat([maxBetAmount], length(bettableSlugsIndex))

    problem = Optimization.OptimizationProblem(profitF, x0, lb=lb, ub=ub)

    sol = solve(problem, BBO_adaptive_de_rand_1_bin_radiuslimited(), maxtime=4.)

    bestSolution = repeat([0.], length(bettableSlugsIndex))
    maxRiskFreeProfit = f(bestSolution, group, markets, limitOrdersBySlug, sortedLimitProbs, noShares, yesShares, bettableSlugsIndex).profitsByEvent |> minimum

    nonZeroIndices = findall(!iszero, sol.u::Vector{Float64})

    for solution in (sol,)
        for indices in powerset(nonZeroIndices)
            betAmount::Vector{Float64} = copy(solution.u)
            betAmount[indices] .= 0

            riskFreeProfit = f(betAmount, group, markets, limitOrdersBySlug, sortedLimitProbs, noShares, yesShares, bettableSlugsIndex).profitsByEvent |> minimum

            if riskFreeProfit > maxRiskFreeProfit
                maxRiskFreeProfit = riskFreeProfit
                bestSolution = betAmount
            end
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

function getMarketsAndBalance!(group, USERNAME)
    markets = Dict{String, Market}()
    myBalance = 0.

    @sync begin 
        for slug in group.slugs
            @async try
                markets[slug] = getMarketBySlug(slug)
            catch err
                bt = catch_backtrace()
                println()
                showerror(stderr, err, bt)
                throw(err)
            end

        end

        myBalance = getUserByUsername(USERNAME).balance
    end

    return (markets=markets, myBalance=myBalance)
end

getSlugs(GROUPS::Dict) = mapreduce(x -> urlToSlug.(x), vcat, keys.(values(GROUPS)))

getSlugs(groups::Vector{Group}) = mapreduce(group -> group.slugs, vcat, groups)

isMarketClosingSoon(market) = market.isResolved || market.closeTime / 1000 < time() + 60 # if resolved or closing in 60 seconds

# Returns all bets since lastBetId exclusive.
function getLiveBets(lastBetId)
    gotAllBets = false
    Bets = Bet[]

    numberOfBetsFetched = 0
    lastBetFetched = nothing

    while !gotAllBets && numberOfBetsFetched <= 100
        toFetch = 3+rand(0:5)
        tmp = getBets(limit=toFetch, before=lastBetFetched)

        for (i, bet) in enumerate(tmp) # probably faster to go in reverse
            if bet.id == lastBetId
                gotAllBets = true
                append!(Bets, tmp[1:i-1])
                break
            end
        end

        numberOfBetsFetched += toFetch

        lastBetFetched = tmp[end].id

        if !gotAllBets
            append!(Bets, tmp)
        end
    end

    return Bets
end

function updateLimitOrders!(MarketData, bet)
    for (i, fill) in enumerate(bet.fills)
        if !isnothing(fill.matchedBetId)
            if i > 1
                delete!(MarketData.limitOrders[bet.outcome], fill.matchedBetId) # Makes sure we delete if user runs out of balance before filling order
            else
                MarketData.limitOrders[bet.outcome][fill.matchedBetId][2] -= fill.amount # we are tracking how much we can trade against the limit order

                if limitOrders[fill.matchedBetId][2] < 1e-2
                    delete!(MarketData.limitOrders[bet.outcome], fill.matchedBetId)
                end
            end
        end
    end

    MarketData.sortedLimitProbs[bet.outcome] = limitOrders[:YES] |> keys |> collect
    if bet.outcome == "YES"
        sort!(MarketData.sortedLimitProbs[bet.outcome])
    elseif bet.outcome == "NO"
        sort!(MarketData.sortedLimitProbs[bet.outcome], rev=true)
    end
end

function addLimitOrder!(limitOrders, bet)
    MarketData[slug].limitOrders[bet.outcome == "YES" ? :NO : :YES][bet.id] = [bet.limitProb, bet.orderAmount - bet.amount]

                limitProb = limitOrders[fill.matchedBetId][1]

                if executedBet.outcome == :YES
                    limitOrders[fill.matchedBetId][2] -=  fill.shares * (1 - limitProb) # The limit order is buying NO
                elseif executedBet.outcome == :NO
                    limitOrders[fill.matchedBetId][2] -= fill.shares * limitProb # The limit order is buying YES
                end

                if limitOrders[fill.matchedBetId][2] < 1e-2
                    delete!(limitOrders, fill.matchedBetId)
                end
    end
end

function arbitrageGroup(group, BotData, MarketData, Arguments)
    markets, botBalance = getMarketsAndBalance!(group, BotData.USERNAME)

    maxBetAmount = botBalance / (3 + 1.5*group.noMarkets)

    # Actually could close in the delay between running and here, or due to reruns. But we don't need to pop the group from groups
    allMarketsClosing = true
    for slug in group.slugs
        if !isMarketClosingSoon(markets[slug])
            allMarketsClosing = false
        end
    end

    if allMarketsClosing
        printstyled("$(group.name) all markets closed at $(Dates.format(now(), "HH:MM:SS.sss"))\n", bold=true, underline=true)
        return false
    end

    plannedBets = PlannedBet[]

    printedGroupName = false

    bettableSlugsIndex = [i for (i, slug) in enumerate(group.slugs) if !isMarketClosingSoon(markets[slug])]

    betAmounts = optimise(group, markets, MarketData, maxBetAmount, MarketData[slug].Shares[:NO], MarketData[slug].Shares[:YES], bettableSlugsIndex)

    oldProb = [markets[slug].probability::Float64 for slug in group.slugs]
    oldNoShares = [MarketData[slug].Shares[:NO] for slug in group.slugs]
    oldYesShares = [MarketData[slug].Shares[:YES] for slug in group.slugs]

    newProfitsByEvent, noShares, yesShares, newProb = f(betAmounts, group, markets, MarketData, oldNoShares, oldYesShares, bettableSlugsIndex)

    oldProfitsByEvent, _, _, _ = f(repeat([0.], length(bettableSlugsIndex)), group, markets, MarketData, oldNoShares, oldYesShares, bettableSlugsIndex)

    profit = minimum(newProfitsByEvent) - minimum(oldProfitsByEvent)
    newYesShares = yesShares .- oldYesShares
    newNoShares = noShares  .- oldNoShares

    if Arguments.printDebug
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
        println(getindex.(Ref(MarketData[slug].Shares[:YES]), group.slugs))
        println(getindex.(Ref(MarketData[slug].Shares[:NO]), group.slugs))
        println(yesShares)
        println(noShares)
        println(group.y_matrix * yesShares)
        println(group.n_matrix * noShares)
        println(sum(abs.(betAmounts)))
        println()
    end

    if profit <= .05 * length(plannedBets)
        return false
    end

    bindingConstraint = false # whether we would bet more if the maxBetAmount was larger
    for (i, j) in enumerate(bettableSlugsIndex)
        amount = betAmounts[i]
        slug = group.slugs[j]

        if isapprox(amount, 0., atol=1e-6)
            continue
        elseif abs(amount) >= 1.
            bet = PlannedBet(abs(amount), newYesShares[j] + newNoShares[j], amount > 0. ? "YES" : "NO", markets[slug]::Market)
            push!(plannedBets, bet)
            
            if abs(amount) >= .98 * maxBetAmount
                printstyled("Bet size is $(100 * abs(amount)/maxBetAmount)% of maxBetAmount\n", color=:orange)
                bindingConstraint = true
            end

            if abs(amount) ≈ 1
                printstyled("Bet size is $(abs(amount))\n", color=:orange)
                bindingConstraint = true
            end
        else
            if !printedGroupName
                printstyled("=== $(group.name) ===\n", color=:bold)
                printedGroupName = true
            end
            println("Bet amount: $(amount) is too small, $(slug), $(markets[slug].question)")
            return false
            break
        end

        if !printedGroupName
            printstyled("=== $(group.name) ===\n", color=:bold)
            printedGroupName = true
        end

        printstyled("\e]8;;$(markets[slug].url)\e\\$(markets[slug].question)\e]8;;\e\\\n", color=:red) # hyperlink
        println("Prior probs:     $(markets[slug].probability * 100)%")
        println("Posterior probs: $(newProb[j]*100)%")
        printstyled("Buy $(bet.shares) $(bet.outcome) shares for $(bet.amount)\n", color=:green)
    end

    if sum(abs.(betAmounts)) >= botBalance - 100
        println("Insufficient Balance $botBalance for $(sum(abs.(betAmounts))) bet")
        return false
    end

    if !printedGroupName
        printstyled("=== $(group.name) ===\n", color=:bold)
        printedGroupName = true
    end
    printstyled("Profits:         $profit\n", color=:yellow)

    if Arguments.confirmBets
        println("Proceed? (y/n)") 
        if readline() !="y"
            return false
        end
    end

    rerun = bindingConstraint

    if Arguments.live
        @sync for bet in plannedBets 
            @async try
                executedBet, ohno = execute(bet, BotData.APIKEY)

                if ohno
                    rerun = true
                end

                # botBalance -= executedBet.amount

                slug =  urlToSlug(bet.market.url)
                updateShares!(MarketData[slug].Shares[:NO], MarketData[slug].Shares[:YES], executedBet)
                updateLimitOrders!(MarketData[slug].limitOrders, executedBet)
            catch err
                bt = catch_backtrace()
                println()
                showerror(stderr, err, bt)
                throw(err)
            end
        end
    end

    return rerun * Arguments.live
end

function arbitrage(GroupData, BotData, MarketData, lastBetId, Arguments)
    bets = getLiveBets(lastBetId)
    println(length(bets))

    seenGroups = Set{Int}()

    for bet in bets
        if bet.contractId in GroupData.contractIdSet && bet.userUsername != BotData.USERNAME 
            slug = contractIdToSlug[bet.contractId]

            if GroupData.contractIdToGroupIndex[bet.contractId] ∉ seenGroups && (isnothing(bet.limitProb) || !(bet.probAfter ≈ bet.probBefore))
                updateLimitOrders!(MarketData[slug].limitOrders, bet)

                rerun = true
                # runs = 0
                
                while rerun# && runs <= 5
                    printstyled("Running $(GroupData.groups[GroupData.contractIdToGroupIndex[bet.contractId]].name)\n", color=:light_cyan)
                    rerun = arbitrageGroup(GroupData.groups[GroupData.contractIdToGroupIndex[bet.contractId]], BotData, MarketData, Arguments)
                end
                push!(seenGroups, GroupData.contractIdToGroupIndex[bet.contractId])

            elseif !isnothing(bet.limitProb)
                MarketData[slug].limitOrders[bet.id] = [bet.limitProb, bet.orderAmount - bet.amount]
            end
        end
    end

    if length(bets) == 0
        return lastBetId
    else
        return bets[1].id
    end
end

function fetchMyShares!(MarketData, markets, USERID)
    @sync for market in markets
        @async try
            MarketData[urlToSlug(market.url)].Shares = stringKeysToSymbol(getPositionsOnMarket(market.id, userId=USERID).totalShares)
        catch err
            bt = catch_backtrace()
            println()
            showerror(stderr, err, bt)
            throw(err)
        end
    end
end

function fetchLimitOrders!(MarketData, markets)
    @sync for market in markets
        @async try
            limitOrders = getLimitOrdersOnMarket(market.id) # need to handle case of more than 1000 limit orders
            
            for bet in limitOrders
                MarketData[urlToSlug(market.url)].limitOrders[bet.id] = [bet.limitProb, bet.orderAmount - bet.amount]
            end
        catch err
            bt = catch_backtrace()
            println()
            showerror(stderr, err, bt)
            throw(err)
        end
    end
end

function readData()
    data = TOML.parsefile("$(@__DIR__)/Arb.toml")
    GROUPS::Dict{String, Dict{String, Vector{String}}} = data["GROUPS"]
    APIKEY = data["APIKEY"]
    USERNAME = data["USERNAME"]

    slugs = getSlugs(GROUPS)
    if !allunique(slugs)
        println("Duplicate slug detected")

        sort!(slugs)
        for i in 1:length(slugs)-1
            if slugs[i] == slugs[i+1]
                println(slugs[i])
            end
        end

        error()
    end

    for slug in slugs
        if '#' in slug
            println("Invalid Slug $slug")
            error()
        end
    end

    return GROUPS, APIKEY, USERNAME
end

function testIndividualGroup(live=false, confirmBets=true, printDebug=true)
    GROUPS, APIKEY, USERNAME = readData()

    groups = Group.(keys(GROUPS), values(GROUPS))
    
    # slugs = getSlugs(GROUPS)
    # markets = getMarkets(slugs)
    # processGroups!(GROUPS, markets)

    printstyled("Fetching my Shares at $(Dates.format(now(), "HH:MM:SS.sss"))\n"; color = :green)
    noSharesBySlug, yesSharesBySlug, amountInvested = fetchMyShares(Dict(groups[1].name => GROUPS[groups[1].name]), USERNAME)
    # noSharesBySlug, yesSharesBySlug, amountInvested = fetchMyShares(GROUPS, USERNAME)
    printstyled("Done fetching my Shares at $(Dates.format(now(), "HH:MM:SS.sss"))\n"; color = :green)

    oldUserBalances = Dict{String, Float64}()

    printstyled("Running at $(Dates.format(now(), "HH:MM:SS.sss"))\n"; color = :blue)
    arbitrageGroup(groups[1], BotData(APIKEY, USERNAME), MarketData(noSharesBySlug, yesSharesBySlug, oldUserBalances, limitOrdersBySlug), Arguments())
    printstyled("Done at $(Dates.format(now(), "HH:MM:SS.sss"))\n"; color = :magenta)
end

function test(groupNames = nothing; live=false, confirmBets=true, printDebug=true)
    GROUPS, APIKEY, USERNAME = readData()  

    if groupNames !== nothing
        GROUPS = Dict(name => GROUPS[name] for name in groupNames)
    end

    groups = Group.(keys(GROUPS), values(GROUPS))

    marketData = MarketData(Dict{Symbol, F}(), Dict(slug => Dict{String, Vector{Float64}}() for slug in getSlugs(groups))) # need to fetch limit orders

    printstyled("Fetching at $(Dates.format(now(), "HH:MM:SS.sss"))\n"; color = :green)

    lastBetId = getBets(limit=1)[1].id
    markets = getMarkets(getSlugs(groups))
    fetchMyShares!(marketData, markets, USERID)
    fetchLimitOrders!(marketData, markets)
    USERID = getUserByUsername(USERNAME).id

    printstyled("Done fetching at $(Dates.format(now(), "HH:MM:SS.sss"))\n"; color = :green)

    contractIdSet = Set(market.id for market in values(markets))
    contractIdToGroupIndex = Dict(markets[slug].id => i for (i, group) in enumerate(groups) for slug in group.slugs)
    contractIdToSlug = Dict(markets[slug].id => slug for group in groups for slug in group.slugs)

    botData = BotData(APIKEY, USERNAME, USERID)
    arguments = Arguments(live, confirmBets, printDebug)
    groupData = GroupData(groups, contractIdSet, contractIdToGroupIndex, contractIdToSlug)

    for group in groups
        rerun = true
        # runs = 0
        
        while rerun# && runs <= 5
            rerun = arbitrageGroup(group, botData, marketData, arguments)
        end
    end

    println("Sleeping")
    sleep(15)

    printstyled("Running at $(Dates.format(now(), "HH:MM:SS.sss"))\n"; color = :blue)

    lastBetId = arbitrage(groupData, botData, marketData, lastBetId, arguments)
    printstyled("Done at $(Dates.format(now(), "HH:MM:SS.sss"))\n"; color = :magenta)
end

function production(groupNames = nothing; live=true, confirmBets=false, printDebug=false, skip=false)
    GROUPS, APIKEY, USERNAME = readData()  

    if groupNames !== nothing
        GROUPS = Dict(name => GROUPS[name] for name in groupNames)
    end

    groups = Group.(keys(GROUPS), values(GROUPS))
    
    printstyled("Fetching my Shares at $(Dates.format(now(), "HH:MM:SS.sss"))\n"; color = :green)
    noSharesBySlug, yesSharesBySlug = fetchMyShares(GROUPS, USERNAME)
    printstyled("Done fetching my Shares at $(Dates.format(now(), "HH:MM:SS.sss"))\n"; color = :green)

    userBalance = Dict{String, Float64}()

    lastBetId = getBets(limit=1)[1].id
    markets = getMarkets(getSlugs(groups))
    contractIdSet = Set(market.id for market in values(markets))

    contractIdToGroupIndex = Dict(markets[slug].id => i for (i, group) in enumerate(groups) for slug in group.slugs)

    botData = BotData(APIKEY, USERNAME)
    marketData = MarketData(noSharesBySlug, yesSharesBySlug, oldUserBalances, limitOrdersBySlug)
    arguments = Arguments(live, confirmBets, printDebug)
    groupData = GroupData(groups, contractIdSet, contractIdToGroupIndex, contractIdToSlug)


    if !skip
        for group in groups
            rerun = true
            # runs = 0
            
            while rerun# && runs <= 5
                rerun = arbitrageGroup(group, botData, marketData, arguments)
            end
        end
    end

    while true
        # oldTime = time()
        printstyled("Running at $(Dates.format(now(), "HH:MM:SS.sss"))\n"; color = :blue)

        lastBetId = arbitrage(groupData, botData, marketData, lastBetId, arguments)
        redeemShares!(noSharesBySlug, yesSharesBySlug) # just needs to be run periodically to prevent overflow

        printstyled("Sleeping at $(Dates.format(now(), "HH:MM:SS.sss"))\n"; color = :magenta)

        # sleep(15)
        sleep(15 + rand())
        # sleep(60 + 2*(rand()-.5) * 5) # - (time() - oldTime) # add some randomness so it can't be exploited based on predicability of betting time.
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

# @printlog "log.txt"
# production()