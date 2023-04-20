using ManifoldMarkets, TOML, Optimization, OptimizationBBO, Dates, Combinatorics, LinearAlgebra, Suppressor, Parameters
using HTTP, JSON3, Dates, OpenSSL
using HTTP.WebSockets

const FEE = 0.03
Base.exit_on_sigint(false)

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
@with_kw mutable struct MarketData 
    Shares::Dict{Symbol, Float64} = Dict{Symbol, Float64}(:NO => 0., :YES => 0.)
    limitOrders::Dict{Symbol, Dict{Float64, Vector{Float64}}} = Dict{Symbol, Dict{Float64, Vector{Float64}}}()
    sortedLimitProbs::Dict{Symbol, Vector{Float64}} = Dict{Symbol, Vector{Float64}}(:NO => [], :YES => [])

    probability::Float64 = -1
    p::Float64 = -1
    pool::Dict{Symbol, Float64} = Dict{Symbol, Float64}()

    id::String = ""
    url::String = ""
    question::String = ""

    isResolved::Bool = false
    closeTime::Int = -1
end

@with_kw mutable struct BotData @deftype String
    const APIKEY
    const Supabase_APIKEY
    const USERNAME
    const USERID
    balance::Float64 = 0
end

struct GroupData
    groups::Vector{Group}
    contractIdSet::Set{String}
    contractIdToGroupIndex::Dict{String, Int}
    contractIdToSlug::Dict{String, String}
end

@with_kw struct Arguments @deftype Bool
    live=false
    confirmBets=true
    printDebug=true
end

struct PlannedBet
    amount::Float64
    shares::Float64
    outcome::String
    redeemedMana::Float64

    id::String
    url::String
    question::String
end

function execute(bet, APIKEY)
    ohno = false
    response = createBet(APIKEY, bet.id, bet.amount, bet.outcome)
    # need to check if returned info matches what we wanted to bet, i.e. if we got less shares than we wanted to. If we got more ig either moved or smth weird with limit orders.

    if response.shares ≉ bet.shares
        io = IOBuffer()

        printstyled("\e]8;;$(bet.url)\e\\$(bet.question)\e]8;;\e\\\n", color=:green) # hyperlink
        println(io, response)
        println(io, response.fills)

        write(stdout, take!(io));

        ohno = true
    end

    return response, ohno
end

function updateShares!(MarketData, newBet, BotData)
    MarketData.Shares[Symbol(newBet.outcome)] += newBet.shares

    if MarketData.Shares[:YES] >= MarketData.Shares[:NO]
        MarketData.Shares[:YES] -= MarketData.Shares[:NO]
        BotData.balance += MarketData.Shares[:NO]
        MarketData.Shares[:NO] = 0.
    elseif MarketData.Shares[:YES] < MarketData.Shares[:NO]
        MarketData.Shares[:NO] -= MarketData.Shares[:YES]
        BotData.balance += MarketData.Shares[:YES]
        MarketData.Shares[:YES] = 0.
    end
end

function f(betAmount, group, MarketData, currentNoShares, currentYesShares, bettableSlugsIndex)
    profitsByEvent = zeros(size(group.y_matrix)[1])

    # newProb = zeros(group.noMarkets) # Makes it obvious which markets we don't bet on. We can print this manually, but this hides errors in fetching market probabilities
    newProb = [MarketData[slug].probability for slug in group.slugs] # So we return the correct results for markets we don't bet on  and closing soon markets

    sharesByEvent = group.y_matrix * currentYesShares + group.n_matrix * currentNoShares

    profitsByEvent = f!(betAmount, group, MarketData, bettableSlugsIndex, sharesByEvent, profitsByEvent)

    noShares = copy(currentNoShares)
    yesShares = copy(currentYesShares)

    for (i, j) in enumerate(bettableSlugsIndex)
        slug = group.slugs[j]
        market = MarketData[slug]

        if abs(betAmount[i]) >= 1.
            shares = 0.
            shares, newProb[j] = betToShares(market.p, market.pool, market.probability, market.limitOrders, market.sortedLimitProbs, betAmount[i])

            if betAmount[i] >= 1.
                yesShares[j] += shares
            elseif betAmount[i] <= -1.
                noShares[j] += shares
            end
        end
    end

    return (profitsByEvent=profitsByEvent, noShares=noShares, yesShares=yesShares, newProbability=newProb)
end

function f!(betAmount, group, MarketData, bettableSlugsIndex, sharesByEvent, profitsByEvent)
    fees = 0.
    totalBetAmount = 0.
    profitsByEvent .= sharesByEvent # we shouldn't count resolved markets.

    @inbounds @fastmath for (i, j) in enumerate(bettableSlugsIndex)
        slug = group.slugs[j]
        market = MarketData[slug]

        if abs(betAmount[i]) >= 1.
            shares = 0.
            @inline shares = betToShares(market.p, market.pool, market.probability, market.limitOrders, market.sortedLimitProbs, betAmount[i]).shares

            fees += FEE

            if betAmount[i] >= 1.
                profitsByEvent .+= @view(group.y_matrix[:, j]) .* shares
            elseif betAmount[i] <= -1.
                profitsByEvent .+= @view(group.n_matrix[:, j]) .* shares
            end
        else
            betAmount[i] = 0.
        end

        totalBetAmount += abs(betAmount[i])
    end

    profitsByEvent .-= totalBetAmount + fees

    return profitsByEvent
end

function optimise(group, MarketData, maxBetAmount, bettableSlugsIndex)
    noShares = [MarketData[slug].Shares[:NO] for slug in group.slugs]
    yesShares = [MarketData[slug].Shares[:YES] for slug in group.slugs]
    sharesByEvent = group.y_matrix * yesShares + group.n_matrix * noShares

    profitsByEvent = similar(sharesByEvent)

    profitF = OptimizationFunction((betAmount, _) -> -minimum( f!(betAmount, group, MarketData, bettableSlugsIndex, sharesByEvent, profitsByEvent) ))

    x0 = repeat([0.], length(bettableSlugsIndex))
    ub = repeat([maxBetAmount], length(bettableSlugsIndex))
    lb = -ub

    problem = Optimization.OptimizationProblem(profitF, x0, lb=lb, ub=ub)

    @time "Adaptive 1" sol = solve(problem, BBO_adaptive_de_rand_1_bin_radiuslimited(), maxtime=1)

    bestSolution = repeat([0.], length(bettableSlugsIndex))
    maxRiskFreeProfit = f(bestSolution, group, MarketData, noShares, yesShares, bettableSlugsIndex).profitsByEvent |> minimum


    nonZeroIndices = findall(!iszero, sol.u::Vector{Float64})

    betAmount::Vector{Float64} = repeat([0.], length(bettableSlugsIndex))

    for indices in powerset(nonZeroIndices)
        betAmount::Vector{Float64} .= sol.u::Vector{Float64}
        betAmount[indices] .= 0

        riskFreeProfit = f(betAmount, group, MarketData, noShares, yesShares, bettableSlugsIndex).profitsByEvent |> minimum

        if riskFreeProfit > maxRiskFreeProfit
            maxRiskFreeProfit = riskFreeProfit
            bestSolution .= betAmount
        end
    end
    

    if bestSolution == repeat([0.], length(bettableSlugsIndex))
        @time "Resampling" sol2 = solve(problem, BBO_resampling_memetic_search(), maxtime=3)

        nonZeroIndices = findall(!iszero, sol2.u::Vector{Float64})

        for indices in powerset(nonZeroIndices)
            betAmount::Vector{Float64} .= sol2.u::Vector{Float64}
            betAmount[indices] .= 0

            riskFreeProfit = f(betAmount, group, MarketData, noShares, yesShares, bettableSlugsIndex).profitsByEvent |> minimum

            if riskFreeProfit > maxRiskFreeProfit
                maxRiskFreeProfit = riskFreeProfit
                bestSolution .= betAmount
            end
        end
    end

    return bestSolution
end

function getMarkets!(MarketData, slugs)
    @sync for slug in slugs
        @async try
            market = getMarketBySlug(slug)

            MarketData[slug].probability = market.probability::Float64
            MarketData[slug].p = market.p::Float64 # can change if a subsidy is given
            MarketData[slug].pool = market.pool
            MarketData[slug].id = market.id
            MarketData[slug].question = market.question
            MarketData[slug].url = market.url

            MarketData[slug].isResolved = market.isResolved
            MarketData[slug].closeTime = market.closeTime
        catch err
            bt = catch_backtrace()
            println()
            showerror(stderr, err, bt)
            throw(err)
        end
    end
end

getSlugs(GROUPS::Dict) = mapreduce(x -> urlToSlug.(x), vcat, keys.(values(GROUPS)))
getSlugs(groups::Vector{Group}) = mapreduce(group -> group.slugs, vcat, groups)

isMarketClosingSoon(market) = market.isResolved || market.closeTime / 1000 < time() + 60 # if resolved or closing in 60 seconds

function arbitrageGroup(group, BotData, MarketData, Arguments)
    rerun = :Success

    # Actually could close in the delay between running and here, or due to reruns. But we don't need to pop the group from groups
    allMarketsClosing = true
    for slug in group.slugs
        if !isMarketClosingSoon(MarketData[slug])
            allMarketsClosing = false
        end
    end

    bettableSlugsIndex = [i for (i, slug) in enumerate(group.slugs) if !isMarketClosingSoon(MarketData[slug])]

    if allMarketsClosing
        printstyled("$(group.name) all markets closed at $(Dates.format(now(), "HH:MM:SS.sss"))\n", bold=true, underline=true)

        rerun = :Success
        return rerun
    end

    plannedBets = PlannedBet[]

    printedGroupName = false

    maxBetAmount = BotData.balance / (2 + 1.5*group.noMarkets)
    redeemManaHack = maximum(slug -> max(MarketData[slug].Shares[:YES] / (1 - MarketData[slug].probability), MarketData[slug].Shares[:NO] / MarketData[slug].probability), group.slugs[bettableSlugsIndex])
    maxBetAmount += min(redeemManaHack, maxBetAmount, 100.)

    betAmounts = optimise(group, MarketData, maxBetAmount, bettableSlugsIndex)
    
    oldProb = [MarketData[slug].probability for slug in group.slugs]
    oldNoShares = [MarketData[slug].Shares[:NO] for slug in group.slugs]
    oldYesShares = [MarketData[slug].Shares[:YES] for slug in group.slugs]

    newProfitsByEvent, noShares, yesShares, newProb = f(betAmounts, group, MarketData, oldNoShares, oldYesShares, bettableSlugsIndex)

    oldProfitsByEvent, _, _, _ = f(repeat([0.], length(bettableSlugsIndex)), group, MarketData, oldNoShares, oldYesShares, bettableSlugsIndex)

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
        println(map(slug -> MarketData[slug].Shares[:YES], group.slugs))
        println(map(slug -> MarketData[slug].Shares[:NO], group.slugs))
        println(yesShares)
        println(noShares)
        println(group.y_matrix * yesShares)
        println(group.n_matrix * noShares)
        println(sum(abs.(betAmounts)))
        println()
    end

    for (i, j) in enumerate(bettableSlugsIndex)
        amount = betAmounts[i]
        slug = group.slugs[j]

        if isapprox(amount, 0., atol=1e-6)
            continue
        elseif abs(amount) >= 1.
            outcome = amount > 0. ? "YES" : "NO"
            shares = newYesShares[j] + newNoShares[j]
            redeemedMana = 0.
            if outcome == "YES"
                redeemedMana = min(MarketData[slug].Shares[:NO], shares)
            elseif outcome == "NO"
                redeemedMana = min(MarketData[slug].Shares[:YES], shares)
            end

            bet = PlannedBet(abs(amount), shares, outcome, redeemedMana, MarketData[slug].id, MarketData[slug].url, MarketData[slug].question)
            push!(plannedBets, bet)
        else
            if !printedGroupName
                printstyled("=== $(group.name) ===\n", color=:bold)
                # printedGroupName = true
                rerun = :BetMore
            end
            println("Bet amount: $(amount) is too small, $(slug)")

            rerun = :Success
            return rerun
        end
    end

    sort!(plannedBets, by = bet -> bet.redeemedMana, rev=true)

    if (profit ≤ 0) && !(profit + FEE * length(plannedBets) ≥ 0 && sum(bet -> bet.redeemedMana, plannedBets, init=0.) > 1.)
        rerun = :Success
        return rerun
    end

    newProbBySlug = Dict(group.slugs[j] => newProb[j] for j in bettableSlugsIndex)
    for (i, bet) in enumerate(plannedBets)
        amount = betAmounts[i]
        slug = urlToSlug(bet.url)
        bet = plannedBets[i]

        if abs(amount) >= .98 * maxBetAmount
            printstyled("Bet size is $(100 * abs(amount)/maxBetAmount)% of maxBetAmount\n", color=:orange)
            # bindingConstraint = true
            rerun = :BetMore
        end

        if abs(amount) ≈ 1
            printstyled("Bet size is $(abs(amount))\n", color=:orange)
            # bindingConstraint = true
            rerun = :BetMore
        end

        if !printedGroupName
            printstyled("=== $(group.name) ===\n", color=:bold)
            printedGroupName = true
        end

        printstyled("\e]8;;$(MarketData[slug].url)\e\\$(MarketData[slug].question)\e]8;;\e\\\n", color=:red) # hyperlink
        println("Prior probs:     $(MarketData[slug].probability * 100)%")
        println("Posterior probs: $(newProbBySlug[slug]*100)%")
        printstyled("Buy $(bet.shares) $(bet.outcome) shares for $(bet.amount), redeeming $(bet.redeemedMana)\n", color=:green)
    end

    if sum(abs.(betAmounts)) >= sum(bet -> bet.redeemedMana, plannedBets) + BotData.balance - 100
        println("Insufficient Balance $(BotData.balance) for $(sum(abs.(betAmounts))) bet redeeming $(sum(bet -> bet.redeemedMana, plannedBets)).")
        rerun = :Success
        return rerun
    end

    if !printedGroupName
        printstyled("=== $(group.name) ===\n", color=:bold)
        printedGroupName = true
    end

    printstyled("Profits:         $(profit + FEE * length(plannedBets))\n", color=:yellow) # no more fee, but we still want to use fee in optimisation

    if Arguments.confirmBets
        println("Proceed? (y/n)") 
        if readline() !="y"
            rerun = :Success
            return rerun
        end
    end

    if Arguments.live
        @sync try 
            for bet in plannedBets 
                @async begin
                    executedBet, ohno = execute(bet, BotData.APIKEY)

                    slug = urlToSlug(bet.url)

                    BotData.balance -= executedBet.amount

                    if ohno
                        rerun = :UnexpectedBet
                    end

                    updateShares!(MarketData[slug], executedBet, BotData)

                    MarketData[slug].probability = executedBet.probAfter # update otherwise we might arb again if bets are duplicated.

                    if !isnothing(executedBet.fills[end].matchedBetId)

                        limitOrder = getBet(executedBet.fills[end].matchedBetId, BotData.Supabase_APIKEY)

                        amountLeft = 0.
                        sharesLeft = 0.
                        if limitOrder.outcome == "NO"
                            sharesLeft = (limitOrder.orderAmount - limitOrder.amount) / (1 - limitOrder.limitProb)
                            amountLeft = sharesLeft * limitOrder.limitProb
                        elseif limitOrder.outcome == "YES"
                            sharesLeft = (limitOrder.orderAmount - limitOrder.amount) / limitOrder.limitProb
                            amountLeft = sharesLeft * (1 - limitOrder.limitProb)
                        end

                        MarketData[slug].limitOrders[Symbol(executedBet.outcome)] = Dict(limitOrder.limitProb => [amountLeft, sharesLeft])

                        if executedBet.outcome == "YES"
                            MarketData[slug].sortedLimitProbs = Dict(:YES=>[limitOrder.limitProb], :NO=>[])
                        elseif executedBet.outcome == "NO"
                            MarketData[slug].sortedLimitProbs = Dict(:YES=>[], :NO=>[limitOrder.limitProb])
                        end
                    end
                end
            end
        catch err
            bt = catch_backtrace()
            println()
            showerror(stderr, err, bt)

            fetchMyShares!(marketDataBySlug, BotData.USERID)

            rerun = :PostFailure
        end

        return rerun
    else
        return :Success
    end
end

function arb(bet, GroupData, BotData, MarketData, Arguments)
    marketId = bet.contractId
    slug = GroupData.contractIdToSlug[marketId]

	if marketId in GroupData.contractIdSet && MarketData[slug].probability ≉ bet.probAfter && bet.userUsername != BotData.USERNAME  
        MarketData[slug].probability = bet.probAfter # prevents this code trigerring if the same bet is fedthrough, doesn't matter that we fetch true prob later as due to async the if check happens before that.

        group = GroupData.groups[GroupData.contractIdToGroupIndex[marketId]]

        rerun = :FirstRun
        runs = 0
        delay = 60
        
        while rerun == :FirstRun || (rerun == :BetMore && runs ≤ 5) || (rerun == :UnexpectedBet && runs ≤ 10) || rerun == :PostFailure
            if rerun == :PostFailure
                sleep(delay)
                delay *= 5
            end 
            printstyled("Running $(group.name)\n", color=:light_cyan)

            try
                if rerun == :FirstRun
                    # @time "Get Market" getMarkets!(MarketData, [slug]) # function takes in a vector of slugs, otherwise iterates over characters of slug
                    @time "Get All Markets" getMarkets!(MarketData, group.slugs)
                else
                    # Needed as we need to update p and pool after betting
                    @time "Get All Markets" getMarkets!(MarketData, group.slugs)
                end
                rerun = arbitrageGroup(group, BotData, MarketData, Arguments)
            catch err
                bt = catch_backtrace()
                println()
                showerror(stderr, err, bt)
                rethrow()
            end
            runs += 1
        end

        for slug in GroupData.groups[GroupData.contractIdToGroupIndex[marketId]].slugs
            MarketData[slug].limitOrders = Dict{Symbol, Dict{Float64, Vector{Float64}}}() # need to reset as we aren't tracking limit orders
            MarketData[slug].sortedLimitProbs = Dict(:YES=>[], :NO=>[])
        end
	end

    return nothing
end

function fetchMyShares!(MarketDataBySlug, USERID)
    @sync for (slug, MarketData) in MarketDataBySlug
        @async try
            positions = getPositionsOnMarket(MarketData.id, userId=USERID)
            if !isempty(positions)
                for (outcome, shares) in positions[1].totalShares
                    MarketData.Shares[Symbol(outcome)] = shares
                end
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
    APIKEY::String = data["APIKEY"]
    Supabase_APIKEY::String = data["SUPABASE_APIKEY"]
    USERNAME::String = data["USERNAME"]

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

    return GROUPS, APIKEY, Supabase_APIKEY, USERNAME
end

function setup(groupNames, live, confirmBets, printDebug)
    GROUPS::Dict{String, Dict{String, Vector{String}}}, APIKEY, Supabase_APIKEY, USERNAME = readData()  

    if groupNames !== nothing
        GROUPS = Dict(name => GROUPS[name] for name in groupNames)
    end

    groups = Group.(keys(GROUPS), values(GROUPS))
    
    marketDataBySlug = Dict(slug => MarketData() for slug in getSlugs(groups))

    printstyled("Fetching at $(Dates.format(now(), "HH:MM:SS.sss"))\n"; color = :green)
    getMarkets!(marketDataBySlug, getSlugs(groups))
    botUser = getUserByUsername(USERNAME)
    USERID = botUser.id
    botBalance = botUser.balance
    printstyled("Fetching Shares at $(Dates.format(now(), "HH:MM:SS.sss"))\n"; color = :green)
    fetchMyShares!(marketDataBySlug, USERID)
    printstyled("Done fetching at $(Dates.format(now(), "HH:MM:SS.sss"))\n"; color = :green)

    contractIdSet = Set(market.id for market in values(marketDataBySlug))
    contractIdToGroupIndex = Dict(marketDataBySlug[slug].id => i for (i, group) in enumerate(groups) for slug in group.slugs)
    contractIdToSlug = Dict(marketDataBySlug[slug].id => slug for group in groups for slug in group.slugs)

    botData = BotData(APIKEY, Supabase_APIKEY, USERNAME, USERID, botBalance)
    arguments = Arguments(live, confirmBets, printDebug)
    groupData = GroupData(groups, contractIdSet, contractIdToGroupIndex, contractIdToSlug)

    return groupData, botData, marketDataBySlug, arguments
end

function production(groupNames = nothing; live=true, confirmBets=false, printDebug=false, skip=false)
    groupData, botData, marketDataBySlug, arguments = setup(groupNames, live, confirmBets, printDebug)

    if !skip
        for group in groupData.groups
            rerun = :FirstRun
            runs = 0
            
            while rerun == :FirstRun || (rerun == :BetMore && runs ≤ 5) || (rerun == :UnexpectedBet && runs ≤ 10) || rerun == :PostFailure 
                @time "Get All Markets" getMarkets!(MarketData, group.slugs)
                rerun = arbitrageGroup(group, botData, marketDataBySlug, arguments)

                runs += 1
            end
        end
    end

    TaskDict = Dict(i => (runAgain=false, task=@async nothing) for i in eachindex(groupData.groups)) # so we don't have to check if there is a task in it or not. Order of tuple matters for parsing

    WebSockets.open(uri(botData.Supabase_APIKEY)) do socket
        println("Opened Socket")
        println(socket)
    
        # send(socket, pushJSON("contracts", "live-contracts"))
        send(socket, pushJSON("contract_bets"))
        println("Sent Intialisation")
    
        try
            @sync begin
                #Reading messages
                @async for msg in socket
                    @async begin
                        msgJSON = JSON3.read(msg)
                        if :payload in keys(msgJSON) && :data in keys(msgJSON.payload)
                            println("Received message: $(msgJSON.payload.data.record.data.userUsername)")

                            try
                                marketId = msgJSON.payload.data.record.data.contractId
                            
                                if marketId in groupData.contractIdSet && !TaskDict[groupData.contractIdToGroupIndex[marketId]].runAgain
                                    TaskDict[groupData.contractIdToGroupIndex[marketId]] = (runAgain = true, task=TaskDict[groupData.contractIdToGroupIndex[marketId]].task)
                                    wait(TaskDict[groupData.contractIdToGroupIndex[marketId]].task)
                                    TaskDict[groupData.contractIdToGroupIndex[marketId]] = (runAgain = false, task=current_task())
                                    arb(msgJSON.payload.data.record.data, groupData, botData, marketDataBySlug, arguments)
                                end
                            catch err
                                bt = catch_backtrace()
                                println()
                                showerror(stderr, err, bt)
                                rethrow(err)
                            end
                        else
                            println(msg)
                        end
                    end
                end
    
                # HeartBeat
                @async while !WebSockets.isclosed(socket)
                    printstyled("Heartbeat at $(Dates.format(now(), "HH:MM:SS.sss"))\n"; color = :blue)
                    send(socket, heartbeatJSON)
                    printstyled("Sleeping at $(Dates.format(now(), "HH:MM:SS.sss"))\n"; color = :blue)
                    sleep(30)
                end

                # Fetch new balance at 8am
                @async while !WebSockets.isclosed(socket)
                    printstyled("Fetcing Balance at $(Dates.format(now(), "HH:MM:SS.sss"))\n"; color = :blue)
                    botData.balance = getUserByUsername(botData.USERNAME).balance
                    printstyled("Sleeping Balance at $(Dates.format(now(), "HH:MM:SS.sss"))\n"; color = :blue)
                    timeTo8 = Second(((today() + Time(8) + Minute(5) - now()) ÷ 1000).value)
                    timeTo8 += timeTo8 > Second(0) ? Second(0) : Second(Day(1))
                    sleep(timeTo8)
                end
            end
        catch err
            if err isa InterruptException
                println("Caught Interrupt")
                return
            end
            bt = catch_backtrace()
            println("Caught")
            showerror(stderr, err, bt)
            rethrow(err)

            println("Finally Caught")
            if !WebSockets.isclosed(socket)
                send(socket, leaveJSON())
                println("Left Channel in catch")
    
                msg = receive(socket)
                println("Received message: $msg")
                # println("Received message: $(JSON.parse(String(msg)))")
    
                close(socket)
            end
        finally
            println("Finally")
            if !WebSockets.isclosed(socket)
                send(socket, leaveJSON("live-contracts"))
                println("Left Channel")
    
                # msg = receive(socket)
                # println("Received message: $msg")
                # println("Received message: $(JSON.parse(String(msg)))")
    
                close(socket)
            end
        end
    
        println("Finished")
        println(WebSockets.isclosed(socket))
    end
end

function test(groupNames = nothing; live=false, confirmBets=true, printDebug=true, skip=false)
    production(groupNames; live=live, confirmBets=confirmBets, printDebug=printDebug, skip=skip)
end

function retryProd(groupNames = nothing; live=true, confirmBets=false, printDebug=false, skip=false)
    delay = 60
    lastRunTime = time()

    while true
        try
            production(groupNames; live=live, confirmBets=confirmBets, printDebug=printDebug, skip=skip)
        catch err
            bt = catch_backtrace()
            println()
            showerror(stderr, err, bt)

            if err isa MethodError || err isa InterruptException
                throw(err)
                break
            end
            
            if time() - lastRunTime > 3600
                delay = 60
            end
            sleep(delay)

            delay *= 5
            lastRunTime = time()
        end
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