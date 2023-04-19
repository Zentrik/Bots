using ManifoldMarkets, TOML, Optimization, OptimizationBBO, Dates, Combinatorics, LinearAlgebra, Suppressor, Parameters

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

@with_kw struct BotData @deftype String
    APIKEY
    Supabase_APIKEY
    USERNAME
    USERID
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

    # response HTTP.post("http://manifold.markets/api/v0/bet", headers = ["Authorization" => "Key " * APIKEY, "Content-Type" => "application/json"], body="{\"amount\":$(bet.amount),\"outcome\":\"$(bet.outcome)\",\"contractId\":\"$(bet.id)\"}")

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

function updateShares!(MarketData, newBet)
    MarketData.Shares[Symbol(newBet.outcome)] += newBet.shares
end

function redeemShares!(MarketData)
    for slug in keys(MarketData)
        if MarketData[slug].Shares[:YES] >= MarketData[slug].Shares[:NO]
            MarketData[slug].Shares[:YES] -= MarketData[slug].Shares[:NO]
            MarketData[slug].Shares[:NO] = 0.
        elseif MarketData[slug].Shares[:YES] < MarketData[slug].Shares[:NO]
            MarketData[slug].Shares[:NO] -= MarketData[slug].Shares[:YES]
            MarketData[slug].Shares[:YES] = 0.
        end
    end

    # @assert mapreduce(slug -> min(noSharesBySlug[slug], yesSharesBySlug[slug]), max, getSlugs(GROUPS)) ≈ 0
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
    # lb = -[maxBetAmount + min(MarketData[slug].Shares[:YES] / (1 - markets[slug].probability::Float64), 150.) for slug in group.slugs[bettableSlugsIndex]] #repeat([-maxBetAmount], length(bettableSlugsIndex))
    # ub = [maxBetAmount + min(MarketData[slug].Shares[:NO] / markets[slug].probability::Float64, 150.) for slug in group.slugs[bettableSlugsIndex]] #repeat([maxBetAmount], length(bettableSlugsIndex))

    
    ub = repeat([maxBetAmount], length(bettableSlugsIndex))
    lb = -ub

    problem = Optimization.OptimizationProblem(profitF, x0, lb=lb, ub=ub)

    @time "Adaptive 1" sol = solve(problem, BBO_adaptive_de_rand_1_bin_radiuslimited(), maxtime=1.5)

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

function getMarketsAndBalance!(MarketData, group, USERNAME)
    myBalance = 0.

    @sync begin 
        for slug in group.slugs
            @async try
                market = getMarketBySlug(slug)

                # if MarketData[slug].probability ≉ market.probability::Float64 || ((MarketData[slug].pool[:YES] ≉ market.pool[:YES] || MarketData[slug].pool[:NO] ≉ market.pool[:NO]) && MarketData[slug].p ≈ market.p::Float64) # || MarketData[slug].p ≉ market.p::Float64
                #     println(MarketData[slug])
                #     println(market)

                #     println(MarketData[slug].pool)
                #     println(market.pool)
                # end

                MarketData[slug].probability = market.probability::Float64
                MarketData[slug].p = market.p::Float64 # can change if a subsidy is given
                MarketData[slug].pool = market.pool

                MarketData[slug].isResolved = market.isResolved
                MarketData[slug].closeTime = market.closeTime
            catch err
                bt = catch_backtrace()
                println()
                showerror(stderr, err, bt)
                throw(err)
            end

        end

        myBalance = getUserByUsername(USERNAME).balance::Float64
    end

    return myBalance
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

    while !gotAllBets && numberOfBetsFetched <= 30
        toFetch = 3+rand(0:8)
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

function arbitrageGroup(group, BotData, MarketData, Arguments)
    rerun = :Success

    @time "Get Markets" botBalance = getMarketsAndBalance!(MarketData, group, BotData.USERNAME)

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

    maxBetAmount = botBalance / (2 + 1.5*group.noMarkets)
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

    # bindingConstraint = false # whether we would bet more if the maxBetAmount was larger

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

    if sum(abs.(betAmounts)) >= sum(bet -> bet.redeemedMana, plannedBets) + botBalance - 100
        println("Insufficient Balance $botBalance for $(sum(abs.(betAmounts))) bet redeeming $(sum(bet -> bet.redeemedMana, plannedBets)).")
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

    # rerun = bindingConstraint ? :BetMore : :Success

    if Arguments.live
        @sync try 
            for bet in plannedBets 
                @async begin
                    executedBet, ohno = execute(bet, BotData.APIKEY)

                    slug = urlToSlug(bet.url)

                    # botBalance -= executedBet.amount

                    if ohno
                        rerun = :UnexpectedBet
                    end

                    updateShares!(MarketData[slug], executedBet)

                    if !isnothing(executedBet.fills)
                        shares = 0.
                        amount = 0.
        
                        for fill in executedBet.fills
                            if isnothing(fill.matchedBetId)
                                shares += fill.shares
                                amount += fill.amount
                            end
                        end
        
                        MarketData[slug].probability = executedBet.probAfter
                        
                        for outcome in (:NO, :YES)
                            MarketData[slug].pool[outcome] += amount
                        end
        
                        MarketData[slug].pool[Symbol(executedBet.outcome)] -= shares
        
                        # MarketData[slug].isResolved = bet.probAfter
                        # MarketData[slug].closeTime = bet.probAfter
                    end

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

function arbitrage(GroupData, BotData, MarketData, lastBetId, Arguments)
    bets = getLiveBets(lastBetId)
    println(length(bets))

    seenGroups = Set{Int}()

    # for bet in Iterators.reverse(bets)
    #     if bet.contractId in GroupData.contractIdSet && bet.userUsername != BotData.USERNAME 
    #         slug = GroupData.contractIdToSlug[bet.contractId]
    #         if !isnothing(bet.isLiquidityProvision) && bet.isLiquidityProvision
    #             # MarketData[slug].pool[:YES] += bet.shares
    #             # MarketData[slug].p = bet.probAfter
    #         elseif !isnothing(bet.fills)
    #             shares = 0.
    #             amount = 0.

    #             for fill in bet.fills
    #                 if isnothing(fill.matchedBetId)
    #                     shares += fill.shares
    #                     amount += fill.amount
    #                 end
    #             end

    #             MarketData[slug].probability = bet.probAfter
                
    #             for outcome in (:NO, :YES)
    #                 MarketData[slug].pool[outcome] += amount
    #             end

    #             MarketData[slug].pool[Symbol(bet.outcome)] -= shares

    #             # MarketData[slug].isResolved = bet.probAfter
    #             # MarketData[slug].closeTime = bet.probAfter
    #         end
    #     end
    # end

    for bet in bets
        if bet.contractId in GroupData.contractIdSet && bet.userUsername != BotData.USERNAME 
            if GroupData.contractIdToGroupIndex[bet.contractId] ∉ seenGroups && (isnothing(bet.limitProb) || !(bet.probAfter ≈ bet.probBefore))
                rerun = :FirstRun
                runs = 0
                delay = 60
                
                while rerun == :FirstRun || (rerun == :BetMore && runs ≤ 5) || (rerun == :UnexpectedBet && runs ≤ 10) || rerun == :PostFailure
                    if rerun == :PostFailure
                        sleep(delay)
                        delay *= 5
                    end 
                    printstyled("Running $(GroupData.groups[GroupData.contractIdToGroupIndex[bet.contractId]].name)\n", color=:light_cyan)
                    rerun = arbitrageGroup(GroupData.groups[GroupData.contractIdToGroupIndex[bet.contractId]], BotData, MarketData, Arguments)

                    runs += 1
                end
                push!(seenGroups, GroupData.contractIdToGroupIndex[bet.contractId])

                for slug in GroupData.groups[GroupData.contractIdToGroupIndex[bet.contractId]].slugs
                    MarketData[slug].limitOrders = Dict{Symbol, Dict{Float64, Vector{Float64}}}() # need to reset as we aren't tracking limit orders
                    MarketData[slug].sortedLimitProbs = Dict(:YES=>[], :NO=>[])
                end
            end
        end
    end

    if length(bets) == 0
        return lastBetId
    else
        return bets[1].id
    end
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

function testIndividualGroup(live=false, confirmBets=true, printDebug=true)
    GROUPS::Dict{String, Dict{String, Vector{String}}}, APIKEY, USERNAME = readData()

    groups = Group.(keys(GROUPS), values(GROUPS))

    groups = [groups[1]]
    
    marketDataBySlug = Dict(slug => MarketData() for slug in getSlugs(groups))

    printstyled("Fetching at $(Dates.format(now(), "HH:MM:SS.sss"))\n"; color = :green)
    getMarkets!(marketDataBySlug, getSlugs(groups))
    USERID = getUserByUsername(USERNAME).id
    fetchMyShares!(marketDataBySlug, USERID)
    printstyled("Done fetching at $(Dates.format(now(), "HH:MM:SS.sss"))\n"; color = :green)

    botData = BotData(APIKEY, USERNAME, USERID)
    arguments = Arguments(live, confirmBets, printDebug)

    printstyled("Running at $(Dates.format(now(), "HH:MM:SS.sss"))\n"; color = :blue)
    arbitrageGroup(groups[1], botData, marketDataBySlug, arguments)
    printstyled("Done at $(Dates.format(now(), "HH:MM:SS.sss"))\n"; color = :magenta)
end

function setup(groupNames, live, confirmBets, printDebug)
    GROUPS::Dict{String, Dict{String, Vector{String}}}, APIKEY, Supabase_APIKEY, USERNAME = readData()  

    if groupNames !== nothing
        GROUPS = Dict(name => GROUPS[name] for name in groupNames)
    end

    groups = Group.(keys(GROUPS), values(GROUPS))
    
    marketDataBySlug = Dict(slug => MarketData() for slug in getSlugs(groups))

    printstyled("Fetching at $(Dates.format(now(), "HH:MM:SS.sss"))\n"; color = :green)
    lastBetId = getBets(limit=1)[1].id
    getMarkets!(marketDataBySlug, getSlugs(groups))
    USERID = getUserByUsername(USERNAME).id
    printstyled("Fetching Shares at $(Dates.format(now(), "HH:MM:SS.sss"))\n"; color = :green)
    fetchMyShares!(marketDataBySlug, USERID)
    printstyled("Done fetching at $(Dates.format(now(), "HH:MM:SS.sss"))\n"; color = :green)

    contractIdSet = Set(market.id for market in values(marketDataBySlug))
    contractIdToGroupIndex = Dict(marketDataBySlug[slug].id => i for (i, group) in enumerate(groups) for slug in group.slugs)
    contractIdToSlug = Dict(marketDataBySlug[slug].id => slug for group in groups for slug in group.slugs)

    botData = BotData(APIKEY, Supabase_APIKEY, USERNAME, USERID)
    arguments = Arguments(live, confirmBets, printDebug)
    groupData = GroupData(groups, contractIdSet, contractIdToGroupIndex, contractIdToSlug)

    return groupData, botData, marketDataBySlug, arguments, lastBetId
end

function production(groupNames = nothing; live=true, confirmBets=false, printDebug=false, skip=false)
    groupData, botData, marketDataBySlug, arguments, lastBetId = setup(groupNames, live, confirmBets, printDebug)

    if !skip
        for group in groupData.groups
            rerun = :FirstRun
            runs = 0
            
            while rerun == :FirstRun || (rerun == :BetMore && runs ≤ 5) || (rerun == :UnexpectedBet && runs ≤ 10) || rerun == :PostFailure 
                rerun = arbitrageGroup(group, botData, marketDataBySlug, arguments)

                runs += 1
            end

            redeemShares!(marketDataBySlug) # just needs to be run periodically to prevent overflow
        end
    end

    while true
        lastRunTime = time()
        printstyled("Running at $(Dates.format(now(), "HH:MM:SS.sss"))\n"; color = :blue)

        lastBetId = arbitrage(groupData, botData, marketDataBySlug, lastBetId, arguments)
        redeemShares!(marketDataBySlug) # just needs to be run periodically to prevent overflow

        printstyled("Sleeping at $(Dates.format(now(), "HH:MM:SS.sss"))\n"; color = :magenta)

        # sleep(15)
        sleep(max(2.5, 7.5 - (time() - lastRunTime)) + rand())
        # sleep(60 + 2*(rand()-.5) * 5) # - (time() - oldTime) # add some randomness so it can't be exploited based on predicability of betting time.
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

# @printlog "log.txt"
# production()