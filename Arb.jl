using ManifoldMarkets, TOML, Optimization, OptimizationBBO, Dates

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
    print(io, "$(plannedBet.market.question)\n  Buy $(plannedBet.shares) $(plannedBet.outcome) shares for $(plannedBet.amount)")
end

function execute(bet, APIKEY)
    createBet(APIKEY, bet.market.id, bet.amount, bet.outcome)
end

function f(betAmount, group, markets, limitOrdersBySlug, noSharesBySlug, yesSharesBySlug, bettableSlugsIndex)
    newProb = [markets[slug].probability for slug in group.slugs]
    noShares = zeros(group.noMarkets)
    yesShares = zeros(group.noMarkets)

    fees = 0.

    for (i, j) in enumerate(bettableSlugsIndex)
        slug = group.slugs[j]
        market = markets[slug]

        shares = betToShares(market, limitOrdersBySlug[slug], betAmount[i]).shares
        newProb[j] = betToShares(market, limitOrdersBySlug[slug], betAmount[i]).probability

        if betAmount[i] >= 1
            yesShares[j] += shares
            fees += 0.1
        elseif betAmount[i] <= -1
            noShares[j] += shares
            fees += 0.1
        else
            betAmount[i] = 0
        end
    end

    profitsByEvent = group.y_matrix * (yesShares + getindex.(Ref(yesSharesBySlug), group.slugs)) + group.n_matrix * (noShares + getindex.(Ref(noSharesBySlug), group.slugs)) .- sum(abs.(betAmount)) .- fees

    return (profitsByEvent=profitsByEvent, noShares=noShares, yesShares=yesShares, newProbability=newProb)
end

function optimise(group, markets, limitOrdersBySlug, maxBetAmount, noSharesBySlug, yesSharesBySlug, bettableSlugsIndex)    
    profitF = OptimizationFunction((betAmount, _) -> -minimum( f(betAmount, group, markets, limitOrdersBySlug, noSharesBySlug, yesSharesBySlug, bettableSlugsIndex).profitsByEvent ))

    x0 = repeat([0.], length(bettableSlugsIndex))
    lb = repeat([-maxBetAmount], length(bettableSlugsIndex))
    ub = repeat([maxBetAmount], length(bettableSlugsIndex))

    problem = Optimization.OptimizationProblem(profitF, x0, lb=lb, ub=ub)

    sol = solve(problem, BBO_adaptive_de_rand_1_bin(), maxtime=30.0)

    return sol
end

function getMarkets(slugs)
    markets = Dict{String, Market}()
    @sync for slug in slugs
        @async markets[slug] = getMarketBySlug(slug)
    end

    return markets
end

function getMarketsAndBets(GROUPS, USERNAME)
    markets = Dict{String, Market}()
    betsBySlug = Dict{String, Vector{Bet}}()
    betsByMe = Dict{String, Vector{Bet}}()

    @sync for group in values(GROUPS), url in keys(group)
        slug = urlToSlug(url)

        @async markets[slug] = getMarketBySlug(slug)
        @async betsBySlug[slug] = getBets(slug=slug)
        @async betsByMe[slug] = getBets(slug=slug, username=USERNAME)
    end

    return (markets=markets, betsBySlug=betsBySlug, betsByMe=betsByMe)
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
            pop!(GROUPS, name)
        end
    end
end

function getLimits(groups, betsBySlug)
    userBalance = Dict{String, Float64}()
    limitOrdersBySlug = Dict()
    splitBets = Dict()

    for group in groups, slug in group.slugs
        limitOrdersByProb, userBalance = sortLimitOrders(betsBySlug[slug], userBalance)
        limitOrdersBySlug[slug], splitBets[slug] = getLimitOrders(limitOrdersByProb, userBalance)
    end

    return (limitOrdersBySlug=limitOrdersBySlug, splitBets=splitBets)
end

function calculateShares(groups, betsByMe)
    netShares = Dict(getSlugs(groups) .=> 0.)
    yesShares = Dict(getSlugs(groups) .=> 0.)
    noShares = Dict(getSlugs(groups) .=> 0.)

    for group in groups, slug in group.slugs
        for bet in betsByMe[slug]
            netShares[slug] += bet.shares * (bet.outcome == "YES" ? 1. : -1.)
        end

        if netShares[slug] > 0.
            yesShares[slug] = netShares[slug]
        else
            noShares[slug] = -netShares[slug]
        end
    end
    return (noShares=noShares, yesShares=yesShares)
end


function arbitrage(GROUPS, APIKEY, USERNAME, live=false, confirmBets=true, printDebug=true)
    markets, betsBySlug, betsByMe = getMarketsAndBets(GROUPS, USERNAME)
    
    processGroups!(GROUPS, markets)

    groups = Group.(keys(GROUPS), values(GROUPS))

    limitOrdersBySlug, splitBets = getLimits(groups, betsBySlug)

    totalNumberOfMarkets = mapreduce(group -> group.noMarkets, +, groups)
    maxBetAmount = getUserByUsername(USERNAME).balance / (3 + totalNumberOfMarkets)

    noSharesBySlug, yesSharesBySlug = calculateShares(groups, betsByMe)

    # sols = map(group -> optimise(group, markets, limitOrdersBySlug, maxBetAmount, noSharesBySlug, yesSharesBySlug, groups))

    for group in groups
        plannedBets = PlannedBet[]
        
        println("=== $(group.name) ===")

        # for market in markets:
        #     skip = skip_market(market)
        #     if skip:
        #         print()
        #         print(skip)
        #         print("Skipping group.\n")
        #         return

        bettableSlugsIndex = [i for (i, slug) in enumerate(group.slugs) if !isMarketClosingSoon(markets[slug])]

        solution = optimise(group, markets, limitOrdersBySlug, maxBetAmount, noSharesBySlug, yesSharesBySlug, bettableSlugsIndex)
        betAmounts = solution.u
        newProfitsByEvent, noShares, yesShares, newProb = f(betAmounts, group, markets, limitOrdersBySlug, noSharesBySlug, yesSharesBySlug, bettableSlugsIndex)

        oldProfitsByEvent, _, _, oldProb = f(repeat([0.], group.noMarkets), group, markets, limitOrdersBySlug, noSharesBySlug, yesSharesBySlug, bettableSlugsIndex)

        profit = minimum(newProfitsByEvent) - minimum(oldProfitsByEvent)

        if printDebug
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
                bet = PlannedBet(abs(amount), yesShares[i] + noShares[i], amount > 0. ? "YES" : "NO", markets[slug])
                push!(plannedBets, bet)
            else
                println("Bet amount: $(amount) is too small, $(slug), $(markets[slug].question)")
                skipMarket = true
                break
            end

            println(slug)
            println("Prior probs:    ", markets[slug].probability)
            println("Posterior probs:", newProb[i])
            println()
        end

        if skipMarket
            continue
        end

        # println()
        # println("Profits:", profit)

                
        if profit <= .01 * length(plannedBets)
            # print()
            println("Insufficient profit $(profit) for $(length(plannedBets)) bets\n")
            continue
        end
        println("Profits: $profit")

        for bet in plannedBets
            println(bet)
        end

        if confirmBets
            println("Proceed? (y/n)") 
            if readline() !="y"
                continue
            end
        end

        # # Make sure markets haven't moved
        # # if any(m.bets[0].createdTime != mf.get_bets(market=m.slug, limit=1)[0].createdTime for m in markets):
        # oldPrice = [markets[slug].probability for slug in group.slugs]
        # newPrice = [client.get_market_by_slug(slug).probability for slug in group.slugs]
        # if not allclose(oldPrice, newPrice) # this won't trigger if limit order gets filled, so need to check if any new bets instead.
        #     print("Markets have moved!\n")
        #     continue
        # end

        if live
            for bet in plannedBets
                execute(bet, APIKEY)
            # print(f"Balance: {my_balance()}\n")
            end
        end
    end
end

function run(groupNames = nothing ;live=false, confirmBets=true, printDebug=true)
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

    arbitrage(GROUPS, APIKEY, USERNAME, live, confirmBets, printDebug)
end

# run()
# run(live=true)
# run(["Next Speaker of the House"])

function prod(groupNames = nothing;live=true, confirmBets=false, printDebug=false)
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

    while true
        # oldTime = time()
        println("Running at $(Dates.format(now(), "HH:MM:SS.sss"))")
        arbitrage(GROUPS, APIKEY, USERNAME, live, confirmBets, printDebug)
        println("Sleeping at $(Dates.format(now(), "HH:MM:SS.sss"))")
        sleep(1.5*60) # - (time() - oldTime)
    end
end