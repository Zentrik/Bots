using Revise, ManifoldMarkets, TOML, Optimization, OptimizationBBO

#%%
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

function f(group, markets, limitOrdersBySlug, currentNoShares, currentYesShares, betAmount)
    newProb = zeros(group.noMarkets)
    noShares = zeros(group.noMarkets)
    yesShares = zeros(group.noMarkets)

    fees = 0.

    for (i, slug) in enumerate(group.slugs)
        market = markets[slug]
        shares = betToShares(market, limitOrdersBySlug[slug], betAmount[i]).shares
        newProb[i] = betToShares(market, limitOrdersBySlug[slug], betAmount[i]).probability

        if betAmount[i] >= 1
            yesShares[i] += shares
            fees += 0.1
        elseif betAmount[i] <= -1
            noShares[i] += shares
            fees += 0.1
        else
            betAmount[i] = 0
        end
    end

    profitsByEvent = group.y_matrix * (yesShares + currentYesShares) + group.n_matrix * (noShares + currentNoShares) .- sum(abs.(betAmount)) .- fees

    return (profitsByEvent=profitsByEvent, yesShares=yesShares, noShares=noShares, newProbability=newProb)
end

function optimise(group, markets, limitOrdersBySlug, maxBetAmount, noShares, yesShares)    
    profitF = OptimizationFunction((betAmount, _) -> -minimum( f(group, markets, limitOrdersBySlug, noShares, yesShares, betAmount).profitsByEvent ))
    x0 = repeat([0], group.noMarkets)
    lb = repeat([-maxBetAmount], group.noMarkets)
    ub = repeat([maxBetAmount], group.noMarkets)
    problem = Optimization.OptimizationProblem(profitF, x0, lb=lb, ub=ub, progress=true)

    sol = solve(problem, BBO_adaptive_de_rand_1_bin_radiuslimited(), maxiters=100000, maxtime=20.0)

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

function processGroups!(GROUPS, markets)
    for (name, group) in GROUPS
        for url in keys(group)
            slug = urlToSlug(url)

            market = markets[slug]
            if market.isResolved || market.closeTime / 1000 < time() + 60 # if resolved or closing in 60 seconds, resolved check should be redundant
                pop!(GROUPS[name], url)
            end
        end
    end

    # Remove any groups that have no/one markets left

    for (name, group) in GROUPS
        if length(group) <= 1
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

    sols = map(group -> optimise(group, markets, limitOrdersBySlug, maxBetAmount, getindex.(Ref(noSharesBySlug), group.slugs), getindex.(Ref(yesSharesBySlug), group.slugs)), groups)

    for (sol, group) in zip(sols, groups)
        plannedBets = PlannedBet[]
        
        println("=== $(group.name) ===")

        # for market in markets:
        #     skip = skip_market(market)
        #     if skip:
        #         print()
        #         print(skip)
        #         print("Skipping group.\n")
        #         return

        betAmounts = sol.u
        newProfitsByEvent, yesShares, noShares, newProb = f(group, markets, limitOrdersBySlug, getindex.(Ref(noSharesBySlug), group.slugs), getindex.(Ref(yesSharesBySlug), group.slugs), betAmounts)

        oldProfitsByEvent = f(group, markets, limitOrdersBySlug, getindex.(Ref(noSharesBySlug), group.slugs), getindex.(Ref(yesSharesBySlug), group.slugs), repeat([0.], group.noMarkets)).profitsByEvent
        profit = minimum(newProfitsByEvent) - minimum(oldProfitsByEvent)

        if printDebug
            println(newProb)
            println(oldProfitsByEvent)
            println(newProfitsByEvent)
            println(profit)
            println()
            println(betAmounts)
            println()
            println(group.y_matrix)
            println(group.n_matrix)
            println()
            println(group.y_matrix * yesShares)
            println(group.n_matrix * noShares)
            println(sum(abs.(betAmounts)))
        end
        
        if profit <= .01 * group.noMarkets
            print()
            println("Insufficient profit $(profit) for $(group.noMarkets) markets\n")
            continue
        end

        for (i, (amount, slug)) in enumerate(zip(betAmounts, group.slugs))
            if isapprox(amount, 0., atol=1e-6)
                continue
            elseif abs(amount) >= 1.
                bet = PlannedBet(abs(amount), yesShares[i] + noShares[i], amount > 0. ? "YES" : "NO", markets[slug])
                push!(plannedBets, bet)
            else
                println("Bet amount: $(amount) is too small, $(slug), $(markets[slug].question)")
            end

            println()
            println(slug)
            println("Prior probs:    ", markets[slug].probability)
            println("Posterior probs:", newProb[i])
        end

        println()
        println("Profits:", profit)

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

function run(;live=false, confirmBets=true, printDebug=true)
    data = TOML.parsefile("ArbBot/Arb.toml")
    GROUPS = data["GROUPS"]
    APIKEY = data["APIKEY"]
    USERNAME = data["USERNAME"]
    
    slugs = getSlugs(GROUPS)
    markets = getMarkets(slugs)
    processGroups!(GROUPS, markets)

    arbitrage(GROUPS, APIKEY, USERNAME, live, confirmBets, printDebug)
end

run()
# run(live=true)