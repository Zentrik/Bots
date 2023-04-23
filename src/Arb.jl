using ManifoldMarkets, TOML, Optimization, OptimizationBBO, Dates, Combinatorics, LinearAlgebra, Parameters
using HTTP, JSON3, Dates, OpenSSL
using HTTP.WebSockets
using SmartAsserts, Logging, LoggingExtras

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

macro smart_assert_showerr(ex, msg=nothing)
    esc(quote
        try
            @smart_assert eval($ex) eval($msg)
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

const FEE = 0.03
Base.exit_on_sigint(false)

struct Group
    name::String
    slugs::Vector{String}
    y_matrix::Matrix{Float32}
    n_matrix::Matrix{Float32}
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

    hasUpdated::typeof(Condition()) = Condition()
    lastOptimisedProb::Float64 = -1
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

function execute(bet, currentProb, APIKEY)
    ohno = false
    response = createBet(APIKEY, bet.id, bet.amount, bet.outcome)
    # need to check if returned info matches what we wanted to bet, i.e. if we got less shares than we wanted to. If we got more ig either moved or smth weird with limit orders.

    @smart_assert_showerr bet.outcome == response.outcome

    if response.shares ≉ bet.shares
        @error "\e]8;;$(bet.url)\e\\$(bet.question)\e]8;;\e\\\n" # hyperlink
        @error response
        @error response.fills

        ohno = true

        @smart_assert_showerr !(length(response.fills) == 1 && isnothing(response.fills[end].matchedBetId) && response.probBefore ≈ currentProb) "$currentProb"
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

    @debug "Running adaptive for $(group.name)"
    @time "Adaptive 1" sol = solve(problem, BBO_adaptive_de_rand_1_bin_radiuslimited(), maxtime=3)

    @debug "Yielding after adaptive, $(group.name)"
    yield()
    @debug "Done yielding after adaptive, $(group.name)"

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
        @debug "Running resampling for $(group.name)"

        @time "Resampling" sol2 = solve(problem, BBO_resampling_memetic_search(), maxtime=3)
        
        @debug "Yielding after resampling, $(group.name)"
        yield()
        @debug "Done yielding after resampling, $(group.name)"

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

function updateMarketData!(MarketData, market)
    # MarketData.probability = market.prob  # this updates slowly?

    MarketData.p = market.p # can change if a subsidy is given
    MarketData.pool = market.pool

    MarketData.probability = poolToProb(MarketData.p, MarketData.pool)

    MarketData.id = market.id
    MarketData.question = market.question
    MarketData.url = "https://manifold.markets/$(market.creatorUsername)/$(market.slug)"

    MarketData.isResolved = market.isResolved
    MarketData.closeTime = market.closeTime
end

function getMarkets!(MarketData, slugs)
    @sync for paritionedSlugs in Iterators.partition(slugs, 20) # 2000 is max characters in uri, so we have about 1900 for slugs, max slug is about 50 : 1900/50 is 38. generous magin to 20
        @async_showerr begin
            contracts_slugs = join(paritionedSlugs, ",")
    
            response = HTTP.get("https://pxidrgkatumlvfqaxcll.supabase.co/rest/v1/contracts?slug=in.($contracts_slugs)", headers= ["apikey" => "eyJhbGciOiJIUzI1NiIsInR5cCI6IkpXVCJ9.eyJpc3MiOiJzdXBhYmFzZSIsInJlZiI6InB4aWRyZ2thdHVtbHZmcWF4Y2xsIiwicm9sZSI6ImFub24iLCJpYXQiOjE2Njg5OTUzOTgsImV4cCI6MTk4NDU3MTM5OH0.d_yYtASLzAoIIGdXUBIgRAGLBnNow7JG2SoaNMQ8ySg", "Content-Type" => "application/json"])
            responseJSON = JSON3.read(response.body)
    
            for contract in responseJSON
                updateMarketData!(MarketData[contract.slug], contract.data)
            end
        end
    end
end

function getMarketsUsingId!(MarketData, slugs) # If we use slugs the request uri is too long
    try
        for slug in slugs
            @smart_assert_showerr !isnothing(MarketData[slug].id)
        end
        
        contracts_ids = join(map(slug -> MarketData[slug].id, slugs), ",") # might be faster to index by id
        # contracts_slugs = join(slugs, ",")

        @smart_assert_showerr length(slugs) < 1000 "too many slugs $(length(slugs)), $slugs when fetching markets"

        response = HTTP.get("https://pxidrgkatumlvfqaxcll.supabase.co/rest/v1/contracts?id=in.($contracts_ids)", headers= ["apikey" => "eyJhbGciOiJIUzI1NiIsInR5cCI6IkpXVCJ9.eyJpc3MiOiJzdXBhYmFzZSIsInJlZiI6InB4aWRyZ2thdHVtbHZmcWF4Y2xsIiwicm9sZSI6ImFub24iLCJpYXQiOjE2Njg5OTUzOTgsImV4cCI6MTk4NDU3MTM5OH0.d_yYtASLzAoIIGdXUBIgRAGLBnNow7JG2SoaNMQ8ySg", "Content-Type" => "application/json"])
        responseJSON = JSON3.read(response.body)

        for contract in responseJSON
            updateMarketData!(MarketData[contract.slug], contract.data)
        end
    catch err
        bt = catch_backtrace()
        println()
        showerror(stderr, err, bt)
        throw(err)
    end
end

getSlugs(GROUPS::Dict) = mapreduce(x -> urlToSlug.(x), vcat, keys.(values(GROUPS)))
getSlugs(groups::Vector{Group}) = mapreduce(group -> group.slugs, vcat, groups)

isMarketClosingSoon(market) = market.isResolved || market.closeTime / 1000 < time() + 60 # if resolved or closing in 60 seconds

function arbitrageGroup(group, BotData, MarketData, Arguments)
    rerun = :Success

    for slug in group.slugs
        @smart_assert_showerr MarketData[slug].probability ≈ poolToProb(MarketData[slug].p, MarketData[slug].pool) "$slug, $(MarketData[slug]), $(MarketData[slug].pool), $(poolToProb(MarketData[slug].p, MarketData[slug].pool))"
    end

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

    maxBetAmount = BotData.balance / (2 + 1.5*group.noMarkets)
    # redeemManaHack = maximum(slug -> max(MarketData[slug].Shares[:YES] / (1 - MarketData[slug].probability), MarketData[slug].Shares[:NO] / MarketData[slug].probability), group.slugs[bettableSlugsIndex])
    # maxBetAmount += min(redeemManaHack, maxBetAmount, 100.)
    # Remove hack for now as it seems to prevent optmiser betting as much as it can which prevents reruns

    betAmounts = optimise(group, MarketData, maxBetAmount, bettableSlugsIndex)

    oldProb = [MarketData[slug].probability for slug in group.slugs]
    oldNoShares = [MarketData[slug].Shares[:NO] for slug in group.slugs]
    oldYesShares = [MarketData[slug].Shares[:YES] for slug in group.slugs]

    newProfitsByEvent, noShares, yesShares, newProb = f(betAmounts, group, MarketData, oldNoShares, oldYesShares, bettableSlugsIndex)

    oldProfitsByEvent, _, _, _ = f(repeat([0.], length(bettableSlugsIndex)), group, MarketData, oldNoShares, oldYesShares, bettableSlugsIndex)

    profit = minimum(newProfitsByEvent) - minimum(oldProfitsByEvent)
    newYesShares = yesShares .- oldYesShares
    newNoShares = noShares  .- oldNoShares

    @debug bettableSlugsIndex
    @debug oldProb
    @debug newProb
    @debug oldProfitsByEvent
    @debug newProfitsByEvent
    @debug profit
    @debug betAmounts
    @debug group.y_matrix
    @debug group.n_matrix
    @debug map(slug -> MarketData[slug].Shares[:YES], group.slugs)
    @debug map(slug -> MarketData[slug].Shares[:NO], group.slugs)
    @debug yesShares
    @debug noShares
    @debug group.y_matrix * yesShares
    @debug group.n_matrix * noShares
    @debug sum(abs.(betAmounts))
    @debug [MarketData[slug].p for slug in group.slugs]
    @debug [MarketData[slug].pool for slug in group.slugs]

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

            @smart_assert_showerr (sign(amount) == sign(newProb[j] - oldProb[j]) || !isempty(MarketData[slug].sortedLimitProbs[Symbol(outcome)])) "$slug, $amount, $(newProb[j]), $(oldProb[j]), $i, $j"
            
            bet = PlannedBet(abs(amount), shares, outcome, redeemedMana, MarketData[slug].id, MarketData[slug].url, MarketData[slug].question)
            push!(plannedBets, bet)
        else
            # rerun = :BetMore
            @warn "Bet amount: $(amount) is too small, $(slug)"

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
            @warn "Bet size is $(100 * abs(amount)/maxBetAmount)% of maxBetAmount"
            # bindingConstraint = true
            rerun = :BetMore
        end

        if abs(amount) ≈ 1
            @warn "Bet size is $(abs(amount))"
            # bindingConstraint = true
            rerun = :BetMore
        end

        @info "\e]8;;$(MarketData[slug].url)\e\\$(MarketData[slug].question)\e]8;;\e\\" # hyperlink
        @info "Prior probs:     $(MarketData[slug].probability * 100)%"
        @info "Posterior probs: $(newProbBySlug[slug]*100)%"
        @info "Buy $(bet.shares) $(bet.outcome) shares for $(bet.amount), redeeming $(bet.redeemedMana)"
    end

    if sum(abs.(betAmounts)) >= sum(bet -> bet.redeemedMana, plannedBets) + BotData.balance - 100
        @error "Insufficient Balance $(BotData.balance) for $(sum(abs.(betAmounts))) bet redeeming $(sum(bet -> bet.redeemedMana, plannedBets))."
        rerun = :InsufficientBalance
        return rerun
    end

    @info "Profits:         $(profit + FEE * length(plannedBets))\n" # no more fee, but we still want to use fee in optimisation

    if Arguments.confirmBets
        println("Proceed? (y/n)") 
        if readline() !="y"
            rerun = :Success
            return rerun
        end
    end

    if Arguments.live
        oldProbBySlug = Dict(group.slugs[j] => oldProb[j] for j in bettableSlugsIndex)

        # I think we are able to fetch messages while we're making bets
        @sync try 
            for bet in plannedBets 
                slug = urlToSlug(bet.url)

                @async_showerr begin
                    # We use oldProb instead of MarketData as MarketData may have updated to a new probability while we bet.
                    executedBet, ohno = execute(bet, oldProbBySlug[slug], BotData.APIKEY)

                    BotData.balance -= executedBet.amount

                    if ohno
                        rerun = :UnexpectedBet
                    end

                    updateShares!(MarketData[slug], executedBet, BotData)

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

                # If data comes in before POST response occurs we don't want to be stuck waiting
                @async_showerr begin
                    @debug "$(Dates.format(now(), "HH:MM:SS.sss")) Waiting begun for $slug"
                    @debug wait(MarketData[slug].hasUpdated)
                    @debug "$(Dates.format(now(), "HH:MM:SS.sss")) Waiting done for $slug"
                end
            end
        catch err
            bt = catch_backtrace()
            println()
            showerror(stderr, err, bt)

            # What if a bet didn't go througth, thus we must fetch
            fetchMyShares!(marketDataBySlug, BotData.USERID)

            rerun = :PostFailure
        end

        
        for (i, slug) in enumerate(group.slugs)
            MarketData[slug].lastOptimisedProb = newProb[i] # not great if we have insufficient balance as later we might have sufficient but think the market is already optimised
        end

        return rerun
    else
        return :Success
    end
end

function fetchMyShares!(MarketDataBySlug, groupData, USERID)
    try
        # https://discourse.julialang.org/t/broadcast-object-property/47104/6
        contracts = join(Base.broadcasted(getproperty, values(MarketDataBySlug), :id), ",")

        @smart_assert_showerr length(MarketDataBySlug) < 1000

        # will fail if we have more than 1000
        response = HTTP.get("https://pxidrgkatumlvfqaxcll.supabase.co/rest/v1/user_contract_metrics?user_id=eq.$USERID&contract_id=in.($contracts)", headers= ["apikey" => "eyJhbGciOiJIUzI1NiIsInR5cCI6IkpXVCJ9.eyJpc3MiOiJzdXBhYmFzZSIsInJlZiI6InB4aWRyZ2thdHVtbHZmcWF4Y2xsIiwicm9sZSI6ImFub24iLCJpYXQiOjE2Njg5OTUzOTgsImV4cCI6MTk4NDU3MTM5OH0.d_yYtASLzAoIIGdXUBIgRAGLBnNow7JG2SoaNMQ8ySg", "Content-Type" => "application/json"])
        responseJSON = JSON3.read(response.body)

        # what if it returns nothing for some contract? presumably that means we never invested so MarketData shoudl already be correct

        for contract_metrics in responseJSON
            slug = groupData.contractIdToSlug[contract_metrics.contract_id]
            for (outcome, shares) in contract_metrics.data.totalShares
                MarketDataBySlug[slug].Shares[outcome] = shares
            end
        end
    catch err
        bt = catch_backtrace()
        println()
        showerror(stderr, err, bt)
        throw(err)
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

function setup(groupNames, live, confirmBets)
    arguments = Arguments(live, confirmBets)

    GROUPS::Dict{String, Dict{String, Vector{String}}}, APIKEY, Supabase_APIKEY, USERNAME = readData()  

    if groupNames !== nothing
        GROUPS = Dict(name => GROUPS[name] for name in groupNames)
    end

    groups = Group.(keys(GROUPS), values(GROUPS))
    
    marketDataBySlug = Dict(slug => MarketData() for slug in getSlugs(groups))

    printstyled("Fetching at $(Dates.format(now(), "HH:MM:SS.sss"))\n"; color = :green)
    @time "Fetching User" botUser = getUserByUsername(USERNAME)
    USERID = botUser.id
    botBalance = botUser.balance
    botData = BotData(APIKEY, Supabase_APIKEY, USERNAME, USERID, botBalance)

    @time "Fetching Markets" getMarkets!(marketDataBySlug, getSlugs(groups))
    contractIdSet = Set(market.id for market in values(marketDataBySlug))
    contractIdToGroupIndex = Dict(marketDataBySlug[slug].id => i for (i, group) in enumerate(groups) for slug in group.slugs)
    contractIdToSlug = Dict(marketDataBySlug[slug].id => slug for group in groups for slug in group.slugs)
    groupData = GroupData(groups, contractIdSet, contractIdToGroupIndex, contractIdToSlug)

    @time "Fetching Shares" fetchMyShares!(marketDataBySlug, groupData, USERID)


    return groupData, botData, marketDataBySlug, arguments
end

function production(groupNames = nothing; live=true, confirmBets=false, skip=false)
    # Should probably move this into the websocket as well
    groupData, botData, marketDataBySlug, arguments = setup(groupNames, live, confirmBets)

    TaskDict = Dict(i => (runAgain=false, task=@async nothing) for i in eachindex(groupData.groups)) # so we don't have to check if there is a task in it or not. Order of tuple matters for parsing

    # currentTask = @async nothing
    WebSockets.open(uri(botData.Supabase_APIKEY)) do socket
        try
            @info "Opened Socket"
            @debug socket
            send(socket, pushJSONContracts())
            # send(socket, pushJSON("contract_bets"))
            @info "Sent Intialisation"
    
            @sync try
                # Reads messages whilst in this loop but blocks arbing due to them, also ensures MarketData updates after we make a bet
                currentTask = @async_showerr if !skip
                    @warn "All: Running all groups at $(Dates.format(now(), "HH:MM:SS.sss"))"
                    for group in groupData.groups            
                        rerun = :FirstRun
                        runs = 0
                        
                        while rerun == :FirstRun || (rerun == :BetMore && runs ≤ 5) || (rerun == :UnexpectedBet && runs ≤ 10) || rerun == :PostFailure
                            @warn "All: Running $(group.name) at $(Dates.format(now(), "HH:MM:SS.sss"))"

                            # Need to wait for current bet to finish before rerunning
                            # this isn't async so surely not a problem
                            rerun = arbitrageGroup(group, botData, marketDataBySlug, arguments)
            
                            runs += 1
                        end
                    end

                    @warn "All: Done all groups at $(Dates.format(now(), "HH:MM:SS.sss"))"
                end
                #Reading messages
                # should switch to @spawn :interactive
                @async_showerr for msg in socket
                    @async_showerr begin 
                        msgJSON = JSON3.read(msg)
                        if !(:payload in keys(msgJSON) && :data in keys(msgJSON.payload))
                            @debug "$(Dates.format(now(), "HH:MM:SS.sss")): $msg"
                            return nothing
                        end

                        market = msgJSON.payload.data.record.data
                        marketId = market.id
                        # @debug market.slug
                    
                        if marketId in groupData.contractIdSet 
                            @smart_assert_showerr market.mechanism == "cpmm-1"

                            slug = groupData.contractIdToSlug[marketId] #market.slug also works
                            oldProb = marketDataBySlug[slug].probability

                            updateMarketData!(marketDataBySlug[slug], market)

                            # println("Received message: $(market.lastBetTime)")
                            # @warn slug
                            if marketDataBySlug[slug].probability ≉ oldProb
                                @debug "$slug at $(marketDataBySlug[slug].probability) at $(Dates.format(now(), "HH:MM:SS.sss"))"
                            # @warn "from prob $(market.prob)"
                            end

                            # Should probably only do this if market data actually changed
                            @debug "$(Dates.format(now(), "HH:MM:SS.sss")) Notified $slug at $(marketDataBySlug[slug].probability)"
                            notify(marketDataBySlug[slug].hasUpdated)

                            # oldProb is to check if market moved?
                            # lastOptimisedProb is to prevent rerunning on already optimised market, accounts for getting new market data due to our own bet
                            if marketDataBySlug[slug].probability ≉ oldProb && marketDataBySlug[slug].probability ≉ marketDataBySlug[slug].lastOptimisedProb && !TaskDict[groupData.contractIdToGroupIndex[marketId]].runAgain
                                
                                TaskDict[groupData.contractIdToGroupIndex[marketId]] = (runAgain = true, task=TaskDict[groupData.contractIdToGroupIndex[marketId]].task)

                                while !istaskdone(currentTask) # say we have 3 tasks A, B, C. A is running when B, C come in so both wait for A to finish. The B runs sets currentTask to itself but C is onlt waiting for A not the new currentTask B so need a while loop.
                                    wait(currentTask)
                                    @debug "Finished waiting for current task, $currentTask, $slug, $(marketDataBySlug[slug].probability)"
                                end

                                while !istaskdone(TaskDict[groupData.contractIdToGroupIndex[marketId]].task)
                                    wait(TaskDict[groupData.contractIdToGroupIndex[marketId]].task)
                                    @debug "Finished waiting for previous run of market, $(TaskDict[groupData.contractIdToGroupIndex[marketId]].task), $slug, $(marketDataBySlug[slug].probability)"
                                end

                                # if we've already optimised on this market data
                                if !TaskDict[groupData.contractIdToGroupIndex[marketId]].runAgain
                                    @debug "No need to rerun $slug"
                                    return nothing
                                end
                                
                                currentTask = current_task()
                                @debug "current task set to $currentTask, $slug"
                                TaskDict[groupData.contractIdToGroupIndex[marketId]] = (runAgain = false, task=current_task())

                                group = groupData.groups[groupData.contractIdToGroupIndex[marketId]]

                                if marketDataBySlug[slug].probability ≈ marketDataBySlug[slug].lastOptimisedProb
                                    @debug "market already optimised $slug"
                                    return nothing
                                end

                                delay = 60
                                runs = 0
                                rerun = :FirstRun

                                while rerun == :FirstRun || (rerun == :BetMore && runs ≤ 5) || (rerun == :UnexpectedBet && runs ≤ 10) || rerun == :PostFailure
                                    if rerun == :PostFailure
                                        sleep(delay)
                                        delay *= 5
                                    end 

                                    # no need to run any previously queued requests.
                                    if TaskDict[groupData.contractIdToGroupIndex[marketId]].runAgain
                                        TaskDict[groupData.contractIdToGroupIndex[marketId]] = (runAgain = false, task=current_task())
                                    end
            
                                    @warn "Running $(group.name) at $(Dates.format(now(), "HH:MM:SS.sss"))"
                                    @debug "current prob $(marketDataBySlug[slug].probability)"
                                    @debug [marketDataBySlug[slug].p for slug in group.slugs]
                                    @debug [marketDataBySlug[slug].pool for slug in group.slugs]
                                    runs += 1
                                    rerun = arbitrageGroup(group, botData, marketDataBySlug, arguments)
                                    @debug rerun
                                end
            
                                for slug in groupData.groups[groupData.contractIdToGroupIndex[marketId]].slugs
                                    marketDataBySlug[slug].limitOrders = Dict{Symbol, Dict{Float64, Vector{Float64}}}() # need to reset as we aren't tracking limit orders
                                    marketDataBySlug[slug].sortedLimitProbs = Dict(:YES=>[], :NO=>[])
                                end
                            end
                        end

                        return nothing
                    end
                end
    
                # HeartBeat
                @async_showerr while !WebSockets.isclosed(socket)
                    # printstyled("Heartbeat at $(Dates.format(now(), "HH:MM:SS.sss"))\n"; color = :blue)
                    send(socket, heartbeatJSON)
                    # printstyled("Sleeping at $(Dates.format(now(), "HH:MM:SS.sss"))\n"; color = :blue)
                    sleep(30)
                end

                # Fetch new balance at 8am
                @async_showerr while !WebSockets.isclosed(socket)
                    # printstyled("Sleeping Balance at $(Dates.format(now(), "HH:MM:SS.sss"))\n"; color = :blue)
                    timeTo8 = Second(((today() + Time(8) + Minute(5) - now()) ÷ 1000).value)
                    timeTo8 += timeTo8 > Second(0) ? Second(0) : Second(Day(1))
                    sleep(timeTo8)

                    printstyled("1: Fetching Balance at $(Dates.format(now(), "HH:MM:SS.sss"))\n"; color = :blue)
                    botData.balance = getUserByUsername(botData.USERNAME).balance

                    # printstyled("Sleeping Balance at $(Dates.format(now(), "HH:MM:SS.sss"))\n"; color = :blue)
                    timeTo8 = Second(((today() + Time(8) + Minute(10) - now()) ÷ 1000).value)
                    timeTo8 += timeTo8 > Second(0) ? Second(0) : Second(Day(1))
                    sleep(timeTo8)

                    if !WebSockets.isclosed(socket)
                        printstyled("2: Fetching Balance at $(Dates.format(now(), "HH:MM:SS.sss"))\n"; color = :blue)
                        botData.balance = getUserByUsername(botData.USERNAME).balance
                    else
                        break
                    end
                end
            catch err
                bt = catch_backtrace()
                println()
                showerror(stderr, err, bt)
                rethrow(err)
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
    end
end

function test(groupNames = nothing; live=false, confirmBets=true, skip=false) 
    production(groupNames; live=live, confirmBets=confirmBets, skip=skip)
end

function retryProd(groupNames = nothing; live=true, confirmBets=false, skip=false)
    delay = 60
    lastRunTime = time()

    while true
        try
            production(groupNames; live=live, confirmBets=confirmBets, skip=skip)
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

function testLogging()
    @info("You won't see this")
    @warn("won't see this either")
    @error("You will only see this")
    @debug "test"
    @smart_assert_showerr 1 != 1
end

# Need to run manually
using Logging, LoggingExtras, Dates
timestamp_logger(logger) = TransformerLogger(logger) do log
    merge(log, (; message = "$(Dates.format(now(), "yyyy-mm-dd HH:MM:SS.sss")) $(log.message)"))
end
const DIR = "Bots\\src\\logs"
global_logger(TeeLogger(
    EarlyFilteredLogger(Bots.not_Bots_message_filter, ConsoleLogger(stderr, Logging.Info)), 
    EarlyFilteredLogger(Bots.not_Bots_message_filter, MinLevelLogger(FileLogger("$DIR\\$(Dates.format(now(), "YYYY-mm-dd HH-MM")).log"), Logging.Info)),
    EarlyFilteredLogger(Bots.not_Bots_message_filter, MinLevelLogger(FileLogger("$DIR\\Debug-$(Dates.format(now(), "YYYY-mm-dd HH-MM")).log"), Logging.Debug)),
    timestamp_logger(MinLevelLogger(FileLogger("$DIR\\Verbose-$(Dates.format(now(), "YYYY-mm-dd HH-mm")).log"), Logging.Debug))
))
# ArbBot.testLogging()