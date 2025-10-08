using ManifoldMarkets, TOML, Optimization, OptimizationBBO, Dates, Combinatorics, LinearAlgebra, Parameters, StructArrays, SciMLBase, LoopVectorization
using HTTP, JSON3, Dates
using HTTP.WebSockets
using SmartAsserts, Logging, LoggingExtras
using TimeZones
import TimeZones: unix2zdt

import ..Bots: @async_showerr, @smart_assert_showerr, display_error, wait_until, MutableOutcomeType, shouldBreak

Base.exit_on_sigint(false)

const F = Float32
const FEE = F(0.25)

struct Group
    name::String
    slugs::Vector{String}
    y_matrix::Matrix{F}
    n_matrix::Matrix{F}
    not_na_matrix::Matrix{F}
    noMarkets::Int64

    function Group(name, GroupDict)
        name = name
        slugs = urlToSlug.(collect(keys(GroupDict)))
        actionVectors = reduce(hcat, values(GroupDict))
        y_matrix = actionVectors .== "YES"
        n_matrix = actionVectors .== "NO"
        not_na_matrix = actionVectors .!= "NA"
        noMarkets = length(slugs)

        return new(name, slugs, y_matrix, n_matrix, not_na_matrix, noMarkets)
    end
end

@with_kw mutable struct MarketData
    shares::MutableOutcomeType{F} = MutableOutcomeType(0, 0)

    limitOrders::MutableOutcomeType{Dict{F, Vector{F}}} = MutableOutcomeType{Dict{F, Vector{F}}}(Dict(), Dict())
    sortedLimitProbs::MutableOutcomeType{Vector{F}} = MutableOutcomeType{Vector{F}}([], [])

    probability::F = -1
    p::F = -1
    pool::MutableOutcomeType{F} = MutableOutcomeType(0, 0)

    id::String = ""
    url::String = ""
    question::String = ""

    isResolved::Bool = false
    closeTime::Int = -1

    hasUpdated::typeof(Condition()) = Condition()
    lastOptimisedProb::F = -1
    lastContractUpdateTimeUTC::ZonedDateTime = ZonedDateTime(0, tz"UTC")
end

@with_kw mutable struct BotData @deftype String
    APIKEY # marking as const seems to break type inference
    Supabase_APIKEY
    USERNAME
    USERID
    balance::F = 0
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
    amount::F
    shares::F
    outcome::String
    redeemedMana::F

    id::String
    url::String
    question::String
end

function execute(bet, currentProb, APIKEY)
    ohno = :Success
    response = createBet(APIKEY, bet.id, bet.amount, bet.outcome)
    # need to check if returned info matches what we wanted to bet, i.e. if we got less shares than we wanted to. If we got more ig either moved or smth weird with limit orders.

    @smart_assert_showerr bet.outcome == response.outcome "$(bet.outcome), $(response.outcome)"

    if !isapprox(response.shares, bet.shares, atol=1e-2)
        @error "\e]8;;$(bet.url)\e\\$(bet.question)\e]8;;\e\\\n" # hyperlink
        @error response
        @error response.fills

        ohno = :LimitOrder

        @smart_assert_showerr !(length(response.fills) == 1 && isnothing(response.fills[end].matchedBetId) && isapprox(response.probBefore, currentProb, atol=1e-4)) "$(response.probBefore) $currentProb $bet"
    end

    if !isapprox(response.probBefore, currentProb, atol=1e-4)
        ohno = :MarketMoved
    end

    return response, ohno
end

function updateShares!(marketData, newBet, botData)
    marketData.shares[Symbol(newBet.outcome)] += newBet.shares

    if marketData.shares.YES >= marketData.shares.NO
        marketData.shares.YES -= marketData.shares.NO
        botData.balance += marketData.shares.NO #/ 2 #hack to deal with repaying loans
        marketData.shares.NO = zero(F)
    elseif marketData.shares.YES < marketData.shares.NO
        marketData.shares.NO -= marketData.shares.YES
        botData.balance += marketData.shares.YES #/ 2#hack to deal with repaying loans
        marketData.shares.YES = zero(F)
    end
end

function f(betAmount, group, marketDataBySlug, currentNoShares, currentYesShares, bettableSlugsIndex)
    profitsByEvent = zeros(eltype(group.y_matrix), size(group.y_matrix)[1])

    # newProb = zeros(group.noMarkets) # Makes it obvious which markets we don't bet on. We can print this manually, but this hides errors in fetching market probabilities
    newProb = [marketDataBySlug[slug].probability for slug in group.slugs] # So we return the correct results for markets we don't bet on  and closing soon markets

    sharesByEvent = group.y_matrix * currentYesShares + group.n_matrix * currentNoShares

    profitsByEvent = f!(betAmount, group, marketDataBySlug, bettableSlugsIndex, sharesByEvent, profitsByEvent)

    noShares = copy(currentNoShares)
    yesShares = copy(currentYesShares)

    for (i, j) in enumerate(bettableSlugsIndex)
        slug = group.slugs[j]
        market = marketDataBySlug[slug]

        if abs(betAmount[i]) >= 1.
            shares = zero(eltype(yesShares))
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

@views @fastmath function f!(betAmount, group, marketDataBySlug, bettableSlugsIndex, sharesByEvent, profitsByEvent)
    fees = zero(eltype(profitsByEvent))
    profitsByEvent .= sharesByEvent # we shouldn't count resolved markets.

    @inbounds for (i, j) in enumerate(bettableSlugsIndex)
        slug = group.slugs[j]
        @inbounds market = marketDataBySlug[slug]

        if abs(betAmount[i]) >= 1.
            shares = zero(eltype(profitsByEvent))
            @inline shares = betToShares(market.p, market.pool, market.probability, market.limitOrders, market.sortedLimitProbs, betAmount[i]).shares

            fees += FEE

            profitsByEvent .-= group.not_na_matrix[:, j] .* abs(betAmount[i])

            if betAmount[i] >= 1.
                profitsByEvent .+= group.y_matrix[:, j] .* shares
            elseif betAmount[i] <= -1.
                profitsByEvent .+= group.n_matrix[:, j] .* shares
            end
        else
            betAmount[i] = zero(eltype(betAmount))
        end
    end

    profitsByEvent .-= fees

    return profitsByEvent
end

fNoLimit!(betAmount, p) = fNoLimit!(betAmount, p.y_matrix, p.n_matrix, p.not_na_matrix, p.pVec, p.poolSOA, p.sharesByEvent, p.sharesYES, p.sharesNO)
function fNoLimit!(betAmount, y_matrix, n_matrix, not_na_matrix, pVec, pool, sharesByEvent, sharesNO, sharesYES)
    sharesYES .= 0
    sharesNO .= 0

    betCount = zero(eltype(sharesByEvent))

    @turbo for i in eachindex(betAmount)
        y = pool.YES[i]
        n = pool.NO[i]

        p = pVec[i]

        sharesYES[i] = ifelse(betAmount[i] >= 1, y + betAmount[i] - y * (n / (n + betAmount[i]))^((1-p)/p), 0)
        sharesNO[i] = ifelse(betAmount[i] <= -1, n - betAmount[i] - n * (y / (y - betAmount[i]))^(p/(1-p)), 0)

        betCount += abs(betAmount[i]) >= 1
    end

    minProfits = F(Inf)

    @turbo for j in axes(n_matrix, 1)
        profit = sharesByEvent[j]

        for i in axes(n_matrix, 2)
            profit += y_matrix[j, i] * sharesYES[i]
            profit += n_matrix[j, i] * sharesNO[i]
            profit -= not_na_matrix[j, i] * abs(betAmount[i])
        end

        minProfits = min(minProfits, profit)
    end

    return -(minProfits - betCount * F(FEE))
end

mutable struct pStruct
    const y_matrix::Matrix{F}
    const n_matrix::Matrix{F}
    const not_na_matrix::Matrix{F}
    const pVec::Vector{F}
    const poolSOA::StructVector{MutableOutcomeType{Float32}, NamedTuple{(:YES, :NO), Tuple{Vector{Float32}, Vector{Float32}}}, Int64}
    const sharesByEvent::Vector{Float32}
    sharesYES::Vector{Float32}
    sharesNO::Vector{Float32}
end
# using FunctionWrappers
# import FunctionWrappers: FunctionWrapper
function optimise(group, marketDataBySlug, maxBetAmount, bettableSlugsIndex)
    noShares = [marketDataBySlug[slug].shares.NO for slug in group.slugs]
    yesShares = [marketDataBySlug[slug].shares.YES for slug in group.slugs]
    sharesByEvent = group.y_matrix * yesShares + group.n_matrix * noShares

    profitsByEvent = similar(sharesByEvent)

    pVec = [marketDataBySlug[slug].p for slug in group.slugs[bettableSlugsIndex]]
    poolSOA = StructArray([marketDataBySlug[slug].pool for slug in group.slugs[bettableSlugsIndex]])

    sharesYES = zeros(F, length(bettableSlugsIndex))
    sharesNO = zeros(F, length(bettableSlugsIndex))
    y_matrix = group.y_matrix[:, bettableSlugsIndex]
    n_matrix = group.n_matrix[:, bettableSlugsIndex]
    not_na_matrix = group.not_na_matrix[:, bettableSlugsIndex]

    if all(slug -> isempty(marketDataBySlug[slug].sortedLimitProbs.YES) && isempty(marketDataBySlug[slug].sortedLimitProbs.NO), group.slugs[bettableSlugsIndex])
        # profitF = OptimizationFunction((betAmount, _) -> fNoLimit!(betAmount, group, pVec, poolSOA, bettableSlugsIndex, sharesByEvent, profitsByEvent))

        profitF = OptimizationFunction(fNoLimit!)

        # profitF = OptimizationFunction((betAmount, _) -> fNoLimit!(betAmount, y_matrix, n_matrix, not_na_matrix, pVec, poolSOA, sharesByEvent, sharesYES, sharesNO))
        # @inline fTest(betAmount, p) = @inline fNoLimit!(betAmount, p.y_matrix, p.n_matrix, p.not_na_matrix, p.pVec, p.poolSOA, p.sharesByEvent, p.sharesYES, p.sharesNO)

        # @benchmark $profitF($betAmounts, $p)

        # profitF = OptimizationFunction((betAmount, p) -> fNoLimit!(betAmount, (; group, pVec, poolSOA, bettableSlugsIndex, sharesByEvent, profitsByEvent)...))
        # profitF = OptimizationFunction(FunctionWrapper{Float64, Tuple{Vector{F}, Tuple{Group, Vector{F}, StructVector{MutableOutcomeType{F}, NamedTuple{(:YES, :NO), Tuple{Vector{F}, Vector{F}}}, Int64}, Vector{Int64}, Vector{F}, Vector{F}}}}((betAmount, p) -> fNoLimit!(betAmount, (; group, pVec, poolSOA, bettableSlugsIndex, sharesByEvent, profitsByEvent)...)))
        # profitF = OptimizationFunction(FunctionWrapper{F, Tuple{Vector{F}}}((betAmount,) -> fNoLimit!(betAmount, group, pVec, poolSOA, bettableSlugsIndex, sharesByEvent, profitsByEvent)))
    else
        profitF = OptimizationFunction((betAmount, _) -> -minimum( f!(betAmount, group, marketDataBySlug, bettableSlugsIndex, sharesByEvent, profitsByEvent)) )
    end

    x0 = zeros(F, length(bettableSlugsIndex))

    ub = F[marketDataBySlug[slug].probability > .95 ? max(0, (marketDataBySlug[slug].shares.NO * (1-marketDataBySlug[slug].probability))) : maxBetAmount for slug in group.slugs[bettableSlugsIndex]]
    lb = F[marketDataBySlug[slug].probability < .05 ? min(0, -(marketDataBySlug[slug].shares.YES * marketDataBySlug[slug].probability)) : -maxBetAmount for slug in group.slugs[bettableSlugsIndex]]

    problem = Optimization.OptimizationProblem(profitF, x0, pStruct(y_matrix, n_matrix, not_na_matrix, pVec, poolSOA, sharesByEvent, sharesYES, sharesNO), lb=lb, ub=ub)#, sense=Optimization.MaxSense)#, abstol=1e-3, MaxStepsWithoutProgress=10^3)

    @debug "Running adaptive for $(group.name)"
    # @descend solve(problem, BBO_de_rand_1_bin_radiuslimited(), maxtime=.15)
    @time sol = solve(problem, BBO_de_rand_1_bin_radiuslimited(), maxtime=.02)
    # @time sol = solve(problem, BBO_de_rand_1_bin_radiuslimited(), maxiters=10^5)

    @debug "Yielding after adaptive, $(group.name)"
    yield()
    @debug "Done yielding after adaptive, $(group.name)"

    bestSolution = zeros(F, length(bettableSlugsIndex))
    maxRiskFreeProfit = f(bestSolution, group, marketDataBySlug, noShares, yesShares, bettableSlugsIndex).profitsByEvent |> minimum

    nonZeroIndices = findall(!iszero, sol.u::Vector{F})

    betAmount::Vector{F} = zeros(length(bettableSlugsIndex))

    for indices in powerset(nonZeroIndices)
        betAmount .= sol.u::Vector{F}
        betAmount[indices] .= 0

        riskFreeProfit = -profitF(betAmount, pStruct(y_matrix, n_matrix, not_na_matrix, pVec, poolSOA, sharesByEvent, sharesYES, sharesNO))
        # riskFreeProfit = f(betAmount, group, marketDataBySlug, noShares, yesShares, bettableSlugsIndex).profitsByEvent |> minimum

        if riskFreeProfit > maxRiskFreeProfit
            maxRiskFreeProfit = riskFreeProfit
            bestSolution .= betAmount
        end
    end

    # if bestSolution == repeat([0.], length(bettableSlugsIndex))
    #     @debug "Running resampling for $(group.name)"

    #     @time "Resampling" sol2 = solve(problem, BBO_resampling_memetic_search(), maxtime=.3)

    #     @debug "Yielding after resampling, $(group.name)"
    #     # sleep(1)
    #     yield()
    #     @debug "Done yielding after resampling, $(group.name)"

    #     nonZeroIndices = findall(!iszero, sol2.u)

    #     for indices in powerset(nonZeroIndices)
    #         betAmount .= sol2.u
    #         betAmount[indices] .= 0

    #         riskFreeProfit = f(betAmount, group, marketDataBySlug, noShares, yesShares, bettableSlugsIndex).profitsByEvent |> minimum

    #         if riskFreeProfit > maxRiskFreeProfit
    #             maxRiskFreeProfit = riskFreeProfit
    #             bestSolution .= betAmount
    #         end
    #     end
    # end

    return bestSolution
end

function updateMarketData!(marketData, market, fs_updated_time)
    # marketData.probability = market.prob  # this updates slowly?

    marketData.p = market.p # can change if a subsidy is given
    marketData.pool.YES = market.pool.YES
    marketData.pool.NO = market.pool.NO

    @debug "Updating market probability from $(marketData.probability) to $(poolToProb(marketData.p, marketData.pool)), slug: $(market.slug)"
    marketData.probability = poolToProb(marketData.p, marketData.pool)

    marketData.id = market.id
    marketData.question = market.question
    marketData.url = "https://manifold.markets/$(market.creatorUsername)/$(market.slug)"

    marketData.isResolved = market.isResolved
    marketData.closeTime = market.closeTime

    marketLastUpdateUTC = ZonedDateTime(DateTime(fs_updated_time, dateformat"yyyy-mm-ddTHH:MM:SS.sss"), tz"UTC")

    @debug "Updating market lastContractUpdateTimeUTC from $(marketData.lastContractUpdateTimeUTC) to $marketLastUpdateUTC, slug: $(market.slug)"
    if marketLastUpdateUTC < marketData.lastContractUpdateTimeUTC
        throw(ErrorException("Updated with old market data?, $(marketData.lastContractUpdateTimeUTC), $marketLastUpdateUTC, $marketData, $market, $fs_updated_time"))
    end

    marketData.lastContractUpdateTimeUTC = marketLastUpdateUTC
end

function updateLimitOrders!(marketData)
    for (i, limitProb) in enumerate(marketData.sortedLimitProbs.YES) # limit orders buying no
        if limitProb + 1e-4 < marketData.probability
            @debug "Deleted limit order buying NO at $limitProb as market probability is now $(marketData.probability)"
            deleteat!(marketData.sortedLimitProbs.YES, i)
            delete!(marketData.limitOrders.YES, limitProb)
            @debug marketData.sortedLimitProbs
            @debug marketData.limitOrders
        end
    end

    for (i, limitProb) in enumerate(marketData.sortedLimitProbs.NO) # limit orders buying yes
        if limitProb - 1e-4 > marketData.probability
            @debug "Deleted limit order buying YES at $limitProb as market probability is now $(marketData.probability)"
            deleteat!(marketData.sortedLimitProbs.NO, i)
            delete!(marketData.limitOrders.NO, limitProb)
            @debug marketData.sortedLimitProbs
            @debug marketData.limitOrders
        end
    end
end

function getMarkets!(marketDataBySlug, slugs, Supabase_APIKEY)
    @sync for paritionedSlugs in Iterators.partition(slugs, 20) # 2000 is max characters in uri, so we have about 1900 for slugs, max slug is about 50 : 1900/50 is 38. generous magin to 20
        @async_showerr begin
            contracts_slugs = join(paritionedSlugs, ",")

            response = HTTP.get("https://pxidrgkatumlvfqaxcll.supabase.co/rest/v1/contracts?slug=in.($contracts_slugs)", headers= ["apikey" => Supabase_APIKEY, "Content-Type" => "application/json"])
            responseJSON = JSON3.read(response.body)

            for contract in responseJSON
                updateMarketData!(marketDataBySlug[contract.slug], contract.data, contract.fs_updated_time)
            end
        end
    end
end

function getMarketsUsingId!(marketDataBySlug, slugs) # If we use slugs the request uri is too long
    for slug in slugs
        @smart_assert_showerr !isnothing(marketDataBySlug[slug].id)
    end

    contracts_ids = join(map(slug -> marketDataBySlug[slug].id, slugs), ",") # might be faster to index by id
    # contracts_slugs = join(slugs, ",")

    @smart_assert_showerr length(slugs) < 1000 "too many slugs $(length(slugs)), $slugs when fetching markets"

    response = HTTP.get("https://pxidrgkatumlvfqaxcll.supabase.co/rest/v1/contracts?id=in.($contracts_ids)", headers= ["apikey" => "eyJhbGciOiJIUzI1NiIsInR5cCI6IkpXVCJ9.eyJpc3MiOiJzdXBhYmFzZSIsInJlZiI6InB4aWRyZ2thdHVtbHZmcWF4Y2xsIiwicm9sZSI6ImFub24iLCJpYXQiOjE2Njg5OTUzOTgsImV4cCI6MTk4NDU3MTM5OH0.d_yYtASLzAoIIGdXUBIgRAGLBnNow7JG2SoaNMQ8ySg", "Content-Type" => "application/json"])
    responseJSON = JSON3.read(response.body)

    for contract in responseJSON
        updateMarketData!(marketDataBySlug[contract.slug], contract.data, contract.fs_updated_time)
    end
end

getSlugs(GROUPS::Dict) = mapreduce(x -> urlToSlug.(x), vcat, keys.(values(GROUPS)))
getSlugs(groups::Vector{Group}) = mapreduce(group -> group.slugs, vcat, groups)

isMarketClosingSoon(market) = market.isResolved || market.closeTime / 1000 < time() + 60 # if resolved or closing in 60 seconds

function arbitrageGroup(group, botData, marketDataBySlug, Arguments, firstSlugToBet=nothing)::Symbol
    rerun = :Success

    for slug in group.slugs
        @smart_assert_showerr marketDataBySlug[slug].probability ≈ poolToProb(marketDataBySlug[slug].p, marketDataBySlug[slug].pool) #"$slug, $(marketDataBySlug[slug]), $(marketDataBySlug[slug].pool), $(poolToProb(marketDataBySlug[slug].p, marketDataBySlug[slug].pool))"
    end

    # Actually could close in the delay between running and here, or due to reruns. But we don't need to pop the group from groups
    allMarketsClosing = true
    for slug in group.slugs
        if !isMarketClosingSoon(marketDataBySlug[slug])
            allMarketsClosing = false
        end
    end

    bettableSlugsIndex = [i for (i, slug) in enumerate(group.slugs) if !isMarketClosingSoon(marketDataBySlug[slug])]

    if allMarketsClosing
        printstyled("$(group.name) all markets closed at $(Dates.format(now(), "HH:MM:SS.sss"))\n", bold=true, underline=true)

        for slug in group.slugs # So we don't rerun these markets
            marketDataBySlug[slug].lastOptimisedProb = marketDataBySlug[slug].probability # not great if we have insufficient balance as later we might have sufficient but think the market is already optimised
        end

        rerun = :Success
        return rerun
    end

    plannedBets = PlannedBet[]

    maxBetAmount = max(50, botData.balance-100) / (2 + 1.5*group.noMarkets)
    redeemManaHack = maximum(slug -> max(marketDataBySlug[slug].shares.YES / (1 - marketDataBySlug[slug].probability), marketDataBySlug[slug].shares.NO / marketDataBySlug[slug].probability), group.slugs[bettableSlugsIndex])
    maxBetAmount += min(redeemManaHack, maxBetAmount, 50.)
    # Remove hack for now as it seems to prevent optmiser betting as much as it can which prevents reruns

    betAmounts = optimise(group, marketDataBySlug, maxBetAmount, bettableSlugsIndex)

    oldProb = [marketDataBySlug[slug].probability for slug in group.slugs]
    oldNoShares = [marketDataBySlug[slug].shares.NO for slug in group.slugs]
    oldYesShares = [marketDataBySlug[slug].shares.YES for slug in group.slugs]

    newProfitsByEvent, noShares, yesShares, newProb = f(betAmounts, group, marketDataBySlug, oldNoShares, oldYesShares, bettableSlugsIndex)

    oldProfitsByEvent, _, _, _ = f(zeros(F, length(bettableSlugsIndex)), group, marketDataBySlug, oldNoShares, oldYesShares, bettableSlugsIndex)

    profit = minimum(newProfitsByEvent) - minimum(oldProfitsByEvent)
    newYesShares = yesShares .- oldYesShares
    newNoShares = noShares  .- oldNoShares

    # @debug bettableSlugsIndex
    @debug oldProb
    @debug newProb
    @debug oldProfitsByEvent
    @debug newProfitsByEvent
    @debug profit
    @debug betAmounts
    @debug group.y_matrix
    @debug group.n_matrix
    @debug group.not_na_matrix
    @debug map(slug -> marketDataBySlug[slug].shares.YES, group.slugs)
    @debug map(slug -> marketDataBySlug[slug].shares.NO, group.slugs)
    @debug oldNoShares
    @debug oldYesShares
    @debug newYesShares
    @debug newNoShares
    @debug yesShares
    @debug noShares
    @debug group.y_matrix * yesShares
    @debug group.n_matrix * noShares
    @debug group.y_matrix * yesShares + group.n_matrix * noShares
    @debug group.not_na_matrix[:, bettableSlugsIndex] * betAmounts
    @debug sum(abs.(betAmounts))
    @debug [marketDataBySlug[slug].p for slug in group.slugs]
    @debug [marketDataBySlug[slug].pool for slug in group.slugs]
    @debug maxBetAmount
    # @debug [marketDataBySlug[slug].limitOrders for slug in group.slugs]
    # @debug [marketDataBySlug[slug].sortedLimitProbs for slug in group.slugs]

    for (i, j) in enumerate(bettableSlugsIndex)
        amount = betAmounts[i]
        slug = group.slugs[j]

        if isapprox(amount, 0., atol=1e-6)
            continue
        elseif abs(amount) >= 1.
            outcome = amount > 0. ? "YES" : "NO"
            shares = newYesShares[j] + newNoShares[j]
            @smart_assert_showerr isapprox(newYesShares[j], 0, atol=1e-2) || isapprox(newNoShares[j], 0, atol=1e-2) "$j"
            redeemedMana = zero(shares)
            if outcome == "YES"
                redeemedMana = min(marketDataBySlug[slug].shares.NO, shares)
            elseif outcome == "NO"
                redeemedMana = min(marketDataBySlug[slug].shares.YES, shares)
            end

            @smart_assert_showerr (sign(amount) == sign(newProb[j] - oldProb[j]) || newProb[j] ≈ marketDataBySlug[slug].sortedLimitProbs[Symbol(outcome)][1]) "$slug, $amount, $(newProb[j]), $(oldProb[j]), $i, $j"

            bet = PlannedBet(abs(amount), shares, outcome, redeemedMana, marketDataBySlug[slug].id, marketDataBySlug[slug].url, marketDataBySlug[slug].question)
            push!(plannedBets, bet)
        else
            # rerun = :BetMore
            @warn "Bet amount: $(amount) is too small, $(slug)"

            rerun = :Success
            return rerun
        end
    end

    if !isnothing(firstSlugToBet)
        for (i, bet) in enumerate(plannedBets)
            if bet.id == marketDataBySlug[firstSlugToBet].id
                plannedBets[i] = plannedBets[1]
                plannedBets[1] = bet
                break
            end
        end
        if length(plannedBets) > 1
            sort!(@views(plannedBets[2:end]), by = bet -> bet.redeemedMana, rev=true)
        end
    else
        sort!(plannedBets, by = bet -> bet.redeemedMana, rev=true)
    end

    # we never return a profit less than 0 so we're not going to redeem bets that have profit less than 0.
    if (profit ≤ 0) && !((profit + FEE * length(plannedBets) ≥ 0) && (sum(bet -> bet.redeemedMana, plannedBets, init=0.) > 1.))
        for (i, slug) in enumerate(group.slugs)
            marketDataBySlug[slug].lastOptimisedProb = newProb[i] # not great if we have insufficient balance as later we might have sufficient but think the market is already optimised
        end

        rerun = :Success
        return rerun
    end

    buffer = IOBuffer()
    newProbBySlug = Dict(group.slugs[j] => newProb[j] for j in bettableSlugsIndex)
    for bet in plannedBets
        amount = bet.amount
        slug = urlToSlug(bet.url)

        if abs(amount) >= .98 * maxBetAmount
            @warn "Bet size is $(string(100 * amount/maxBetAmount))% of maxBetAmount"
            # bindingConstraint = true
            rerun = :BetMore
        end

        if abs(amount) ≈ 1
            @warn "Bet size is $(string(amount))"
            # bindingConstraint = true
            rerun = :BetMore
        end

        println(buffer, "\e]8;;$(string(marketDataBySlug[slug].url))\e\\$(string(marketDataBySlug[slug].question))\e]8;;\e\\
        Prior probs:     $(string(marketDataBySlug[slug].probability * 100))%
        Posterior probs: $(string(newProbBySlug[slug]*100))%
        Buy $(string(bet.shares)) $(string(bet.outcome)) shares for $(string(amount)), redeeming $(string(bet.redeemedMana))")
    end

    if (sum(abs.(betAmounts)) >= sum(bet -> bet.redeemedMana, plannedBets)/2 + botData.balance - 100) && (sum(abs.(betAmounts)) >= sum(bet -> bet.redeemedMana, plannedBets)/2) || plannedBets[1].amount >= botData.balance # /2 is hack to deal with repaying loans
        println(buffer, "Expected Profits:         $profit")
        @info String(take!(buffer))

        @error "Insufficient Balance $(botData.balance) for $(sum(abs.(betAmounts))) bet redeeming $(sum(bet -> bet.redeemedMana, plannedBets))."
        rerun = :InsufficientBalance
        return rerun
    end

    if Arguments.confirmBets
        println(buffer, "Expected Profits:         $profit")
        @info String(take!(buffer))

        println("Proceed? (y/n)")
        if readline() !="y"
            rerun = :Success
            return rerun
        end
    end

    if Arguments.live
        oldProbBySlug = Dict(group.slugs[j] => oldProb[j] for j in bettableSlugsIndex)
        timers = Vector{Timer}(undef, length(plannedBets))

        for (i, slug) in enumerate(group.slugs)
            marketDataBySlug[slug].lastOptimisedProb = newProb[i] # not great if we have insufficient balance as later we might have sufficient but think the market is already optimised
        end

        movedMarkets = 0
        skip = false
        waitTimedOut = false

        if any(oldProb .!= [marketDataBySlug[slug].probability for slug in group.slugs])
            @info "Market moved after optimisation"
            return :UnexpectedBet
        end

        try
            # need to add check to see if we're actually receiving messages
            @sync for (i, bet) in enumerate(plannedBets)
                @async_showerr begin
                    if skip
                        return nothing
                    end

                    slug = urlToSlug(bet.url)

                    if oldProbBySlug[slug] != marketDataBySlug[slug].probability
                        movedMarkets += 1
                        rerun = :UnexpectedBet
                        @info "$slug moved from $(oldProbBySlug[slug]) to $(marketDataBySlug[slug].probability)"
                        skip = true

                        return nothing
                    end

                    # We use oldProb instead of MarketData as MarketData may have updated to a new probability while we bet.
                    start_time = time()
                    executedBet, ohno = execute(bet, oldProbBySlug[slug], botData.APIKEY)
                    println("\e]8;;$(string(marketDataBySlug[slug].url))\e\\$(string(marketDataBySlug[slug].question))\e]8;;\e\\: $(time() - start_time)s")

                    botData.balance -= executedBet.amount

                    if ohno != :Success
                        if executedBet.outcome == "YES"
                            yesShares[group.slugs .== slug] .= oldYesShares[group.slugs .== slug] .+ executedBet.shares
                        elseif executedBet.outcome == "NO"
                            noShares[group.slugs .== slug] .= oldNoShares[group.slugs .== slug] .+ executedBet.shares
                        end

                        rerun = :UnexpectedBet
                    end

                    if ohno == :MarketMoved
                        movedMarkets += 1
                    end

                    updateShares!(marketDataBySlug[slug], executedBet, botData)

                    if !isnothing(executedBet.fills[end].matchedBetId)
                        # Limit order might not update in time?
                        limitOrder = getBet(executedBet.fills[end].matchedBetId, botData.Supabase_APIKEY)
                        @debug limitOrder

                        amountLeft = zero(F)
                        sharesLeft = zero(F)
                        if limitOrder.outcome == "NO"
                            sharesLeft = (limitOrder.orderAmount - limitOrder.amount) / (1 - limitOrder.limitProb)
                            amountLeft = sharesLeft * limitOrder.limitProb
                        elseif limitOrder.outcome == "YES"
                            sharesLeft = (limitOrder.orderAmount - limitOrder.amount) / limitOrder.limitProb
                            amountLeft = sharesLeft * (1 - limitOrder.limitProb)
                        end

                        marketDataBySlug[slug].limitOrders[Symbol(executedBet.outcome)] = Dict(F(limitOrder.limitProb) => F[amountLeft, sharesLeft])

                        @debug marketDataBySlug[slug].limitOrders

                        if executedBet.outcome == "YES"
                            marketDataBySlug[slug].sortedLimitProbs.YES = [limitOrder.limitProb]
                        elseif executedBet.outcome == "NO"
                            marketDataBySlug[slug].sortedLimitProbs.NO = [limitOrder.limitProb]
                        end

                        @debug marketDataBySlug[slug].sortedLimitProbs
                    else # Should only happen if theres a limit order that we expected to hit but ended up not, e.g. user runs out of balance.
                        @debug "Removing limit orders on $slug, $(executedBet.outcome), $(marketDataBySlug[slug].limitOrders), $(marketDataBySlug[slug].sortedLimitProbs)"
                        marketDataBySlug[slug].limitOrders[Symbol(executedBet.outcome)] = Dict()
                        marketDataBySlug[slug].sortedLimitProbs[Symbol(executedBet.outcome)] = []
                        @debug "Removed limit orders $slug, $(marketDataBySlug[slug].limitOrders), $(marketDataBySlug[slug].sortedLimitProbs)"
                    end

                    @debug "$(Dates.format(now(), "HH:MM:SS.sss")) Waiting begun for $slug"
                    # Don't place @debug on wait_until otherwise error is not thrown.
                    timers[i] = Timer(60) do _
                        waitTimedOut = true
                        notify(marketDataBySlug[slug].hasUpdated, "Wait timed out waiting for new bet, UTC: $(unix2zdt(executedBet.createdTime/10^3)), $(marketDataBySlug[slug].lastContractUpdateTimeUTC) on $slug", error=true)
                    end
                    try
                        while marketDataBySlug[slug].lastContractUpdateTimeUTC < unix2zdt(executedBet.createdTime/10^3)
                            @debug "$(Dates.format(now(), "HH:MM:SS.sss")) Waiting for new bet, UTC: $(unix2zdt(executedBet.createdTime/10^3)), $(marketDataBySlug[slug].lastContractUpdateTimeUTC)"
                            wait(marketDataBySlug[slug].hasUpdated)
                        end
                    finally
                        close(timers[i])
                    end

                    @debug "$(Dates.format(now(), "HH:MM:SS.sss")) Waiting done for $slug"
                end

                sleep(0.4) # so manifold processes requests in desired order. doesn't seem to help much with a 100ms sleep
            end
        catch err
            println(buffer, "Expected Profits:         $profit") # no more fee, but we still want to use fee in optimisation
            println(buffer, "Actual Profits:           $(@turbo minimum(group.y_matrix * yesShares + group.n_matrix * noShares - group.not_na_matrix[:, bettableSlugsIndex] * abs.(betAmounts)) - minimum(oldProfitsByEvent) - FEE * length(plannedBets))")
            @info String(take!(buffer))

            for i in eachindex(timers)
                if isdefined(timers, i)
                    close(timers[i])
                end
            end

            if waitTimedOut
                rethrow()
            end

            rerun = :PostFailure

            return rerun
        end

        println(buffer, "Expected Profits:         $profit") # no more fee, but we still want to use fee in optimisation
        println(buffer, "Actual Profits:           $(@turbo minimum(group.y_matrix * yesShares + group.n_matrix * noShares - group.not_na_matrix[:, bettableSlugsIndex] * abs.(betAmounts)) - minimum(oldProfitsByEvent) - FEE * length(plannedBets))")

        # Doesn't account for redemptions
        # println(buffer, "Actual Profits:         $(minimum(group.y_matrix * [marketDataBySlug[slug].shares.YES for slug in group.slugs] + group.n_matrix * [marketDataBySlug[slug].shares.NO for slug in group.slugs] - group.not_na_matrix[:, bettableSlugsIndex] * abs.(betAmounts)) - minimum(oldProfitsByEvent))") # no more fee, but we still want to use fee in optimisation
        @info String(take!(buffer))

        if movedMarkets >= 2
            @debug "Sleeping for two seconds as $movedMarkets markets moved"
            sleep(2)
        end
        return rerun
    else
        println(buffer, "Expected Profits:         $profit")
        @info String(take!(buffer))

        return :Success
    end
end

fetchMyShares!(marketDataBySlug, groupData, USERID) = fetchMyShares!(marketDataBySlug,  groupData, USERID, Base.broadcasted(getproperty, values(marketDataBySlug), :id))

function fetchMyShares!(marketDataBySlug, groupData, USERID, slugs)
    # https://discourse.julialang.org/t/broadcast-object-property/47104/6
    contracts = join(slugs, ",")

    @smart_assert_showerr length(marketDataBySlug) < 1000

    # will fail if we have more than 1000
    response = HTTP.get("https://pxidrgkatumlvfqaxcll.supabase.co/rest/v1/user_contract_metrics?user_id=eq.$USERID&contract_id=in.($contracts)", headers= ["apikey" => "eyJhbGciOiJIUzI1NiIsInR5cCI6IkpXVCJ9.eyJpc3MiOiJzdXBhYmFzZSIsInJlZiI6InB4aWRyZ2thdHVtbHZmcWF4Y2xsIiwicm9sZSI6ImFub24iLCJpYXQiOjE2Njg5OTUzOTgsImV4cCI6MTk4NDU3MTM5OH0.d_yYtASLzAoIIGdXUBIgRAGLBnNow7JG2SoaNMQ8ySg", "Content-Type" => "application/json"])
    responseJSON = JSON3.read(response.body)

    # what if it returns nothing for some contract? presumably that means we never invested so MarketData should already be correct

    updated_time_per_contract = Dict{String, ZonedDateTime}()

    for contract_metrics in responseJSON
        slug = groupData.contractIdToSlug[contract_metrics.contract_id]
        time = ZonedDateTime(DateTime(contract_metrics.fs_updated_time, dateformat"yyyy-mm-ddTHH:MM:SS.sss"), tz"UTC")
        if haskey(updated_time_per_contract, slug) && updated_time_per_contract[slug] >= time
            continue
        end
        updated_time_per_contract[slug] = time
        for (outcome, shares) in contract_metrics.data.totalShares
            marketDataBySlug[slug].shares[outcome] = shares
        end
    end
end

function readData()
    data = TOML.parsefile("$(@__DIR__)/Arb.toml")
    GROUPS::Dict{String, Dict{String, Vector{String}}} = data["GROUPS"]
    APIKEY::String = data["APIKEY"]
    Supabase_APIKEY::String = data["SUPABASE_APIKEY"]
    USERNAME::String = data["USERNAME"]

    function checkURL(url)
        if !occursin(r"https://manifold.markets/[\w\d]+/[\w\d]+", url)
            error("$url is not a valid url")
        end
    end
    foreach(x -> checkURL.(x), keys.(values(GROUPS)))

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
        if '#' in slug || '?' in slug
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

    @time "Fetching Markets" getMarkets!(marketDataBySlug, getSlugs(groups), Supabase_APIKEY)
    contractIdSet = Set(market.id for market in values(marketDataBySlug))
    contractIdToGroupIndex = Dict(marketDataBySlug[slug].id => i for (i, group) in enumerate(groups) for slug in group.slugs)
    contractIdToSlug = Dict(marketDataBySlug[slug].id => slug for group in groups for slug in group.slugs)
    groupData = GroupData(groups, contractIdSet, contractIdToGroupIndex, contractIdToSlug)

    @time "Fetching Shares" fetchMyShares!(marketDataBySlug, groupData, USERID)


    return groupData, botData, marketDataBySlug, arguments
end

function runGroup(group, groupData, botData, marketDataBySlug, arguments, taskValue=nothing, slug=nothing)
    delay = 60
    runs = 0
    rerun = :FirstRun

    numberOfPostFailures = 0

    while rerun == :FirstRun || (rerun == :BetMore && runs ≤ 5) || (rerun == :UnexpectedBet && runs ≤ 1) || rerun == :PostFailure
        if rerun == :PostFailure
            numberOfPostFailures += 1
            if numberOfPostFailures >= 2
                throw(ErrorException("Multiple failures to create bets."))
            end

            sleep(delay)
            delay *= 5

            # What if a bet didn't go througth, thus we must fetch
            @info "Fetching shares and balance"
            fetchMyShares!(marketDataBySlug, groupData, botData.USERID, group.slugs)
            botData.balance = getUserByUsername(botData.USERNAME).balance
            getMarkets!(marketDataBySlug, group.slugs, botData.Supabase_APIKEY)
        end

        # if rerun != :FirstRun && rerun != :numberOfPostFailures
        #     @info "Fetching markets"
        #     getMarkets!(marketDataBySlug, group.slugs, botData.Supabase_APIKEY)
        # end

        @warn "Running $(group.name) at $(Dates.format(now(), "HH:MM:SS.sss"))"
        @debug [marketDataBySlug[slug].p for slug in group.slugs]
        @debug [marketDataBySlug[slug].pool for slug in group.slugs]
        runs += 1

        if !isnothing(taskValue)
            # no need to run any previously queued requests.
            if taskValue.runAgain
                taskValue = (runAgain = false, task=current_task())
            end
        end

        rerun = arbitrageGroup(group, botData, marketDataBySlug, arguments, slug)
        @debug rerun
    end
end

function production(groupNames = nothing; live=true, confirmBets=false, skip=false)
    # Should probably move this into the websocket as well
    groupData, botData, marketDataBySlug, arguments = setup(groupNames, live, confirmBets)

    TaskDict = Dict(i => (runAgain=false, task=@async nothing) for i in eachindex(groupData.groups)) # so we don't have to check if there is a task in it or not. Order of tuple matters for parsing

    # currentTask = @async nothing

    errorOccurred = false

    # # Fetch new balance every hour, hack to work around repaying loans when redeeming
    # @async_showerr while true
    #     @info "3: Fetching Balance at $(Dates.format(now(), "HH:MM:SS.sss"))"
    #     botData.balance = getUserByUsername(botData.USERNAME).balance

    #     sleep(Hour(1))
    # end

    # while true
    #     @time "Fetching Markets" getMarkets!(marketDataBySlug, getSlugs(groupData.groups), botData.Supabase_APIKEY)

    #     @warn "All: Running all groups at $(Dates.format(now(), "HH:MM:SS.sss"))"
    #     for group in groupData.groups
    #         if errorOccurred
    #             break
    #         end
    #         runGroup(group, groupData, botData, marketDataBySlug, arguments)

    #         for slug in group.slugs
    #             marketDataBySlug[slug].limitOrders = MutableOutcomeType(Dict{F, Vector{F}}(), Dict{F, Vector{F}}()) # need to reset as we aren't tracking limit orders
    #             marketDataBySlug[slug].sortedLimitProbs = MutableOutcomeType(F[], F[])
    #         end
    #     end

    #     @warn "All: Done all groups at $(Dates.format(now(), "HH:MM:SS.sss"))"
    #     sleep(60+rand()*10-5)
    # end

    WebSockets.open(uri(botData.Supabase_APIKEY), suppress_close_error=true, headers=["User-Agent" => "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/70.0.3538.102 Safari/537.36 Edge/18.19582"]) do socket
        try
            @info "Opened Socket"
            @debug socket
            # send(socket, pushJSONContracts("contracts", "live-contracts-"))
            send(socket, pushJSONContracts())
            # send(socket, pushJSONBets())
            @info "Sent Intialisation"

            Base.Experimental.@sync begin # So that we exit on error straight away and close socket
                # Reads messages whilst in this loop but blocks arbing due to them, also ensures MarketData updates after we make a bet
                currentTask = @async_showerr if !skip && !errorOccurred
                    # sleep(5) # idk might work to fix bets not coming through on first run
                    @warn "All: Running all groups at $(Dates.format(now(), "HH:MM:SS.sss"))"
                    for group in groupData.groups
                        if errorOccurred
                            break
                        end
                        runGroup(group, groupData, botData, marketDataBySlug, arguments)

                        for slug in group.slugs
                            marketDataBySlug[slug].limitOrders = MutableOutcomeType(Dict{F, Vector{F}}(), Dict{F, Vector{F}}()) # need to reset as we aren't tracking limit orders
                            marketDataBySlug[slug].sortedLimitProbs = MutableOutcomeType(F[], F[])
                        end
                    end

                    @warn "All: Done all groups at $(Dates.format(now(), "HH:MM:SS.sss"))"
                end
                #Reading messages
                # should switch to @spawn :interactive
                @async_showerr for msg in socket
                    @async_showerr begin
                        # GC.enable(false)
                        msgJSON = JSON3.read(msg)
                        if !(:data in keys(msgJSON.payload))
                            if :status in keys(msgJSON.payload) && msgJSON.payload.status == "error"
                                # @error "$(Dates.format(now(), "HH:MM:SS.sss")): $msg"
                                throw(ErrorException("$(Dates.format(now(), "HH:MM:SS.sss")): $msg"))
                            else
                                @debug "$(Dates.format(now(), "HH:MM:SS.sss")): $msg"
                                return nothing
                            end
                        end

                        if !(:data in keys(msgJSON.payload.data.record)) # message doesn't seem to include probabilities or anything useful, only data about relevancy of market.
                            return nothing
                        end

                        market = msgJSON.payload.data.record.data
                        marketId = market.id
                        # @debug market.slug

                        if marketId in groupData.contractIdSet
                            @smart_assert_showerr market.mechanism == "cpmm-1"

                            slug = groupData.contractIdToSlug[marketId] #market.slug also works
                            oldProb = marketDataBySlug[slug].probability

                            marketLastUpdateUTC = ZonedDateTime(DateTime(msgJSON.payload.data.record.fs_updated_time, dateformat"yyyy-mm-ddTHH:MM:SS.sss"), tz"UTC")
                            marketLastUpdateLocal = astimezone(marketLastUpdateUTC, localzone())

                            # println("Received message: $(marketLastUpdate)")
                            # @warn slug
                            # if marketDataBySlug[slug].probability ≉ oldProb
                                @debug "Received message $slug at $(poolToProb(market.p, market.pool)) with lastContractUpdateTimeBritish $(marketLastUpdateLocal) at $(now())"
                            # end
                            # @warn "from prob $(market.prob)"

                            if marketLastUpdateUTC < marketDataBySlug[slug].lastContractUpdateTimeUTC
                                @debug "Old bet? $marketLastUpdateLocal, $(astimezone(marketDataBySlug[slug].lastContractUpdateTimeUTC, localzone())) for $slug at $(marketDataBySlug[slug].probability) at $(now())"
                                return nothing
                            end
                            updateMarketData!(marketDataBySlug[slug], market, msgJSON.payload.data.record.fs_updated_time)
                            updateLimitOrders!(marketDataBySlug[slug]) # need to delete any limit orders that we know have been filled

                            # Should probably only do this if market data actually changed
                            # probability won't change if we are hitting a limit order

                            # We do this before checking if bet is 10s old, so that we can stop waiting on bet information to come through after betting
                            @debug "$(now()) Notified $slug at $(marketDataBySlug[slug].probability) with $(marketLastUpdateLocal)"
                            notify(marketDataBySlug[slug].hasUpdated)

                            # Big problem if we receive bets delayed, if bet delayed by 10s just kill program
                            # Problem, when a bet is made we receive ~3 messages and the first couple may not have updated lastBetTime
                            if marketLastUpdateLocal < now(localzone()) - Second(6)
                                @debug "6s old bet $(marketLastUpdateLocal), $(now()), $(now(localzone())), $(now(localzone()) - Second(6)) for $slug at $(marketDataBySlug[slug].probability)"
                                return nothing
                                # throw(ErrorException("Got very old bet on $slug, $market at $(time())"))
                            end

                            groupIndex = groupData.contractIdToGroupIndex[marketId]

                            # oldProb is to check if market moved?
                            # lastOptimisedProb is to prevent rerunning on already optimised market, accounts for getting new market data due to our own bet
                            @debug oldProb
                            @debug marketDataBySlug[slug].lastOptimisedProb
                            @debug TaskDict[groupIndex]
                            # if marketDataBySlug[slug].probability ≉ oldProb && marketDataBySlug[slug].probability ≉ marketDataBySlug[slug].lastOptimisedProb && !TaskDict[groupIndex].runAgain

                            # When a bet occurs, two messages are sent, one with the old lastBetTime and one with the new but both have the correct probability so on second run oldProb doesn't update
                            if marketDataBySlug[slug].probability ≉ marketDataBySlug[slug].lastOptimisedProb && !TaskDict[groupIndex].runAgain

                                TaskDict[groupIndex] = (runAgain = true, task=TaskDict[groupIndex].task)

                                while !istaskdone(currentTask) # say we have 3 tasks A, B, C. A is running when B, C come in so both wait for A to finish. The B runs sets currentTask to itself but C is onlt waiting for A not the new currentTask B so need a while loop.
                                    wait(currentTask)
                                    @debug "Finished waiting for current task, $currentTask, $slug, $(marketDataBySlug[slug].probability)"
                                end

                                while !istaskdone(TaskDict[groupIndex].task)
                                    wait(TaskDict[groupIndex].task)
                                    @debug "Finished waiting for previous run of market, $(TaskDict[groupIndex].task), $slug, $(marketDataBySlug[slug].probability)"
                                end

                                # if we've already optimised on this market data
                                if !TaskDict[groupIndex].runAgain
                                    @debug "No need to rerun $slug"
                                    return nothing
                                end

                                currentTask = current_task()
                                @debug "current task set to $currentTask, $slug"
                                TaskDict[groupIndex] = (runAgain = false, task=current_task())

                                if marketDataBySlug[slug].probability ≈ marketDataBySlug[slug].lastOptimisedProb
                                    @debug "market already optimised $slug"
                                    return nothing
                                end

                                group = groupData.groups[groupIndex]

                                @debug "current prob $(marketDataBySlug[slug].probability)"
                                runGroup(group, groupData, botData, marketDataBySlug, arguments, TaskDict[groupIndex], slug)

                                @debug "Removing limit orders after running on $(group.name), $(marketDataBySlug[slug].limitOrders), $(marketDataBySlug[slug].sortedLimitProbs)"
                                for slug in group.slugs
                                    marketDataBySlug[slug].limitOrders = MutableOutcomeType(Dict{F, Vector{F}}(), Dict{F, Vector{F}}()) # need to reset as we aren't tracking limit orders
                                    marketDataBySlug[slug].sortedLimitProbs = MutableOutcomeType(F[], F[])
                                end
                                @debug "Removed limit orders after running on $(group.name), $(marketDataBySlug[slug].limitOrders), $(marketDataBySlug[slug].sortedLimitProbs)"
                            end
                        end

                        # GC.enable(true)
                        return nothing
                    end
                end

                # HeartBeat
                @async_showerr begin
                    sleep(30)

                    while !WebSockets.isclosed(socket)
                        # printstyled("Heartbeat at $(Dates.format(now(), "HH:MM:SS.sss"))\n"; color = :blue)
                        @debug "Heartbeat at $(Dates.format(now(), "HH:MM:SS.sss"))"
                        try
                            send(socket, heartbeatJSON)
                        catch
                            @error "Heartbeat Failed"
                            rethrow()
                        end
                        # printstyled("Sleeping at $(Dates.format(now(), "HH:MM:SS.sss"))\n"; color = :blue)
                        # @debug "Sleeping at $(Dates.format(now(), "HH:MM:SS.sss"))"
                        sleep(30)
                    end
                end

                # Fetch new balance every hour, hack to work around repaying loans when redeeming
                @async_showerr while !WebSockets.isclosed(socket)
                    sleep(Hour(1))

                    if !WebSockets.isclosed(socket)
                        @info "3: Fetching Balance at $(Dates.format(now(), "HH:MM:SS.sss"))"
                        botData.balance = getUserByUsername(botData.USERNAME).balance
                    else
                        break
                    end
                end

                # Fetch new balance every hour, hack to work around repaying loans when redeeming
                @async_showerr while !WebSockets.isclosed(socket)
                    if readline() == "R"
                        if WebSockets.isclosed(socket) break end
                        @info "4: Fetching Balance at $(Dates.format(now(), "HH:MM:SS.sss"))"
                        botData.balance = getUserByUsername(botData.USERNAME).balance
                        @info "4: New Balance is $(botData.balance)"
                    end
                end

                # Restart every 5 hours, hack to work around some markets not being arbed/ websocket stop sending messages
                @async_showerr while !WebSockets.isclosed(socket)
                    sleep(Hour(5))

                    if !WebSockets.isclosed(socket)
                        error("Restarting")
                    else
                        break
                    end
                end

                # Fetch new balance at 8am
                @async_showerr while !WebSockets.isclosed(socket)
                    # printstyled("Sleeping Balance at $(Dates.format(now(), "HH:MM:SS.sss"))\n"; color = :blue)
                    timeTo8 = Second(((today() + Time(8) + Minute(5) - now()) ÷ 1000).value)
                    timeTo8 += timeTo8 > Second(0) ? Second(0) : Second(Day(1))
                    sleep(timeTo8)

                    if !WebSockets.isclosed(socket)
                        printstyled("1: Fetching Balance at $(Dates.format(now(), "HH:MM:SS.sss"))\n"; color = :blue)
                        botData.balance = getUserByUsername(botData.USERNAME).balance
                    else
                        break
                    end

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
            end
        finally
            errorOccurred = true
            println("Finally")
            if !WebSockets.isclosed(socket)
                send(socket, leaveJSON())
                println("Left Channel")

                close(socket)
            end
        end
    end
end

function test(groupNames = nothing; live=false, confirmBets=true, skip=false)
    production(groupNames; live=live, confirmBets=confirmBets, skip=skip)
end

function retryProd(runs=1, groupNames = nothing; live=true, confirmBets=false, skip=false)
    delay = 60
    lastRunTime = time()
    exit = false

    for run in 1:runs
        if exit
            break
        end
        try
            println("Running at ", Dates.format(now(), "HH:MM:SS.sss"))

            production(groupNames; live=live, confirmBets=confirmBets, skip=skip)
        catch
            println("Caught at ", Dates.format(now(), "HH:MM:SS.sss"))

            for (exception, _) in current_exceptions()
                @show exception
                if shouldBreak(exception)
                    exit = true
                    break
                end
            end

            if run == runs
                exit = true
            end

            if exit
                break
            end

            if time() - lastRunTime > 60^2
                delay = 60
            end
            delay = min(delay, 60*30)
            sleep(delay)

            delay *= 5
            lastRunTime = time()
        end
    end
    println("Exited at ", Dates.format(now(), "HH:MM:SS.sss"))
end