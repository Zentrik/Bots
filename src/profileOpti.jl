using Pkg, Revise
Pkg.activate("Bots")
using Bots
# includet("$(@__DIR__)/Bots/src/Arb.jl")
# includet("Bots/src/Arb.jl")
includet("Arb.jl")

function setup()
    groupNames = ["Nuclear Weapons Detonation 2023"]
    arguments = Arguments(live=false, confirmBets=false)

    GROUPS::Dict{String, Dict{String, Vector{String}}}, APIKEY, Supabase_APIKEY, USERNAME = readData()

    if groupNames !== nothing
        GROUPS = Dict(name => GROUPS[name] for name in groupNames)
    end

    delete!(GROUPS["Nuclear Weapons Detonation 2023"], "https://manifold.markets/levifinkelstein/will-a-nuclear-weapon-be-detonated-241e2078f53d")

    groups = Group.(keys(GROUPS), values(GROUPS))

    marketDataBySlug = Dict(slug => MarketData() for slug in getSlugs(groups))

    USERID = ""
    botBalance = 1e6
    botData = BotData(APIKEY, Supabase_APIKEY, USERNAME, USERID, botBalance)

    for (i, slug) in pairs(first(groups).slugs)
        marketDataBySlug[slug].p = [0.2064787513236006, 0.47838962384408346, 0.2775101996372839, 0.2746828557228065, 0.4719656613593342, 0.35674236340711246, 0.3888933925736188, 0.19860222200238398, 0.13090253629326612, 0.22641172770985657, 0.3772018543113427][i]

        tmp = [Dict(:NO => 747.0842215717668, :YES => 627.2281129211327), Dict(:NO => 1256.8017456523282, :YES => 3349.4668405130005), Dict(:NO => 1775.7429973184794, :YES => 2007.8599621286187), Dict(:NO => 651.5863191298266, :YES => 13290.223096339212), Dict(:NO => 1275.3749538446677, :YES => 3312.527407917233), Dict(:NO => 1656.1560922691876, :YES => 2668.974243088883), Dict(:NO => 1123.5804751778835, :YES => 2077.7396315911797), Dict(:NO => 601.1182589140856, :YES => 1630.4899648820406), Dict(:NO => 438.71196606081827, :YES => 8178.8869199350875), Dict(:NO => 392.4544324342122, :YES => 1278.6356354844572), Dict(:NO => 1593.0915854731445, :YES => 2803.758652539501)][i]

        marketDataBySlug[slug].pool.YES = tmp[:YES]
        marketDataBySlug[slug].pool.NO = tmp[:NO]

        marketDataBySlug[slug].probability = [0.23659927985789392, 0.2560261355742694, 0.2535635219835542, 0.01822861476055391, 0.2560261432969136, 0.2560261339933592, 0.2560261456758842, 0.08371589271111027, 0.008014385513421455, 0.08242745421311462, 0.25602613697979976][i]

        if i in (4, 9)
            marketDataBySlug[slug].isResolved = true
        else
            marketDataBySlug[slug].isResolved = false
        end

        marketDataBySlug[slug].closeTime = round(time()*1000*2)

        marketDataBySlug[slug].url = "https://manifold.markets/test/$slug"
    end

    contractIdSet = Set(market.id for market in values(marketDataBySlug))
    contractIdToGroupIndex = Dict(marketDataBySlug[slug].id => i for (i, group) in enumerate(groups) for slug in group.slugs)
    contractIdToSlug = Dict(marketDataBySlug[slug].id => slug for group in groups for slug in group.slugs)
    groupData = GroupData(groups, contractIdSet, contractIdToGroupIndex, contractIdToSlug)

    for (i, slug) in pairs(first(groupData.groups).slugs)
        marketDataBySlug[slug].shares.YES = [91.11652268745661, 0.0, 1243.686397515011, 7.844391802791506e-12, 0.0, 24143.332383301884, 0.0, 9.919176591210999e-11, 0.0, 0.0, 6993.464526210724][i]
        marketDataBySlug[slug].shares.NO = [0.0, 836.6293889312217, 0.0, 0.0, 26630.748861967342, 0.0, 5004.22157881651, 0.0, 1.0174971976084635e-11, 9.424638847121969e-11, 0.0][i]
    end

    return groupData, botData, marketDataBySlug, arguments
end

groupData, botData, marketDataBySlug, arguments = setup()
@assert arguments.live == false

arbitrageGroup(first(groupData.groups), botData, marketDataBySlug, arguments)

@profview arbitrageGroup(first(groupData.groups), botData, marketDataBySlug, arguments)

using BenchmarkTools, Logging
global_logger(ConsoleLogger(Logging.Error))
@benchmark arbitrageGroup($(first(groupData.groups)), $botData, $marketDataBySlug, $arguments)

using StaticArrays
struct GroupStatic
    name::String
    slugs::Vector{String}
    y_matrix::SMatrix{8, 11, Bool, 88}
    n_matrix::SMatrix{8, 11, Bool, 88}
    not_na_matrix::SMatrix{8, 11, Bool, 88}
    noMarkets::Int64
end

GroupStatic((; name, slugs, y_matrix, n_matrix, not_na_matrix, noMarkets)) = GroupStatic(name, slugs, y_matrix, n_matrix, not_na_matrix, noMarkets)

@benchmark arbitrageGroup($(GroupStatic(first(groupData.groups))), $botData, $marketDataBySlug, $arguments)
@profview arbitrageGroup(GroupStatic(first(groupData.groups)), botData, marketDataBySlug, arguments)

group = first(groupData.groups)
betAmounts = F[14.742921088001445, -17.597491561818767, 0.0, 0.0, -17.30302644505758, -12.794756348500918, 0.0, 0.0, 1.0000000001389506]
bettableSlugsIndex = [i for (i, slug) in enumerate(group.slugs) if !marketDataBySlug[slug].isResolved]
noShares = [marketDataBySlug[slug].shares[:NO] for slug in group.slugs]
yesShares = [marketDataBySlug[slug].shares[:YES] for slug in group.slugs]
sharesByEvent = group.y_matrix * yesShares + group.n_matrix * noShares
profitsByEvent = similar(sharesByEvent)

f!(betAmounts, first(groupData.groups), marketDataBySlug, bettableSlugsIndex, sharesByEvent, profitsByEvent)
fNoLimit!(betAmounts, first(groupData.groups), [marketDataBySlug[slug].p for slug in group.slugs[bettableSlugsIndex]], StructArray([marketDataBySlug[slug].pool for slug in group.slugs[bettableSlugsIndex]]), bettableSlugsIndex, sharesByEvent, profitsByEvent)

@benchmark fNoLimit!($betAmounts, $(first(groupData.groups)), $([marketDataBySlug[slug].p for slug in group.slugs[bettableSlugsIndex]]), $(StructArray([marketDataBySlug[slug].pool for slug in group.slugs[bettableSlugsIndex]])), $bettableSlugsIndex, $sharesByEvent, $profitsByEvent)

@benchmark fNoLimit!($betAmounts, $(GroupStatic(first(groupData.groups))), $([marketDataBySlug[slug].p for slug in group.slugs[bettableSlugsIndex]]), $(StructArray([marketDataBySlug[slug].pool for slug in group.slugs[bettableSlugsIndex]])), $bettableSlugsIndex, $sharesByEvent, $profitsByEvent)

@benchmark f!($betAmounts, $(GroupStatic(first(groupData.groups))), $marketDataBySlug, $bettableSlugsIndex, $sharesByEvent, $profitsByEvent)

@benchmark f!($betAmounts, $((first(groupData.groups))), $marketDataBySlug, $bettableSlugsIndex, $sharesByEvent, $profitsByEvent)

using Cthulhu
@descend arbitrageGroup(first(groupData.groups), botData, marketDataBySlug, arguments)
pVec = [marketDataBySlug[slug].p for slug in group.slugs[bettableSlugsIndex]]
poolSOA = StructArray([marketDataBySlug[slug].pool for slug in group.slugs[bettableSlugsIndex]])
@descend fNoLimit!(betAmounts, group, pVec, poolSOA, bettableSlugsIndex, sharesByEvent, profitsByEvent)
# @descend minimum( fNoLimit!(betAmounts, group, pVec, poolSOA, bettableSlugsIndex, sharesByEvent, profitsByEvent))

profitF = OptimizationFunction((betAmount, _) -> fNoLimit!(betAmount, group, pVec, poolSOA, bettableSlugsIndex, sharesByEvent, profitsByEvent))
x0 = repeat([zero(F)], length(bettableSlugsIndex))
maxBetAmount = 1000.
ub = F[marketDataBySlug[slug].probability > .98 ? 0. : maxBetAmount for slug in group.slugs[bettableSlugsIndex]]
lb = F[marketDataBySlug[slug].probability < .03 ? 0. : -maxBetAmount for slug in group.slugs[bettableSlugsIndex]]

problem = Optimization.OptimizationProblem(profitF, x0, lb=lb, ub=ub)
@descend solve(problem, BBO_de_rand_1_bin_radiuslimited(), maxtime=.25)

fClosure = (betAmount, _) -> fNoLimit!(betAmount, group, pVec, poolSOA, bettableSlugsIndex, sharesByEvent, profitsByEvent)
# fClosure = (betAmount, _) -> Float64(-minimum( fNoLimit!(betAmount, group, pVec, poolSOA, bettableSlugsIndex, sharesByEvent, profitsByEvent)))

fClosure(betAmounts, "")
@descend fClosure(betAmounts, "")

fExperimental = Base.Experimental.@opaque (betAmount, _) -> fNoLimit!(betAmount, group, pVec, poolSOA, bettableSlugsIndex, sharesByEvent, profitsByEvent)

fExperimental(betAmounts, "")
@benchmark fExperimental($betAmounts, "")

@benchmark fClosure($betAmounts, "")
# @benchmark Float64(-minimum( fNoLimit!(betAmounts, group, pVec, poolSOA, bettableSlugsIndex, sharesByEvent, profitsByEvent)))
@benchmark fNoLimit!($betAmounts, $group, $pVec, $poolSOA, $bettableSlugsIndex, $sharesByEvent, $profitsByEvent)

@descend f!(betAmounts, (GroupStatic(first(groupData.groups))), marketDataBySlug, bettableSlugsIndex, sharesByEvent, profitsByEvent)

using JET

@views function f2!(betAmount, group, MarketData, bettableSlugsIndex, sharesByEvent, profitsByEvent)
    fees = 0.
    profitsByEvent .= sharesByEvent # we shouldn't count resolved markets.

    @inbounds @fastmath for (i, j) in enumerate(bettableSlugsIndex)
        slug = group.slugs[j]
        market = MarketData[slug]

        if abs(betAmount[i]) >= 1.
            shares = 0.
            @inline shares = betToShares(market.p, market.pool, market.probability, market.limitOrders, market.sortedLimitProbs, betAmount[i]).shares

            fees += FEE

            profitsByEvent .-= group.not_na_matrix[:, j] .* abs(betAmount[i])

            if betAmount[i] >= 1.
                profitsByEvent .+= group.y_matrix[:, j] .* shares
            elseif betAmount[i] <= -1.
                profitsByEvent .+= group.n_matrix[:, j] .* shares
            end
        else
            betAmount[i] = 0.
        end
    end

    profitsByEvent .-= fees

    return profitsByEvent
end
@benchmark f2!($betAmounts, $((first(groupData.groups))), $marketDataBySlug, $bettableSlugsIndex, $sharesByEvent, $profitsByEvent)

@benchmark f2!($betAmounts, $(GroupStatic(first(groupData.groups))), $marketDataBySlug, $bettableSlugsIndex, $sharesByEvent, $profitsByEvent)

using ManifoldMarkets

function f3!(betAmount, group, MarketData, bettableSlugsIndex, sharesByEvent, profitsByEvent)
    fees = 0.
    profitsByEvent .= sharesByEvent # we shouldn't count resolved markets.

    @inbounds @fastmath for (i, j) in enumerate(bettableSlugsIndex)
        slug = group.slugs[j]
        market = MarketData[slug]

        @inline shares = ifelse(abs(betAmount[i]) >= 1., ManifoldMarkets.betCPPM(market.p, market.pool, abs(betAmount[i]), betAmount[i] > 0 ? :YES : :NO).shares, 0.)

        fees += ifelse(abs(betAmount[i]) >= 1., FEE, 0.)

        # profitsByEvent .-= group.not_na_matrix[:, j] .* ifelse(abs(betAmount[i]) >= 1., abs(betAmount[i]), 0.)

        if betAmount[i] >= 1.
            profitsByEvent .+= group.y_matrix[:, j] .* shares
            profitsByEvent .-= group.not_na_matrix[:, j] .* betAmount[i]
        elseif betAmount[i] <= -1.
            profitsByEvent .+= group.n_matrix[:, j] .* shares
            profitsByEvent .-= group.not_na_matrix[:, j] .* -betAmount[i]
        end
    end

    profitsByEvent .-= fees

    return profitsByEvent
end

f3!(betAmounts, (GroupStatic(first(groupData.groups))), marketDataBySlug, bettableSlugsIndex, sharesByEvent, profitsByEvent)

@benchmark f3!($betAmounts, $((first(groupData.groups))), $marketDataBySlug, $bettableSlugsIndex, $sharesByEvent, $profitsByEvent)
@benchmark f3!($betAmounts, $(GroupStatic(first(groupData.groups))), $marketDataBySlug, $bettableSlugsIndex, $sharesByEvent, $profitsByEvent)

@descend f3!(betAmounts, (GroupStatic(first(groupData.groups))), marketDataBySlug, bettableSlugsIndex, sharesByEvent, profitsByEvent)

function f4!(betAmount, group, MarketData, bettableSlugsIndex, sharesByEvent, profitsByEvent)
    fees = 0.
    profitsByEvent .= sharesByEvent # we shouldn't count resolved markets.

    @inbounds @fastmath for (i, j) in enumerate(bettableSlugsIndex)
        market = MarketData[i]

        y = market.pool[:YES]
        n = market.pool[:NO]

        fees += ifelse(abs(betAmount[i]) >= 1., FEE, 0.)

        p = market.p
        # implement Maniswap
        if betAmount[i] >= 1
            shares = y + betAmount[i] - y * (n / (n + abs(betAmount[i])))^((1-p)/p)
        elseif betAmount[i] <= -1
            shares = n + betAmount[i] - n * (y / (y + abs(betAmount[i])))^(p/(1-p))
        end

        # profitsByEvent .-= group.not_na_matrix[:, j] .* ifelse(abs(betAmount[i]) >= 1., abs(betAmount[i]), 0.)

        if betAmount[i] >= 1.
            profitsByEvent .+= group.y_matrix[:, j] .* shares
            profitsByEvent .-= group.not_na_matrix[:, j] .* betAmount[i]
        elseif betAmount[i] <= -1.
            profitsByEvent .+= group.n_matrix[:, j] .* shares
            profitsByEvent .-= group.not_na_matrix[:, j] .* -betAmount[i]
        end
    end

    profitsByEvent .-= fees

    return profitsByEvent
end
f4!(betAmounts, (GroupStatic(first(groupData.groups))), [marketDataBySlug[slug] for slug in group.slugs[bettableSlugsIndex]], bettableSlugsIndex, sharesByEvent, profitsByEvent)

@benchmark f4!($betAmounts, $(GroupStatic(first(groupData.groups))), $([marketDataBySlug[slug] for slug in group.slugs[bettableSlugsIndex]]), $bettableSlugsIndex, $sharesByEvent, $profitsByEvent)

function f5!(betAmount, group, MarketData, bettableSlugsIndex, sharesByEvent, profitsByEvent)
    fees = 0.
    profitsByEvent .= sharesByEvent # we shouldn't count resolved markets.

    @inbounds @fastmath for (i, j) in enumerate(bettableSlugsIndex)
        market = MarketData[i]

        y = market.pool[1]
        n = market.pool[2]

        fees += ifelse(abs(betAmount[i]) >= 1., FEE, 0.)

        p = market.p
        # implement Maniswap
        if betAmount[i] >= 1
            shares = y - y * (n / (n + betAmount[i]))^((1-p)/p)
        elseif betAmount[i] <= -1
            shares = n - n * (y / (y - betAmount[i]))^(p/(1-p))
        end

        # profitsByEvent .-= group.not_na_matrix[:, j] .* ifelse(abs(betAmount[i]) >= 1., abs(betAmount[i]), 0.)

        if betAmount[i] >= 1.
            profitsByEvent .+= group.y_matrix[:, j] .* shares
            profitsByEvent .-= group.not_na_matrix[:, j] .* betAmount[i]
        elseif betAmount[i] <= -1.
            profitsByEvent .+= group.n_matrix[:, j] .* shares
            profitsByEvent .-= group.not_na_matrix[:, j] .* -betAmount[i]
        end
    end

    profitsByEvent .-= fees

    return profitsByEvent
end
mutable struct MarketDataTest
    p::Float64
    pool::SVector{2, Float64}
end
@benchmark f5!($betAmounts, $(GroupStatic(first(groupData.groups))), $([MarketDataTest(marketDataBySlug[slug].p, SVector{2, Float64}(marketDataBySlug[slug].pool[:YES], marketDataBySlug[slug].pool[:NO])) for slug in group.slugs[bettableSlugsIndex]]), $bettableSlugsIndex, $sharesByEvent, $profitsByEvent)
@descend f5!(betAmounts, (GroupStatic(first(groupData.groups))), [MarketDataTest(marketDataBySlug[slug].p, SVector{2, Float64}(marketDataBySlug[slug].pool[:YES], marketDataBySlug[slug].pool[:NO])) for slug in group.slugs[bettableSlugsIndex]], bettableSlugsIndex, sharesByEvent, profitsByEvent)

function f6!(betAmount, group, MarketData, bettableSlugsIndex, profitsByEvent)
    fees = zero(eltype(profitsByEvent))

    @inbounds @fastmath @simd for i in eachindex(bettableSlugsIndex)
        market = MarketData[i]
        j = bettableSlugsIndex[i]

        y = market.pool[1]
        n = market.pool[2]

        fees += ifelse(abs(betAmount[i]) >= 1., typeof(fees)(FEE), zero(fees))

        p = market.p

        if betAmount[i] >= 1
            shares = y + betAmount[i] - y * (n / (n + betAmount[i]))^((1-p)/p)
            profitsByEvent += group.y_matrix[:, j] * shares
            profitsByEvent -= group.not_na_matrix[:, j] * betAmount[i]
        elseif betAmount[i] <= -1
            shares = n - betAmount[i] - n * (y / (y - betAmount[i]))^(p/(1-p))
            profitsByEvent += group.n_matrix[:, j] * shares
            profitsByEvent -= group.not_na_matrix[:, j] * -betAmount[i]
        end
    end

    profitsByEvent = profitsByEvent .- fees

    return profitsByEvent
end

f6!(betAmounts, (GroupStatic(first(groupData.groups))), [MarketDataTest(marketDataBySlug[slug].p, SVector{2, Float64}(marketDataBySlug[slug].pool[:YES], marketDataBySlug[slug].pool[:NO])) for slug in group.slugs[bettableSlugsIndex]], bettableSlugsIndex, SVector{8, Float64}(sharesByEvent))
@benchmark f6!($betAmounts, $(GroupStatic(first(groupData.groups))), $([MarketDataTest(marketDataBySlug[slug].p, SVector{2, Float64}(marketDataBySlug[slug].pool[:YES], marketDataBySlug[slug].pool[:NO])) for slug in group.slugs[bettableSlugsIndex]]), $bettableSlugsIndex, $(SVector{8, Float64}(sharesByEvent)))
@code_llvm debuginfo=:none f6!(betAmounts, (GroupStatic(first(groupData.groups))), [MarketDataTest(marketDataBySlug[slug].p, SVector{2, Float64}(marketDataBySlug[slug].pool[:YES], marketDataBySlug[slug].pool[:NO])) for slug in group.slugs[bettableSlugsIndex]], bettableSlugsIndex, SVector{8, Float64}(sharesByEvent))

mutable struct MarketDataTest32
    p::Float32
    pool::SVector{2, Float32}
end
f6!(Float32.(betAmounts), (GroupStatic(first(groupData.groups))), [MarketDataTest32(marketDataBySlug[slug].p, SVector{2, Float32}(marketDataBySlug[slug].pool[:YES], marketDataBySlug[slug].pool[:NO])) for slug in group.slugs[bettableSlugsIndex]], bettableSlugsIndex, SVector{8, Float32}(sharesByEvent))
@benchmark f6!($(Float32.(betAmounts)), $(GroupStatic(first(groupData.groups))), $([MarketDataTest32(marketDataBySlug[slug].p, SVector{2, Float32}(marketDataBySlug[slug].pool[:YES], marketDataBySlug[slug].pool[:NO])) for slug in group.slugs[bettableSlugsIndex]]), $bettableSlugsIndex, $(SVector{8, Float32}(sharesByEvent)))
@code_llvm debuginfo=:none f6!(Float32.(betAmounts), (GroupStatic(first(groupData.groups))), [MarketDataTest32(marketDataBySlug[slug].p, SVector{2, Float32}(marketDataBySlug[slug].pool[:YES], marketDataBySlug[slug].pool[:NO])) for slug in group.slugs[bettableSlugsIndex]], bettableSlugsIndex, SVector{8, Float32}(sharesByEvent))
@descend f6!(Float32.(betAmounts), (GroupStatic(first(groupData.groups))), [MarketDataTest32(marketDataBySlug[slug].p, SVector{2, Float32}(marketDataBySlug[slug].pool[:YES], marketDataBySlug[slug].pool[:NO])) for slug in group.slugs[bettableSlugsIndex]], bettableSlugsIndex, SVector{8, Float32}(sharesByEvent))

@code_native syntax=:intel debuginfo=:none f6!(Float32.(betAmounts), (GroupStatic(first(groupData.groups))), [MarketDataTest32(marketDataBySlug[slug].p, SVector{2, Float32}(marketDataBySlug[slug].pool[:YES], marketDataBySlug[slug].pool[:NO])) for slug in group.slugs[bettableSlugsIndex]], bettableSlugsIndex, SVector{8, Float32}(sharesByEvent))

using StructArrays
f6!(Float32.(betAmounts), (GroupStatic(first(groupData.groups))), StructArray([MarketDataTest32(marketDataBySlug[slug].p, SVector{2, Float32}(marketDataBySlug[slug].pool[:YES], marketDataBySlug[slug].pool[:NO])) for slug in group.slugs[bettableSlugsIndex]], unwrap = T -> (T<:AbstractVector)), bettableSlugsIndex, SVector{8, Float32}(sharesByEvent))
@benchmark f6!($(Float32.(betAmounts)), $(GroupStatic(first(groupData.groups))), $(StructArray([MarketDataTest32(marketDataBySlug[slug].p, SVector{2, Float32}(marketDataBySlug[slug].pool[:YES], marketDataBySlug[slug].pool[:NO])) for slug in group.slugs[bettableSlugsIndex]], unwrap = T -> (T<:AbstractVector))), $bettableSlugsIndex, $(SVector{8, Float32}(sharesByEvent)))
@code_llvm debuginfo=:none f6!(Float32.(betAmounts), (GroupStatic(first(groupData.groups))), StructArray([MarketDataTest32(marketDataBySlug[slug].p, SVector{2, Float32}(marketDataBySlug[slug].pool[:YES], marketDataBySlug[slug].pool[:NO])) for slug in group.slugs[bettableSlugsIndex]], unwrap = T -> (T<:AbstractVector)), bettableSlugsIndex, SVector{8, Float32}(sharesByEvent))
@code_native syntax=:intel debuginfo=:none f6!(Float32.(betAmounts), (GroupStatic(first(groupData.groups))), StructArray([MarketDataTest32(marketDataBySlug[slug].p, SVector{2, Float32}(marketDataBySlug[slug].pool[:YES], marketDataBySlug[slug].pool[:NO])) for slug in group.slugs[bettableSlugsIndex]], unwrap = T -> (T<:AbstractVector)), bettableSlugsIndex, SVector{8, Float32}(sharesByEvent))

@benchmark f6!($(Float32.(betAmounts)), $(GroupStatic(first(groupData.groups))), $(StructArray([MarketDataTest32(marketDataBySlug[slug].p, SVector{2, Float32}(marketDataBySlug[slug].pool[:YES], marketDataBySlug[slug].pool[:NO])) for slug in group.slugs[bettableSlugsIndex]], unwrap = T -> (T<:AbstractVector))), $(SVector{9, Int}(bettableSlugsIndex)), $(SVector{8, Float32}(sharesByEvent)))

struct GroupStatic32
    name::String
    slugs::Vector{String}
    y_matrix::SMatrix{8, 11, Float32, 88}
    n_matrix::SMatrix{8, 11, Float32, 88}
    not_na_matrix::SMatrix{8, 11, Float32, 88}
    noMarkets::Int64
end
GroupStatic32((; name, slugs, y_matrix, n_matrix, not_na_matrix, noMarkets)) = GroupStatic32(name, slugs, y_matrix, n_matrix, not_na_matrix, noMarkets)

f6!(Float32.(betAmounts), GroupStatic32(first(groupData.groups)), StructArray([MarketDataTest32(marketDataBySlug[slug].p, SVector{2, Float32}(marketDataBySlug[slug].pool[:YES], marketDataBySlug[slug].pool[:NO])) for slug in group.slugs[bettableSlugsIndex]], unwrap = T -> (T<:AbstractVector)), bettableSlugsIndex, SVector{8, Float32}(sharesByEvent))
#BEST?
@benchmark f6!($(Float32.(betAmounts)), $(GroupStatic32(first(groupData.groups))), $(StructArray([MarketDataTest32(marketDataBySlug[slug].p, SVector{2, Float32}(marketDataBySlug[slug].pool[:YES], marketDataBySlug[slug].pool[:NO])) for slug in group.slugs[bettableSlugsIndex]], unwrap = T -> (T<:AbstractVector))), $bettableSlugsIndex, $(SVector{8, Float32}(sharesByEvent)))
@code_llvm debuginfo=:none f6!(Float32.(betAmounts), GroupStatic32(first(groupData.groups)), StructArray([MarketDataTest32(marketDataBySlug[slug].p, SVector{2, Float32}(marketDataBySlug[slug].pool[:YES], marketDataBySlug[slug].pool[:NO])) for slug in group.slugs[bettableSlugsIndex]], unwrap = T -> (T<:AbstractVector)), bettableSlugsIndex, SVector{8, Float32}(sharesByEvent))
@code_native syntax=:intel debuginfo=:none f6!(Float32.(betAmounts), GroupStatic32(first(groupData.groups)), StructArray([MarketDataTest32(marketDataBySlug[slug].p, SVector{2, Float32}(marketDataBySlug[slug].pool[:YES], marketDataBySlug[slug].pool[:NO])) for slug in group.slugs[bettableSlugsIndex]], unwrap = T -> (T<:AbstractVector)), bettableSlugsIndex, SVector{8, Float32}(sharesByEvent))

function f7!(betAmount, group, MarketData, bettableSlugsIndex, profitsByEvent)
    fees = zero(eltype(profitsByEvent))
    shares = @MVector zeros(eltype(betAmount), 11)
    bets = @MVector zeros(eltype(betAmount), 11)

    @inbounds @fastmath for (i, j) in enumerate(bettableSlugsIndex)
        market = MarketData[i]

        y = market.pool[1]
        n = market.pool[2]

        fees += ifelse(abs(betAmount[i]) >= 1., typeof(fees)(FEE), zero(fees))

        p = market.p
        # implement Maniswap
        if betAmount[i] >= 1
            shares[j] = y - y * (n / (n + betAmount[i]))^((1-p)/p)
            bets[j] = betAmount[i]
        elseif betAmount[i] <= -1
            shares[j] = n - n * (y / (y - betAmount[i]))^(p/(1-p))
            bets[j] = -betAmount[i]
        end
    end

    profitsByEvent = profitsByEvent + group.y_matrix * shares + group.n_matrix * shares + group.not_na_matrix * bets .- fees

    return profitsByEvent
end

f7!(betAmounts, (GroupStatic(first(groupData.groups))), [MarketDataTest(marketDataBySlug[slug].p, SVector{2, Float64}(marketDataBySlug[slug].pool[:YES], marketDataBySlug[slug].pool[:NO])) for slug in group.slugs[bettableSlugsIndex]], bettableSlugsIndex, SVector{8, Float64}(sharesByEvent))
@benchmark f7!($betAmounts, $(GroupStatic(first(groupData.groups))), $([MarketDataTest(marketDataBySlug[slug].p, SVector{2, Float64}(marketDataBySlug[slug].pool[:YES], marketDataBySlug[slug].pool[:NO])) for slug in group.slugs[bettableSlugsIndex]]), $bettableSlugsIndex, $(SVector{8, Float64}(sharesByEvent)))
@code_llvm debuginfo=:none f7!(betAmounts, (GroupStatic(first(groupData.groups))), [MarketDataTest(marketDataBySlug[slug].p, SVector{2, Float64}(marketDataBySlug[slug].pool[:YES], marketDataBySlug[slug].pool[:NO])) for slug in group.slugs[bettableSlugsIndex]], bettableSlugsIndex, SVector{8, Float64}(sharesByEvent))

f7!(Float32.(betAmounts), (GroupStatic(first(groupData.groups))), StructArray([MarketDataTest32(marketDataBySlug[slug].p, SVector{2, Float32}(marketDataBySlug[slug].pool[:YES], marketDataBySlug[slug].pool[:NO])) for slug in group.slugs[bettableSlugsIndex]], unwrap = T -> (T<:AbstractVector)), bettableSlugsIndex, SVector{8, Float32}(sharesByEvent))
@benchmark f7!($(Float32.(betAmounts)), $(GroupStatic(first(groupData.groups))), $([MarketDataTest(marketDataBySlug[slug].p, SVector{2, Float64}(marketDataBySlug[slug].pool[:YES], marketDataBySlug[slug].pool[:NO])) for slug in group.slugs[bettableSlugsIndex]]), $bettableSlugsIndex, $(SVector{8, Float32}(sharesByEvent)))

f7!(SVector{9, Float32}(betAmounts), (GroupStatic(first(groupData.groups))), StructArray([MarketDataTest32(marketDataBySlug[slug].p, SVector{2, Float32}(marketDataBySlug[slug].pool[:YES], marketDataBySlug[slug].pool[:NO])) for slug in group.slugs[bettableSlugsIndex]], unwrap = T -> (T<:AbstractVector)), bettableSlugsIndex, SVector{8, Float32}(sharesByEvent))
@benchmark f7!($(SVector{9, Float32}(betAmounts)), $(GroupStatic(first(groupData.groups))), $([MarketDataTest(marketDataBySlug[slug].p, SVector{2, Float64}(marketDataBySlug[slug].pool[:YES], marketDataBySlug[slug].pool[:NO])) for slug in group.slugs[bettableSlugsIndex]]), $bettableSlugsIndex, $(SVector{8, Float32}(sharesByEvent)))


function f8!(betAmount, group, MarketData, bettableSlugsIndex, profitsByEvent)
    fees = zero(eltype(profitsByEvent))

    @inbounds @fastmath for (i, j) in enumerate(bettableSlugsIndex)
        market = MarketData[i]

        y = market.pool[1]
        n = market.pool[2]

        p = market.p

        fees += ifelse(abs(betAmount[i]) >= 1., typeof(fees)(FEE), zero(fees))

        profitsByEvent += ifelse(betAmount[i] >= 1, group.y_matrix[:, j] * (y - y * (n / (n + betAmount[i]))^((1-p)/p)), ifelse(betAmount[i] <= 1, -group.n_matrix[:, j] * (n - n * (y / (y - betAmount[i]))^(p/(1-p))), zero(profitsByEvent)))

        profitsByEvent -= group.not_na_matrix[:, j] * abs(betAmount[i])
    end

    return profitsByEvent .- fees
end

@descend f8!(Float32.(betAmounts), GroupStatic32(first(groupData.groups)), StructArray([MarketDataTest32(marketDataBySlug[slug].p, SVector{2, Float32}(marketDataBySlug[slug].pool[:YES], marketDataBySlug[slug].pool[:NO])) for slug in group.slugs[bettableSlugsIndex]], unwrap = T -> (T<:AbstractVector)), bettableSlugsIndex, SVector{8, Float32}(sharesByEvent))
@benchmark f8!($(Float32.(betAmounts)), $(GroupStatic32(first(groupData.groups))), $(StructArray([MarketDataTest32(marketDataBySlug[slug].p, SVector{2, Float32}(marketDataBySlug[slug].pool[:YES], marketDataBySlug[slug].pool[:NO])) for slug in group.slugs[bettableSlugsIndex]], unwrap = T -> (T<:AbstractVector))), $bettableSlugsIndex, $(SVector{8, Float32}(sharesByEvent)))
@code_llvm debuginfo=:none f8!(Float32.(betAmounts), GroupStatic32(first(groupData.groups)), StructArray([MarketDataTest32(marketDataBySlug[slug].p, SVector{2, Float32}(marketDataBySlug[slug].pool[:YES], marketDataBySlug[slug].pool[:NO])) for slug in group.slugs[bettableSlugsIndex]], unwrap = T -> (T<:AbstractVector)), bettableSlugsIndex, SVector{8, Float32}(sharesByEvent))


@benchmark f8!($(SVector{9, Float32}(betAmounts)), $(GroupStatic32(first(groupData.groups))), $(StructArray([MarketDataTest32(marketDataBySlug[slug].p, SVector{2, Float32}(marketDataBySlug[slug].pool[:YES], marketDataBySlug[slug].pool[:NO])) for slug in group.slugs[bettableSlugsIndex]], unwrap = T -> (T<:AbstractVector))), $bettableSlugsIndex, $(SVector{8, Float32}(sharesByEvent)))

using LoopVectorization
function f9!(betAmount, group, MarketData, bettableSlugsIndex, profitsByEvent)
    for i in eachindex(bettableSlugsIndex)
        j = bettableSlugsIndex[i]
        market = MarketData[i]

        y = market.pool[1]
        n = market.pool[2]

        p = market.p

        profitsByEvent += ifelse(betAmount[i] >= 1, group.y_matrix[:, j] * (y - y * (n / (n + betAmount[i]))^((1-p)/p)), ifelse(betAmount[i] <= 1, -group.n_matrix[:, j] * (n - n * (y / (y - betAmount[i]))^(p/(1-p))), zero(betAmount[i])))

        profitsByEvent -= group.not_na_matrix[:, j] * abs(betAmount[i])
    end

    return profitsByEvent
end

f9!(Float32.(betAmounts), GroupStatic32(first(groupData.groups)), StructArray([MarketDataTest32(marketDataBySlug[slug].p, SVector{2, Float32}(marketDataBySlug[slug].pool[:YES], marketDataBySlug[slug].pool[:NO])) for slug in group.slugs[bettableSlugsIndex]], unwrap = T -> (T<:AbstractVector)), bettableSlugsIndex, SVector{8, Float32}(sharesByEvent))
@benchmark f9!($(Float32.(betAmounts)), $(GroupStatic32(first(groupData.groups))), $(StructArray([MarketDataTest32(marketDataBySlug[slug].p, SVector{2, Float32}(marketDataBySlug[slug].pool[:YES], marketDataBySlug[slug].pool[:NO])) for slug in group.slugs[bettableSlugsIndex]], unwrap = T -> (T<:AbstractVector))), $bettableSlugsIndex, $(SVector{8, Float32}(sharesByEvent)))
@code_llvm debuginfo=:none f8!(Float32.(betAmounts), GroupStatic32(first(groupData.groups)), StructArray([MarketDataTest32(marketDataBySlug[slug].p, SVector{2, Float32}(marketDataBySlug[slug].pool[:YES], marketDataBySlug[slug].pool[:NO])) for slug in group.slugs[bettableSlugsIndex]], unwrap = T -> (T<:AbstractVector)), bettableSlugsIndex, SVector{8, Float32}(sharesByEvent))


@benchmark f9!($(SVector{9, Float32}(betAmounts)), $(GroupStatic32(first(groupData.groups))), $(StructArray([MarketDataTest32(marketDataBySlug[slug].p, SVector{2, Float32}(marketDataBySlug[slug].pool[:YES], marketDataBySlug[slug].pool[:NO])) for slug in group.slugs[bettableSlugsIndex]], unwrap = T -> (T<:AbstractVector))), $bettableSlugsIndex, $(SVector{8, Float32}(sharesByEvent)))

function f10!(betAmount, group, MarketData, bettableSlugsIndex, profitsByEvent)
    @inbounds @fastmath for (i, j) in enumerate(bettableSlugsIndex)
        market = MarketData[i]

        y = market.pool[1]
        n = market.pool[2]

        p = market.p

        profitsByEvent += ifelse(betAmount[i] >= 1, group.y_matrix[:, j] * (y - y * (n / (n + betAmount[i]))^((1-p)/p)) .- eltype(profitsByEvent)(FEE), ifelse(betAmount[i] <= 1, -group.n_matrix[:, j] * (n - n * (y / (y - betAmount[i]))^(p/(1-p))) .- eltype(profitsByEvent)(FEE), zero(profitsByEvent)))
    end

    return profitsByEvent
end

f10!(Float32.(betAmounts), GroupStatic32(first(groupData.groups)), StructArray([MarketDataTest32(marketDataBySlug[slug].p, SVector{2, Float32}(marketDataBySlug[slug].pool[:YES], marketDataBySlug[slug].pool[:NO])) for slug in group.slugs[bettableSlugsIndex]], unwrap = T -> (T<:AbstractVector)), bettableSlugsIndex, SVector{8, Float32}(sharesByEvent))
@benchmark f10!($(Float32.(betAmounts)), $(GroupStatic32(first(groupData.groups))), $(StructArray([MarketDataTest32(marketDataBySlug[slug].p, SVector{2, Float32}(marketDataBySlug[slug].pool[:YES], marketDataBySlug[slug].pool[:NO])) for slug in group.slugs[bettableSlugsIndex]], unwrap = T -> (T<:AbstractVector))), $bettableSlugsIndex, $(SVector{8, Float32}(sharesByEvent)))
@code_llvm debuginfo=:none f10!(Float32.(betAmounts), GroupStatic32(first(groupData.groups)), StructArray([MarketDataTest32(marketDataBySlug[slug].p, SVector{2, Float32}(marketDataBySlug[slug].pool[:YES], marketDataBySlug[slug].pool[:NO])) for slug in group.slugs[bettableSlugsIndex]], unwrap = T -> (T<:AbstractVector)), bettableSlugsIndex, SVector{8, Float32}(sharesByEvent))

@benchmark f10!($(SVector{9, Float32}(betAmounts)), $(GroupStatic32(first(groupData.groups))), $(StructArray([MarketDataTest32(marketDataBySlug[slug].p, SVector{2, Float32}(marketDataBySlug[slug].pool[:YES], marketDataBySlug[slug].pool[:NO])) for slug in group.slugs[bettableSlugsIndex]], unwrap = T -> (T<:AbstractVector))), $bettableSlugsIndex, $(SVector{8, Float32}(sharesByEvent)))

function f6Test!(betAmount, y_matrix, n_matrix, not_na_matrix, pVec, poolYES, poolNO, profitsByEvent)
    fees = zero(eltype(profitsByEvent))

    @inbounds @fastmath @simd for i in eachindex(betAmount)
        y = poolYES[i]
        n = poolNO[i]

        fees += ifelse(abs(betAmount[i]) >= 1., typeof(fees)(FEE), zero(fees))

        p = pVec[i]

        if betAmount[i] >= 1
            shares = y + betAmount[i] - y * (n / (n + betAmount[i]))^((1-p)/p)
            profitsByEvent += y_matrix[:, i] * shares
            profitsByEvent -= not_na_matrix[:, i] * betAmount[i]
        elseif betAmount[i] <= -1
            shares = n + betAmount[i] - n * (y / (y - betAmount[i]))^(p/(1-p))
            profitsByEvent += n_matrix[:, i] * shares
            profitsByEvent -= not_na_matrix[:, i] * -betAmount[i]
        end
    end

    profitsByEvent = profitsByEvent .- fees

    return profitsByEvent
end

tmp = StructArray([MarketDataTest32(marketDataBySlug[slug].p, (marketDataBySlug[slug].pool[:YES], marketDataBySlug[slug].pool[:NO])) for slug in group.slugs[bettableSlugsIndex]], unwrap = T -> (T<:AbstractVector))

f6Test!(Float32.(betAmounts), GroupStatic(first(groupData.groups)).y_matrix[:, bettableSlugsIndex], GroupStatic(first(groupData.groups)).n_matrix[:, bettableSlugsIndex], GroupStatic(first(groupData.groups)).not_na_matrix[:, bettableSlugsIndex], tmp.p, StructArrays.components(tmp.pool)[1], StructArrays.components(tmp.pool)[2], SVector{8, Float32}(sharesByEvent))
@benchmark f6Test!($(Float32.(betAmounts)), $(SMatrix{8, 9, Bool}(GroupStatic(first(groupData.groups)).y_matrix[:, bettableSlugsIndex])), $(SMatrix{8, 9, Bool}(GroupStatic(first(groupData.groups)).n_matrix[:, bettableSlugsIndex])), $(SMatrix{8, 9, Bool}(GroupStatic(first(groupData.groups)).not_na_matrix[:, bettableSlugsIndex])), $tmp.p, $(StructArrays.components(tmp.pool)[1]), $(StructArrays.components(tmp.pool)[2]), $(SVector{8, Float32}(sharesByEvent)))

@code_warntype f6Test!(Float32.(betAmounts), GroupStatic(first(groupData.groups)).y_matrix[:, bettableSlugsIndex], GroupStatic(first(groupData.groups)).n_matrix[:, bettableSlugsIndex], GroupStatic(first(groupData.groups)).not_na_matrix[:, bettableSlugsIndex], tmp.p, StructArrays.components(tmp.pool)[1], StructArrays.components(tmp.pool)[2], SVector{8, Float32}(sharesByEvent))


using SIMD
function f6SIMD!(betAmount, y_matrix, n_matrix, not_na_matrix, pVec, poolYES, poolNO, profitsByEvent)
    @inbounds @fastmath for i in LoopVecRange{8}(betAmount, unsafe=true)

        y = poolYES[i]
        n = poolNO[i]

        p = pVec[i]

        if any(betAmount[i] >= 1)
            shares = vifelse(betAmount[i] >= 1, y + betAmount[i] - y * (n / (n + betAmount[i]))^((1-p)/p), 0)
            for j in 1:8
                profitsByEvent += y_matrix[:, i.i + j - 1] * shares[j]
                profitsByEvent -= ifelse(betAmount[i.i+j-1] <= -1, not_na_matrix[:, i.i+j-1] * betAmount[i.i+j-1], zeros(Float32, 8))
            end
        elseif any(betAmount[i] <= -1)
            shares = vifelse(betAmount[i] <= -1, n + betAmount[i] - n * (y / (y - betAmount[i]))^(p/(1-p)), 0)
            for j in 1:8
                profitsByEvent += n_matrix[:, i.i+j-1] * shares[j]
                profitsByEvent -= ifelse(betAmount[i.i+j-1] <= -1, not_na_matrix[:, i.i+j-1] * -betAmount[i.i+j-1], zeros(Float32, 8))
            end
        end
    end

    return profitsByEvent
end

tmp = StructArray([MarketDataTest32(marketDataBySlug[slug].p, (marketDataBySlug[slug].pool[:YES], marketDataBySlug[slug].pool[:NO])) for slug in group.slugs[bettableSlugsIndex]], unwrap = T -> (T<:AbstractVector))

f6SIMD!(Float32.(betAmounts), first(groupData.groups).y_matrix, first(groupData.groups).n_matrix, first(groupData.groups).not_na_matrix, tmp.p, StructArrays.components(tmp.pool)..., SVector{8, Float32}(sharesByEvent))
@benchmark f6SIMD!($(Float32.(betAmounts)), $(first(groupData.groups).y_matrix), $(first(groupData.groups).n_matrix), $(first(groupData.groups).not_na_matrix), $tmp.p, $(StructArrays.components(tmp.pool)[1]), $(StructArrays.components(tmp.pool)[2]), $(SVector{8, Float32}(sharesByEvent)))

@code_warntype f6SIMD!(Float32.(betAmounts), first(groupData.groups).y_matrix, first(groupData.groups).n_matrix, first(groupData.groups).not_na_matrix, tmp.p, StructArrays.components(tmp.pool)..., SVector{8, Float32}(sharesByEvent))

@benchmark f6SIMD!($(Float32.(betAmounts)), $(GroupStatic(first(groupData.groups)).y_matrix), $(GroupStatic(first(groupData.groups)).n_matrix), $(GroupStatic(first(groupData.groups)).not_na_matrix), $(StructArray([MarketDataTest32(marketDataBySlug[slug].p, SVector{2, Float32}(marketDataBySlug[slug].pool[:YES], marketDataBySlug[slug].pool[:NO])) for slug in group.slugs[bettableSlugsIndex]], unwrap = T -> (T<:AbstractVector))), $bettableSlugsIndex, $(SVector{8, Float32}(sharesByEvent)))
@code_llvm debuginfo=:none f6!(Float32.(betAmounts), (GroupStatic(first(groupData.groups))), StructArray([MarketDataTest32(marketDataBySlug[slug].p, SVector{2, Float32}(marketDataBySlug[slug].pool[:YES], marketDataBySlug[slug].pool[:NO])) for slug in group.slugs[bettableSlugsIndex]], unwrap = T -> (T<:AbstractVector)), bettableSlugsIndex, SVector{8, Float32}(sharesByEvent))
@code_native syntax=:intel debuginfo=:none f6!(Float32.(betAmounts), (GroupStatic(first(groupData.groups))), StructArray([MarketDataTest32(marketDataBySlug[slug].p, SVector{2, Float32}(marketDataBySlug[slug].pool[:YES], marketDataBySlug[slug].pool[:NO])) for slug in group.slugs[bettableSlugsIndex]], unwrap = T -> (T<:AbstractVector)), bettableSlugsIndex, SVector{8, Float32}(sharesByEvent))

function f11!(betAmount, y_matrix, n_matrix, not_na_matrix, pVec, poolYES, poolNO, profitsByEvent)
    fees = zero(eltype(profitsByEvent))

    @turbo for i in eachindex(betAmount)
        y = poolYES[i]
        n = poolNO[i]

        # fees += ifelse(abs(betAmount[i]) >= 1., typeof(fees)(FEE), zero(fees))

        p = pVec[i]

        yesShares = y + betAmount[i] - y * (n / (n + betAmount[i]))^((1-p)/p)
        noShares = n + betAmount[i] - n * (y / (y - betAmount[i]))^(p/(1-p))

        profitsByEvent += ifelse(betAmount[i] >= 1,
                            y_matrix[:, i] * yesShares,
                        ifelse(betAmount[i] <= 1,
                            n_matrix[:, i] * noShares,
                        zero(betAmount[i])))

        profitsByEvent -= not_na_matrix[:, i] * abs(betAmount[i])
    end

    profitsByEvent = profitsByEvent .- fees

    return profitsByEvent
end

tmp = StructArray([MarketDataTest32(marketDataBySlug[slug].p, (marketDataBySlug[slug].pool[:YES], marketDataBySlug[slug].pool[:NO])) for slug in group.slugs[bettableSlugsIndex]], unwrap = T -> (T<:AbstractVector))

f11!(Float32.(betAmounts), GroupStatic(first(groupData.groups)).y_matrix[:, bettableSlugsIndex], GroupStatic(first(groupData.groups)).n_matrix[:, bettableSlugsIndex], GroupStatic(first(groupData.groups)).not_na_matrix[:, bettableSlugsIndex], tmp.p, StructArrays.components(tmp.pool)[1], StructArrays.components(tmp.pool)[2], SVector{8, Float32}(sharesByEvent))
@benchmark f11!($(Float32.(betAmounts)), $(SMatrix{8, 9, Bool}(GroupStatic(first(groupData.groups)).y_matrix[:, bettableSlugsIndex])), $(SMatrix{8, 9, Bool}(GroupStatic(first(groupData.groups)).n_matrix[:, bettableSlugsIndex])), $(SMatrix{8, 9, Bool}(GroupStatic(first(groupData.groups)).not_na_matrix[:, bettableSlugsIndex])), $tmp.p, $(StructArrays.components(tmp.pool)[1]), $(StructArrays.components(tmp.pool)[2]), $(SVector{8, Float32}(sharesByEvent)))

function f11Test!(betAmount, y_matrix, n_matrix, not_na_matrix, pVec, poolYES, poolNO, profitsByEvent)
    @turbo for i in eachindex(betAmount)
        y = poolYES[i]
        n = poolNO[i]
        p = pVec[i]

        yesShares = y + betAmount[i] - y * (n / (n + betAmount[i]))^((1-p)/p)
        noShares = n + betAmount[i] - n * (y / (y - betAmount[i]))^(p/(1-p))

        maskedYesShares = ifelse(betAmount[i] >= 1, yesShares, zero(yesShares))
        maskedNoShares = ifelse(betAmount[i] <= -1, noShares, zero(noShares))

        for j in eachindex(profitsByEvent)
            profitsByEvent[j] += y_matrix[j, i] * maskedYesShares
            profitsByEvent[j] += n_matrix[j, i] * maskedNoShares
            profitsByEvent[j] -= not_na_matrix[j, i] * abs(betAmount[i])
        end
    end

    return profitsByEvent
end

tmp2 = GroupStatic32(first(groupData.groups)).y_matrix[:, bettableSlugsIndex]
tmp3 = GroupStatic32(first(groupData.groups)).n_matrix[:, bettableSlugsIndex]
tmp4 = GroupStatic32(first(groupData.groups)).not_na_matrix[:, bettableSlugsIndex]
f11Test!(Float32.(betAmounts), tmp2, tmp3, tmp4, tmp.p, StructArrays.components(tmp.pool)[1], StructArrays.components(tmp.pool)[2], zeros(Float32, 8))

@benchmark f11Test!($(Float32.(betAmounts)), $tmp2, $tmp3, $tmp4, $tmp.p, $(StructArrays.components(tmp.pool)[1]), $(StructArrays.components(tmp.pool)[2]), $(zeros(Float32, 8)))

@views @fastmath function f6Simple!(betAmount, group, MarketData, bettableSlugsIndex, sharesByEvent, profitsByEvent)
    fees = zero(eltype(profitsByEvent))
    profitsByEvent .= sharesByEvent

    @inbounds for i in eachindex(bettableSlugsIndex)
        market = MarketData[i]
        j = bettableSlugsIndex[i]

        y = market.YES[i]
        n = market.NO[i]

        fees += ifelse(abs(betAmount[i]) >= 1., typeof(fees)(FEE), zero(fees))

        p = market.p

        if betAmount[i] >= 1
            shares = y + betAmount[i] - y * (n / (n + betAmount[i]))^((1-p)/p)
            @. profitsByEvent += group.y_matrix[:, j] * shares
            @. profitsByEvent -= group.not_na_matrix[:, j] * betAmount[i]
        elseif betAmount[i] <= -1
            shares = n + betAmount[i] - n * (y / (y - betAmount[i]))^(p/(1-p))
            @. profitsByEvent += group.n_matrix[:, j] * shares
            @. profitsByEvent -= group.not_na_matrix[:, j] * -betAmount[i]
        end
    end

    profitsByEvent .-= fees

    return profitsByEvent
end
mutable struct MarketDataTest323
    p::Float32
    YES::Float32
    NO::Float32
end

f6Simple!(Float32.(betAmounts), first(groupData.groups), StructArray([MarketDataTest323(marketDataBySlug[slug].p, marketDataBySlug[slug].pool[:YES], marketDataBySlug[slug].pool[:NO]) for slug in group.slugs[bettableSlugsIndex]]), bettableSlugsIndex, sharesByEvent, similar(sharesByEvent))

@benchmark f6Simple!($(Float32.(betAmounts)), $(first(groupData.groups)), $(StructArray([MarketDataTest323(marketDataBySlug[slug].p, marketDataBySlug[slug].pool[:YES], marketDataBySlug[slug].pool[:NO]) for slug in group.slugs[bettableSlugsIndex]])), $bettableSlugsIndex, $sharesByEvent, $(similar(sharesByEvent)))

@benchmark f6Simple!($(Float32.(betAmounts)), $(first(groupData.groups)), $([MarketDataTest323(marketDataBySlug[slug].p, marketDataBySlug[slug].pool[:YES], marketDataBySlug[slug].pool[:NO]) for slug in group.slugs[bettableSlugsIndex]]), $bettableSlugsIndex, $sharesByEvent, $(similar(sharesByEvent)))

@benchmark f6Simple!($(Float32.(betAmounts)), $(GroupStatic32(first(groupData.groups))), $(StructArray([MarketDataTest323(marketDataBySlug[slug].p, marketDataBySlug[slug].pool[:YES], marketDataBySlug[slug].pool[:NO]) for slug in group.slugs[bettableSlugsIndex]])), $bettableSlugsIndex, $sharesByEvent, $(similar(sharesByEvent)))

@benchmark f6Simple!($(Float32.(betAmounts)), $(GroupStatic32(first(groupData.groups))), $([MarketDataTest323(marketDataBySlug[slug].p, marketDataBySlug[slug].pool[:YES], marketDataBySlug[slug].pool[:NO]) for slug in group.slugs[bettableSlugsIndex]]), $bettableSlugsIndex, $sharesByEvent, $(similar(sharesByEvent)))

@benchmark f6Simple!($(Float32.(betAmounts)), $(first(groupData.groups)), $([MarketDataTest323(marketDataBySlug[slug].p, marketDataBySlug[slug].pool[:YES], marketDataBySlug[slug].pool[:NO]) for slug in group.slugs[bettableSlugsIndex]]), $bettableSlugsIndex, $sharesByEvent, $(similar(sharesByEvent)))

@code_warntype f6Simple!(Float32.(betAmounts), first(groupData.groups), StructArray([MarketDataTest323(marketDataBySlug[slug].p, marketDataBySlug[slug].pool[:YES], marketDataBySlug[slug].pool[:NO]) for slug in group.slugs[bettableSlugsIndex]]), bettableSlugsIndex, sharesByEvent, similar(sharesByEvent))


function f6SimpleStatic!(betAmount, group, MarketData, bettableSlugsIndex, profitsByEvent)
    fees = zero(eltype(profitsByEvent))

    @views @inbounds @fastmath for i in eachindex(bettableSlugsIndex)
        market = MarketData[i]
        j = bettableSlugsIndex[i]

        y = market.YES[i]
        n = market.NO[i]

        fees += ifelse(abs(betAmount[i]) >= 1., typeof(fees)(FEE), zero(fees))

        p = market.p

        if betAmount[i] >= 1
            shares = y + betAmount[i] - y * (n / (n + betAmount[i]))^((1-p)/p)
            profitsByEvent += group.y_matrix[:, j] * shares
            profitsByEvent -= group.not_na_matrix[:, j] * betAmount[i]
        elseif betAmount[i] <= -1
            shares = n + betAmount[i] - n * (y / (y - betAmount[i]))^(p/(1-p))
            profitsByEvent += group.n_matrix[:, j] * shares
            profitsByEvent -= group.not_na_matrix[:, j] * -betAmount[i]
        end
    end

    profitsByEvent = profitsByEvent .- fees

    return profitsByEvent
end

@benchmark f6SimpleStatic!($(Float32.(betAmounts)), $(GroupStatic32(first(groupData.groups))), $(StructArray([MarketDataTest323(marketDataBySlug[slug].p, marketDataBySlug[slug].pool[:YES], marketDataBySlug[slug].pool[:NO]) for slug in group.slugs[bettableSlugsIndex]])), $bettableSlugsIndex, $(SVector{8, Float32}(sharesByEvent)))

@benchmark f6SimpleStatic!($(Float32.(betAmounts)), $(first(groupData.groups)), $(StructArray([MarketDataTest323(marketDataBySlug[slug].p, marketDataBySlug[slug].pool[:YES], marketDataBySlug[slug].pool[:NO]) for slug in group.slugs[bettableSlugsIndex]])), $bettableSlugsIndex, $(SVector{8, Float32}(sharesByEvent)))

@benchmark f6SimpleStatic!($(Float32.(betAmounts)), $(GroupStatic32(first(groupData.groups))), $(StructArray([MarketDataTest323(marketDataBySlug[slug].p, marketDataBySlug[slug].pool[:YES], marketDataBySlug[slug].pool[:NO]) for slug in group.slugs[bettableSlugsIndex]])), $bettableSlugsIndex, $sharesByEvent)

@benchmark f6SimpleStatic!($(Float32.(betAmounts)), $(first(groupData.groups)), $(StructArray([MarketDataTest323(marketDataBySlug[slug].p, marketDataBySlug[slug].pool[:YES], marketDataBySlug[slug].pool[:NO]) for slug in group.slugs[bettableSlugsIndex]])), $bettableSlugsIndex, $sharesByEvent)

using SIMD
tmp = StructArray([MarketDataTest32(marketDataBySlug[slug].p, (marketDataBySlug[slug].pool[:YES], marketDataBySlug[slug].pool[:NO])) for slug in group.slugs[bettableSlugsIndex]], unwrap = T -> (T<:AbstractVector))

function f6SIMD!(betAmount, y_matrix, n_matrix, not_na_matrix, pVec, poolYES, poolNO, sharesByEvent, profitsByEvent)
    profitsByEvent .= sharesByEvent

    @inbounds @fastmath for i in LoopVecRange{8}(betAmount, unsafe=true)
        y = poolYES[i]
        n = poolNO[i]

        p = pVec[i]

        sharesYES = vifelse(betAmount[i] >= 1, y + betAmount[i] - y * (n / (n + betAmount[i]))^((1-p)/p), 0)

        sharesNO = vifelse(betAmount[i] <= -1, n - betAmount[i] - n * (y / (y - betAmount[i]))^(p/(1-p)), 0)

        @inbounds @fastmath for j in eachindex(sharesByEvent)
            profitsByEvent[i] += y_matrix[j*size(y_matrix, 1) + i] * sharesYES
            profitsByEvent[i] += n_matrix[j*size(y_matrix, 1) + i] * sharesNO
        end
    end

    return profitsByEvent
end

f6SIMD!(Float32.(betAmounts[1:8]), vec(first(groupData.groups).y_matrix), vec(first(groupData.groups).n_matrix), first(groupData.groups).not_na_matrix, tmp.p, StructArrays.components(tmp.pool)..., SVector{8, Float32}(sharesByEvent), zero(sharesByEvent))
@benchmark f6SIMD!($(Float32.(betAmounts[1:8])), $(first(groupData.groups).y_matrix), $(first(groupData.groups).n_matrix), $(first(groupData.groups).not_na_matrix), $tmp.p, $(StructArrays.components(tmp.pool)[1]), $(StructArrays.components(tmp.pool)[2]), $(sharesByEvent), $(zero(sharesByEvent)))

@code_native syntax=:intel debuginfo=:none f6SIMD!(Float32.(betAmounts[1:8]), vec(first(groupData.groups).y_matrix), first(groupData.groups).n_matrix, first(groupData.groups).not_na_matrix, tmp.p, StructArrays.components(tmp.pool)..., sharesByEvent, zero(sharesByEvent))

using LoopVectorization
function f91!(betAmount, y_matrix, n_matrix, not_na_matrix, pVec, poolYES, poolNO, sharesByEvent, profitsByEvent, sharesYES, sharesNO)
    profitsByEvent .= sharesByEvent
    sharesYES .= 0
    sharesNO .= 0

    fees = zero(eltype(profitsByEvent))

    @turbo for i in eachindex(betAmount)
        y = poolYES[i]
        n = poolNO[i]

        p = pVec[i]

        sharesYES[i] = ifelse(betAmount[i] >= 1, y + betAmount[i] - y * (n / (n + betAmount[i]))^((1-p)/p), 0)
        sharesNO[i] = ifelse(betAmount[i] <= -1, n - betAmount[i] - n * (y / (y - betAmount[i]))^(p/(1-p)), 0)

        for j in axes(n_matrix, 1)
            profitsByEvent[j] += y_matrix[j, i] * sharesYES[i]
            profitsByEvent[j] += n_matrix[j, i] * sharesNO[i]
            profitsByEvent[j] -= not_na_matrix[j, i] * abs(betAmount[i])
            # profitsByEvent[j] -= abs(betAmount[i])
        end

        fees += sum(abs(betAmount[i]) >= 1)
    end

    return @turbo minimum(profitsByEvent) - fees  * FEE
end

f91!(Float32.(betAmounts), first(groupData.groups).y_matrix[:, bettableSlugsIndex], first(groupData.groups).n_matrix[:, bettableSlugsIndex], first(groupData.groups).not_na_matrix[:, bettableSlugsIndex], tmp.p, StructArrays.components(tmp.pool)..., sharesByEvent, zero(sharesByEvent), zero(betAmounts), zero(betAmounts))
@benchmark f91!($(Float32.(betAmounts)), $(first(groupData.groups).y_matrix[:, bettableSlugsIndex]), $(first(groupData.groups).n_matrix[:, bettableSlugsIndex]), $(first(groupData.groups).not_na_matrix[:, bettableSlugsIndex]), $tmp.p, $(StructArrays.components(tmp.pool)[1]), $(StructArrays.components(tmp.pool)[2]), $(sharesByEvent), $(zero(sharesByEvent)), $(zero(betAmounts)), $(zero(betAmounts)))

@benchmark f91!($(Float32.(betAmounts[1:8])), $(first(groupData.groups).y_matrix), $(first(groupData.groups).n_matrix), $(SparseMatrixCSC{Bool, Int32}(1 .- first(groupData.groups).not_na_matrix)), $tmp.p, $(StructArrays.components(tmp.pool)[1]), $(StructArrays.components(tmp.pool)[2]), $(sharesByEvent), $(zero(sharesByEvent)), $(zero(betAmounts)), $(zero(betAmounts)))
@benchmark f91!($(Float32.(betAmounts[1:8])), $(first(groupData.groups).y_matrix[:, bettableSlugsIndex]), $(first(groupData.groups).n_matrix[:, bettableSlugsIndex]), $(first(groupData.groups).not_na_matrix[:, bettableSlugsIndex]), $tmp.p, $(StructArrays.components(tmp.pool)[1]), $(StructArrays.components(tmp.pool)[2]), $(sharesByEvent), $(zero(sharesByEvent)), $(zero(betAmounts)), $(zero(betAmounts)))

@benchmark f91!($(Float32.(betAmounts[1:8])), $(first(groupData.groups).y_matrix), $(first(groupData.groups).n_matrix), $(BitMatrix(first(groupData.groups).not_na_matrix)), $tmp.p, $(StructArrays.components(tmp.pool)[1]), $(StructArrays.components(tmp.pool)[2]), $(sharesByEvent), $(zero(sharesByEvent)), $(zero(betAmounts)), $(zero(betAmounts)))

@benchmark f91!($(Float32.(betAmounts[1:8])), $(BitMatrix(first(groupData.groups).y_matrix)), $(BitMatrix(first(groupData.groups).n_matrix)), $(BitMatrix(first(groupData.groups).not_na_matrix)), $tmp.p, $(StructArrays.components(tmp.pool)[1]), $(StructArrays.components(tmp.pool)[2]), $(sharesByEvent), $(zero(sharesByEvent)), $(zero(betAmounts)), $(zero(betAmounts)))


@code_native debuginfo=:none syntax=:intel f9!(Float32.(betAmounts), first(groupData.groups).y_matrix, first(groupData.groups).n_matrix, first(groupData.groups).not_na_matrix, tmp.p, StructArrays.components(tmp.pool)..., sharesByEvent, zero(sharesByEvent), zero(betAmounts), zero(betAmounts))

function powTest(betAmount, y_matrix, n_matrix, not_na_matrix, pVec, poolYES, poolNO, sharesByEvent, profitsByEvent)
    shares = zero(betAmount)
    @turbo for i in eachindex(betAmount)
        shares[i] = ifelse(betAmount[i] >= 1, poolYES[i] + betAmount[i] - poolYES[i] * (poolNO[i]/(poolNO[i] + betAmount[i]))^((1-pVec[i])/pVec[i]), 0)

        for j in axes(y_matrix, 1)
            profitsByEvent[j] += shares[i] # y_matrix[j, i] * shares[i]
        end
    end

    return profitsByEvent
end

powTest(Float32.(betAmounts), first(groupData.groups).y_matrix, first(groupData.groups).n_matrix, first(groupData.groups).not_na_matrix, tmp.p, StructArrays.components(tmp.pool)..., SVector{8, Float32}(sharesByEvent), zero(sharesByEvent))
@benchmark powTest($(Float32.(betAmounts)), $(first(groupData.groups).y_matrix), $(first(groupData.groups).n_matrix), $(first(groupData.groups).not_na_matrix), $tmp.p, $(StructArrays.components(tmp.pool)[1]), $(StructArrays.components(tmp.pool)[2]), $(sharesByEvent), $(zero(sharesByEvent)))
@code_native debuginfo=:none syntax=:intel powTest(Float32.(betAmounts), first(groupData.groups).y_matrix, first(groupData.groups).n_matrix, first(groupData.groups).not_na_matrix, tmp.p, StructArrays.components(tmp.pool)..., SVector{8, Float32}(sharesByEvent), zero(sharesByEvent))

using LoopVectorization
function f92!(betAmount, y_matrix, n_matrix, not_na_matrix, pVec, poolYES, poolNO, sharesByEvent, sharesYES, sharesNO)
    sharesYES .= 0
    sharesNO .= 0

    fees = zero(eltype(sharesByEvent))

    @turbo for i in eachindex(betAmount)
        y = poolYES[i]
        n = poolNO[i]

        p = pVec[i]

        sharesYES[i] = ifelse(betAmount[i] >= 1, y + betAmount[i] - y * (n / (n + betAmount[i]))^((1-p)/p), 0)
        sharesNO[i] = ifelse(betAmount[i] <= -1, n - betAmount[i] - n * (y / (y - betAmount[i]))^(p/(1-p)), 0)

        fees += sum(abs(betAmount[i]) >= 1)
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

    return minProfits - fees * F(FEE)
end

f92!(Float32.(betAmounts), first(groupData.groups).y_matrix[:, bettableSlugsIndex], first(groupData.groups).n_matrix[:, bettableSlugsIndex], first(groupData.groups).not_na_matrix[:, bettableSlugsIndex], tmp.p, StructArrays.components(tmp.pool)..., sharesByEvent, zero(betAmounts), zero(betAmounts))
@benchmark f92!($(Float32.(betAmounts)), $(first(groupData.groups).y_matrix[:, bettableSlugsIndex]), $(first(groupData.groups).n_matrix[:, bettableSlugsIndex]), $(first(groupData.groups).not_na_matrix[:, bettableSlugsIndex]), $tmp.p, $(StructArrays.components(tmp.pool)[1]), $(StructArrays.components(tmp.pool)[2]), $(sharesByEvent), $(zero(betAmounts)), $(zero(betAmounts)))
@benchmark f92!($(Float32.(betAmounts[1:8])), $(first(groupData.groups).y_matrix[:, bettableSlugsIndex]), $(first(groupData.groups).n_matrix[:, bettableSlugsIndex]), $(first(groupData.groups).not_na_matrix[:, bettableSlugsIndex]), $tmp.p, $(StructArrays.components(tmp.pool)[1]), $(StructArrays.components(tmp.pool)[2]), $(sharesByEvent), $(zero(betAmounts)), $(zero(betAmounts)))

pVec = [marketDataBySlug[slug].p for slug in group.slugs[bettableSlugsIndex]]
poolSOA = StructArray([marketDataBySlug[slug].pool for slug in group.slugs[bettableSlugsIndex]])

function f93!(betAmount, y_matrix, n_matrix, not_na_matrix, pVec, pool, sharesByEvent, sharesNO, sharesYES)
    sharesYES .= 0
    sharesNO .= 0

    fees = zero(eltype(sharesByEvent))

    @turbo for i in eachindex(betAmount)
        y = pool.YES[i]
        n = pool.NO[i]

        p = pVec[i]

        sharesYES[i] = ifelse(betAmount[i] >= 1, y + betAmount[i] - y * (n / (n + betAmount[i]))^((1-p)/p), 0)
        sharesNO[i] = ifelse(betAmount[i] <= -1, n - betAmount[i] - n * (y / (y - betAmount[i]))^(p/(1-p)), 0)

        fees += sum(abs(betAmount[i]) >= 1)
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

    return minProfits - fees * F(FEE)
end

f93!(Float32.(betAmounts), first(groupData.groups).y_matrix[:, bettableSlugsIndex], first(groupData.groups).n_matrix[:, bettableSlugsIndex], first(groupData.groups).not_na_matrix[:, bettableSlugsIndex], pVec, poolSOA, sharesByEvent, zero(betAmounts), zero(betAmounts))
@benchmark f93!($(Float32.(betAmounts)), $(first(groupData.groups).y_matrix[:, bettableSlugsIndex]), $(first(groupData.groups).n_matrix[:, bettableSlugsIndex]), $(first(groupData.groups).not_na_matrix[:, bettableSlugsIndex]), $pVec, $poolSOA, $sharesByEvent, $(zero(betAmounts)), $(zero(betAmounts)))
@benchmark f93!($(Float32.(betAmounts[1:8])), $(first(groupData.groups).y_matrix[:, bettableSlugsIndex]), $(first(groupData.groups).n_matrix[:, bettableSlugsIndex]), $(first(groupData.groups).not_na_matrix[:, bettableSlugsIndex]), $pVec, $poolSOA, $sharesByEvent, $(zero(betAmounts)), $(zero(betAmounts)))


@code_native debuginfo=:none syntax=:intel f93!(Float32.(betAmounts), first(groupData.groups).y_matrix[:, bettableSlugsIndex], first(groupData.groups).n_matrix[:, bettableSlugsIndex], first(groupData.groups).not_na_matrix[:, bettableSlugsIndex], pVec, poolSOA, sharesByEvent, zero(betAmounts), zero(betAmounts))

function f9Test2!(betAmount, y_matrix, n_matrix, not_na_matrix, pVec, poolYES, poolNO, sharesByEvent, sharesYES, sharesNO)
    sharesYES .= 0
    sharesNO .= 0

    @turbo for i in eachindex(betAmount)
        y = poolYES[i]
        n = poolNO[i]

        p = pVec[i]

        sharesYES[i] = ifelse(betAmount[i] >= 1, y + betAmount[i] - y * (n / (n + betAmount[i]))^((1-p)/p), 0)
        sharesNO[i] = ifelse(betAmount[i] <= -1, n - betAmount[i] - n * (y / (y - betAmount[i]))^(p/(1-p)), 0)
    end

    minProfits = F(Inf)

    @turbo for j in axes(n_matrix, 2)
        profit = sharesByEvent[j]

        for i in axes(n_matrix, 1)
            profit += y_matrix[i, j] * sharesYES[i]
            profit += n_matrix[i, j] * sharesNO[i]
            profit -= not_na_matrix[i, j] * abs(betAmount[i])
        end

        minProfits = min(minProfits, profit)
    end

    return minProfits
end

f9Test2!(Float32.(betAmounts), first(groupData.groups).y_matrix[:, bettableSlugsIndex]', first(groupData.groups).n_matrix[:, bettableSlugsIndex]', first(groupData.groups).not_na_matrix[:, bettableSlugsIndex]', tmp.p, StructArrays.components(tmp.pool)..., sharesByEvent, zero(betAmounts), zero(betAmounts))
@benchmark f9Test2!($(Float32.(betAmounts)), $(first(groupData.groups).y_matrix[:, bettableSlugsIndex]' |> collect), $(first(groupData.groups).n_matrix[:, bettableSlugsIndex]' |> collect), $(first(groupData.groups).not_na_matrix[:, bettableSlugsIndex]' |> collect), $tmp.p, $(StructArrays.components(tmp.pool)[1]), $(StructArrays.components(tmp.pool)[2]), $(sharesByEvent), $(zero(betAmounts)), $(zero(betAmounts)))
@benchmark f9Test2!($(Float32.(betAmounts[1:8])), $(first(groupData.groups).y_matrix[:, bettableSlugsIndex]' |> collect), $(first(groupData.groups).n_matrix[:, bettableSlugsIndex]' |> collect), $(first(groupData.groups).not_na_matrix[:, bettableSlugsIndex]' |> collect), $tmp.p, $(StructArrays.components(tmp.pool)[1]), $(StructArrays.components(tmp.pool)[2]), $(sharesByEvent), $(zero(betAmounts)), $(zero(betAmounts)))