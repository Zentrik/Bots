using ManifoldMarkets

# @code_warntype getBets()
bets = getBets()

betsDict = Dict{String, Vector{Bet}}()
@time for url in ("https://manifold.markets/Lifejacker/will-there-be-83500-confirmed-cases", "https://manifold.markets/Lifejacker/will-there-be-83500-confirmed-cases-e402502ec82e", "https://manifold.markets/IsaacKing/will-this-markets-probability-be-at", "https://manifold.markets/BionicD0LPH1N/will-isaac-kings-69-market-resolve", "https://manifold.markets/NcyRocks/will-joe-biden-run-for-president-in", "https://manifold.markets/LivInTheLookingGlass/will-joe-biden-seek-a-second-term-i", "https://manifold.markets/ahalekelly/will-vladimir-putin-still-be-the-le", "https://manifold.markets/njmkw/if-ukraine-recaptures-cherson-in-20-abc72a4594fc")
    begin
        slug = urlToSlug(url)
        betsDict[slug] = getBets(slug=slug)
    end
end

betsDict = Dict{String, Vector{Bet}}()
@time @sync for url in ("https://manifold.markets/Lifejacker/will-there-be-83500-confirmed-cases", "https://manifold.markets/Lifejacker/will-there-be-83500-confirmed-cases-e402502ec82e", "https://manifold.markets/IsaacKing/will-this-markets-probability-be-at", "https://manifold.markets/BionicD0LPH1N/will-isaac-kings-69-market-resolve", "https://manifold.markets/NcyRocks/will-joe-biden-run-for-president-in", "https://manifold.markets/LivInTheLookingGlass/will-joe-biden-seek-a-second-term-i", "https://manifold.markets/ahalekelly/will-vladimir-putin-still-be-the-le", "https://manifold.markets/njmkw/if-ukraine-recaptures-cherson-in-20-abc72a4594fc")
    @async begin
        slug = urlToSlug(url)
        betsDict[slug] = getBets(slug=slug)
    end
end


###

slug = "china-officially-abandons-covid-zer"

market = getMarketBySlug(slug)
bets = getAllBets(slug=slug)
limitOrdersAmountsShares, splitBets = getLimitOrders(bets) 
betToShares(market, limitOrdersAmountsShares, -2780)

### This market has loads of limit orders

slug = "dan-stock-permanent"

market = getMarketBySlug(slug)
bets = getAllBets(slug=slug)
bets = getBets(slug=slug)
sortedLimits, userBalances = sortLimitOrders(bets)
limitOrdersAmountsShares, splitBets = getLimitOrders(bets)
betToShares(market, limitOrdersAmountsShares, -2100)