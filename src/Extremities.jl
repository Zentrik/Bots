using ManifoldMarkets, TOML, Optimization, OptimizationBBO, Dates, Combinatorics, LinearAlgebra, Parameters
using HTTP, JSON3, Dates, OpenSSL
using HTTP.WebSockets
using SmartAsserts, Logging, LoggingExtras

import ..Bots: @spawn_showerr, @smart_assert_showerr, display_error, wait_until, MutableOutcomeType, shouldBreak

Base.exit_on_sigint(false)

@with_kw mutable struct MarketData 
    shares::MutableOutcomeType{Float64} = MutableOutcomeType(0, 0)
end

@with_kw mutable struct BotData @deftype String
    APIKEY
    Supabase_APIKEY
    FIREBASE_APIKEY
    USERNAME
    USERID
    balance::Float64 = 0
end

@with_kw struct Arguments @deftype Bool
    live=false
    confirmBets=true
end

struct PlannedBet
    amount::Float64
    shares::Float64
    outcome::Symbol

    id::String
end

function execute!(botData, bet)
    ohno = false
    response = createBet(botData.APIKEY, bet.id, bet.amount, bet.outcome)
    # need to check if returned info matches what we wanted to bet, i.e. if we got less shares than we wanted to. If we got more ig either moved or smth weird with limit orders.

    botData.balance -= bet.amount

    @smart_assert_showerr String(bet.outcome) == response.outcome "$(bet.outcome), $(response.outcome)"

    if response.shares / bet.shares < 0.1 && (response.probAfter > 0.02 || response.probAfter < .98)
        @error bet.id # hyperlink
        @error response
        @error response.fills

        ohno = true
    end

    return response, ohno
end

function updateShares!(marketData, newBet, BotData)
    marketData.shares[Symbol(newBet.outcome)] += newBet.shares

    if marketData.shares[:YES] >= marketData.shares[:NO]
        marketData.shares[:YES] -= marketData.shares[:NO]
        BotData.balance += marketData.shares[:NO]
        marketData.shares[:NO] = 0.
    elseif marketData.shares[:YES] < marketData.shares[:NO]
        marketData.shares[:NO] -= marketData.shares[:YES]
        BotData.balance += marketData.shares[:YES]
        marketData.shares[:YES] = 0.
    end
end

function stringKeysToSymbol(dict)
    return Dict(Symbol(key) => value for (key, value) in dict)
end

function fetchMyPortfolio(botData)
    response = HTTP.get("https://pxidrgkatumlvfqaxcll.supabase.co/rest/v1/user_contract_metrics?user_id=eq.$(botData.USERID)", headers= ["apikey" => botData.Supabase_APIKEY, "Content-Type" => "application/json", "User-Agent" => "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/70.0.3538.102 Safari/537.36 Edge/18.19582"])
    responseJSON = JSON3.read(response.body)

    marketById = typeof(Dict("sad" => MarketData()))()
    for user_contract_metrics in responseJSON
        marketById[user_contract_metrics.contract_id] = MarketData()
        for (outcome, shares) in user_contract_metrics.data.totalShares
            marketById[user_contract_metrics.contract_id].shares[outcome] = shares
        end
    end

    return marketById
end

function readData()
    data = TOML.parsefile("$(@__DIR__)/Extremities.toml")
    APIKEY::String = data["APIKEY"]
    Supabase_APIKEY::String = data["SUPABASE_APIKEY"]
    FIREBASE_APIKEY::String = data["FIREBASE_APIKEY"]
    USERNAME::String = data["USERNAME"]
    SMARTUSERS = data["SMARTUSERS"]
    dumbUsers = data["dumbUsers"]

    return APIKEY, Supabase_APIKEY, FIREBASE_APIKEY, USERNAME, SMARTUSERS, dumbUsers
end

function setup(live, confirmBets)
    arguments = Arguments(live, confirmBets)

    APIKEY, Supabase_APIKEY, FIREBASE_APIKEY, USERNAME, SMARTUSERS, dumbUsers = readData()  
    
    @time "Fetching User" botUser = getUserByUsername(USERNAME)
    USERID = botUser.id
    botBalance = botUser.balance
    botData = BotData(APIKEY, Supabase_APIKEY, FIREBASE_APIKEY, USERNAME, USERID, botBalance)

    @time "Fetching Portfolio" marketById = fetchMyPortfolio(botData)

    return marketById, botData, arguments, Set(reduce(vcat, values(SMARTUSERS))), Set(reduce(vcat, values(dumbUsers)))
end

const pushJSON = JSON3.write(
	Dict("topic"=> "realtime:live-bets-", 
		"event"=>"phx_join", 
		"payload"=>Dict(
			"config"=>Dict("broadcast"=>Dict("ack"=> false, "self"=> false), 
						"presence"=>Dict("key"=>""), 
						"postgres_changes"=>[Dict("event"=>"INSERT", "schema"=>"public", "table"=>"contract_bets", "filter"=>"prob_after=lt.0.005"), Dict("event"=>"INSERT", "schema"=>"public", "table"=>"contract_bets", "filter"=>"prob_after=gt.0.995")]
			),
		),
		"ref"=>"1"
    )
)

isMarketClosed(market) = market.isResolved || market.closeTime / 1000 < time() # if resolved or closed


function production(; live=true, confirmBets=false)
    marketById, botData, arguments, smartUsers, dumbUsers = setup(live, confirmBets)

    WebSockets.open(uri(botData.Supabase_APIKEY), suppress_close_error=true, headers= ["User-Agent" => "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/70.0.3538.102 Safari/537.36 Edge/18.19582"]) do socket
        try
            @info "Opened Socket"
            @debug socket

            # We could fetch contracts instead which saves us making a HTTP request but we really want the bettors id
            # send(socket, pushJSON)
            send(socket, pushJSONBets())
            @info "Sent Intialisation"
    
            Base.Experimental.@sync begin
                #Reading messages
                # should switch to @spawn
                @spawn_showerr for msg in socket
                    msgJSON = JSON3.read(msg)
                    if !(:data in keys(msgJSON.payload))
                        if :status in keys(msgJSON.payload) && msgJSON.payload.status == "error"
                            # @error "$(Dates.format(now(), "HH:MM:SS.sss")): $msg"
                            throw(ErrorException("$(Dates.format(now(), "HH:MM:SS.sss")): $msg"))
                        else
                            @debug "$(Dates.format(now(), "HH:MM:SS.sss")): $msg"
                            continue
                        end
                    end

                    # FILTER out updates, only want inserts

                    bet = msgJSON.payload.data.record.data

                    if .005 < bet.probAfter < .995
                        continue
                    end

                    # @debug msgJSON
                    @debug bet

                    if bet.isChallenge
                        @debug "Bet was challenge"
                        continue
                    end

                    if bet.isAnte
                        @debug "Bet is ante"
                        continue
                    end

                    # # On creation of dpm-2 market the subsidy is placed as a bet, https://manifold.markets/JuJumper/who-will-be-the-next-patriarch-of-m
                    # if bet.probBefore == 0 && bet.probAfter == 1
                    #     @debug "Is this a dpm-2 market $(bet.probBefore), $(bet.probAfter)"
                    #     continue
                    # end

                    if bet.outcome ∉ ("YES", "NO")
                        @debug "Not a binary market"
                        continue
                    end
                    
                    if bet.createdTime/1000 < time() - 10
                        @debug "Bet happened over 10s ago, $(bet.createdTime) in unix time"
                        continue
                    end

                    marketId = bet.contractId
                    @debug marketId

                    sign = bet.probAfter < 0.005 ? 1 : -1

                    outcome = sign == 1 ? :YES : :NO

                    oppositeOutcome = outcome == :YES ? :NO : :YES

                    # if marketId in keys(marketById) && marketById[marketId].shares[oppositeOutcome] > 0
                    #     bet = PlannedBet(20, marketById[marketId].shares[oppositeOutcome], outcome, marketId)
                    #     response, ohno = execute!(botData, bet)

                    #     # Too annoying to get token every time it expires
                    #     # response = HTTP.post("https://api-nggbo3neva-uc.a.run.app/sellshares", headers=["Authorization" => "Bearer " * botData.FIREBASE_APIKEY, "Content-Type" => "application/json"], body=Dict("contractId"=>marketId, "outcome"=>oppositeOutcome, "shares" => marketById[marketId].shares[oppositeOutcome]), socket_type_tls=OpenSSL.SSLStream)

                    #     updateShares!(marketById[marketId], response, botData)
                    # end
                    # continue on so we bet more if big prob change


                    # Bigger probem is big prob changes based on news, we will incorrectly bet in this case 

                    # should account for anti-bot market: bet to 90% then big move. As long as we only bet 10 and don't bet on a market multiple times should be fine
                    # bet from 1 to 10, then 10 to 1 is an easy way to exploit this bot. rn using market.prob might help as that is some sort of lagged value. Could look at probChanges?

                    # if market creator is the bettor/ market new ? how to test this
                    
                    if bet.userId in smartUsers 
                        @debug "$(bet.userId) is a smart bettor"
                        continue
                    end

                    if bet.userId == botData.USERID
                        @debug "I made this bet"
                        continue
                    end

                    # either not a liquid market or small prob change/ what matter is how much an idiot bet not how much the market slipped
                    # odd number to make it harder to work out, could add rand() ppl start working it out
                    # might need different logic if selling, but given we assuming market prob before is 10-90 probs didn't matter
                    if abs(bet.shares) < 218.1567984
                        @debug "$(abs(bet.shares)) is too small"
                        continue
                    end

                    if !((bet.probBefore > 0.1 && bet.probAfter < 0.005) || (bet.probBefore < 0.9 && bet.probAfter > 0.995))
                        @debug "Bet Probabilities not in desired range, $(bet.probBefore), $(bet.probAfter)"
                        continue
                    end

                    @debug "market id: $marketId from $(bet.probBefore) to $(bet.probAfter) at $(Dates.format(now(), "HH:MM:SS.sss"))"

                    if marketId in keys(marketById) 
                        if marketById[marketId].shares[outcome] > 1 
                            @debug "$(marketById[marketId].shares[outcome]) already own shares, $marketId, $(marketById[marketId].shares), $outcome"
                            continue
                        elseif marketById[marketId].shares[outcome == :YES ? :NO : :YES] > 500
                            amount = 40.
                        end
                    else
                        # Could cause race
                        marketById[marketId] = MarketData()
                        amount = 20.
                    end

                    if bet.userId in dumbUsers
                        amount *= 2
                    end

                    timenow = time_ns()
                    response = HTTP.get("https://pxidrgkatumlvfqaxcll.supabase.co/rest/v1/contracts?id=eq.$marketId", headers = ["apikey" => botData.Supabase_APIKEY, "Content-Type" => "application/json", "User-Agent" => "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/70.0.3538.102 Safari/537.36 Edge/18.19582"], socket_type_tls=OpenSSL.SSLStream)
                    @info "Getting Market took $((time_ns() - timenow)/10^6) ms"
                    responseJSON = JSON3.read(response.body)[1]

                    @debug responseJSON

                    market = responseJSON.data

                    if market.mechanism != "cpmm-1"
                        @debug "$(market.mechanism) not cpmm-1"
                        continue
                    end

                    if responseJSON.creator_id == bet.userId 
                        @debug "$(bet.userId) is the market creator"
                        continue
                    end

                    if market.closeTime - time() * 1000 < 24*60*60*1000 + 6.234*60*1000 # closes in a day
                        @debug "$marketId closes in $(market.closeTime/1000 - time())s"
                        continue
                    end

                    if time() * 1000 - market.createdTime < 1.12654*24*60*60*1000 && market.outcomeType != "STONK" # created a day ago
                        @debug "$marketId created $(time() - market.createdTime/1000)s ago"
                        continue
                    end

                    if isMarketClosed(market)
                        @debug "$marketId closed, $(market.isResolved), close time is $(market.closeTime)"
                        continue
                    end

                    newProb = poolToProb(market.p, market.pool)

                    if !((market.prob > 0.1 && newProb < 0.005) || (bet.probBefore < 0.9 && newProb > 0.995))
                        @debug "Market Probabilities not in desired range, $(market.prob), $(newProb), $(market.p), $(market.pool)"
                        continue
                    end
                    shares = amount / (bet.probAfter + sign * .25)
                    bet = PlannedBet(amount, shares, outcome, marketId)

                    @info "Making bet on $marketId for $amount buying $outcome at $(Dates.format(now(), "HH:MM:SS.sss"))"
                    
                    if arguments.live
                        response, ohno = execute!(botData, bet)
                        updateShares!(marketById[marketId], response, botData)
                        # handle non cppm-1 markets

                        # sell shares if probability moved too far from extreme

                        # only works for cpmm-1
                        if ohno 
                            bet = PlannedBet(amount, marketById[marketId].shares[outcome], oppositeOutcome, marketId)
                            response, ohno = execute!(botData, bet)
                            # response = HTTP.post("https://api-nggbo3neva-uc.a.run.app/sellshares", headers=["Authorization" => "Bearer " * botData.FIREBASE_APIKEY, "Content-Type" => "application/json"], body=Dict("contractId"=>marketId, "outcome"=> response.outcome, "shares" => respone.shares), socket_type_tls=OpenSSL.SSLStream)

                            # Could cause data race
                            updateShares!(marketById[marketId], response, botData)
                        end
                    end

                    continue
                end
    
                # HeartBeat
                @spawn_showerr while !WebSockets.isclosed(socket)
                    try
                        send(socket, heartbeatJSON)
                    catch
                        @error "Heartbeat Failed"
                        rethrow()
                    end
                    sleep(30)
                end

                # Fetch new balance at 8am
                @spawn_showerr while !WebSockets.isclosed(socket)
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
            println("Finally")
            if !WebSockets.isclosed(socket)
                send(socket, leaveJSON("live-bets"))
                println("Left Channel")
    
                close(socket)
            end
        end
    end
end

function test(; live=false, confirmBets=true) 
    production(; live=live, confirmBets=confirmBets)
end

function retryProd(runs=1; live=true, confirmBets=false)
    delay = 60
    lastRunTime = time()

    for run in 1:runs
        try
            production(; live=live, confirmBets=confirmBets)
        catch 
            tmp = false
            tmp = for (exception, _) in current_exceptions()
                if shouldBreak(exception)
                    return true
                end
            end

            if tmp == true
                break
            end

            println(Dates.format(now(), "HH:MM:SS.sss"))

            if run == runs
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