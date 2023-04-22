using ManifoldMarkets, TOML, Optimization, OptimizationBBO, Dates, Combinatorics, LinearAlgebra, Parameters
using HTTP, JSON3, Dates, OpenSSL
using HTTP.WebSockets
using SmartAsserts, Logging, LoggingExtras

macro async_showerr(ex)
    esc(quote
        @async try
            eval($ex)
        catch err
            bt = catch_backtrace()
            # println()
            # showerror(stderr, err, bt)
            @error "Something went wrong" exception = (err, bt)
            rethrow()
        end
    end)
end

Base.exit_on_sigint(false)

@with_kw mutable struct MarketData 
    shares::Dict{Symbol, Float64} = Dict{Symbol, Float64}(:NO => 0., :YES => 0.)
end

@with_kw mutable struct BotData @deftype String
    const APIKEY
    const Supabase_APIKEY
    const FIREBASE_APIKEY
    const USERNAME
    const USERID
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

    botData -= bet.amount

    try 
        @smart_assert bet.outcome == response.outcome
    catch err
        @error err.msg
        rethrow(err)
    end

    if abs(response.shares - bet.shares) / bet.shares < 0.05
        @error bet.id # hyperlink
        @error response
        @error response.fills

        ohno = true
    end

    return response, ohno
end

function updateShares!(MarketData, newBet, BotData)
    MarketData.shares[Symbol(newBet.outcome)] += newBet.shares

    if MarketData.shares[:YES] >= MarketData.shares[:NO]
        MarketData.shares[:YES] -= MarketData.shares[:NO]
        BotData.balance += MarketData.shares[:NO]
        MarketData.shares[:NO] = 0.
    elseif MarketData.shares[:YES] < MarketData.shares[:NO]
        MarketData.shares[:NO] -= MarketData.shares[:YES]
        BotData.balance += MarketData.shares[:YES]
        MarketData.shares[:YES] = 0.
    end
end

function stringKeysToSymbol(dict)
    return Dict(Symbol(key) => value for (key, value) in dict)
end

function fetchMyPortfolio(botData)
    try
        response = HTTP.get("https://pxidrgkatumlvfqaxcll.supabase.co/rest/v1/user_contract_metrics?user_id=eq.$(botData.USERID)", headers= ["apikey" => botData.Supabase_APIKEY, "Content-Type" => "application/json"])
        responseJSON = JSON3.read(response.body)

        marketById = typeof(Dict("sad" => Extremities.MarketData()))()
        for user_contract_metrics in responseJSON
            marketById[user_contract_metrics.contract_id] = MarketData(stringKeysToSymbol(user_contract_metrics.data.totalShares))
        end

        return marketById
    catch err
        bt = catch_backtrace()
        println()
        showerror(stderr, err, bt)
        throw(err)
    end
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

pushJSON = JSON3.write(
	Dict("topic"=> "realtime:live-bets", 
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

    WebSockets.open(uri(botData.Supabase_APIKEY)) do socket
        try
            @info "Opened Socket"
            @debug socket

            # We could fetch contracts instead which saves us making a HTTP request but we really want the bettors id
            send(socket, pushJSON)
            @info "Sent Intialisation"
    
            @sync try
                #Reading messages
                # should switch to @spawn
                @async_showerr for msg in socket
                    @async_showerr begin 
                        msgJSON = JSON3.read(msg)
                        if !(:payload in keys(msgJSON) && :data in keys(msgJSON.payload))
                            @debug "$(Dates.format(now(), "HH:MM:SS.sss")): $msg"
                            return nothing
                        end

                        # FILTER out updates, only want inserts

                        bet = msgJSON.payload.data.record.data

                        # @debug msgJSON
                        @debug bet

                        if bet.isChallenge
                            @debug "Bet was challenge"
                            return nothing
                        end

                        if bet.outcome ∉ ("YES", "NO")
                            @debug "Not a binary market"
                            return nothing
                        end
                        
                        if bet.fills[end].timestamp/1000 < time() - 10
                            @debug "Bet happened over 10s ago, $(bet.fills[end].timestamp) in unix time"
                            return nothing
                        end

                        # should account for anti-bot market: bet to 90% then big move. As long as we only bet 10 and don't bet on a market multiple times should be fine
                        # bet from 1 to 10, then 10 to 1 is an easy way to exploit this bot. rn using market.prob might help as that is some sort of lagged value. Could look at probChanges?

                        # if market creator is the bettor/ market new ? how to test this
                        
                        if bet.userId in smartUsers 
                            @debug "$(bet.userId) is a smart bettor"
                            return
                        end

                        marketId = bet.contractId
                        @debug marketId

                        # either not a liquid market or small prob change/ what matter is how much an idiot bet not how much the market slipped
                        # odd number to make it harder to work out, could add rand() ppl start working it out
                        # might need different logic if selling, but given we assuming market prob before is 10-90 probs didn't matter
                        if abs(bet.shares) < 218.1567984
                            @debug "$(abs(bet.shares)) is too small"
                            return
                        end

                        if bet.isAnte
                            @debug "Bet is ante"
                            return
                        end

                        # On creation of dpm-2 market the subsidy is placed as a bet, https://manifold.markets/JuJumper/who-will-be-the-next-patriarch-of-m
                        if bet.probBefore == 0 && bet.probAfter == 1
                            @debug "Is this a dpm-2 market $(bet.probBefore), $(bet.probAfter)"
                            return
                        end

                        if !(bet.probBefore > 0.1 && bet.probAfter < 0.005) || (bet.probBefore < 0.9 && bet.probAfter > 0.995)
                            @debug "Bet Probabilities not in desired range, $(bet.probBefore), $(bet.probAfter)"
                            return
                        end

                        @debug "market id: $marketId from $(bet.probBefore) to $(bet.probAfter) at $(Dates.format(now(), "HH:MM:SS.sss"))"

                        sign = bet.probAfter < 0.005 ? 1 : -1

                        outcome = sign == 1 ? :YES : :NO

                        if marketId in keys(marketById) 
                            if marketById[marketId].shares[outcome] > 1
                                @debug "$(marketById[marketId].shares[outcome]) already own shares, $marketId, $(marketById[marketId].shares), $outcome"
                                return
                            elseif marketById[marketId].shares[outcome == :YES ? :NO : :YES] > 50
                                amount = 10.
                            end
                        else
                            marketById[marketId] = MarketData()
                            amount = 5.
                        end

                        if bet.userId in dumbUsers
                            amount *= 2
                        end

                        response = HTTP.get("https://pxidrgkatumlvfqaxcll.supabase.co/rest/v1/contracts?id=eq.$marketId", headers = ["apikey" => botData.Supabase_APIKEY, "Content-Type" => "application/json"], socket_type_tls=OpenSSL.SSLStream)
                        responseJSON = JSON3.read(response.body)[1]

                        @debug responseJSON

                        if responseJSON.data.mechanism != "cpmm-1"
                            @debug "$(responseJSON.data.mechanism) not cpmm-1"
                            return
                        end

                        if responseJSON.creator_id == bet.userId 
                            @debug "$(bet.userId) is the market creator"
                            return
                        end

                        if responseJSON.data.closeTime - time() * 1000 < 60*60*1000 + 6.234*60*1000 # closes in an hour
                            @debug "$marketId closes in $(responseJSON.data.closeTime/1000 - time())s"
                            return
                        end

                        if time() * 1000 - responseJSON.data.createdTime < 1.12654*24*60*60*1000# created a day ago
                            @debug "$marketId created $(time() - responseJSON.data.createdTime/1000)s ago"
                            return
                        end

                        if isMarketClosed(responseJSON.data)
                            @debug "$marketId closed, $(responseJSON.data.isResolved), close time is $(responseJSON.data.closeTime)"
                            return
                        end

                        newProb = poolToProb(market.p, market.pool)

                        if !(market.prob > 0.1 && newProb < 0.005) || (bet.probBefore < 0.9 && newProb > 0.995)
                            @debug "Market Probabilities not in desired range, $(market.prob), $(newProb), $(market.p), $(market.pool)"
                            return
                        end
                        shares = amount / (bet.probAfter + sign * .25)
                        bet = PlannedBet(amount, shares, outcome, marketId)

                        @info "Making bet on $marketId for $amount buying $outcome at $(Dates.format(now(), "HH:MM:SS.sss"))"
                        
                        if arguments.live
                            response, ohno = execute!(botData, bet)
                            # handle non cppm-1 markets

                            # sell shares is probability moved too far from extreme

                            # only works for cpmm-1
                            if ohno 
                                response = HTTP.post("https://api-nggbo3neva-uc.a.run.app/sellshares", headers=["Authorization" => "Bearer " * botData.FIREBASE_APIKEY, "Content-Type" => "application/json"], body=Dict("contractId"=>marketId, "outcome"=> response.outcome, "shares" => respone.shares), socket_type_tls=OpenSSL.SSLStream)
                            else
                                updateShares!(MarketData, response, botData)
                            end
                        end

                        return nothing
                    end
                end
    
                # HeartBeat
                @async_showerr while !WebSockets.isclosed(socket)
                    send(socket, heartbeatJSON)
                    sleep(30)
                end

                # Fetch new balance at 8am
                @async_showerr while !WebSockets.isclosed(socket)
                    printstyled("Fetching Balance at $(Dates.format(now(), "HH:MM:SS.sss"))\n"; color = :blue)
                    botData.balance = getUserByUsername(botData.USERNAME).balance
                    timeTo8 = Second(((today() + Time(8) + Minute(5) - now()) ÷ 1000).value)
                    timeTo8 += timeTo8 > Second(0) ? Second(0) : Second(Day(1))
                    sleep(timeTo8)
                    printstyled("Fetcing Balance at $(Dates.format(now(), "HH:MM:SS.sss"))\n"; color = :blue)
                    botData.balance = getUserByUsername(botData.USERNAME).balance
                    timeTo8 = Second(((today() + Time(8) + Minute(10) - now()) ÷ 1000).value)
                    timeTo8 += timeTo8 > Second(0) ? Second(0) : Second(Day(1))
                    sleep(timeTo8)
                    printstyled("Fetcing Balance at $(Dates.format(now(), "HH:MM:SS.sss"))\n"; color = :blue)
                    botData.balance = getUserByUsername(botData.USERNAME).balance
                    timeTo8 = Second(((today() + Time(8) + Minute(30) - now()) ÷ 1000).value)
                    timeTo8 += timeTo8 > Second(0) ? Second(0) : Second(Day(1))
                    sleep(timeTo8)
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

function retryProd(; live=true, confirmBets=false)
    delay = 60
    lastRunTime = time()

    while true
        try
            production(; live=live, confirmBets=confirmBets)
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
    try 
        @smart_assert 1 != 1
    catch err
        @error err.msg
        # rethrow(err)
    end
end

# Need to run manually
# using Logging, LoggingExtras, Dates

# timestamp_logger(logger) = TransformerLogger(logger) do log
#     merge(log, (; message = "$(Dates.format(now(), "yyyy-mm-dd HH:MM:SS.sss")) $(log.message)"))
# end
# const DIR = "H:\\Code\\ManifoldMarkets.jl\\Bots\\src\\logs-Extremities"
# global_logger(TeeLogger(
#     EarlyFilteredLogger(Bots.not_Bots_message_filter, ConsoleLogger(stderr, Logging.Info)), 
#     EarlyFilteredLogger(Bots.not_Bots_message_filter, MinLevelLogger(FileLogger("$DIR\\$(Dates.format(now(), "YYYY-mm-dd")).log"), Logging.Info)),
#     EarlyFilteredLogger(Bots.not_Bots_message_filter, MinLevelLogger(FileLogger("$DIR\\Debug-$(Dates.format(now(), "YYYY-mm-dd")).log"), Logging.Debug)),
#     timestamp_logger(MinLevelLogger(FileLogger("$DIR\\Verbose-$(Dates.format(now(), "YYYY-mm-dd")).log"), Logging.Debug))
# ))
# Extremities.testLogging()