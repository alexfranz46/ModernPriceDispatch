using Printf
using Serialization

""" Prints the optimal offers and bids in the format of the forms
"""
function print_policy_trades(stage, yIn, uv, b)
    tp = ceil(Int, stage/6)

    uv = uv[stage, yIn + 1, :]
    b = round.(b[:, tp], digits=2)

    pivot = findfirst(x -> x >= 0, uv)
    
    if pivot !== nothing
        offers = [uv[pivot]; diff(uv[pivot:end])]

        idx = pivot
        first = true
        for val in offers
            if val != 0
                if first
                    println("OFFERS")
                    print("From 0MW to ")
                    first = false
                else
                    print("       plus ")
                end
                # println("$(lpad(val, 2))MW @ \$$(lpad(b[idx], 6))/MWh")
                @printf("%2dMW   @   \$%6.2f/MWh\n", val, b[idx])
            end
            idx = idx + 1
        end
    end

    if pivot != 1
        bids = [uv[pivot - 1]; diff(reverse(uv[1:pivot-1]))]

        idx = pivot - 1
        first = true
        for val in bids
            if val != 0
                if first
                    println("BIDS")
                    print("From 0MW to ")
                    first = false
                else
                    print("       plus ")
                end
                # println("$(lpad(-val, 2))MW below \$$(b[idx+1])/MWh")
                @printf("%2dMW below \$%6.2f/MWh\n", -val, b[idx+1])
            end
            idx = idx - 1
        end
    end
end

function print_SoC_dependent_policy_trades(stage, E, uv, b)
    yPrev = uv[stage, 1, :]
    previousE = 0
    for e in 1:E
        yNow = uv[stage, e + 1, :]
        
        if yNow != yPrev
        println("SoC ∈ [$(previousE), $(e-1)] MWh:")
        print_policy_trades(stage, e-1, uv, b)
        println()
        
        previousE = e
        yPrev = yNow
        end
    end
    # Print the final block
    println("SoC ∈ [$(previousE), $(E)] MWh:")
    print_policy_trades(stage, E, uv, b)
    println()
end

function performance_from_starting_y()
    E = 240
    TP = 48
    T = 288
    uvDecision = deserialize(joinpath("Serials", "uvDecision2.jls"))
    PriceVals = deserialize(joinpath("Serials", "PriceVals.jls"))
    PriceBounds = deserialize(joinpath("Serials", "PriceBounds.jls"))
    
    pActual =   [
        110, 105, 104, 101,  82,  82, 104, 104, 104, 104, 104, 104,
        106, 104, 104, 104, 104, 103, 104, 104, 104, 103, 103, 103,
        103, 103, 103, 103, 100, 100, 102, 102, 101, 102, 101, 101,
        100, 101, 100, 101, 100, 100, 100, 100, 100, 100, 100, 100,
        100, 100, 103, 103, 103, 103, 100, 104, 104, 104, 106, 106,
        83, 100, 103, 103, 108, 109,  83,  84, 101, 105, 117, 129,
        83,  86, 109, 132, 133, 135, 111, 132, 136, 137, 140, 145,
        95, 130, 130, 147, 149, 151, 145, 160, 175, 190, 190, 157,
        200, 151, 149, 149, 139, 130, 149, 149, 149, 131, 113, 113,
        149, 113, 113, 113, 113, 103, 149, 111,  92,  92, 109, 111,
        110,  92,  91,  89,  89,  89,  90,  90,  89,  89,  88,  88,
        90,  90,  90,  89,  89,  87,  77,  80,  77,  81,  84,  74,
        86,  87,  87,  85,  78,  78,  80,  80,  80,  79,  79,  78,
        82,  82,  82,  81,  81,  79,  82,  81,  68,  68,  68,  24,
        67,  59,  59,  58,  58,  58,  45,  45,  45,  45,  35,  34,
        77,  76,  78,  80,  82,  84,  82,  89,  90,  89,  90,  92,
        90,  92,  92,  92,  95, 108,  92,  94,  98, 100, 110, 110,
        92,  92,  94,  97, 108, 110,  95,  98,  99,  99,  98, 104,
        105, 100, 149, 100, 109, 140, 151, 152, 105, 141, 104, 105,
        151, 109, 151, 105,  99,  98, 103,  96,  96,  95,  95,  94,
        113,  94,  93,  91,  91,  90,  95,  95,  95,  92,  92,  91,
        92, 119, 119, 112,  94,  92, 123, 111,  94,  91,  91,  90,
        94,  92,  90,  88,  88,  87, 109,  93,  89,  87,  87,  85,
        93, 109,  93,  85,  84,  84,  86,  92,  85,  84,  84,  83
    ]

    iActual = []
    for (t, p) in enumerate(pActual)
        push!(iActual, find_band(Float64(p), PriceBounds[:,ceil(Int, t/6)]))
    end

    dayRevenue = zeros(E+1)
    i = find_band(110.0, PriceBounds[:,TP])  # starting band
    for y0 in 0:E
        caseObjective = 0
        yHistory = [y0]  # starting with an empty battery
        pHistory = [PriceVals[i, TP]]
        for stage in 1:T
            yIn = yHistory[stage]
            
            # Update trading period
            tp = ceil(Int, stage/6)

            # fetch and record optimal decision for stage/state
            yNext = yIn - uvDecision[stage, yIn + 1, i] 
            push!(yHistory, yNext)
            
            # fetch and record price for stage/state
            pNow = PriceVals[i, tp]
            push!(pHistory, pNow)

            # calculate objective
            caseObjective += pActual[stage]*(yIn - yNext)
            # if y0 == 0 || y0 == 240
            #     println(caseObjective)
            # end

            # use random distribution to find j, which becomes i for next iteration
            j = iActual[stage]
            i = j
        end

        dayRevenue[y0+1] = caseObjective

        if y0 == 0 || y0 == 240
            t = 1:T

            ticks = 0:0.5:24
            labels = [mod(x, 1) == 0 ? string(Int(x)) : "" for x in ticks]

            p = plot(t/12, yHistory[2:end], 
                xlabel="time of day (HH)",
                ylabel="Battery charge (MWh)", 
                legend=:topleft, 
                label="y", 
                color=:blue,
                linewidth=2,
                xticks=(ticks, labels))

            plot!(twinx(), t/12, pActual,
                ylabel="Price (\$/MWh)", 
                legend=:topright, 
                label="p", 
                color=:red,
                linestyle=:solid,
                linewidth=2, 
                xticks=:none) # Hides overlapping x-ticks from the second axis

            display(p)
        end
    end
end