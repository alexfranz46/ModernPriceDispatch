function performance_with_historic_data(case::BESTCase; y0=0, plim=500)

    #TODO Reformat in same style as everything else...

    policy = case.policy
    data = case.data.testCaseData
    bess = case.bess
    PriceBounds = case.price.PriceBounds

    # create individual visual for each day
    for day in groupby(data, :TradingDate)

        # Set initial conditions
        dailyRevenue = 0
        yHistory = [y0]

        # Follow decisions stage by stage
        for (stage, period) in enumerate(eachrow(day))

            # Variables
            tp = period.TradingPeriod
            yIn = yHistory[end]
            pNow = period.DollarsPerMegawattHour
            i = find_band(pNow, PriceBounds[:,tp])

            # Update storage
            yNext = policy.storageDecisions[stage, yIn+1, i]
            push!(yHistory, yNext)

            # Upddate revenue
            dailyRevenue += pNow*(yIn - yNext)
        end

        # Plot results 

        # format x-axis
        vals = range(0, 24, length=nrow(day)+1)
        ticks = 0:0.5:24
        labels = [mod(x, 1) == 0 ? string(Int(x)) : "" for x in ticks]

        # Plot storage on left axis
        p=plot(vals, yHistory,
            title="$(day.TradingDate[1]); revenue=\$$(round(dailyRevenue, digits=2))",
            xlabel="time of day (HH)",
            ylabel="Battery charge (MWh)", 
            ylims=(-0.03*bess.storageCapacityMWh,1.03*bess.storageCapacityMWh),
            legend=false, 
            color=:blue,
            linewidth=2,
            ytickfontcolor=:blue,
            y_guidefontcolor=:blue,
            xticks=(ticks, labels))

        # Plot price on right axis
        p2 = twinx()
        plot!(p2, vals[1:end-1], day.DollarsPerMegawattHour,
            ylabel="Price (\$/MWh)", 
            legend=false, 
            color=:red,
            linestyle=:solid,
            linewidth=2, 
            ytickfontcolor=:red,
            y_guidefontcolor=:red,
            xticks=:none) # Hides overlapping x-ticks from the second axis
        
        pmax = maximum(day.DollarsPerMegawattHour)
        if plim < pmax
            ylims!(p2, (-0.03*pmax, 1.03*pmax))
            display(p)
            annotate!(p2, [(12, pmax/2, text("INCOMPLETE PRICE SCALE"))])
        end

        ylims!(p2, (-0.03*plim,1.03*plim))
        display(p)
    end
end

#TODO: performance with simulated data?