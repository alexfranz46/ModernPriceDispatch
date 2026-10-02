#TODO Reformat in same style as everything else...

function historic_performance(case::BESTCase; y0=0, plim=500)

    policies = case.policies
    data = case.data.testCaseData
    bess = case.bess
    PriceBounds = case.price.PriceBounds

    colors = cgrad(:Blues, length(policies) + 1, categorical=true)

    # create individual visual for each day
    for day in groupby(data, :TradingDate)

        # Initialize plot
        p = plot(xlabel="time of day (HH)", ylabel="Battery charge (MWh)", legend=:topleft)

        # format x-axis
        vals = range(0, 24, length=nrow(day)+1)
        ticks = 0:0.5:24
        labels = [mod(x, 1) == 0 ? string(Int(x)) : "" for x in ticks]

        # Plot each all solvers on same canvas
        for (idx, (solver, policy))  in enumerate(case.policies)

            # Set initial conditions
            dailyRevenue = 0
            yHistory = [y0]

            # Follow decisions stage by stage
            for (stage, interval) in enumerate(eachrow(day))

                # Variables
                yIn = yHistory[stage]
                i = find_band(interval.DollarsPerMegawattHour, PriceBounds[:,interval.TradingPeriod])

                # Update storage
                yNext = policy.storageDecisions[stage, yIn+1, i]
                push!(yHistory, yNext)

                # Upddate revenue
                dailyRevenue += interval.DollarsPerMegawattHour*(yIn - yNext)
            end

            # Plot storage on left axis
            plot!(p, vals, yHistory,
                title="$(day.TradingDate[1]); revenue=\$$(round(dailyRevenue, digits=2))",
                ylims=(-0.03*bess.storageCapacityMWh,1.03*bess.storageCapacityMWh),
                label=String(solver),
                color=colors[idx+1],
                linewidth=2,
                ytickfontcolor=:blue,
                y_guidefontcolor=:blue,
                xticks=(ticks, labels)
            )
        end

        # plot price on right axis
        p2 = twinx(p)
        plot!(p2, vals[1:end-1], day.DollarsPerMegawattHour,
        ylabel="Price (\$/MWh)", 
        legend=false, 
        color=:red,
        linestyle=:solid,
        linewidth=2, 
        ytickfontcolor=:red,
        y_guidefontcolor=:red,
        xticks=:none) # Hides overlapping x-ticks from the second axis

        # format right y-axis
        pmax = maximum(day.DollarsPerMegawattHour)
        if plim < pmax
            ylims!(p2, (-0.03*pmax, 1.03*pmax))
            annotate!(p2, [(12, pmax/2, text("ADJUSTED PRICE SCALE"))])
        else
            ylims!(p2, (-0.03*plim,1.03*plim))
        end

        display(p)
    end
end


function historic_performance_1hr(case::BESTCase; f=conservative, y0=zeros(Int, 12), U=zeros(Int, 12 - 1, 10), plim=500)

    policies = case.policies
    data = case.data.testCaseData
    bess = case.bess
    PriceBounds = case.price.PriceBounds

    colors = cgrad(:Blues, length(policies) + 1, categorical=true)

    # create individual visual for each day
    for day in groupby(data, :TradingDate)

        # Initialize plot
        p = plot(xlabel="time of day (HH)", ylabel="Battery charge (MWh)", legend=:topleft)

        # format x-axis
        vals = range(0, 24, length=nrow(day)+1)
        ticks = 0:0.5:24
        labels = [mod(x, 1) == 0 ? string(Int(x)) : "" for x in ticks]

        # Plot each all solvers on same canvas
        for (idx, (solver, policy))  in enumerate(case.policies)

            # Set initial conditions
            dailyRevenue = 0
            yHistory = [yi for yi in y0]

            # Follow decisions stage by stage
            for (stage, interval) in enumerate(eachrow(day))

                # Variables
                yIn = f(case, yHistory[stage], U)
                i = find_band(interval.DollarsPerMegawattHour, PriceBounds[:,interval.TradingPeriod])

                if yIn < 0  
                    # @warn "infeasible f(⋅)"
                    yIn = 0
                elseif yIn > case.bess.storageCapacityMWh
                    # @warn "infeasible f(⋅)"
                    yIn = case.bess.storageCapacityMWh
                end
                # Update storage
                yNext = policy.storageDecisions[stage, yIn+1, i]
                push!(yHistory, yNext)

                # Upddate revenue
                dailyRevenue += interval.DollarsPerMegawattHour*(yIn - yNext)

                # Update unrealised trades
                uNext = vec(policy.tradeDecisions[stage, yIn+1, :])
                # println("$stage: $yNext")
                U = [U[2:end, :]; uNext']
            end

            # Plot storage on left axis
            plot!(p, vals, yHistory[12:end],
                title="$(day.TradingDate[1]); revenue=\$$(round(dailyRevenue, digits=2))",
                ylims=(-0.03*bess.storageCapacityMWh,1.03*bess.storageCapacityMWh),
                label=String(solver),
                color=colors[idx+1],
                linewidth=2,
                ytickfontcolor=:blue,
                y_guidefontcolor=:blue,
                xticks=(ticks, labels)
            )
        end

        # plot price on right axis
        p2 = twinx(p)
        plot!(p2, vals[1:end-1], day.DollarsPerMegawattHour,
        ylabel="Price (\$/MWh)", 
        legend=false, 
        color=:red,
        linestyle=:solid,
        linewidth=2, 
        ytickfontcolor=:red,
        y_guidefontcolor=:red,
        xticks=:none) # Hides overlapping x-ticks from the second axis

        # format right y-axis
        pmax = maximum(day.DollarsPerMegawattHour)
        if plim < pmax
            ylims!(p2, (-0.03*pmax, 1.03*pmax))
            annotate!(p2, [(12, pmax/2, text("ADJUSTED PRICE SCALE"))])
        else
            ylims!(p2, (-0.03*plim,1.03*plim))
        end

        display(p)
    end
end

#TODO: performance with simulated data?






function bellman_visual_check_grid(case::BESTCase, policy::String, xdim::Int, seriesdim::Int, subplotdim::Int; highlight_peaks::Bool=true, plot_derivatives::Bool=false, save_image::Bool=false)
    # Dimension convention
    # 1 = time
    # 2 = storage
    # 3 = price

    # Reject incorrect dim inputs
    if Set((xdim, seriesdim, subplotdim)) != Set((1, 2, 3))
        throw(ArgumentError("Inputs must be a permutation of 1, 2, and 3 with no repetition."))
    end

    # Get BellmanVals
    V = case.policies[Symbol(policy)].BellmanVals

                
    # Plot derivative dVdxdim instead
    if plot_derivatives
        V=diff(V, dims=xdim)
    end

    # Set all axes limits
    global_ylims = (
        minimum(V),
        maximum(V)
    )

    # Dimension
    dim_names = ["t", "y", "π"]

    # Get representative values for longer sets
    dim_vals = [
        # round.(Int, range(1, case.T, length=case.PBANDS)),
        collect(1:30:case.T),
        round.(Int, range(0, case.bess.storageCapacityMWh, length=case.PBANDS)),
        collect(1:case.PBANDS)
    ]

    peak_times = Dict(
        91  => (label = "Morning Peak", color = colorant"gold"),
        211 => (label = "Evening Peak", color = colorant"dodgerblue"),
    )

    # Set color gradient
    colors = cgrad([:red, :green], case.PBANDS, categorical=true)
    get_series_color(seriesval, j) = highlight_peaks && seriesdim == 1 && haskey(peak_times, seriesval) ?
        peak_times[seriesval].color :
        colors[j]

    # Get vals
    xvals       = dim_vals[xdim]
    seriesvals  = dim_vals[seriesdim]
    subplotvals = dim_vals[subplotdim]

    # Figure layout
    l = @layout [
        ylab{0.002w} grid(2,5) leg{0.05w}
        _ xlab{0.002h} _
    ]

    # define plot
    p = plot(layout=l, size=(1850,700))

    # Empty left panel to fit y-axis titles
    plot!(p[1], framestyle=:none, ticks=nothing)

    # plot V
    for (i, subplotval) in enumerate(subplotvals)

        # Logic for whether to display axes ticks and titles 
        show_x = (i - 1) ÷ 5 == 1
        show_y = (i - 1) % 5 == 0

        # Append peak label to subplot title if subplot is t
        subplot_title = "$(dim_names[subplotdim]) = $subplotval"
        if highlight_peaks && subplotdim == 1 && haskey(peak_times, subplotval)
            subplot_title *= " ($(peak_times[subplotval].label))"
        end

        # Format subplot i
        plot!(
            p[i+1],
            title = subplot_title,
            xlabel = show_x ? dim_names[xdim] : "",
            ylabel = show_y ? "V" : "",
            xformatter = show_x ? identity : _ -> "",
            yformatter = show_y ? identity : _ -> "",
            legend = false,
        )

        # plot V for suplot i
        for (j, seriesval) in enumerate(seriesvals)
            
            # Draw series in peak colour if series is t
            line_color = get_series_color(seriesval, j)

            # Format V for subplot and series
            inds = Any[Colon(), Colon(), Colon()]
            storage_shift(dim, val) = dim == 2 ? val + 1 : val
            inds[seriesdim]  = storage_shift(seriesdim,  seriesval)
            inds[subplotdim] = storage_shift(subplotdim, subplotval)
            y = vec(V[inds...])

            # draw series
            plot!(
                p[i+1],
                y,
                color = line_color,
                label = false,
                ylims=global_ylims,
            )

            # Draw verctical lines at peak times if x axis is time
            if highlight_peaks && xdim == 1
                for (t_peak, (_, peak_color)) in peak_times
                    vline!(
                        p[i+1],
                        [t_peak],
                        color = peak_color,
                        linewidth = 1,
                        linestyle = :dash,
                        label = false,
                        linealpha = 0.2,
                    )
                end
            end
        end
    end

    # Legend panel

    # Draw fake series for legend to recognise
    for (j, seriesval) in enumerate(seriesvals)
        # Draw series in peak colour if series is t
        line_color = get_series_color(seriesval, j)
        
        plot!(
            p[12],
            [NaN],
            [NaN],
            color = line_color,
            label = string(seriesval),
        )
    end

    # display legend
    plot!(
        p[12],
        framestyle = :none,
        grid = false,
        legend = :left,
        legend_title = dim_names[seriesdim],
        ticks = nothing,
    )

    # Empty left panel to fit x-axis titles
    plot!(p[13], framestyle=:none, ticks=nothing)

    # return complete plot
    if save_image
        filename = "plots/" * policy * "_" * dim_names[xdim] * "_" * dim_names[seriesdim] * "_" * dim_names[subplotdim] * ".pdf"
        savefig(p, filename)
    else
        display(p)
    end
end


function bellman_visual_check_gif( case::BESTCase, policy::String, xdim::Int, seriesdim::Int, animdim::Int)

    # Dimension convention
    # 1 = time
    # 2 = storage
    # 3 = price

    if Set((xdim, seriesdim, animdim)) != Set((1,2,3))
        throw(ArgumentError("Inputs must be a permutation of 1, 2, and 3."))
    end

    @warn "bellman_visual_check_gif() can take a few minutes to run."

    V = case.policies[Symbol(policy)].BellmanVals

    global_ylims = (
        minimum(V),
        maximum(V)
    )

    dim_names = ["t", "y", "π"]

    # Full dimensions for animation
    full_vals = [
        collect(1:case.T),
        collect(0:case.bess.storageCapacityMWh),
        collect(1:case.PBANDS)
    ]

    # Representative values for series only
    rep_vals = [
        round.(Int, range(1, case.T, length=10)),
        round.(Int, range(0, case.bess.storageCapacityMWh, length=10)),
        collect(1:case.PBANDS)
    ]

    xvals      = full_vals[xdim]
    animvals   = full_vals[animdim]
    seriesvals = rep_vals[seriesdim]

    colors = cgrad( [:red, :green], length(seriesvals), categorical=true)

    anim = @animate for animval in animvals

        println("$animval")

        p = plot(
            xlabel = dim_names[xdim],
            ylabel = "V",
            title = "$(dim_names[animdim]) = $animval",
            legend = :outerright,
            legend_title = dim_names[seriesdim],
            ylims = global_ylims,
        )

        for (j, seriesval) in enumerate(seriesvals)

            inds = Any[Colon(), Colon(), Colon()]
            storage_shift(dim, val) = dim == 2 ? val + 1 : val
            inds[seriesdim]  = storage_shift(seriesdim, seriesval)
            inds[animdim] = storage_shift(animdim, animvals)
            y = vec(V[inds...])

            plot!(
                p,
                y,
                color = colors[j],
                label = string(seriesval),
            )
        end
    end

    gif(anim, "output.gif", fps=20)
end