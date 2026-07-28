using DataFrames
using StatsBase
using Dates 
using TimeZones
using Serialization


function interpolate_missing_data(g)
    missingCount = 6 - nrow(g)

    recordedMins = [minute(zdt) + second(zdt)/60 for zdt in g.PublishDateTime]
    
    # identify which timestamps are missing
    if isodd(g.TradingPeriod[1])  # odd trading periods
        expectedMins = [0, 5, 10, 15, 20, 25]
    else                          # even trading periods
        expectedMins = [30, 35, 40, 45, 50, 55]
    end

    # Circular distance on a 60-minute clock (aka 8:59.55s is close to 09:00.00s)
    circ_dist(a, b) = begin
        d = abs(a - b)
        min(d, 60 - d)
    end

    # return time difference to closet expected for each recorded
    distances = [minimum(circ_dist.(recordedMins, expectedMins[i])) for i in 1:6]
    # MAYBE TODO: it is possible under this method for two observations to be assigned to the same expected time: 
    # eg: 00:04.59s, 00:05.10s, 00:10.20s -> 1st element closest, 2nd element second closest
    #                                     -> BUT: 2nd element is close to 00:05.00 slot, already "taken" by 1st element
    #                                     -> THEREFORE: 3rd element, while further from the expected it is closest to, is
    #                                        in actuality closer to 00:10.00 (which is not "taken") than the second.

    # find missing indices
    idxs = partialsortperm(distances, 1:missingCount, rev=true)

    # convert missing indices to missing times
    date = g.TradingDate[1]
    tp = g.TradingPeriod[1]
    h = (tp - 1) ÷ 2
    s = 0
    tz = timezone(g.PublishDateTime[1])
    times = [
        ZonedDateTime(year(date), month(date), day(date), h, expectedMins[i], s, tz) 
        for i in idxs
    ]

    # append those timestamps with mean
    val = mean(g.DollarsPerMegawattHour)

    # Build rows to append
    newRows = DataFrame(
        TradingDate = fill(date, missingCount),
        TradingPeriod = fill(tp, missingCount),
        PublishDateTime = times,
        DollarsPerMegawattHour = fill(val, missingCount),
    )

    newG = vcat(g, newRows)

    # sort
    sort!(newG, :PublishDateTime)
    return newG
end


function remove_extra_data(g)    
    # set ideal interval between observations
    target = Millisecond(Minute(5)).value

    # sort list
    sort!(g.PublishDateTime)

    # calculate the interavals of time between each observation and find worst
    intervals = [g.PublishDateTime[i+1] - g.PublishDateTime[i] for i in 1:nrow(g)-1]
    errors = [abs(inv.value - target) for inv in intervals]
    worstIdx = argmax(errors)

    # The anomaly is caused by either observation i or observation i+1
    candA = worstIdx       # Left observation
    CandB = worstIdx + 1   # Right observation
    
    # remove observations from copy
    filteredA = deleteat!(copy(g), candA)
    filteredB = deleteat!(copy(g), CandB)

    # calculate new intervals
    intervalsA = [filteredA.PublishDateTime[i+1] - filteredA.PublishDateTime[i] for i in 1:nrow(g)-2]
    intervalsB = [filteredB.PublishDateTime[i+1] - filteredB.PublishDateTime[i] for i in 1:nrow(g)-2]

    # calculate new errors
    errA = [abs(inv.value - target) for inv in intervalsA]
    errB = [abs(inv.value - target) for inv in intervalsB]

    # trim group and retrun
    if maximum(errA) < maximum(errB)
        return filteredA
    else
        return filteredB
    end
end


function remove_mislabelled_data(g)
            # Check all entries are from the same trading period
        date = g.TradingDate[1]
        tp = g.TradingPeriod[1]
        h = (tp - 1) ÷ 2
        tz = timezone(g.PublishDateTime[1])

        # Compute expected minute range
        if isodd(tp)
            expectedMins = (0, 25)
        else
            expectedMins = (30, 55)
        end

        # Build lower and upper expected timestamps
        lb = ZonedDateTime(year(date), month(date), day(date), h, expectedMins[1], 0, tz)
        ub = ZonedDateTime(year(date), month(date), day(date), h, expectedMins[2], 0, tz)

        # Tolerance band (±2.5 minutes)
        tol = Second(150)

        # Keep timestamps within [lb - tol, ub + tol]
        keep = map(g.PublishDateTime) do zdt
            (lb - tol) <= zdt <= (ub + tol)
        end

        # Filter the group
        return g[keep, :]
end


function clean_DEP_data(df)
    # initialize returned df
    result = DataFrame()

    # handle daylight savings (total tp !=48)
    # TODO

    for g in groupby(df, [:TradingDate, :TradingPeriod])

        # identify and remove data published at the wrong time
        g = remove_mislabelled_data(g)

        # guarantee every day has exactly 288 observations
        n = nrow(g)

        if n == 6
            append!(result, g)
        elseif n >= 7
            while nrow(g) != 6
                g = remove_extra_data(g)
            end
            append!(result, g)
        elseif n <= 5
            g = interpolate_missing_data(g)
            append!(result, g)
        end
    end

    return result
end

df= deserialize("df.jls")
clean_DEP_data(df)