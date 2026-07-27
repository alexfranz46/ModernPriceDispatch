using DataFrames
using Dates 
using TimeZones
using Serialization


function interpolate_missing_data(g)
    # identify which timestamps are missing
    if mod(g.TradingPeriod, 2) == 0
        #0,5,...25
        mins = [0.0, 5.0, 10.0, 15.0, 20.0, 25.0]
    else
        mins = [30.0, 35.0, 40.0, 45.0, 50.0, 55.0]
    end

    #TODO: need to get the correct vals from zdt values in g

    # s = 0
    # h = hour(zdt)
    # d = day(zdt)
    # m = month(zdt)
    # y = year(zdt)
    # tz = timezone(zdt)

    times = [ZonedDateTime(y,m,d,h,mins[i],tz) for i in 1:6]

    # append those timestamps with mean
    val = mean(g.DollarsPerMegawattHour)
    for pdt in times
        push!(g, [g.TradingDate, g.TradingPeriod, pdt, val])
    end

    # sort
    sort!(g)
    return g
end


function remove_extra_data(g)    
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


function clean_DEP_data(df)
    # initialize returned df
    result = DataFrame()

    # iterate over every group
    for g in groupby(df, [:TradingDate, :TradingPeriod])

        n = nrow(g)

        if n == 6
            append!(result, g)
        elseif n >= 7
            while nrow(g) != 6
                g = remove_extra_data(g)
            end
            append!(result, g)
        elseif n <= 5
            append!(result, g)
        end
    end

    return result
end
