using DataFrames
using CSV
using StatsBase
using Dates 
using TimeZones

#=  Stochastic dynamic method to optimize BESS usage in 5min increments 
    while ignoring gate closure.
=#


""" Returns the band a price falls into
"""
function find_band(p::Float64, bounds::Vector{Float64})
    for (i, ub) in enumerate(bounds)
        if p <= ub
            return i - 1
        end
    end
    return length(bounds) + 1   # the final band (upper bound = Inf)
end

# set fineness of temporal mesh
TP = 48

# collect all CSV files
csv_files = String[]
root = joinpath(@__DIR__, "DispatchEnergyPrices")
for (path, _, files) in walkdir(root)
    for f in files
        endswith(f, ".csv") && push!(csv_files, joinpath(path, f))
    end
end
println("Importing data from $(length(csv_files)) files.")

# Combine all CSVs into one DataFrame
df = vcat([CSV.read(file, DataFrame) for file in csv_files]...)

# Keep only OTA2201 connection point data
df = filter(:PointOfConnection => p -> p == "OTA2201", df)

# Keep only date, trading period and price data
select!(df, [:TradingDate, :TradingPeriod, :DollarsPerMegawattHour])

# TODO: Need to make sure that all rows are sorted by time... 
# PublishDateTime not a good way to do it, as they are sometimes all the same (see first file)

println("Finished tidying data.")

# Define price bands
# For band i: lb = bounds[i], ub = bounds[i+1]
# since prices tend to cluster for low values, p is in a band if: lb < p <= ub
numBands = 5
deciles = range(0, 1, length = numBands+1)[2:end-1]
PriceBounds = Matrix{Float64}(undef, numBands+1, TP)
PriceBounds[1,:] .= -Inf
PriceBounds[end,:] .= Inf
PriceVals = Matrix{Float64}(undef, numBands, TP)

for tp in 1:TP
    # Prices in trading period
    prices = Vector(df[df.TradingPeriod .== tp, :DollarsPerMegawattHour])

    # price boundaries
    PriceBounds[2:end-1, tp] = quantile(prices, deciles)

    # check for bands where lb == ub, happens if not enough data
    badBands = findall(i -> PriceBounds[i, tp] == PriceBounds[i+1, tp], 1:numBands)
    if !isempty(badBands)
        @warn "lb==ub in period $tp at indices: $badBands"
    end

    # price values
    for i in 1:numBands
        lb = PriceBounds[i, tp]
        ub = PriceBounds[i+1, tp]
        PriceVals[i, tp] = mean(prices[(prices .> lb) .& (prices .<= ub)])
    end
end

# Define transition Matrix (i -> j for each t)
TransitionMatrix = zeros(Float64, numBands, numBands, TP)
#                                 band i,   band j,   t

# TODO: NEW LOOP FORMAT
# for row in ordered_df, what is current transition, classify it,
# end loop, divide by totals for each row/column (whichever makes sense) to get probs
rowI = nothing
firstRow = true
for rowJ in eachrow(df)
    # DEBUG ONLY TODO remove
    global firstRow, rowI # Tells the loop to use the outer variables

    # ignore first row, because we don't have the data on what we transitioned from
    if firstRow
        rowI = rowJ
        firstRow = false
        continue  # skip rest of iteration
    end

    # identify the indices of this occurence
    tpI = rowI.TradingPeriod
    tpJ = rowJ.TradingPeriod
    i = find_band(rowI.DollarsPerMegawattHour, PriceBounds[:, tpI])
    j = find_band(rowJ.DollarsPerMegawattHour, PriceBounds[:, tpJ])

    # Store result
    TransitionMatrix[i, j, tpI] += 1

    # transition to next t: r(t+1) becomes new r(t)
    rowI = rowJ
end

for tp in 1:TP
    for i in 1:numBands
        # calculate the total sum of rows
        totalI = sum(TransitionMatrix[i, :, tp])
        # TODO: account for totalI==0?

        # transform count into probabilties
        TransitionMatrix[i, :, tp] ./= totalI
    end
end

println("Markovian price process fully defined.")

# Code below is mostly copied form BESS_optimal_stoch_dyn.jl

# With out Markovian process, instead of iterating over noise, we iterate of each band

# Time horizon
T = 12*24
t = 1:T

# battery
E = 240  # max bettery storage (MWh)
r = 10  # max discharge (MWh/5min)
s = 10  # max charge (MWh/5min)
y0 = 0  # initial battery charge (MWh)

# Matrix to store expected bellman function values
BellmanVals = fill(Inf, T + 1, E+1, numBands)
BellmanVals[end,:,:] .= 0  # (termination condition)

# Create matrices to store the optimal y(t) for each stage, y(t-1) and p
yDecision = fill(-1, T, numBands, E+1, numBands)  # first numBands is pNext, second is stage var pNow

for stage in reverse(t)
    tp = ceil(Int, stage/6)

    for i in 1:numBands  # p(t)
        # price stage variable
        pNow = PriceVals[i, tp]

        for yNow in 0:E  # y(t)
            # Lower bound for y based on y stage varaible
            y⁻ = max(0, yNow - r)
            y⁺ = min(E, yNow + s)
            yRange = y⁻:y⁺

            # initialize solution
            tempObj = fill(Inf, numBands)
            yOptimal = fill(Inf, numBands)

            for j in 1:numBands  # p(t+1) [NOTE: Independent of pNow]
                # TODO this price could be in a different tp...
                # TODO Should this be a here and now decision?
                # TODO 
                
                # matrix to record cost to go before expectation
                objSim = fill(Inf, y⁺-y⁻+1)

                for (yIdx, yNext) in enumerate(yRange)  # y(t+1)
                    # find discharge
                    yDiff = yNext - yNow
                    if yDiff == 0
                        u = 0
                        v = 0
                    elseif yDiff >= 0
                        u = 0
                        v = yDiff
                    elseif yDiff <= 0
                        u = -yDiff
                        v = 0
                    end

                    # lookup new state expected bellman value
                    bellmanVal = BellmanVals[stage + 1, yNext + 1, j]

                    # record 
                    objSim[yIdx] = pNow*(u - v) + bellmanVal
                end

                # TODO: i think this section doesnt make sense for this setup
                # this maximization is essentially looking at which future price is best for us...
                # which is not something we can control... does that mean we should simply do the
                # expectation below without finding yOptimal?

                # find optimal solution for this future price
                tempObj[j], optIdxPair = findmax(objSim)
                yOptimal[j] = yRange[optIdxPair[1]]
            end

            # Find expected bellman value
            optObj = sum(tempObj .* TransitionMatrix[i, :, tp])

            # DEBUG:
            println("$(findmax(abs.(yOptimal.-mean(yOptimal)))[1])")

            # Store results
            BellmanVals[stage, yNow + 1, i] = optObj
            yDecision[stage, :, yNow + 1, i] = yOptimal
        end
    end 
end

i0 = ceil(Int, numBands/2)
finalObj = sum(BellmanVals[1, y0 + 1, :] .* TransitionMatrix[i0, :, TP])
println(finalObj)

yPrev = y0
println("t=0: y=$y0")
for i in t
    # vscodedisplay(BellmanVals[i,:,:], "t=$i")
    # @show(yDecision[i, 3, yPrev+1, 3])
    println("t=$i: y=$(yDecision[i, 1, yPrev+1, 1])")
    yPrev = yDecision[i, 1, yPrev+1, 1]
    # vscodedisplay(yDecision[i, 3, :, 3], "t=$i")
end