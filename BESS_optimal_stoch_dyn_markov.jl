using DataFrames
using CSV
using StatsBase
using Dates 
using TimeZones
using Random
using Distributions
using Serialization
using Plots

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

# TODO: CHANGE LOGIC TO INCORPORATE THEM instead
# remove daylight saving trading period 49 and 50
df = filter(:TradingPeriod => tp -> tp <= 48, df)

# Keep only date, trading period and price data
select!(df, [:TradingDate, :TradingPeriod, :DollarsPerMegawattHour])

# TODO: Need to make sure that all rows are sorted by time... 
# PublishDateTime not a good way to do it, as they are sometimes all the same (see first file)

println("Finished tidying data.")

# Define price bands
# For band i: lb = bounds[i], ub = bounds[i+1]
# since prices tend to cluster for low values, p is in a band if: lb < p <= ub
numBands = 10
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

# serialize("TransitionMatrix.jls", TransitionMatrix)
# serialize("PriceVals.jls", PriceVals)

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
yDecision = fill(-1, T, E+1, numBands)  # first numBands is pNext, second is stage var pNow

for stage in reverse(t)
    tp = ceil(Int, stage/6)

    for i in 1:numBands  # p(t)
        # price stage variable
        pNow = PriceVals[i, tp]

        for yIn in 0:E  # y(t-1)
            # Lower bound for y based on y stage varaible
            y⁻ = max(0, yIn - r)
            y⁺ = min(E, yIn + s)
            yRange = y⁻:y⁺

            # initialize solution
            objSim = fill(Inf, y⁺-y⁻+1)

            for (yIdx, yNext) in enumerate(yRange)  # y(t+1)
                # calculate expected bellman value (here-and-now over j)
                expBellmanVal = sum(BellmanVals[stage + 1, yNext + 1, :] .* TransitionMatrix[i, :, tp])  # This is NaN

                # record 
                objSim[yIdx] = pNow*(yIn - yNext) + expBellmanVal
            end

            # find and record best decision
            BellmanVals[stage, yIn + 1, i], optIdxPair = findmax(objSim)
            yDecision[stage, yIn + 1, i] = yRange[optIdxPair[1]]
        end
    end 
end

# Run and print policy

# set common random numbers
Seeds = 1:10

for seedNum in Seeds
    Random.seed!(seedNum)

    # Define probability distributions
    dists = Matrix{DiscreteNonParametric}(undef, numBands, TP)
    for tp in 1:TP
        for i in 1:numBands
            # Set values and probabilities for price band and trading period
            dists[i, tp] = DiscreteNonParametric(1:numBands, TransitionMatrix[i, :, tp])
        end
    end

    # Simulate
    objective = 0
    i = ceil(Int, numBands/2)  # starting in central band
    yHistory = [0]  # starting with an empty battery
    pHistory = [PriceVals[i, TP]]
    for stage in t
        yIn = yHistory[stage]
        
        # Update trading period
        tp = ceil(Int, stage/6)

        # fetch and record optimal decision for stage/state
        yNext = yDecision[stage, yIn + 1, i] 
        push!(yHistory, yNext)
        
        # fetch and record price for stage/state
        pNow = PriceVals[i, tp]
        push!(pHistory, pNow)

        # calculate objective
        objective += pNow*(yIn - yNext)

        # use random distribution to find j, which becomes i for next iteration
        j = rand(dists[i, tp])
        i = j
    end

    # Plot results 
    # Plot charge on left axis
    p = plot(t, yHistory[2:end], 
        ylabel="Battery charge (MWh)", 
        legend=:topleft, 
        label="y", 
        color=:blue,
        linewidth=2,
        title="seed=$seedNum")

    # Plot price on right axis
    plot!(twinx(), t, pHistory[2:end],
        ylabel="Price (\$/MWh)", 
        legend=:topright, 
        label="p", 
        color=:red,
        linestyle=:solid,
        linewidth=2, 
        xticks=:none) # Hides overlapping x-ticks from the second axis

    display(p)
end