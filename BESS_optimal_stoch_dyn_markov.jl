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

# Ensure times are sorted by published time
# df.PublishDateTime= ZonedDateTime.(df.PublishDateTime)

# Keep only date, trading period and price data
select!(df, [:TradingDate, :TradingPeriod, :DollarsPerMegawattHour])

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
        @warn "lb==ub in period $p at indices: $badBands"
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
    # ignore first row, because we don't have the data on what we transitioned from
    if firstRow
        rowI = rowJ
        firstRow = false
        continue
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

vscodedisplay(TransitionMatrix[:,:,1])

for tp in 1:TP
    for i in 1:numBands
        # calculate the total sum of rows
        totalI = sum(TransitionMatrix[i, :, tp])
        # TODO: account for totalI==0?

        # transform count into probabilties
        TransitionMatrix[i, :, tp] ./= totalI
    end
end

vscodedisplay(TransitionMatrix[:,:,1])