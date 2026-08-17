using DataFrames
using CSV
using StatsBase
using Dates 
using TimeZones
using Random
using Distributions
using Serialization
using Plots
include("clean_DEP_data.jl")

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

# fetches data from serial instead of running if data.jls exists
file_path = joinpath("Serials", "df.jls")

if isfile(file_path)
    # Load the existing file
    df = deserialize(file_path)
    println("Loaded clean data from cache.")
else
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

    # format datetime
    fmt = dateformat"yyyy-mm-ddTHH:MM:SS.sssz"
    df.PublishDateTime = ZonedDateTime.(String.(df.PublishDateTime),fmt)

    mkpath(dirname(file_path))
    serialize(file_path, df)
    println("Saved data to cache.")
end

# TODO: CHANGE LOGIC TO INCORPORATE THEM instead
# remove daylight saving trading period 49 and 50
df = filter(:TradingPeriod => tp -> tp <= 48, df)
# df = filter(:DollarsPerMegawattHour => p -> p <= 500, df)

# Keep only date, trading period and price data
select!(df, [:TradingDate, :TradingPeriod, :PublishDateTime, :DollarsPerMegawattHour])

# Clean data
df_clean = clean_DEP_data(df)
# df = clean_DEP_data(df)

if false
    CSV.write("raw_data.csv", df)
    CSV.write("clean_data.csv", df_clean)
end

if true
    count = combine(groupby(df, [:TradingDate, :TradingPeriod]), nrow => :count)   
    count_clean = combine(groupby(df_clean, [:TradingDate, :TradingPeriod]), nrow => :count)   
end

# TODO: Need to make sure that all rows are sorted by time... 
# PublishDateTime not a good way to do it, as they are sometimes all the same (see first file)

println("Finished tidying data.")

# # EXPLORE DATA

if false
    p = plot(df.PublishDateTime, df.DollarsPerMegawattHour, 
        label="original", 
        color=:blue,
        linewidth=0.5,
        m = :diamond, 
        markersize = 6, 
        markercolor = :blue
    )

    # Plot price on right axis
    plot!(df_clean.PublishDateTime, df_clean.DollarsPerMegawattHour,
        xlabel="time",
        ylabel="Price (\$/MWh)", 
        legend=:topright, 
        label="clean", 
        color=:red,
        linestyle=:solid,
        linewidth=0.5,
        m = :diamond, 
        markersize = 3, 
        markercolor = :red,
        xlim=(DateTime(df.PublishDateTime[259]),DateTime(df.PublishDateTime[265])),
        ylim=(0,0.1)
    ) # Hides overlapping x-ticks from the second axis

    display(p)
end

if false
    h = histogram(df.DollarsPerMegawattHour)
    display(h)

    arr = [0.01,0.05,0.1,0.5,1,5,10,50,100,250,500,1000,5000,10000]
    for (i,p) in enumerate(arr)
        c = nrow(filter(:DollarsPerMegawattHour => pi -> pi <= p, df))
        println("$(round(Int,100*c/nrow(df)))% of prices <= $p \$/MWh")
    end
end

# Define price bands
# For band i: lb = bounds[i], ub = bounds[i+1]
# since prices tend to cluster for low values, p is in a band if: lb < p <= ub
numBands = 10
Vals_file_path = joinpath("Serials", "PriceVals.jls")
Bounds_file_path = joinpath("Serials", "PriceBounds.jls")
Trans_file_path = joinpath("Serials", "TransitionMatrix.jls")

if (isfile(Vals_file_path) && isfile(Vals_file_path) && isfile(Trans_file_path))
    # Load the existing file
    PriceVals = deserialize(Vals_file_path)
    PriceBounds = deserialize(Bounds_file_path)
    TransitionMatrix = deserialize(Trans_file_path)
    println("Loaded Markov data from cache.")
else
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

    # go through data row-by-row and keep tally of each transition type
    rowI = nothing
    firstRow = true
    for rowJ in eachrow(df)

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

    # convert tally to probabilities
    for tp in 1:TP
        for i in 1:numBands
            # calculate the total sum of rows
            totalI = sum(TransitionMatrix[i, :, tp])
            # TODO: account for totalI==0?

            # transform count into probabilties
            TransitionMatrix[i, :, tp] ./= totalI
        end
    end

    # save data as serials

    mkpath(dirname(Vals_file_path))
    serialize(Vals_file_path, PriceVals)
    
    mkpath(dirname(Bounds_file_path))
    serialize(Bounds_file_path, PriceBounds)
    
    mkpath(dirname(Trans_file_path))
    serialize(Trans_file_path, TransitionMatrix)
    
    println("Saved Markov data to cache.")
end

println("Markovian price process fully defined.")

# check serials match dimensions
size(PriceVals) == (numBands, TP) || error("PriceVals.jls has incorrect dimensions!")
size(PriceBounds) == (numBands + 1, TP) || error("PriceBounds.jls has incorrect dimensions!")
size(TransitionMatrix) == (numBands, numBands, TP) || error("TransitionMatrix.jls has incorrect dimensions!")

if false
    plt = heatmap(TransitionMatrix[:,:,1], title="tp=1", 
        legend=false, widen=false,
        aspect_ratio=:equal, clim=(0.0 ,1.0), 
        xticks=1:10, yticks=1:10, 
        xlabel="i", ylabel="j"
    )

    # Loop through all 48 matrices and build the animation frames
    anim = @animate for i in 1:48
        heatmap!(plt, TransitionMatrix[:,:,i], title="tp=$i")
    end

    # Save the animation to a file
    gif(anim, "TransitionMatrix.gif", fps=5)
end

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
yDecision = fill(-1, T, E+1, numBands)
uvDecision = copy(yDecision)

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
                expBellmanVal = sum(BellmanVals[stage + 1, yNext + 1, :] .* TransitionMatrix[i, :, tp])

                # record 
                objSim[yIdx] = pNow*(yIn - yNext) + expBellmanVal
            end

            # find and record best decision
            BellmanVals[stage, yIn + 1, i], optIdxPair = findmax(objSim)
            yDecision[stage, yIn + 1, i] = yRange[optIdxPair[1]]
            uvDecision[stage, yIn + 1, i] = yIn - yRange[optIdxPair[1]]
        end
    end 
end

if false
    Bellman_file_path = joinpath("Serials", "BellmanVals.jls")
    uvDecision_file_path = joinpath("Serials", "uvDecision.jls")
    serialize(Bellman_file_path, BellmanVals)
    serialize(uvDecision_file_path, uvDecision)
end

if false
    plt = heatmap(uvDecision[1,:,:], title="t=1", 
        legend=false, widen=false,
        aspect_ratio=:equal, clim=(-10.0, 10.0),
        xlim=(1,numBands), ylim=(0,T),
        xticks=1:10
    )

    # Loop through all 48 matrices and build the animation frames
    anim = @animate for i in 1:T
        heatmap!(plt, uvDecision[i,:,:], title="t=$i", aspect_ratio=:equal)
    end

    # Save the animation to a file
    gif(anim, "uv.gif", fps=5)
end

i0 = find_band(110.0, PriceBounds[:,TP])
finalObj = sum(TransitionMatrix[i0, :, TP].*BellmanVals[1, y0 + 1, :])
println("expected=$(round(Int, finalObj))")

# investigate ADR monotonicity
if true
    p = plot(xlim=(1,10), xticks=1:10, ylim=(-10,10), legend=false)
    
    monotonousADR = Matrix{Bool}(undef, T, E+1)
    nonMonotonousCount = Matrix{Int}(undef, T, E+1)
    nonMonotonousBands = Matrix{Vector{Int}}(undef, T, E+1)
    for t in 1:T
        for e in 1:E+1
            monotonousADR[t, e] = issorted(uvDecision[t, e, :])
            slice = uvDecision[t, e, :]
            
            # Count how many adjacent pairs violate the sorted rule
            count = 0
            idxs = []
            for i in 1:(length(slice) - 1)
                if slice[i] > slice[i+1]
                    count += 1
                    push!(idxs, i+1)
                end
            end
                    
            nonMonotonousCount[t,e] = count
            nonMonotonousBands[t,e] = idxs

            if count != 0
                # add non monot to plot
                plot!(p, 1:10, slice)
            end
        end
    end

    display(p)

    hm = heatmap(Int.(monotonousADR), c = [:red, :green], 
        legend = false, 
        aspect_ratio = :equal, 
        xlabel="Battery Charge", 
        ylabel="5-min period", 
        xlim=(0,E), 
        ylim=(1,T)
    )

    display(hm)
end



# Simulate 1/07/2026
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

caseObjective = 0
i = find_band(110.0, PriceBounds[:,TP])  # starting band
yHistory = [0]  # starting with an empty battery
pHistory = [PriceVals[i, TP]]
for stage in t
    yIn = yHistory[stage]
    
    # Update trading period
    tp = ceil(Int, stage/6)

    # fetch and record optimal decision for stage/state
    yNext = yIn - uvDecisionOG[stage, yIn + 1, i] 
    push!(yHistory, yNext)
    
    # fetch and record price for stage/state
    pNow = PriceVals[i, tp]
    push!(pHistory, pNow)

    # calculate objective
    caseObjective += pActual[stage]*(yIn - yNext)

    # use random distribution to find j, which becomes i for next iteration
    j = iActual[stage]
    i = j
end

# Plot results 
ticks = 0:0.5:24
labels = [mod(x, 1) == 0 ? string(Int(x)) : "" for x in ticks]

# Plot storage on left axis
p = plot(t/12, yHistory[2:end], 
    xlabel="time of day (HH)",
    ylabel="Battery charge (MWh)", 
    legend=:topleft, 
    label="y", 
    color=:blue,
    linewidth=2,
    xticks=(ticks, labels))

# Plot price on right axis
plot!(twinx(), t/12, pActual,
    ylabel="Price (\$/MWh)", 
    legend=:topright, 
    label="p", 
    color=:red,
    linestyle=:solid,
    linewidth=2, 
    xticks=:none) # Hides overlapping x-ticks from the second axis

display(p)

println("1Jul=$caseObjective")


# # MASS SIMULATION

# Run and print policy

# set common random numbers
Seeds = 1:1250

# format x axis plot
ticks = 0:0.5:24
labels = [mod(x, 1) == 0 ? string(Int(x)) : "" for x in ticks]

simObjectives = []
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
    
    push!(simObjectives, objective)

    # # Plot results 
    # # Plot storage on left axis
    # p = plot(t/12, yHistory[2:end], 
    #     xlabel="time of day (HH)",
    #     ylabel="Battery charge (MWh)", 
    #     legend=:topleft, 
    #     label="y", 
    #     color=:blue,
    #     linewidth=2,
    #     xticks=(ticks, labels),
    #     title="seed=$seedNum")

    # # Plot price on right axis
    # plot!(twinx(), t/12, pHistory[2:end],
    #     ylabel="Price (\$/MWh)", 
    #     legend=:topright, 
    #     label="p", 
    #     color=:red,
    #     linestyle=:solid,
    #     linewidth=2, 
    #     xticks=:none) # Hides overlapping x-ticks from the second axis

    # display(p)
end

println("sim=$(round(Int, mean(simObjectives)))")
h = histogram(simObjectives, legend=:none)
display(h)

# # MASS SIMULATION (AGAINST REAL DATA)

# get dates which have exactly 288 datapoints
valid_dates = combine(groupby(df_clean, :TradingDate), nrow => :count)
valid_dates = Set(valid_dates.TradingDate[valid_dates.count .== 288])
df_valid = filter(row -> row.TradingDate in valid_dates, df_clean)
# ---

histObjective = Float64[]
for day in unique(df_valid.TradingDate)
    # get day data
    df_day = filter(:TradingDate => d -> d == day, df_valid)

    
    pHistory = df_day.DollarsPerMegawattHour
    
    objective = 0
    i = ceil(Int, numBands/2)  # starting in central band
    yHistory = [0]  # starting with an empty battery
    for stage in t
        yIn = yHistory[stage]
        
        # Update trading period
        tp = ceil(Int, stage/6)

        # fetch and record optimal decision for stage/state
        yNext = yDecision[stage, yIn + 1, i] 
        push!(yHistory, yNext)
        
        # fetch and record price for stage/state
        pNow = pHistory[stage]

        # calculate objective
        objective += pNow*(yIn - yNext)

        # use random distribution to find j, which becomes i for next iteration
        j = find_band(pNow, PriceBounds[:,tp])
        i = j
    end
    
    push!(histObjective, objective)
end

println("historic=$(round(Int, mean(histObjective)))")
xlimMatch = xlims(h)
h = histogram(histObjective, legend=:none, xlims=xlimMatch)  # same xlims as previous
display(h)
println("historic observations outside histogram xlims: 
    $(sum(x -> x < (xlimMatch[1]), histObjective)) below; 
    $(sum(x -> x > (xlimMatch[2]), histObjective)) above."
)

