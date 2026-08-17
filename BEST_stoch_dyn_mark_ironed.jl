using JuMP
import HiGHS
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

# # GET RAW DATA

# # fetches data from serial instead of running if data.jls exists
# df_file_path = joinpath("Serials", "df.jls")

# if isfile(df_file_path)
#     # Load the existing file
#     df = deserialize(df_file_path)
#     println("Loaded clean data from cache.")
# else
#     # collect all CSV files
#     csv_files = String[]
#     root = joinpath(@__DIR__, "DispatchEnergyPrices")
#     for (path, _, files) in walkdir(root)
#         for f in files
#             endswith(f, ".csv") && push!(csv_files, joinpath(path, f))
#         end
#     end
#     println("Importing data from $(length(csv_files)) files.")

#     # Combine all CSVs into one DataFrame
#     df = vcat([CSV.read(file, DataFrame) for file in csv_files]...)

#     # Keep only OTA2201 connection point data
#     df = filter(:PointOfConnection => p -> p == "OTA2201", df)

#     # format datetime
#     fmt = dateformat"yyyy-mm-ddTHH:MM:SS.sssz"
#     df.PublishDateTime = ZonedDateTime.(String.(df.PublishDateTime),fmt)

#     # save data as serials
#     mkpath(dirname(df_file_path))
#     serialize(df_file_path, df)
#     println("Saved raw data to cache.")
# end

# # TODO: CHANGE LOGIC TO INCORPORATE THEM instead
# # remove daylight saving trading period 49 and 50
# df = filter(:TradingPeriod => tp -> tp <= 48, df)
# # df = filter(:DollarsPerMegawattHour => p -> p <= 500, df)

# # Keep only date, trading period and price data
# select!(df, [:TradingDate, :TradingPeriod, :PublishDateTime, :DollarsPerMegawattHour])

# println("Finished tidying data.")

# # DEFINE MARKOV PRICE PROCESS

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

# # Dynamic Program

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
    println("Solving stage $stage")

    tp = ceil(Int, stage/6)

    # constants
    p = PriceVals[:, tp]  # price for each band
    m = TransitionMatrix[:, :, tp]  # transition probabilities for each bands [i -> j]
    V = BellmanVals[stage + 1, :, :]  # Bellman vals for each battery SoC and band
    # eV = V * m'  # expected bellman vals for each battery SoC and band



    eV = similar(V)
    for e in 0:E, i in 1:numBands
        eV[e+1,i] = sum(V[e+1,j] * m[i,j] for j in 1:numBands)
    end

    # if stage == 275
    #     vscodedisplay(V, "V")
    #     vscodedisplay(eV, "eV")
    #     vscodedisplay(m, "m")
    #     vscodedisplay(p, "p")
    # end

    # @show any(isnan, V)
    # @show any(isinf, V)
    # @show any(isnan, eV)
    # @show any(isnan, m)

    # ...

    startTime = time()


    for yIn in 0:E  # y(t-1)
        
        
        # bound for y based on y stage varaible
        y⁻ = max(0, yIn - r)
        y⁺ = min(E, yIn + s)
        yRange = y⁻:y⁺

        # LP
        model = Model(HiGHS.Optimizer)
        set_silent(model) # disable solver printing

        # Variables
        @variable(model, b[1:numBands, yRange], Bin)

        # keep track of objective value for each band
        @expression(model, objBand[i in 1:numBands], sum(b[i,e] * (p[i] * (yIn - e) + eV[e+1, i]) for e in yRange))

        # keep track of optimal SoC determined from binary variable
        @expression(model, y[i=1:numBands], sum(e*b[i,e] for e in yRange))

        # constraints
        @constraint(model, monotonous[i in 2:numBands], y[i - 1] >= y[i])  # >= because decreasing SoC at the end of period means more discharge at high price and more charge at low prices
        @constraint(model, unique[i in 1:numBands], sum(b[i,e] for e in yRange) == 1)  # only one b=1 all others =0 for each i

        # objective
        @objective(model, Max, sum(objBand))
        
        if stage == 275 && yIn == 54
            print(model) 
        end
        

        # SOLVE
        optimize!(model)  # NOTE: 90% of solve time is in this optimization step...

        if !is_solved_and_feasible(model)
            error("Solver did not find an optimal solution for stage $stage and SoC $yIn")
        end

        # update bellmans and decision
        BellmanVals[stage, yIn + 1, :] = value.(objBand)
        yDecision[stage, yIn + 1, :] = round.(Int, value.(y))
        uvDecision[stage, yIn + 1, :] = yIn .- round.(Int, value.(y))

        if !issorted(yIn .- round.(Int, value.(y)))
            println("t=$stage, SoC=$yIn, y=$(yIn .- round.(Int, value.(y)))")
        end

    end 
    
    endTime = time()
    println("average LP solve time is $((endTime-startTime)/(E+1))")
end

if false
    Bellman_file_path = joinpath("Serials", "BellmanVals2.jls")
    uvDecision_file_path = joinpath("Serials", "uvDecision2.jls")
    serialize(Bellman_file_path, BellmanVals)
    serialize(uvDecision_file_path, uvDecision)
end

i0 = find_band(110.0, PriceBounds[:,TP])
finalObj = sum(TransitionMatrix[i0, :, TP].*BellmanVals[1, y0 + 1, :])
println("expected=$(round(Int, finalObj))")



# # # Compare

# # load from original markovian
# Bellman_file_path = joinpath("Serials", "BellmanVals.jls")
# uvDecision_file_path = joinpath("Serials", "uvDecision.jls")
# BellmanValsOG = deserialize(Bellman_file_path)
# uvDecisionOG = deserialize(uvDecision_file_path)

# # count differnce decision
# countIroned = 0
# countChanged = 0
# for ti in t
#     for e in 1:E+1
#         slice = uvDecision[ti,e,:]
#         slice2 = uvDecisionOG[ti,e,:]
#         if !issorted(slice2)
#             countIroned = countFix + 1
#         elseif slice != slice2
#             countChanged = countChanged + 1
#         end
#     end
# end

# println("uvDecision changes because of monotonous requirement")
# println("$(countIroned) ironed")
# println("$(countChanged) otherwise changed")

