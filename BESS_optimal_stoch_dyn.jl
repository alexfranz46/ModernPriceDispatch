using Plots
using Random
using Distributions

#=  Stochastic dynamic method to optimize BESS usage in 5min increments 
    while ignoring gate closure.
=#

# v1: discharging-only, optimize over 1 hour

T = 12*24
t = 1:T

# battery
E = 240  # max bettery storage (MWh)
r = 10  # max discharge    # TODO Modify to be every 5 minutes (MWh/5min)
s = 10  # max charge       # TODO Modify to be every 5 minutes (MWh/5min)
y0 = 0  # initial battery charge (MWh)

# Prices 1 hour 
# based on RTP 1/07/2026, OTA2201, tp=38&39, rounded to nearest 10
# p0 = 140  # last RTP from tp=37
# p = [150, 150, 110, 140, 100, 110, 150, 110, 150, 110, 100, 100]
#      1    2    3    4    5    6    7    8    9    10   11   12

# Prices full day
# based on RTP 1/07/2026, OTA2201, rounded to nearest 10
# Price means calculated by trading period, noise realizations every 5 minutes
p0 = 110
# pTP = [ 97, 100, 102, 100, 100, 100, 100, 100, 100, 103, 100, 103, 
#        113, 133, 132, 172, 157, 133, 115, 110,  93,  90,  90,  78, 
#         85,  80,  80,  65,  62,  37,  80,  88,  95, 100,  98, 100, 
#        117, 127, 120,  95,  93,  90, 103,  98,  90,  92,  88,  83]
pActual = [110, 105, 104, 101,  82,  82, 104, 104, 104, 104, 104, 104,
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
            93, 109,  93,  85,  84,  84,  86,  92,  85,  84,  84,  83]


# Noise
# ϵ = 5
# pNoise = [-20, -10, 0, 10, 20]
# probNoise = fill(0.2, 5)

# (deterministic version)
# ϵ = 1
# pNoise = [0]
# probNoise = [1]

ϵ = 5
pNoise = [-20, -10, 0, 10, 20]
probNoise = [0.2, 0.2, 0.2, 0.2, 0.2]

# Matrix to store expected bellman function values
BellmanVals = fill(Inf, T + 1, E+1, ϵ)
BellmanVals[end,:,:] .= 0  # (termination condition)

# Create matrices to store the optimal Q(t), y(t) for each stage, Q(t-1), y(t-1)
yDecision = fill(-1, T, ϵ, E+1) 

# MAIN LOOP

for stage in reverse(t)
    TP = ceil(Int, stage/6)

    for θ in 1:ϵ  # p(t)
        # price stage variable
        # pNow = pTP[TP] + pNoise[θ]  # version for 30-minute pricing
        pNow = pActual[stage] + pNoise[θ]

        for yNow in 0:E  # y(t)
            # Lower bound for y based on y stage varaible
            y⁻ = max(0, yNow - r)
            y⁺ = min(E, yNow + s)
            yRange = y⁻:y⁺

            # initialize solution
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
                expBellmanVal = sum(BellmanVals[stage + 1, yNext + 1, :] .* probNoise)

                # record 
                objSim[yIdx] = pNow*(u - v) + expBellmanVal
            end

            # find and record best decision
            BellmanVals[stage, yNow + 1,  θ], optIdxPair = findmax(objSim)
            yDecision[stage, θ, yNow + 1] = yRange[optIdxPair[1]]
        end
    end 
end

finalObj = sum(probNoise.*BellmanVals[1, y0 + 1, :])
println("expected=$(round(Int, finalObj))")

# no noise simulation
caseObjective = 0
yHistory = [y0]
for i in t
    push!(yHistory, yDecision[i, 3, yHistory[end]+1])
    # caseObjective += pTP[ceil(Int, i/6)]*(yHistory[end-1] - yHistory[end])  # verion for 30-minute pricing
    caseObjective += pActual[i]*(yHistory[end-1] - yHistory[end])
end
println("1Jul=$caseObjective")


# format x axis plot
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
# plot!(twinx(), t/12, repeat(pTP, inner=6),  # version for 30-minute pricing
plot!(twinx(), t/12, pActual,
    ylabel="Price (\$/MWh)", 
    legend=:topright, 
    label="p", 
    color=:red,
    linestyle=:solid,
    linewidth=2, 
    xticks=:none) # Hides overlapping x-ticks from the second axis

display(p)


# # mass simulation

# # set common random numbers
# Random.seed!(1)

# # simulation samples
# totalSamples = 10000

# # Set noise array
# noiseDist = DiscreteNonParametric(pNoise, probNoise)

# # initialize return arrays
# simObj = fill(Inf, totalSamples)

# for i in 1:totalSamples
#     # determine demand for sample
#     sampleNoise = rand(noiseDist, T)
#     pReal = pActual + sampleNoise
    
#     # @show(sampleNoise)

#     # initialize 
#     yPrior = y0
#     tempObj = 0
#     tempQCurt = 0

#     # Solve per stage
#     for stage in t
#         # @show(stage)
        
#         # extract optimal decision from policy
#         noiseIdx = findfirst(==(sampleNoise[stage]), dNoise)

#         yCurrent = yDecision[stage, noiseIdx, yPrior + 1, noiseIdx]

#         # add obj to counter
#         tempObj += 0
#     end

#     simObj[i] = tempObj
# end








## NOTE ##
# confident this implementation is correct because, regardless of where it starts, answers converge and follow the same 'path'.