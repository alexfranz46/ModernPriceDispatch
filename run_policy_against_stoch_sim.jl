# include("first_try_ADR_deterministic.jl")
include("first_try_ADR_stoch_dyn.jl")
using Random
using UnPack
using Distributions
using Plots

#=
use dynamic recursion to sovle stochastic demand problem:
    d +- ε = {-4, -2, 0, 2, 4}, with equal probabilities.
=#

function solve_using_stochastic_policy(totalSamples::Int, noiseDist::Distribution, QPolicy::Array{}, yPolicy::Array{}, params::modelParameters)
    # unpack parameters
    @unpack T, t, d, η, L, Q0, y0, dNoise = params
    
    # initialize return arrays
    simObj = fill(Inf, totalSamples)
    QCurt = fill(0, totalSamples)

    for i in 1:totalSamples
        # determine demand for sample
        sampleNoise = rand(noiseDist, T)
        dReal = d + sampleNoise
        
        # @show(sampleNoise)

        # initialize 
        QPrior = Q0
        yPrior = y0
        tempObj = 0
        tempQCurt = 0

        # Solve per stage
        for stage in t
            # @show(stage)
            
            # extract optimal decision from policy
            noiseIdx = findfirst(==(sampleNoise[stage]), dNoise)
            QCurrent = QPolicy[stage, noiseIdx, QPrior + 1, yPrior + 1]
            yCurrent = yPolicy[stage, noiseIdx, QPrior + 1, yPrior + 1]

            # @show(QCurrent)
            # @show(yCurrent)

            # decompose y
            yDiff = yCurrent - yPrior
            u, ηv = decompose_y(yDiff)

            # find lost load and curtailed load
            zEql = dReal[stage] - QCurrent - u + ηv/η

            # @show(zEql)

            # assume excess generation is curtailled
            # TODO: This is likely the part which is causing the error...
            if zEql < 0
                tempQCurt -= zEql
                zEql = 0
            end

            # add obj to counter
            tempObj += gen_cost(QCurrent, params) + L*zEql

            # @show(tempObj)

            # update y(t-1)
            QPrior = QCurrent
            yPrior = yCurrent
        end

        simObj[i] = tempObj
        QCurt[i] = tempQCurt
    end

    return simObj, QCurt
end;


function run()
    # set common random numbers
    Random.seed!(1)

    # Set noise array
    noise = [-4, -2, 0, 2, 4]
    prob = [0.2, 0.2, 0.2, 0.2, 0.2]
    noise = [-4, -2, 0, 8, 16]
    prob = [0.2, 0.2, 0.5, 0.05, 0.05]
    noiseDist = DiscreteNonParametric(noise, prob)
    
    # Define model parameters
    params = modelParameters(
        24,                                             # T
        1:24,                                           # t
        [40, 41, 42, 43, 35, 40, 40, 25, 10, 8, 6, 5,
        5, 6, 8, 10, 20, 30, 55, 72, 75, 70, 64, 60],   # d

        10,                                             # B
        [10, 20, 30, 40, 50, 70, 90, 110, 150, 200],    # mc
        [5, 5, 5, 5, 5, 5, 10, 10, 10, 10],             # bmax

        70,     # Qmax
        10,     # rho

        30,     # E
        1,      # η
        15,     # r
        15,     # s
        
        1000,   # L
        
        35,     # Q0
        0,      # y0

        5,      # noiseScenarios
        noise,  # dNoise
        prob    # probNoise
    )

    # simulation samples
    totalSamples = 10000

    # initialize results storing array
    simObj = fill(undef, totalSamples)

    # POLICIES

    # # Deterministic policy (NOTE: Requires include("filename.jl") statement)
    # objDeterm, QPolicy, yPolicy = solve_deterministic_DP(params);

    # Stochastic
    objExpect, QPolicy, yPolicy = solve_stochastic_DP(params)

    # HARDCODED VERSIONS

    # # Deterministic
    # objDeterm = 48470.0
    # QPolicy = [40, 41, 42, 43, 35, 40, 40, 25, 10, 8, 6, 5, 5, 7, 10, 19, 29, 39, 49, 59, 64, 70, 64, 60]
    # yPolicy = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 3, 12, 21, 30, 24, 11, 0, 0, 0, 0]

    # # Stochastic demand +-4
    # objDeterm = 57510.0
    # QPolicy = [44, 45, 46, 47, 40, 43, 44, 29, 14, 12, 10, 10, 12, 15, 15, 20, 30, 40, 50, 60, 70, 70, 66, 62]
    # yPolicy = [0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 4, 9, 12, 18, 24, 30, 23, 9, 2, 0, 0, 0]

    # # Stochastic demand -4 +16
    # objDeterm =71330.0
    # QPolicy = [40, 41, 42, 43, 40, 40, 40, 33, 20, 20, 20, 20, 20, 25, 25, 25, 30, 40, 50, 60, 62, 70, 64, 60]
    # yPolicy = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 11, 19, 26, 28, 30, 25, 13, 0, 0, 0, 0]

    # Solve
    simObj, QCurt = solve_using_stochastic_policy(totalSamples, noiseDist, QPolicy, yPolicy, params);

    println("optimal: $objExpect")

    # collect and print expected value
    μ = mean(simObj)
    sErr = sqrt(mean((simObj.-μ).^2)/(totalSamples))
    println("E[x]=$μ +- $sErr")
    h = histogram(simObj)
    display(h)
    # println("$(simObj)")
end;

# run script
run();