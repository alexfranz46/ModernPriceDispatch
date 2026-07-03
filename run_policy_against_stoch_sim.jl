include("first_try_ADR_dynamic.jl");
using Random
using UnPack
using Distributions
using Plots

#=
use dynamic recursion to sovle stochastic demand problem:
    d +- ε = {-4, -2, 0, 2, 4}, with equal probabilities.
=#

function solve_stochastic_using_policy(totalSamples::Int, noiseDist::Distribution, QPolicy::Vector{}, yPolicy::Vector{}, params::modelParameters)
    # unpack parameters
    @unpack T, t, d, η, L, y0 = params
    
    # initialize return arrays
    simObj = fill(0, totalSamples)
    QCurt = fill(0, totalSamples)

    for i in 1:totalSamples
        # determine demand for sample
        dRealized = d + rand(noiseDist, T)
        
        # initialize 
        yPrior = y0
        tempObj = 0
        tempQCurt = 0

        # Solve per stage
        for stage in t
            # decompose y
            yCurrent = yPolicy[stage]
            yDiff = yCurrent - yPrior
            u, ηv = decompose_y(yDiff)

            # find lost load and curtailed load
            zEql = dRealized[stage] - QPolicy[stage] - u + ηv/η

            # assume excess generation is curtailled
            # TODO: This is likely the part which is causing the error...
            if zEql <= 0  # TODO: CHANGE TO < 0
                tempQCurt -= zEql
                zEql = 0
            end

            # add obj to counter
            tempObj += gen_cost(QPolicy[stage], params) + L*zEql

            # update y(t-1)
            yPrior = yPolicy[stage]
        end

        simObj[i] = tempObj
        QCurt[i] = tempQCurt
    end

    return simObj, QCurt
end;


function main2()
    # set common random numbers
    Random.seed!(1)

    # Set noise array
    noise = [-4, -2, 0, 2, 4]
    prob = [0.2, 0.2, 0.2, 0.2, 0.2]
    # noise = [-4, -2, 0, 8, 16]
    # prob = [0.2, 0.2, 0.5, 0.05, 0.05]
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
        0       # y0
    )

    # simulation samples
    totalSamples = 1000

    # initialize results storing array
    simObj = fill(undef, totalSamples)

    # Deterministic policy
    # objDeterm, QPolicy, yPolicy = solve_deterministic_DP(params);

    # Hardcoded version
    objDeterm = 48470.0
    QPolicy = [40, 41, 42, 43, 35, 40, 40, 25, 10, 8, 6, 5, 5, 7, 10, 19, 29, 39, 49, 59, 64, 70, 64, 60]
    yPolicy = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 3, 12, 21, 30, 24, 11, 0, 0, 0, 0]

    # Solve
    simObj, QCurt = solve_stochastic_using_policy(totalSamples, noiseDist, QPolicy, yPolicy, params);
    

    println("deterministic: $objDeterm")

    # collect and print expected value
    μ = mean(simObj)
    println("E[x]=$μ")
    h = histogram(simObj)
    display(h)
end;

# run script
main2();