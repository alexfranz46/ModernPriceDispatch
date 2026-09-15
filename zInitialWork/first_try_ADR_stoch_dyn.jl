using UnPack

#=
Use dynamic recursion to solve stochastic problem
THIS WAS ORGINALLY A COPY OF first_try_ADR_dynamic.jl
=#

""" Define model

    # Fields
    - `T`: quantity of observed trading periods
    - `t`: list of trading period indices
    - `d`: demand for each trading period (MWh)
    - `B`: qunatity of generation price bands/tranches
    - `mc`: generation costs per price band/tranche (\$/MWh)
    - `bmax`: generation capacity per price band/tranche (MWh)
    - `Qmax`: maximum generation capcity (MWh)
    - `rho`: generator ramp-up constraint (MWh/h)
    - `E`: maximum storage capacity (MWh)
    - `η`: battery round-trip losses (dimensionless)
    - `r`: battery discharge constraint (MWh/h)
    - `s`: battery charge constraint (MWh/h)
    - `L`: Value of Lost Loas (\$/MWh)
    - `Q0`: initial generation
    - `y0`: initial storage
    - `noiseScenarios`: count of noise scenarios considered
    - `dNoise`: array of possible demand noise values (MWh)
    - `probNoise`: respective probabilities of noise above

    # Constructors
    - `modelParameters`: initializes parameters and checks Vector dimensions
"""
struct modelParameters
    # trading periods and demand
    T::Int
    t::Vector{Int}
    d::Vector{Int}
    
    # Generation tranches
    B::Int
    mc::Vector{Int}
    bmax::Vector{Int}

    # Generator specs
    Qmax::Int
    rho::Int

    # Battery specs
    E::Int
    η::Int
    r::Int
    s::Int

    # Value of Lost Load
    L::Int

    # Initial state conditions
    Q0::Int
    y0::Int

    # Noise
    noiseScenarios::Int
    dNoise::Vector{Int} 
    probNoise::Vector{Float64} 

    # initialize and check dimensions
    function modelParameters(T, t, d, B, mc, bmax, Qmax, rho, E, η, r, s, L, Q0, y0, noiseScenarios, dNoise, probNoise)
        T == length(t) ||
        throw(ArgumentError("t must have length T"))

        T == length(d) ||
        throw(ArgumentError("d must have length T"))

        B == length(mc) ||
        throw(ArgumentError("mc must have length B"))

        B == length(bmax) ||
        throw(ArgumentError("bmax must have length B"))

        noiseScenarios == length(dNoise) ||
        throw(ArgumentError("dNoise must have same length noiseScenarios"))
        
        noiseScenarios == length(probNoise) ||
        throw(ArgumentError("probNoise must have same length noiseScenarios"))

        sum(probNoise) == 1 ||
        throw(ArgumentError("Sum of probabilities must be 1"))

        new(T, t, d, B, mc, bmax, Qmax, rho, E, η, r, s, L, Q0, y0, noiseScenarios, dNoise, probNoise)
    end
end;


""" Function to solve system for one input pair of Q, y

    # Arguments 
    - `stage`: current stage
    - `QIn`: generation in previous stage (state variable)
    - `yIn`: storage in previous stage (state variable)
    - `BellmanVals`: bellman cost-to-go values for each [stage, Q, y] triplet
    - `params`: struct of model parameters

    # Returns
    - `optObj`: best possible bellman cost-to-go value for iteration
    - `QOptimal`: best Q state decision for iteration
    - `yOptimal`: best y state decision for iteration
"""
function solve_single_iteration(stage::Int, QIn::Int, yIn::Int, BellmanVals::Array{Float64}, params::modelParameters)
    # unpack parameters
    @unpack d, Qmax, rho, E, η, r, s, L, noiseScenarios, dNoise, probNoise = params

    # determine bounds for feasible Q(t) and y(t) given Q(t-1) and y(t-1) iterations
    # Note: the bounds are inclusive, i.e. Q⁻ <= Q <= Q⁺
    Q⁻ = 0
    Q⁺ = min(Qmax, QIn + rho)
    y⁻ = max(0, yIn - r)
    y⁺ = min(E, yIn + η*s)

    # Define 
    QFeasRange = Q⁻:Q⁺
    yFeasRange = y⁻:y⁺

    # initialize return arrays
    tempObj = fill(Inf, noiseScenarios)
    QOptimal = fill(Inf, noiseScenarios)
    yOptimal = fill(Inf, noiseScenarios)
    
    # Loop over noise scenarios
    for θ in eachindex(dNoise)
        # Find actual demand
        dReal = d[stage] + dNoise[θ]

        # initialize objective
        # additional matrices required to record current+future costs for each feasible Q(t),y(t).
        objSim = fill(Inf, Q⁺-Q⁻+1, y⁺-y⁻+1)  # THESE NEED TO BE WIPED AFTER LOOP...

        # Loop over all the possible values of Q(t), y(t) for the given QIn=Q(t-1), yIn=y(t-1).
        for (QFeasIdx, QFeas) in enumerate(QFeasRange)  # Q(t)
            for (yFeasIdx, yFeas) in enumerate(yFeasRange)  # y(t)
                # decompose y into u and ηv 
                yDiff = yFeas - yIn
                u, ηv = decompose_y(yDiff)
                
                # Find lost load required to balance supply and demand (equilibrium)
                zEql = dReal - QFeas - u + ηv/η
                
                # Extract logical meaning from equilibrium
                if zEql > dReal     # this only occurs if Q < v, meaning that the battery is charging with non-existant generation
                    continue        # hence, this is an infeasible scenario, with a bellman function of Inf
                elseif zEql < 0     # this means that supply exceeds demand, and some energy is curtailed
                    zEql = 0        # hence, this is a wasteful scenario, but not directly punished and with no Lost Load
                    # continue        # this version means spilling/curtailment of electricity is not allowed. THIS IS WHAT ANDY DOES!
                end

                # Lookup bellman function value at t, Q(t), y(t)
                bellmanVal = BellmanVals[stage + 1, QFeas + 1, yFeas + 1]  # works because indices of _Grid = _ value at indices - 1

                # Find objective
                objSim[QFeasIdx,yFeasIdx] = gen_cost(QFeas, params) + L*zEql + bellmanVal
            end
        end

        # return optimized answer for noise scenario
        tempObj[θ], optIdxPair = findmin(objSim)
        OptQIdx, OptyIdx = Tuple(optIdxPair)
        QOptimal[θ] = QFeasRange[OptQIdx]
        yOptimal[θ] = yFeasRange[OptyIdx]

    end

    # find expected objective using dot product with probablities
    optObj = sum(tempObj .* probNoise)

    return optObj, QOptimal, yOptimal
    
    # # BELOW IS A HERE-AND-NOW SOLVER

    # # Loop over all the possible values of Q(t), y(t) for the given QIn=Q(t-1), yIn=y(t-1).
    # for (QFeasIdx, QFeas) in enumerate(QFeasRange)  # Q(t)
    #     for (yFeasIdx, yFeas) in enumerate(yFeasRange)  # y(t)
    #         # decompose y into u and ηv 
    #         yDiff = yFeas - yIn
    #         u, ηv = decompose_y(yDiff)

    #         # Lookup bellman function value at t, Q(t), y(t)
    #         bellmanVal = BellmanVals[stage, QFeas + 1, yFeas + 1]  # works because indices of _Grid = _ value at indices - 1

    #         # initialize expected objective
    #         tempObj = 0
            
    #         # Loop over noise possibilities
    #         for θ in eachindex(dNoise)
    #             # Find actal demand
    #             dReal = d[stage] + dNoise[θ]
                
    #             # Find lost load required to balance supply and demand (equilibrium)
    #             zEql = dReal - QFeas - u + ηv/η
                
    #             # Extract logical meaning from equilibrium
    #             if zEql > dReal     # this only occurs if Q < v, meaning that the battery is charging with non-existant generation
    #                 tempObj = Inf   # hence, this is an infeasible scenario, with a bellman function of Inf
    #                 break
    #             elseif zEql < 0     # this means that supply exceeds demand, and some energy is curtailed
    #                 zEql = 0        # hence, this is a wasteful scenario, but not directly punished and with no Lost Load
    #             end

    #             # Update expected objective
    #             tempObj += probNoise[θ] * (gen_cost(QFeas, params) + L*zEql + bellmanVal)
    #         end    
            
    #         # Set expected objective
    #         objSim[QFeasIdx,yFeasIdx] = tempObj
    #     end
    # end

    # # vscodedisplay(objSim)

    # # return optimized answer
    # optObj, optIdxPair = findmin(objSim)
    # OptQIdx, OptyIdx = Tuple(optIdxPair)
    # QOptimal = QFeasRange[OptQIdx]
    # yOptimal = yFeasRange[OptyIdx]
    
    # return optObj, QOptimal, yOptimal
end;


""" Decomposes y into components u and ηv

    # Arguments
    - `yDiff`: y(t) - y(t-1)

    # Returns
    - `u`: discharging
    - `ηv`: charging accounting for round-trip losses
"""
function decompose_y(yDiff::Int)
    if yDiff == 0
        u = 0
        ηv = 0
    elseif yDiff >= 0
        u = 0
        ηv = yDiff
    elseif yDiff <= 0
        u = -yDiff
        ηv = 0
    else
        throw(ErrorException("Failed to decompose y into u and ηv."))
    end

    return u, ηv
end


""" Calculates generation cost for given total generation.

    # Arguments 
    - `Q`: Generation amount (MWh)
    - `params`: struct of model parameters

    # Returns
    - `totalCost`: Total cost of generation (\$)
"""
function gen_cost(Q::Int, params::modelParameters)
    # unpack parameters
    @unpack mc, bmax = params
    
    # initialize counter
    idx = 1
    totalCost = 0

    # cost of full non-final capacity tranches
    while Q > bmax[idx]
        totalCost += mc[idx]*bmax[idx]
        Q -= bmax[idx]
        idx += 1
    end 

    # cost of final (sometimes non-full) capacity tranche
    totalCost += mc[idx]*Q

    return totalCost
end;


""" Main function:
    1- Define model Parameters.
    2- Solve deterministic integer lookahead problem using dynamic recursion.
    3- Print optimal objetive and actions.
"""
function solve_stochastic_DP(params::modelParameters)

    @unpack T, t, Qmax, E, Q0, y0, noiseScenarios, dNoise = params  # remove unnecessary params

    # crate a discretized grid for the two state variables: q and y
    QStateRange = 0:Qmax
    yStateRange = 0:E

    # Create matrix to store expected bellman function values for each stage, Q, y
    BellmanVals = fill(Inf, T + 1, Qmax+1, E+1)

    # set all t=24 entries to 0 (termination condition)
    BellmanVals[end,:,:] .= 0

    # Create matrices to store the optimal Q(t), y(t) for each stage, Q(t-1), y(t-1)
    QDecision = fill(-1, T, noiseScenarios, Qmax+1, E+1)
    yDecision = fill(-1, T, noiseScenarios, Qmax+1, E+1)

    # objective
    finalObj = Inf
    
    # MAIN LOOP    
    for stage in reverse(t)
        # println("Solve stage #$stage")

        # special case for t=1, with known initial values
        # if stage == 1
        #     # call solver
        #     optObj, QOptimal, yOptimal = solve_single_iteration(stage, Q0, y0, BellmanVals, params)
            
        #     # Store optimal 0-stage/1-state decisions
        #     BellmanVals[stage, Q0 + 1, y0 + 1] = optObj
        #     QDecision[stage, :, Q0 + 1, y0 + 1] = QOptimal
        #     yDecision[stage, :, Q0 + 1, y0 + 1] = yOptimal

        #     # find objective for starting conditions
        #     finalObj = BellmanVals[1, Q0 + 1, y0 + 1]

        #     break  # End solver
        # end

        # Loop over all possible values of Q(t-1) and y(t-1)
        for QState in QStateRange  # Q(t-1)
            for yState in yStateRange  # y(t-1)            
                # Call solver
                optObj, QOptimal, yOptimal = solve_single_iteration(stage, QState, yState, BellmanVals, params)

                # Store optimal state/stage decisions
                BellmanVals[stage, QState + 1, yState + 1] = optObj
                QDecision[stage, :, QState + 1, yState + 1] = QOptimal
                yDecision[stage, :, QState + 1, yState + 1] = yOptimal
            end 
        end
    end

    # find objective for starting conditions
    finalObj = BellmanVals[1, Q0 + 1, y0 + 1]

    # TODO: Remove? Policy cannot be defined for stochastic, just retrun optimal action array.
    # # Define policy
    # QPolicy = fill(-1, T, noiseScenarios)
    # yPolicy = fill(-1, T, noiseScenarios)
    # for ti in t
    #     # Store values in temporary vars
    #     if ti == 1
    #         QTemp = Q0
    #         yTemp = y0
    #     else
    #         QTemp = QPolicy[ti - 1, :]
    #         yTemp = yPolicy[ti - 1, :]
    #     end

    #     # Save optimal decision to policy 
    #     QPolicy[ti, :] = QDecision[ti, :, QTemp + 1, yTemp + 1]
    #     yPolicy[ti, :] = yDecision[ti, :, QTemp + 1, yTemp + 1]
    #     # push!(QPolicy, QDecision[ti, :, QTemp + 1, yTemp + 1])
    #     # push!(yPolicy, yDecision[ti, :, QTemp + 1, yTemp + 1])
    # end

    # Displays
    if false
        tDisp = 18
        # vscodedisplay(QDecision[tDisp,1,:,:], "Q[t=$tDisp]")
        # vscodedisplay(yDecision[tDisp,1,:,:], "y[t=$tDisp]")
        for i in 1:noiseScenarios
            vscodedisplay(QDecision[tDisp,i,:,:], "Q[t=$tDisp,θ=$(dNoise[i])]")
            vscodedisplay(yDecision[tDisp,i,:,:], "y[t=$tDisp,θ=$(dNoise[i])]")
        end
    end

    if false
        tDisp = 1
        vscodedisplay(BellmanVals[tDisp,:,:], "V(t=$tDisp)")
        global tempV = BellmanVals[tDisp,:,:]
        # for ti in 1:T
        #     vscodedisplay(BellmanVals[ti,:,:], "V(t=$ti)")
        # end
    end

    # for debugging
    global BellmanVals
    global QDecision
    global yDecision

    return finalObj, QDecision, yDecision
end;


""" old code from above function to print results in pretty format.
"""
# function printPolicy()
#     # Print results
#     println("OBJECTIVE: $finalObj")
#     println("@ t = 0: Q=$Q0, y=$y0")
#     QRes = Q0
#     yRes = y0
#     for i in t
#         println("@ t = $i: Q=$(QDecision[i, QRes + 1, yRes + 1]), y=$(yDecision[i, QRes + 1, yRes + 1])")
#         QTemp = QRes
#         yTemp = yRes
#         QRes = QDecision[i, QTemp + 1, yTemp + 1]
#         yRes = yDecision[i, QTemp + 1, yTemp + 1]
#     end
# end;


function main()
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

        5,                          # noiseScenarios
        # [-4, -2, 0, 2, 4],          # dNoise
        # [0.2, 0.2, 0.2, 0.2, 0.2]   # probNoise
        [-4, -2, 0, 8, 16],         # dNoise
        [0.2, 0.2, 0.5, 0.05, 0.05] # probNoise
    )

    # Solve
    obj, Q, y = solve_stochastic_DP(params)
    
    # display answer
    # println("paper got: 52,377")
    println(obj)
    # println(Q)
    # println(y)
end;

# RUN SCRIPT
main();