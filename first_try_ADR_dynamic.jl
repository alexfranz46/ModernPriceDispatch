using UnPack

#=
Use the JuMP tutorials to solve deterministic model using dynamic recursion
=#

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

    # initialize and check dimensions
    function modelParameters(T, t, d, B, mc, bmax, Qmax, rho, E, η, r, s, L, Q0, y0)
        
        println("T=$T, len(t)=$(length(t))")

        T == length(t) ||
        throw(ArgumentError("t must have length T"))

        T == length(d) ||
        throw(ArgumentError("d must have length T"))

        B == length(mc) ||
        throw(ArgumentError("mc must have length B"))

        B == length(bmax) ||
        throw(ArgumentError("bmax must have length B"))

        new(T, t, d, B, mc, bmax, Qmax, rho, E, η, r, s, L, Q0, y0)
    end

end;

""" function to solve system for one input pair of Q, y
"""
function solve_single_iteration(stage::Int, QIn::Int, yIn::Int, BellmanVals::Array{Float64}, params::modelParameters)#Qmax::Int, E::Int, rho::Int, r::Int, s::Int, η::Int, d::Vector{Int}, mc::Vector{Int}, bmax::Vector{Int}, L::Int)
    # unpack parameters
    @unpack T, t, B, mc, bmax, d, Qmax, E, η, r, s, rho, L, Q0, y0 = params  # remove unnecessary params
    
    # determine bounds for feasible Q(t) and y(t) given Q(t-1) and y(t-1) iterations
    # Note: the bounds are inclusive, i.e. Q⁻ <= Q <= Q⁺
    Q⁻ = 0
    Q⁺ = min(Qmax, QIn + rho)
    y⁻ = max(0, yIn - r)
    y⁺ = min(E, yIn + η*s)

    # Define 
    QFeasRange = Q⁻:Q⁺
    yFeasRange = y⁻:y⁺

    # additional matrices required to record current+future costs for each feasible Q(t),y(t).
    objSim = fill(Inf, Q⁺-Q⁻+1, y⁺-y⁻+1)  # THESE NEED TO BE WIPED AFTER LOOP...
    
    # Loop over all the possible values of Q(t), y(t) for the given QIn=Q(t-1), yIn=y(t-1).
    for (QFeasIdx, QFeas) in enumerate(QFeasRange)  # Q(t)
        for (yFeasIdx, yFeas) in enumerate(yFeasRange)  # y(t)
            # decompose y into u and ηv 
            yDiff = yFeas - yIn
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

            # Domain check u and ηv
            if u < 0
                throw(DomainError(u, "u must be positive."))
            elseif ηv <0
                throw(DomainError(ηv, "ηv must be positive."))
            end
            
            # Find range of lost load required to balance supply and demand (equilibrium)
            zEql = d[stage] - QFeas - u + ηv/η
            
            # if equilibrium is infeasible, reject solution (obj=Inf)
            if !(0 <= zEql <= d[stage])
                continue
            end

            # Lookup bellman function value at t, Q(t), y(t)
            bellmanVal = BellmanVals[stage, QFeas + 1, yFeas + 1]  # works because indices of _Grid = _ value at indices - 1

            # Find objective
            objSim[QFeasIdx,yFeasIdx] = gen_cost(QFeas, mc, bmax) + L*zEql + bellmanVal
            
            # temp = gen_cost(QFeas, mc, bmax) + L*zEql + bellmanVal
            # if temp < 0
            #     println("@ (Q,yFeas)=($QFeas,$yFeas), obj=$temp, cᵗ(Q)=$(gen_cost(QFeas, mc, bmax)), Lz=$(L*zEql), F=$bellmanVal")
            # end

            # println("Q=$QFeas, u=$u, v=$(ηv/η), y=$yFeas, z=$zMin")
        end
    end

    # return optimized answer
    optObj, optIdxPair = findmin(objSim)
    OptQIdx, OptyIdx = Tuple(optIdxPair)
    QOptimal = QFeasRange[OptQIdx]
    yOptimal = yFeasRange[OptyIdx]
    
    return optObj, QOptimal, yOptimal
end;

""" calculates generation cost for given total generation
"""
function gen_cost(Q::Int, mc::Vector, bmax::Vector)
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

""" Main function
    Runs code to solve deterministic integer lookahead problem using dynamic methods.
"""
function main()
    # Define model parameters
    params = modelParameters(
        24,                                             # T
        1:24,                                  # t
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

    @unpack T, t, d, B, mc, bmax, Qmax, rho, E, η, r, s, L, Q0, y0 = params  # remove unnecessary params

    # crate a discretized grid for the two state variables: q and y
    QStateRange = 0:Qmax
    yStateRange = 0:E

    # Create matrix of Inf to store bellman function values for each stage, Q, y
    BellmanVals = fill(Inf, T, Qmax+1, E+1)

    # set all t=24 entries to 0
    BellmanVals[end,:,:] .= 0

    # Create matrices to store the optimal Q(t), y(t) for each stage, Q(t-1), y(t-1)
    QDecision = fill(-1, T, Qmax+1, E+1)
    yDecision = fill(-1, T, Qmax+1, E+1)

    # objective
    finalObj = Inf
    
    # mininize single-period equation to solve for each state/stage pair, starting at t=24

    # MAIN LOOP    
    for stage in reverse(t)
    println("Start t=$stage")

        # special case for t=1, with known initial values
        if stage == 1
            optObj, QOptimal, yOptimal = solve_single_iteration(stage, Q0, y0, BellmanVals, params)
            
            # TODO: post solve updates to bellman and _decisions
            finalObj = optObj
            QDecision[stage, Q0 + 1, y0 + 1] = QOptimal
            yDecision[stage, Q0 + 1, y0 + 1] = yOptimal

            break
        end

        # Loop over all possible values of Q(t-1) and y(t-1)
        for QState in QStateRange  # Q(t-1)
            for yState in yStateRange  # y(t-1)            
                optObj, QOptimal, yOptimal = solve_single_iteration(stage, QState, yState, BellmanVals, params)

                # Store optimal state/stage decisions
                BellmanVals[stage - 1, QState + 1, yState + 1] = optObj
                QDecision[stage, QState + 1, yState + 1] = QOptimal
                yDecision[stage, QState + 1, yState + 1] = yOptimal

                # println("stage=$(stage - 1); optCost=$optObj; Q(t)=$(QFeasRange[OptQIdx]), y(t)=$(yFeasRange[OptyIdx])")
            end 
        end

        # TODO: find best something? Don't think there is anything
        # ...

    end

    # final stuff

    # println(findall(==(-1),yDecision))
    # println(count(==(-1),yDecision))


    println("OBJECTIVE: $finalObj")
    println("@ t = 0: Q=$Q0, y=$y0")
    QRes = Q0
    yRes = y0
    for i in t
        println("@ t = $i: Q=$(QDecision[i, QRes + 1, yRes + 1]), y=$(yDecision[i, QRes + 1, yRes + 1])")
        QTemp = QRes
        yTemp = yRes
        QRes = QDecision[i, QTemp + 1, yTemp + 1]
        yRes = yDecision[i, QTemp + 1, yTemp + 1]
    end
    # println(QDecision)
    # println(yDecision)


end

# RUN SCRIPT
main()