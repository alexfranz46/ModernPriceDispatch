# INSERT STATIC POLICY OPTIMIZERS IN THIS FILE.

""" Original BEST algorithm. 
    Note: Does not enforce monotonicity.
"""
function _BEST_legacy(case::BESTCase)
    
    # Extract relavent inputs
    T = case.T
    PBANDS = case.PBANDS
    bess = case.bess
    PriceVals = case.price.PriceVals
    TransitionMatrix = case.price.TransitionMatrix

    # Initialize returns
    sol = StaticPolicy(case)
    BellmanVals = sol.BellmanVals
    storageDecisions = sol.storageDecisions
    tradeDecisions = sol.tradeDecisions

    # Termination condition
    BellmanVals[end,:,:] .= 0

    # Iterative solver
    for stage in reverse(1:T)
        tp = ceil(Int, stage/6)

        for yIn in 0:bess.storageCapacityMWh  # y(t-1)
            # Lower bound for y based on y stage varaible
            y⁻ = max(0, yIn - bess.dischargeCapacityInterval) 
            y⁺ = min(bess.storageCapacityMWh, yIn + bess.chargeCapacityInterval)
            yRange = y⁻:y⁺

            for i in 1:PBANDS  # p(t)
                # Price stage variable
                pNow = PriceVals[i, tp]
                
                # Initialize solution
                objSim = fill(Inf, y⁺-y⁻+1)

                for (yIdx, yNext) in enumerate(yRange)  # y(t)
                    # Calculate expected bellman value (here-and-now over j)
                    expBellmanVal = sum(BellmanVals[stage + 1, yNext + 1, :] .* TransitionMatrix[i, :, tp])

                    # Record objective 
                    objSim[yIdx] = pNow*(yIn - yNext) + expBellmanVal
                end

                # Find and keep best decision
                BellmanVals[stage, yIn + 1, i], optIdxPair = findmax(objSim)
                storageDecisions[stage, yIn + 1, i] = yRange[optIdxPair[1]]
                tradeDecisions[stage, yIn + 1, i] = yIn - yRange[optIdxPair[1]]
            end
        end 
    end

    return sol
end


""" Optimized BEST algorithm which guarantees monotoic solutions. 
"""
function _BEST_monotonic(case)

    # Extract relavent inputs
    T = case.T
    PBANDS = case.PBANDS
    bess = case.bess
    PriceVals = case.price.PriceVals
    TransitionMatrix = case.price.TransitionMatrix

    # Initialize returns
    sol = StaticPolicy(case)
    BellmanVals = sol.BellmanVals
    storageDecisions = sol.storageDecisions
    tradeDecisions = sol.tradeDecisions

    # Termination condition
    BellmanVals[end,:,:] .= 0

    # Iterative solver
    for stage in reverse(1:T)
        tp = ceil(Int, stage/6)

        for yIn in 0:bess.storageCapacityMWh  # y(t-1)
            # Lower bound for y based on y stage varaible
            y⁻ = max(0, yIn - bess.dischargeCapacityInterval) 
            y⁺ = min(bess.storageCapacityMWh, yIn + bess.chargeCapacityInterval)
            yRange = y⁻:y⁺

            for i in 1:PBANDS  # p(t)
                # Price stage variable
                pNow = PriceVals[i, tp]
                
                # Initialize solution
                objSim = fill(Inf, y⁺-y⁻+1)

                for (yIdx, yNext) in enumerate(yRange)  # y(t)
                    # Calculate expected bellman value (here-and-now over j)
                    expBellmanVal = sum(BellmanVals[stage + 1, yNext + 1, :] .* TransitionMatrix[i, :, tp])

                    # Record objective 
                    objSim[yIdx] = pNow*(yIn - yNext) + expBellmanVal
                end

                # Find and keep best decision
                BellmanVals[stage, yIn + 1, i], optIdxPair = findmax(objSim)
                storageDecisions[stage, yIn + 1, i] = yRange[optIdxPair[1]]
                tradeDecisions[stage, yIn + 1, i] = yIn - yRange[optIdxPair[1]]
            end

            # Check saved solution is monotonic
            if !issorted(tradeDecisions[stage, yIn + 1, :])

                # Calclate belman values
                expBellmanVals = fill(Inf, bess.storageCapacityMWh+1, PBANDS)
                for e in yRange, i in 1:PBANDS
                    expBellmanVals[e+1,i] = sum(BellmanVals[stage+1,e+1,j] * TransitionMatrix[i,j,tp] for j in 1:PBANDS)
                end

                # Run ironed LP method
                model = Model(HiGHS.Optimizer)
                set_silent(model) # disable solver printing

                # Variables
                @variable(model, b[1:PBANDS, yRange], Bin)

                # keep track of objective value for each band
                @expression(model, objBand[i in 1:PBANDS], sum(b[i,e] * (PriceVals[i, tp] * (yIn - e) + expBellmanVals[e+1, i]) for e in yRange))

                # keep track of optimal SoC determined from binary variable
                @expression(model, y[i in 1:PBANDS], sum(e*b[i,e] for e in yRange))

                # constraints
                @constraint(model, monotonic[i in 2:PBANDS], y[i - 1] >= y[i])  # >= because decreasing SoC at the end of period means more discharge at high price and more charge at low prices
                @constraint(model, unique[i in 1:PBANDS], sum(b[i,e] for e in yRange) == 1)  # only one b=1 all others =0 for each i

                # objective
                @objective(model, Max, sum(objBand))  # TODO: Modify to have transition matrix

                # SOLVE
                optimize!(model)  # NOTE: 90% of solve time is in this optimization step...

                if !is_solved_and_feasible(model)
                    error("Solver did not find an optimal solution for stage $stage and SoC $yIn")
                end

                # update bellmans and decision
                BellmanVals[stage, yIn + 1, :] = value.(objBand)
                storageDecisions[stage, yIn + 1, :] = round.(Int, value.(y))
                tradeDecisions[stage, yIn + 1, :] = yIn .- round.(Int, value.(y))
            end
        end 
    end

    return sol
end
