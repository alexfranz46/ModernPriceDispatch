# INSERT STATIC POLICY OPTIMIZERS IN THIS FILE.

""" monotonic method without lp WIP!!
"""
function _BEST_temp(case::BESTCase)
    
    # Extract relavent inputs
    T = case.T
    PERPERIOD = case.PERPERIOD
    PBANDS = case.PBANDS
    bess = case.bess
    PriceVals = case.price.PriceVals
    TransitionMatrix = case.price.TransitionMatrix

    # Initialize returns
    sol = StaticPolicy(case, false)
    BellmanVals = sol.BellmanVals
    storageDecisions = sol.storageDecisions
    tradeDecisions = sol.tradeDecisions

    # Termination condition
    BellmanVals[end,:,:] .= 0

    # Iterative solver
    for stage in reverse(1:T)
        println("Solving stage $stage")

        tp = ceil(Int, stage/PERPERIOD)

        V = BellmanVals[stage+1,:,:]
        P = PriceVals[:,tp]

        for yIn in 0:bess.storageCapacityMWh  # y(t-1)

            # feasible bounds of trades
            u⁻ = max(yIn - bess.storageCapacityMWh, -bess.chargeCapacityInterval) 
            u⁺ = min(yIn, bess.dischargeCapacityInterval)

            for i in 1:PBANDS  # p(t-1)
                
                # Initialize variables for enumerative optimizer 
                uBest = nothing
                vBest = -Inf

                # Probability of each price band for next stage
                M = TransitionMatrix[i, :, tp]
                
                for uCandIdx in case.stacks.feasibleIdxs[(u⁻, u⁺)]

                    uCandidate = case.stacks[uCandIdx]

                    # define indices for relavent bellman values of next stage
                    # indices = CartesianIndex.(yIn .- uIn .+ 1, uCandIdx, 1:PBANDS)
                    
                    # Calculate bellman value for this stage and candidate trade stack
                    # obj = sum(M .* (P.*uCandidate .+ V[indices]) )
                    obj = 0
                    for j in 1:PBANDS
                        obj += M[j] * (P[j]* uCandidate[j] + V[yIn - uCandidate[j] + 1, j])
                    end

                    # save details if best decision yet
                    if obj > vBest
                        vBest = obj
                        uBest = uCandidate
                    end
                end

                # while true

                #     # define indices for relavent bellman values of next stage
                #     # indices = CartesianIndex.(yIn .- uCandidate .+ 1, 1:PBANDS)
                #     #                                ^ SHOULD BE uIn CHECK FORMULA!!
                    
                #     # Calculate bellman value for this stage and candidate trade stack
                #     # obj = sum(M .* (PriceVals[:,tp].*uCandidate .+ V[indices]) )

                #     obj = 0
                #     for j in 1:PBANDS
                #         obj += M[j] * (P[j]* uCandidate[j] + V[yIn + 1, uCandIdx, j])
                #     end  # CHECK THIS IS CORRECT!!!!!!!

                #     # save details if best decision yet
                #     if obj > vBest
                #        uBest = copy(uCandidate)
                #        vBest = copy(obj) 
                #     end

                #     # iterate candidate trade stack
                #     next_valid_stack!(uCandidate, u⁺) || break
                # end

                # record optimized decision
                BellmanVals[stage, yIn + 1, i] = vBest
                tradeDecisions[stage, yIn + 1, i, :] .= uBest 
                storageDecisions[stage, yIn + 1, i, :] .= yIn .- uBest
            end
        end 
    end

    return sol
end


""" gc method WIP!!
"""
function _BEST_gc(case::BESTCase)
    
    # Extract relavent inputs
    T = case.T
    PERPERIOD = case.PERPERIOD
    PBANDS = case.PBANDS
    bess = case.bess
    PriceVals = case.price.PriceVals
    TransitionMatrix = case.price.TransitionMatrix

    # Initialize returns
    sol = StaticPolicy(case, true)
    BellmanVals = sol.BellmanVals
    storageDecisions = sol.storageDecisions
    tradeDecisions = sol.tradeDecisions

    # Termination condition
    BellmanVals[end,:,:,:] .= 0

    # Iterative solver
    for stage in reverse(1:T)
        println("Solving stage $stage")

        tp = ceil(Int, stage/PERPERIOD)

        V = BellmanVals[stage+1,:,:,:]
        P = PriceVals[:,tp]

        for yIn in 0:bess.storageCapacityMWh  # y(t-h_GC)

            # feasible bounds of unrealised trades
            uIn⁻ = max(yIn - bess.storageCapacityMWh, -bess.chargeCapacityInterval) 
            uIn⁺ = min(yIn, bess.dischargeCapacityInterval)

            for i in 1:PBANDS  # p(t-h_GC)

                # Probability of each price band for next stage
                M = TransitionMatrix[i, :, tp]

                for uInIdx in case.stacks.feasibleIdxs[(uIn⁻, uIn⁺)]# u(t-h_GC)
                    
                    # @show(yIn, uIn)
                                    
                    uIn = case.stacks[uInIdx]

                    # feasible bounds of trades
                    uCand⁻ = max(yIn - uIn[1] - bess.storageCapacityMWh, -bess.chargeCapacityInterval) 
                    uCand⁺ = min(yIn - uIn[end], bess.dischargeCapacityInterval)

                    # @show(uCand⁻, uCand⁺)
                        
                    # Initialize variables for enumerative optimizer
                    vBest = -Inf
                    uBest = case.stacks.idleStack  # (0, 0, ..., 0)

                    for uCandIdx in case.stacks.feasibleIdxs[(uCand⁻, uCand⁺)]

                        uCandidate = case.stacks[uCandIdx]

                        # define indices for relavent bellman values of next stage
                        # indices = CartesianIndex.(yIn .- uIn .+ 1, uCandIdx, 1:PBANDS)
                        
                        # Calculate bellman value for this stage and candidate trade stack
                        # obj = sum(M .* (P.*uCandidate .+ V[indices]) )
                        obj = 0
                        for j in 1:PBANDS
                            obj += M[j] * (P[j]* uIn[j] + V[yIn - uIn[j] + 1, uCandIdx, j])
                            #              ^ HERE TO  ^ doesn't change in for uCandIdx loop...
                        end

                        # save details if best decision yet
                        if obj > vBest
                            vBest = obj
                            uBest = uCandidate
                        end
                    end

                    # for (uCandIdx, uCandidate) in enumerate(case.stacks)  # DECISION

                    #     # calculate bellman value if feasible
                    #     if (uCand⁻ > uCandidate[1] && uCandidate[end] > uCand⁺)  # Decision

                    #         # define indices for relavent bellman values of next stage
                    #         indices = CartesianIndex.(yIn .- uCandidate .+ 1, uCandIdx, 1:PBANDS)
                            
                    #         # Calculate bellman value for this stage and candidate trade stack
                    #         obj = sum(M .* (PriceVals[:,tp].*uCandidate .+ V[indices]) )

                    #         # save details if best decision yet
                    #         if obj > vBest
                    #             vBest = obj
                    #             uBest = uCandidate
                    #         end
                    #     end
                        
                    # record optimized decision
                    BellmanVals[stage, yIn + 1, uInIdx, i] = vBest
                    tradeDecisions[stage, yIn + 1, uInIdx, i, :] .= uBest 
                    storageDecisions[stage, yIn + 1, uInIdx, i, :] .= yIn .- uBest
                    # end
                end

                # for (uCandidx, uCandidate) in enumerate(case.stacks)
                #     is_feasible(uCandidate, yIn, uIn) || continue
                #     ...
                # end

                # while true  # u(t-h_GC)

                #     # feasible bounds of trades
                #     uCand⁻ = max(yIn - uIn[1] - bess.storageCapacityMWh, -bess.dischargeCapacityInterval) 
                #     uCand⁺ = min(yIn - uIn[end], bess.chargeCapacityInterval)
                #     uCandidate = fill(uCand⁻, PBANDS)
                #     uCandIdx = 1

                #     # Initialize variables for enumerative optimizer 
                #     uBest = nothing
                #     vBest = -Inf

                #     while true  # Decision

                #         # define indices for relavent bellman values of next stage
                #         indices = CartesianIndex.(yIn .- uIn .+ 1, uCandIdx, 1:PBANDS)
                        
                #         # Calculate bellman value for this stage and candidate trade stack
                #         obj = sum(M .* (PriceVals[:,tp].*uCandidate .+ V[indices]) )

                #         # save details if best decision yet
                #         if obj > vBest
                #             uBest = copy(uCandidate)
                #             vBest = copy(obj) 
                #         end

                #         # iterate candidate trade stack
                #         next_valid_stack!(uCandidate, uCand⁺) || break
                #         uCandIdx += 1
                #     end

                    
                #     # record optimized decision
                #     BellmanVals[stage, yIn + 1, uInIdx, i] = vBest
                #     tradeDecisions[stage, yIn + 1, uInIdx, i, :] .= uBest 
                #     storageDecisions[stage, yIn + 1, uInIdx, i, :] .= yIn .- uBest
                    
                #     # iterate unrealised trade stack
                #     next_valid_stack!(uIn, uIn⁺) || break
                #     uInIdx += 1
                # end
            end
        end 
    end

    return sol
end



# OLDER MODELS... (Static policy N=3)

""" Original BEST algorithm. 
    Note: Does not enforce monotonicity.
"""
function _BEST_legacy(case::BESTCase)
    
    # Extract relavent inputs
    T = case.T
    PERPERIOD = case.PERPERIOD
    PBANDS = case.PBANDS
    bess = case.bess
    PriceVals = case.price.PriceVals
    TransitionMatrix = case.price.TransitionMatrix

    # Initialize returns
    sol = StaticPolicy(case, 3)
    BellmanVals = sol.BellmanVals
    storageDecisions = sol.storageDecisions
    tradeDecisions = sol.tradeDecisions

    # Termination condition
    BellmanVals[end,:,:] .= 0

    # Iterative solver
    for stage in reverse(1:T)
        tp = ceil(Int, stage/PERPERIOD)

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
                tradeDecisions[stage, yIn + 1, i] = yRange[optIdxPair[1]] - yIn
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
    sol = StaticPolicy(case, 3)
    BellmanVals = sol.BellmanVals
    storageDecisions = sol.storageDecisions
    tradeDecisions = sol.tradeDecisions

    # Termination condition
    BellmanVals[end,:,:] .= 0

    # Iterative solver
    for stage in reverse(1:T)
        tp = ceil(Int, stage/PERPERIOD)

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
                tradeDecisions[stage, yIn + 1, i] = yRange[optIdxPair[1]] - yIn
            end

            # Check saved solution is monotonic
            if !issorted(tradeDecisions[stage, yIn + 1, :], rev=true)

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
                tradeDecisions[stage, yIn + 1, :] = round.(Int, value.(y)) .- yIn
            end
        end 
    end

    return sol
end
