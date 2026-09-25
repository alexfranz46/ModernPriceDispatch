function first_order_markov_price_process(case::BESTCase)
    
    # Extract relavent inputs
    TP = case.TP
    PBANDS = case.PBANDS
    trainingData = case.data.trainingData

    # Initialize retruns
    PriceVals = Matrix{Float64}(undef, PBANDS, TP)
    PriceBounds = Matrix{Float64}(undef, PBANDS+1, TP)
    TransitionMatrix = zeros(Float64, PBANDS, PBANDS, TP)

    # Unbounded lowest and highest bounds
    PriceBounds[1, :] .= -Inf  # TODO: 0?
    PriceBounds[end, :] .= Inf

    # Calculate cutoff values for qunatiles
    cutoffs = range(0, 1, length = PBANDS+1)[2:end-1]

    # Define bands
    for tp in 1:TP

        # Prices in trading period
        prices = Vector(trainingData[trainingData.TradingPeriod .== tp, :DollarsPerMegawattHour])

        # Price boundaries
        PriceBounds[2:end-1, tp] = quantile(prices, cutoffs)

        # Check for bands where lb == ub, happens if not enough data
        badBands = findall(i -> PriceBounds[i, tp] == PriceBounds[i+1, tp], 1:PBANDS)
        if !isempty(badBands)
            @warn "Upper and lower bounds of price bands $badBands in trading period $tp are equal. Consider decreasing the number of bands."
        end

        # Price values
        for i in 1:PBANDS
            lb = PriceBounds[i, tp]
            ub = PriceBounds[i+1, tp]
            PriceVals[i, tp] = mean(prices[(prices .> lb) .& (prices .<= ub)])
        end
    end

    # Define transition matrix

    # Initialize
    rowI = nothing
    firstRow = true

    # Iterate over every row ordered by time
    for rowJ in eachrow(trainingData)

        # Ignore first row
        if firstRow
            rowI = rowJ
            firstRow = false
            continue  # skip rest of iteration
        end

        # identify the indices of ttransition
        tpI = rowI.TradingPeriod
        tpJ = rowJ.TradingPeriod
        i = find_band(rowI.DollarsPerMegawattHour, PriceBounds[:, tpI])
        j = find_band(rowJ.DollarsPerMegawattHour, PriceBounds[:, tpJ])

        # Update tally
        TransitionMatrix[i, j, tpI] += 1

        # Update for next transition
        rowI = rowJ
    end

    # Convert tallies to probabilities
    for tp in 1:TP
        for i in 1:PBANDS
            
            # calculate the total sum of rows
            totalI = sum(TransitionMatrix[i, :, tp])
            # TODO: account for totalI==0?

            # transform count into probabilties
            TransitionMatrix[i, :, tp] ./= totalI
        end
    end

    return PriceProcess(PriceVals, PriceBounds, TransitionMatrix)
end


function new_price_process(case::BESTCase)
    # Extract relavent inputs
    TP = case.TP
    #...

    # Initialize retruns
    PriceVals = nothing
    PriceBounds = nothing
    TransitionMatrix = nothing

    # INSERT PRICE PROCESS METHOD HERE

    return PriceProcess(PriceVals, PriceBounds, TransitionMatrix)
end