""" Returns the band a price falls into
"""
function find_band(p::Float64, bounds::Vector{Float64})
    for (i, ub) in enumerate(bounds)
        if p <= ub
            return i - 1
        end
    end
end


# """ Checks if a certain variable exists in the serials cache
# """
# function serial_path(variable)
#     return joinpath("Serials", variable*".jls")
# end




####

""" Stores data for node
"""
struct NodeData
    node::String
    trainingData::DataFrame
    testCaseData::DataFrame
end


""" Stores BESS parameters
"""
struct BESS
    storageCapacityMWh::Int
    dischargeCapacityInterval::Int
    chargeCapacityInterval::Int
end


""" Stores price behaviour
"""
struct PriceProcess
    PriceVals::Matrix{Float64}
    PriceBounds::Matrix{Float64}
    TransitionMatrix::Array{Float64,3}
end


""" Stores static policy for optimal BESS operation
"""
mutable struct StaticPolicy
    BellmanVals::Array{Float64,3}
    storageDecisions::Array{Int,3}
    tradeDecisions::Array{Int,3}
end


""" Case object...
"""
mutable struct BESTCase
    
    # Principal parameters
    TP::Int
    PERPERIOD::Int
    T::Int
    PBANDS::Int

    # Sub-structures
    data::Union{Nothing,NodeData}
    bess::Union{Nothing,BESS}
    price::Union{Nothing,PriceProcess}
    policy::Union{Nothing,StaticPolicy}

    # Initialize with Constructor
    function BESTCase(TP::Int, PERPERIOD::Int, PBANDS::Int)
        
        new(
            # Principal parameters
            TP,
            PERPERIOD,
            TP*PERPERIOD,
            PBANDS,

            # Sub-structures
            nothing,
            nothing,
            nothing,
            nothing,
        )
    end
end


""" Outer constructor for StaticPolicy
"""
function StaticPolicy(case::BESTCase)
        
    T = case.T
    PBANDS = case.PBANDS
    storage = case.bess.storageCapacityMWh
    
    StaticPolicy(
        Array{Float64, 3}(undef, T+1, storage + 1, PBANDS),  # TODO: implement switch-case for termination condition?
        Array{Int, 3}(undef, T, storage + 1, PBANDS),
        Array{Int, 3}(undef, T, storage + 1, PBANDS),
    )
end


""" Updates case object with training and test case data
"""
function set_data!(case::BESTCase, node::String, clear_cache::Bool; include_partial_clean_data::Bool=true)

    # Clear cache 
    serialDir = joinpath(@__DIR__, "Serials", "NodeData")
    filepath = joinpath(serialDir, node*".jls")
    if clear_cache

        # Iterate over every file in directory
        for serial in readdir(serialDir, join=true)
            
            # Identify file for current node
            if serial == filepath
                
                # Remove file
                rm(serial, recursive=true, force=true)
                println("Cleared $node.jls from cache.")
            end
        end
    end

    # Load full dataset
    if isfile(filepath)

        # Data is already in cache
        trainingData = deserialize(filepath)
        println("Loaded $node.jls from cache.")
    else

        # Load data from EA database
        trainingData = load_and_process_data(node)

        # Save data to cache
        serialize(filepath, trainingData)
        println("Saved $node.jls to cache.")
    end

    # Identify clean data for test cases 
    testCaseData = get_clean_test_case_data(case.TP, case.PERPERIOD, trainingData, include_partial_clean_data)

    # Initialize
    case.data = NodeData(
        node,
        trainingData,
        testCaseData
    )

    return nothing
end


function set_BESS!(case::BESTCase, storageCapacityMWh::Int, dischargeCapacityMW::Int, chargeCapacityMW::Int)
        
    # Validate inputs
    storageCapacityMWh > 0  ||  throw(ArgumentError("storageCapacityMWh must be positive."))
    dischargeCapacityMW > 0 ||  throw(ArgumentError("dischargeCapacityMW must be positive."))
    chargeCapacityMW > 0    ||  throw(ArgumentError("chargeCapacityMW must be positive."))

    # Convert MW power capacities in MWh/time-step capacities
    exactDischargeCapacityInterval = dischargeCapacityMW * 24 / case.T
    exactChargeCapacityInterval = chargeCapacityMW * 24 / case.T
    
    # Integer approximation
    dischargeCapacityInterval = round(Int, exactDischargeCapacityInterval, RoundDown)
    chargeCapacityInterval = round(Int, exactChargeCapacityInterval, RoundDown)
    
    # Warn if approximation is significant
    if abs(exactDischargeCapacityInterval - dischargeCapacityInterval) > 1e-6 || abs(exactChargeCapacityInterval - chargeCapacityInterval) > 1e-6
        @warn """
        This solver approximates charge and discharge limits using integer energy transfers per time step.
        
        Requested:
        discharge = $exactDischargeCapacityInterval MWh/interval
        charge = $exactChargeCapacityInterval MWh/interval
        
        Approximated:
        discharge = $dischargeCapacityInterval MWh/interval
        charge = $chargeCapacityInterval MWh/interval
        """
    end

    # Initialize
    case.bess = BESS(
        storageCapacityMWh,
        dischargeCapacityInterval,
        chargeCapacityInterval,
    )
    
    return nothing
end


function set_price!(case::BESTCase, price_model)
    

    #TODO?: Load from serials?

    # MarkovDataFiles = [
    #     joinpath(serialDir, "PriceVals.jls"),
    #     joinpath(serialDir, "PriceBounds.jls"),
    #     joinpath(serialDir, "TransitionMatrix.jls")
    # ]
    # if all(f -> isfile(f), MarkovDataFiles)

    #     # Data is already in cache
    #     PriceVals = deserialize(MarkovDataFiles[1])
    #     println("Loaded $(MarkovDataFiles[1]) from cache.")
    
    #     PriceBounds = deserialize(MarkovDataFiles[2])
    #     println("Loaded $(MarkovDataFiles[2]) from cache.")
    
    #     TransitionMatrix = deserialize(MarkovDataFiles[3])
    #     println("Loaded $(MarkovDataFiles[3]) from cache.")
    # else
        
    #     # Load data from EA database
    #     PriceVals, PriceBounds, TransitionMatrix = first_order_markov_price_process(rawData, TP, PBANDS)
    
    #     # Save data to cache
    #     serialize(MarkovDataFiles[1], PriceVals)
    #     println("Saved $(MarkovDataFiles[1]) to cache.")
    
    #     serialize(MarkovDataFiles[2], PriceBounds)
    #     println("Saved $(MarkovDataFiles[2]) to cache.")
        
    #     serialize(MarkovDataFiles[3], TransitionMatrix)
    #     println("Saved $(MarkovDataFiles[3]) to cache.")
    # end



    # Initialize using inputted method
    case.price = price_model(case)

    return nothing
end

function solve_policy!(case::BESTCase, policy_optimizer)
    
    case.policy = policy_optimizer(case)

    return nothing
end
