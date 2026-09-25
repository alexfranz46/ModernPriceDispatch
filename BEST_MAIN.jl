# # # # # # # # # # #
 # # #PACKAGES # # # 
# # # # # # # # # # #

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


# # # # # # # # # # #
 # #SUBDIRECTORIES # 
# # # # # # # # # # #

mkpath("Plots")
mkpath("Serials/NodeData")
mkpath("Serials/PriceProcess")
#TODO: CHECK DispatchEnergyPrices DIRECTORY EXISTS IN CORRECT FORMAT!


# # # # # # # # # # #
 # # # MODULES # # # 
# # # # # # # # # # #

include("BEST_helper_functions.jl")
include("BEST_data_managers.jl")
include("BEST_price_models.jl")
include("BEST_static_policy_optimizers.jl")
include("BEST_performance_tests.jl")


# # # # # # # # # # #
 # # # TOGGLES # # # 
# # # # # # # # # # #

clearNodeData = false


# # # # # # # # # # #
 # # #PARAMETERS # # 
# # # # # # # # # # #

TP = 48  # number of trading periods in a day
PERPERIOD = 6  # number of sub-trading period time steps (i.e.: 6x 5-min RTD in a 30-min tp)
PBANDS = 10  # number of price bands
node = "OTA2201"  # Prices from this GIP/GXP
storageCapacityMWh = 240  # BESS Energy Capacity (MWh)
dischargeCapacityMW = 120  # BESS discharge power capacity (MW)
chargeCapacityMW = 120  # BESS charge power capacity (MW)


# # # # # # # # # # #
 # # # PROBLEM # # # 
# # # # # # # # # # #

# Defire problem case
case = BESTCase(TP, PERPERIOD, PBANDS)

# Handle data
set_data!(case, node, clearNodeData)#, include_partial_clean_data=false)

# Define price process
set_price!(case, first_order_markov_price_process)

# Define BESS
set_BESS!(case, storageCapacityMWh, dischargeCapacityMW, chargeCapacityMW)

# Optimize
solve_policy!(case, _BEST_monotonic)

# Plots
performance_with_historic_data(case)
