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
PERPERIOD = 6  # number of trading intervals per trading period (i.e.: 6x 5-min in a 30-min tp)
GCCOUNT = 0  # number of trading intervals between Gate Closure and real time (i.e.: 12x 5-min in a 1-hr GC)
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
# set_data!(case, node, clearNodeData)
set_data!(case, node, clearNodeData, include_partial_clean_data=false)

# Define price process
set_price!(case, first_order_markov_price_process)

# Define BESS
set_BESS!(case, storageCapacityMWh, dischargeCapacityMW, chargeCapacityMW)

# Optimize
solve_policy!(case, _BEST_legacy)
solve_policy!(case, _BEST_monotonic)
# solve_policy!(case, [_BEST_legacy, _BEST_monotonic])

# Plots
historic_performance(case)
# historic_performance_1hr(case)  # TODO WIP

# bellman_visual_check_grid(case, "_BEST_legacy", 1, 2, 3)
# bellman_visual_check_grid(case, "_BEST_legacy", 1, 3, 2)
# bellman_visual_check_grid(case, "_BEST_legacy", 2, 1, 3)
# bellman_visual_check_grid(case, "_BEST_legacy", 3, 1, 2)
# bellman_visual_check_grid(case, "_BEST_legacy", 2, 3, 1)
# bellman_visual_check_grid(case, "_BEST_legacy", 3, 2, 1)

# bellman_visual_check_gif(case, "_BEST_legacy", 2, 3, 1)





# # WORKSPACE


# tallies = sum(case.price.TransitionTally, dims=2)[:,1,:]
# plot(tallies)

# least_common_band = [a[1] for a in argmin(tallies, dims=1)]
# # check how far from average/outliers...?




#investigate bellman values

# V = case.policies[:_BEST_monotonic].BellmanVals

# # y-SERIES
# colors = cgrad([:red, :green], case.bess.storageCapacityMWh÷20+1, categorical=true)
 
# for t in 1:case.T+1
#     p=plot()
#     for e in 0:20:case.bess.storageCapacityMWh
#         plot!(p, V[t,e+1,:], title="t=$t", xlabel="π", ylabel="V_t(y,π)", legend_title = "y", label=string(e), color = colors[e÷20+1], legend=:outerright, ylims=(0,1.5*10^5))
#     end
#     display(p)
# end

# for i in 1:PBANDS
#     p=plot()
#     for e in 0:20:case.bess.storageCapacityMWh
#         plot!(p, V[:,e+1,i], title="π=$i", xlabel="t", ylabel="V_t(y,π)", legend_title = "y", label=string(e), color = colors[e÷20+1], legend=:outerright, ylims=(0,1.5*10^5))
#     end
#     display(p)
# end


# # π-SERIES
# colors = cgrad([:red, :green], PBANDS, categorical=true)

# for e in 0:case.bess.storageCapacityMWh
#     p=plot()
#     for i in 1:PBANDS
#         plot!(p, V[:,e+1,i], title="y=$e", xlabel="t", ylabel="V_t(y,π)", legend_title = "π", label=string(i), color = colors[i], legend=:outerright, ylims=(0,1.5*10^5))
#     end
#     display(p)
# end

# # anim = @animate 
# for t in 1:case.T+1
#     p=plot()
#     for i in 1:PBANDS
#         plot!(p, V[t,:,i], title="t=$t", xlabel="y", ylabel="V_t(y,π)", legend_title = "π", label=string(i), color = colors[i], legend=:outerright, ylims=(0,1.5*10^5))
#     end
#     display(p)
# end
# # mp4(anim, "Enveloppe.mp4", fps=20)

# # NOTE: Devrivatives
# # for t in 1:case.T+1
# #     p=plot()
# #     for i in 1:PBANDS
# #         dVdy = diff(V[t,:,i])
# #         plot!(p, dVdy, title="t=$t", xlabel="y", ylabel="∂V_t(y,π)/∂y", legend_title = "π", label=string(i), color = colors[i], legend=:outerright, ylims=(0,700))
# #     end
# #     display(p)
# # end


# # t-SERIES
# colors = cgrad([:red, :green], case.T÷12+1, categorical=true)

# for i in 1:PBANDS
#     p=plot()
#     for t in 1:12:case.T+1
#         plot!(p, V[t,:,i], title="π=$i", xlabel="e", ylabel="V_t(y,π)", legend_title = "t", label=string(t), color = colors[t÷12+1], legend=:outerright, ylims=(0,1.5*10^5))
#     end
#     display(p)
# end

# for e in 0:case.bess.storageCapacityMWh
#     p=plot()
#     for t in 1:12:case.T+1
#         plot!(p, V[t,e+1,:], title="y=$e", xlabel="π", ylabel="V_t(y,π)", legend_title = "t", label=string(t), color = colors[t÷12+1], legend=:outerright, ylims=(0,1.5*10^5))
#     end
#     display(p)
# end




# monotonicADR = Matrix{Bool}(undef, case.T, case.bess.storageCapacityMWh+1)

# for t in 1:case.T
#     for e in 1:case.bess.storageCapacityMWh+1
#         slice = vec(case.policies[:_BEST_monotonic].BellmanVals[t, e, :])
        
#         monotonicADR[t, e] = issorted(slice, rev=true)
#     end
# end

# hm = heatmap(Int.(monotonicADR'), c = [:red, :green], 
#     legend = false, 
#     aspect_ratio = :equal, 
#     xlabel="5-min period", 
#     ylabel="State of Charge (MWh)", 
#     xlim=(1,case.T), 
#     ylim=(1,case.bess.storageCapacityMWh+1),
#     size=(528,448)
# )


# monotonicP = Matrix{Bool}(undef, case.T, case.PBANDS)

# for t in 1:case.T
#     for i in 1:case.PBANDS
#         slice = vec(case.policies[:_BEST_monotonic].BellmanVals[t, :, i])
        
#         monotonicP[t, i] = issorted(slice)
#     end
# end

# hm = heatmap(Int.(monotonicP'), c = [:red, :green], 
#     legend = false, 
#     aspect_ratio = :equal, 
#     xlabel="5-min period", 
#     ylabel="State of Charge (MWh)", 
#     xlim=(1,case.T), 
#     ylim=(1,case.PBANDS),
#     size=(528,448)
# )





# p3d=surface(1:case.bess.storageCapacityMWh+1,1:case.T+1,V[:,:,1], camera = (210, 30))
# display(p3d)