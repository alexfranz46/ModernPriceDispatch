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
using StatsPlots
include("clean_DEP_data.jl")
include("BEST_LP_validation.jl")

mkpath("Plots")

# get BEST policies
uvDecision = deserialize(joinpath("Serials", "uvDecision.jls"))
uvDecisionIron = deserialize(joinpath("Serials", "uvDecisionIron.jls"))
PriceBounds = deserialize(joinpath("Serials", "PriceBounds.jls"))

""" Returns the band a price falls into
"""
function find_band(p::Float64, bounds::Vector{Float64})
    for (i, ub) in enumerate(bounds)
        if p <= ub
            return i - 1
        end
    end
end

# get and clean test data (50 days)
csv_files = String[]
root = joinpath(@__DIR__, "DispatchEnergyPrices", "2026_test_set")
for (path, _, files) in walkdir(root)
    for f in files
        endswith(f, ".csv") && push!(csv_files, joinpath(path, f))
    end
end
println("Importing data from $(length(csv_files)) files.")

# Combine all CSVs into one DataFrame
df_test = vcat([CSV.read(file, DataFrame) for file in csv_files]...)

# Keep only OTA2201 connection point data
df_test = filter(:PointOfConnection => p -> p == "OTA2201", df_test)

# format datetime
fmt = dateformat"yyyy-mm-ddTHH:MM:SS.sssz"
df_test.PublishDateTime = ZonedDateTime.(String.(df_test.PublishDateTime),fmt)

df_test = filter(:TradingPeriod => tp -> tp <= 48, df_test)

# Keep only date, trading period and price data
select!(df_test, [:TradingDate, :TradingPeriod, :PublishDateTime, :DollarsPerMegawattHour])

# make sure every tp has only 6 entries
df_test = clean_DEP_data(df_test)

# check all days in test set are available
valid_dates = combine(groupby(df_test, :TradingDate), nrow => :count)
valid_dates = Set(valid_dates.TradingDate[valid_dates.count .== 288])
df_valid = filter(row -> row.TradingDate in valid_dates, df_test)

days = unique(df_valid.TradingDate)
focusDays = [Date("2026-07-23"), Date("2026-08-09"), Date("2026-07-27")]
dc = 0

if false
    df = deserialize(joinpath("Serials", "df.jls"))
    df = filter(:TradingPeriod => tp -> tp <= 48, df)

    select!(df, [:TradingDate, :TradingPeriod, :PublishDateTime, :DollarsPerMegawattHour])

    df_clean = clean_DEP_data(df)
    valid_dates = combine(groupby(df_clean, :TradingDate), nrow => :count)
    valid_dates = Set(valid_dates.TradingDate[valid_dates.count .== 288])
    df_valid = filter(row -> row.TradingDate in valid_dates, df_clean)

    days = unique(df_valid.TradingDate)
    focusDays = ["2024-08-06"]
end

# # Solve noGC, BEST, ironBEST for each day
noGC = Float64[]  # 50x1 array of revenues
BEST = Float64[]
ironBEST = Float64[]
for day in days
    
    # get day data
    df_day = filter(:TradingDate => d -> d == day, df_valid)
    
    pHistory = df_day.DollarsPerMegawattHour

    if length(pHistory) == 288
        objectiveNoGC, yHistoryNoGC = noGC_LP(pHistory)
        yHistoryNoGC = [0 ; yHistoryNoGC]
        push!(noGC, objectiveNoGC)
    else
        error("wrong number of prices for $(Day.TradingDate)")
    end

    # p=plot(1/12:1/12:24, pHistory, title="$(day)")
    # display(p)

    objective = 0
    objectiveIron = 0
    i = 5  # starting in central band
    yHistory = [0]  # starting with an empty battery
    yHistoryIron = [0]
    for stage in 1:288
        yIn = yHistory[stage]
        yInIron = yHistoryIron[stage]
        
        # Update trading period
        tp = ceil(Int, stage/6)

        # fetch and record optimal decision for stage/state
        yNext = yIn - uvDecision[stage, yIn + 1, i] 
        push!(yHistory, yNext)
        yNextIron = yIn - uvDecisionIron[stage, yIn + 1, i] 
        push!(yHistoryIron, yNextIron)
        
        # fetch and record price for stage/state
        pNow = pHistory[stage]

        # calculate objective
        objective += pNow*(yIn - yNext)
        objectiveIron += pNow*(yInIron - yNextIron)

        # use random distribution to find j, which becomes i for next iteration
        j = find_band(pNow, PriceBounds[:,tp])
        i = j
    end
    
    push!(BEST, objective)
    push!(ironBEST, objectiveIron)

    if day in focusDays
        dc += 1

        # Plot results 
        ticks = 0:0.5:24
        labels = [mod(x, 1) == 0 ? string(Int(x)) : "" for x in ticks]

        # Plot storage on left axis
        p = plot(0:1/12:24, yHistoryNoGC, 
            xlabel="time of day (HH)",
            ylabel="State of charge (MWh)", 
            title="$day",
            legend=:topright,
            legendfontsize=5, 
            label="Perf. f-sight", 
            color=:silver,
            linewidth=1,
            xticks=(ticks, labels)
        )

        plot!(0:1/12:24, yHistory,
            label="BEST", 
            color=:blue,
            linewidth=1
        )

        plot!(0:1/12:24, yHistoryIron,
            label="monoBEST", 
            color=:hotpink,
            linewidth=1
        )

        # Plot price on right axis
        plot!(twinx(), 1/12:1/12:24, pHistory,
            ylabel="Price (\$/MWh)", 
            # ylims=(0, 349),
            legend=false, 
            color=:red,
            linestyle=:solid,
            linewidth=1, 
            xticks=:none,
            ytickfontcolor=:red,
            y_guidefontcolor=:red,
            yformatter = :plain
        ) 

        # display(p)
        savefig(p, "plots/ResultsDay$dc.pdf")

        if dc == 2
            ylims!(p, (0, 349)) 
            # display(p)
            savefig(p, "plots/ResultsDay$(dc)Scaled.pdf")
        end
    end
end

RevenuePerMethod=[noGC BEST ironBEST]
RevenuePerMethod ./= 1000

# plot bar chart to compare the three
labels = ["Perfect foresight" "BEST" "monoBEST"]

# sort from highest to lowest revenues
sort_indices = sortperm(vec(RevenuePerMethod[:,1]), rev=true)
sorted_days = days[sort_indices]
sorted_RevenuePerMethod = RevenuePerMethod[sort_indices, :]

# number of visual outliers
N = 5

# Highest revenue items
days_top = sorted_days[1:N]
rev_top = sorted_RevenuePerMethod[1:N, :]

# Include #5 in the lower plot as reference
days_bottom = sorted_days[N:end]
rev_bottom = sorted_RevenuePerMethod[N:end, :]

p1 = groupedbar(
    days_top,
    rev_top,
    label=labels,
    legend=false,
    ylabel="Daily revenue/loss from arbitrage (000's of \$)",
    xrotation=90,
    color = [:silver :blue :hotpink],
    bottom_margin = 10Plots.mm,
    left_margin = 10Plots.mm, 
    yguidefont = font(14),
    ytickfont = font(10),
    yformatter = :plain
)

p2 = groupedbar(
    days_bottom,
    rev_bottom,
    label=labels,
    ymirror=true,
    legend=:topright,
    ylabel="Daily revenue/loss from arbitrage (000's of \$)",
    xrotation=90,
    annotations=[
        ((0.012,0.17), text("same day", 9)),
        ((0.01,0.14), text("⟷", 30, color=:darkgreen))
    ],
    color = [:silver :blue :hotpink],
    bottom_margin = 10Plots.mm,
    right_margin = 10Plots.mm, 
    yguidefont = font(14),
    ytickfont = font(10),
    legendfont = font(14),
    yformatter = :plain
)

# Align y-axis = 0 for both subplots
ymin2 = minimum(rev_bottom)
ymax2 = maximum(rev_bottom)
r = abs(ymin2) / (ymax2 - ymin2)
ymax1 = maximum(rev_top)
ymin1 = -r/(1-r) * ymax1
ylims!(p1, (ymin1, ymax1))
ylims!(p2, (ymin2, ymax2))

pjoin = plot(
    p1, p2,
    layout=@layout([a{0.1w} b{0.9w}]),
    size=(1200,800)
)

display(pjoin)
savefig(pjoin, "plots/ResultsMain.pdf")

# for (idx, day) in enumerate(days)
#     if day in focusDays
#         println(RevenuePerMethod[idx, :])
#     end
# end