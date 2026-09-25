function load_and_process_data(node::String)
    
    # collect all CSV files
    csv_files = String[]
    root = joinpath(@__DIR__, "DispatchEnergyPrices") 
    for (path, _, files) in walkdir(root)        
        for f in files
            endswith(f, ".csv") && push!(csv_files, joinpath(path, f))
        end
    end
    # TODO: if no csv found Throw error asking to run dowload_EA_datasets.py!

    println("Importing data from $(length(csv_files)) files.")

    # Combine all CSVs into one DataFrame
    rawData = vcat([CSV.read(file, DataFrame) for file in csv_files]...)

    # Keep only OTA2201 connection point data
    rawData = filter(:PointOfConnection => p -> p == node, rawData) # TODO:REMOVE

    # Format datetime
    format = dateformat"yyyy-mm-ddTHH:MM:SS.sssz"
    rawData.PublishDateTime = ZonedDateTime.(String.(rawData.PublishDateTime), format)

    # account for daylight saving days with 46 or 50 trading periods
    daylightSavings = Set(first(g.TradingDate) for g in groupby(rawData, :TradingDate) if maximum(g.TradingPeriod) != 48)
    rawData = filter(row -> !(row.TradingDate in daylightSavings), rawData)  # remove daylight savings data
    
    # TODO: alternative method which changes dst to match nzt...
    # shortDays = Set(first(g.TradingDate) for g in groupby(intervalsByTP, :TradingDate) if nrow(g) == 46)
    # longDays = Set(first(g.TradingDate) for g in groupby(intervalsByTP, :TradingDate) if nrow(g) == 50)
    # transform!(rawData, [:TradingDate, :TradingPeriod] => 
    #     ByRow((d, tp) -> begin 
    #         if d in shortDays && tp >= 5 tp + 2 
    #         elseif d in longDays && tp >= 7 tp - 2 
    #         else tp 
    #         end 
    #     end) => :TradingPeriod 
    # )

    # Keep only date, trading period, timestamp and price data
    select!(rawData, [:TradingDate, :TradingPeriod, :PublishDateTime, :DollarsPerMegawattHour])

    return rawData
end


function get_clean_test_case_data(TP::Int, PERPERIOD::Int, rawData::DataFrame, include_partial_clean_data::Bool)
    
    # Count number of observations per TP and Day
    intervalsByDay = combine(groupby(rawData, :TradingDate), nrow => :Count)
    intervalsByTP = combine(groupby(rawData, [:TradingDate, :TradingPeriod]), nrow => :Count)

    # Find days with exactly T data points
    cleanData = Set(first(g.TradingDate) for g in groupby(intervalsByTP, :TradingDate) if all(g.Count .== PERPERIOD))
    partialCleanData = Set(intervalsByDay.TradingDate[intervalsByDay.Count .== TP*PERPERIOD])
    
    if !isempty(cleanData) && !include_partial_clean_data
        return filter(row -> row.TradingDate in cleanData, rawData)
    elseif !isempty(partialCleanData)
        @warn """
        Some of the returned data is only "partially complete":

        Complete data: days with $(TP*PERPERIOD) price points, and exactly $(TP) periods with $(PERPERIOD) price points each.
        Partially complete data: days with $(TP*PERPERIOD) price points, but not necessarily $(TP) periods with $(PERPERIOD) price points each.

        To exclude this partially complete data, set kwarg 'include_partial_clean_data' to false.
        """
        return filter(row -> row.TradingDate in partialCleanData, rawData)
    else
        error("The provided data has no days with complete data. Consider running clean_DEP_data() for approximated complete data.")
        # TODO: include clean_DEP_data.jl
    end
end

