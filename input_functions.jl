using CSV
using DataFrames
using StatsBase
csv_path = "../../../Data/ODE_Paper_2"

"""
    function add_future_frequency!(data) 

Modifies temperature data file to incorporate season, year, and RCP-specific future heatwave frequencies as projected by CMIP5 (delta = days).
"""
function add_future_frequency!(data)
    annual_count = collect(1:1:365)
    data.Time .= annual_count

    x = data.futureDeltaMAX
    data[!, "HeatWaveAdditional"] = zeros(length(x))
    grouped_data = groupby(data, :futureDeltaFrequency)
    exceptions = collect(1:365)[data.frequencyBanned .== 1]

    for group in grouped_data
        x = group.futureDeltaMAX
        max_heat_wave_periods = max(1, ceil(Int, group.futureDeltaFrequency[1] / 3.0))
        t = Vector{Float64}(undef, length(x) - 3)
        for i in 1:(length(x) - 3)
            t[i] = mean(x[i:(i + 2)])
            if !isfinite(t[i])
                @error "invalid inputs $(x[i:i+2]) at index $i"
            end
        end
        ixs = Vector{Int}(indexin(t, sort(t, rev=true)))
        indexes = group.Time[ixs]
        valid_indexes = []
        for i in indexes
            if i ∈ exceptions
                continue
            end
            if !isempty(valid_indexes) && minimum(abs.(valid_indexes .- i)) < 4
                continue
            end
            push!(valid_indexes, i)
        end
        for i in valid_indexes[1:max_heat_wave_periods]
            if i > 3
                data[(i - 2):i, "HeatWaveAdditional"] = data[(i - 2):i, "futureDeltaPeak"]
                data[(i - 2):i, "futureHeatwave"] =
                    data[(i - 2):i, "futureDeltaPeak"] + data[(i - 2):i, "futureTAVG"]
            else
                data[1:i, "HeatWaveAdditional"] = data[1:i, "futureDeltaPeak"]
                data[1:i, "futureHeatwave"] =
                    data[1:i, "futureDeltaPeak"] + data[1:i, "futureTAVG"]
            end
        end
    end
    return data
end

# Run frequency function on all applicable directory files 
heatwave_dict = Dict()
for file in readdir(csv_path)
    if occursin(r"^\s*(?:heatwave_2|$)", file)
        myfile = "$(file)"
        dfname = split(myfile, '.')[1]
        data = sort!(CSV.read(joinpath(csv_path, "$(file)"), DataFrame), [:DoY])
        freq_data = add_future_frequency!(data)
        updated_heatwave = freq_data
        heatwave_dict[dfname] = updated_heatwave
    end
end

# Write results 
for (key, value) in heatwave_dict
    df = value
    CSV.write(joinpath(csv_path, "freq_$key.csv"), df)
end
