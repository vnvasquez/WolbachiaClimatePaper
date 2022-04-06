using GeneDrive
using StatsBase

mutable struct WolbachiaExperiment
    node::Node
    smoothing_type::Union{Nothing, String, Tuple{String, Float64}}
    smoothing_window::Union{Nothing, Int64}
    smoothed_temperature_data::Union{Nothing, Vector{Float64}}
end

function WolbachiaExperiment(node, smoothing_type::Tuple, smoothing_window)
    smoothed_temperature_data = [
        i < smoothing_window ? StatsBase.mean(node.temperature.values[begin:i]) :
        StatsBase.mean(
            node.temperature.values[(i - smoothing_window + 1):i],
            StatsBase.eweights(smoothing_window, smoothing_type[2]),
        ) for i in 1:length(node.temperature.values)
    ]

    WolbachiaExperiment(node, smoothing_type, smoothing_window, smoothed_temperature_data)
end

function WolbachiaExperiment(node, smoothing_type::String, smoothing_window)
    smoothed_temperature_data = [
        i < smoothing_window ? StatsBase.mean(node.temperature.values[begin:i]) :
        StatsBase.mean(node.temperature.values[(i - smoothing_window + 1):i]) for
        i in 1:length(node.temperature.values)
    ]

    WolbachiaExperiment(node, smoothing_type, smoothing_window, smoothed_temperature_data)
end

function WolbachiaExperiment(node)
    WolbachiaExperiment(node, nothing, nothing, nothing)
end
