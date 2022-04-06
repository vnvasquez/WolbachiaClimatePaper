function get_smoothed_temperature_value(
    smoothed_temperature_data::Vector{Float64},
    t,
    temperature_model::TimeSeriesTemperature,
    temp_value_from_inputs::Float64,
)
    return smoothed_temperature_data[Int(floor(t))]
end

function get_smoothed_temperature_value(
    smoothed_temperature_data::Nothing,
    t,
    temperature_model::TimeSeriesTemperature,
    temp_value_from_inputs::Float64,
)
    return GeneDrive.get_temperature_value(temperature_model, temp_value_from_inputs, t)
end

function format_dynamic_model_results(data::WolbachiaExperiment, sol)
    return format_dynamic_model_results(data.node, sol)
end

function format_dynamic_model_results(data::WolbachiaExperiment)
    return format_dynamic_model_results(data.node)
end

# Because these two accessor functions are not currently exported 
import GeneDrive: get_organisms, count_genotypes

function get_organisms(data::WolbachiaExperiment)
    return get_organisms(data.node)
end

function count_genotypes(data::WolbachiaExperiment, species::Type{AedesAegypti})
    return count_genotypes(data.node, species)
end
