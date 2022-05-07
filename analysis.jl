using DataStructures

function compare_decadal_trends(avg_info::Dict)
    df = DataStructures.SortedDict()
    for (name, dataset) in avg_info
        timesteps = [v["timesteps"][1:(end - 1)] for v in values(dataset)]
        avg_temp = [v["avg_temp"] for v in values(dataset)]
        rep_score = [v["replacement_efficacy_score"] for v in values(dataset)]
        sup_score = [v["suppression_efficacy_score"] for v in values(dataset)]
        days_to_achieve_FREQthresholds =
            [v["days_to_achieve_FREQthresholds"] for v in values(dataset)]
        FREQthresholds_in_percent =
            [v["FREQthresholds_in_percent"] for v in values(dataset)]

        dfnames = collect(keys(avg_info[name]))
        df[name] = DataFrame(
            name=dfnames,
            timesteps=timesteps,
            avg_temp=avg_temp,
            rep_score=rep_score,
            sup_score=sup_score,
            days_to_achieve_FREQthreshold=days_to_achieve_FREQthresholds,
            FREQthresholds_in_percent=FREQthresholds_in_percent,
        )
    end
    return df
end
