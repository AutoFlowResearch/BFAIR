import numpy as np
import pandas as pd
from pyopenms import FeatureMap, FeatureXMLFile


def extractNamesAndIntensities(feature_dir, sample_names, database):
    metabolites_unique = database[0].unique()
    extracted_data_dict = {}
    cnt = 0
    for name in sample_names:
        features = FeatureMap()
        FeatureXMLFile().load(
            feature_dir + "/" + name + ".featureXML", features
        )
        for f in features:
            if (
                f.getMetaValue("PeptideRef").decode("utf-8")
                in metabolites_unique
            ):
                formula = database[
                    database[0] == f.getMetaValue("PeptideRef").decode("utf-8")
                ][1]
                extracted_data_dict[cnt] = {
                    "sample_group_name": name,
                    "Metabolite": f.getMetaValue("PeptideRef").decode("utf-8"),
                    "Formula": list(formula)[0],
                    "Intensity": f.getMetaValue("peak_apex_int"),
                }
                cnt = cnt + 1
    return pd.DataFrame.from_dict(extracted_data_dict, "index")


def calculateMeanVarRSD(
    extracted_data_all, sample_name_2_replicate_groups, min_reps=3
):
    stats_all_dict = {}
    cnt = 0
    for replicate_group in sample_name_2_replicate_groups[
        "replicate_group_name"
    ].unique():
        list_of_compound_intensities = pd.DataFrame(
            columns=["Metabolite", "Formula", "Intensities"]
        )
        for sample_name_index, sample_name in sample_name_2_replicate_groups[
            sample_name_2_replicate_groups["replicate_group_name"]
            == replicate_group
        ].iterrows():
            for met_index, met in extracted_data_all[
                extracted_data_all["sample_group_name"]
                == sample_name["sample_group_name"]
            ].iterrows():
                list_of_compound_intensities = (
                    list_of_compound_intensities.append(
                        pd.DataFrame(
                            [
                                {
                                    "Metabolite": met["Metabolite"],
                                    "Formula": met["Formula"],
                                    "Intensities": met["Intensity"],
                                }
                            ]
                        ),
                        ignore_index=True,
                    )
                )
        for met_name in list_of_compound_intensities["Metabolite"].unique():
            intensities = list_of_compound_intensities[
                list_of_compound_intensities["Metabolite"] == met_name
            ]
            if len(intensities) >= min_reps:
                mean = np.mean(intensities["Intensities"])
                var = np.var(intensities["Intensities"])
                rsd = np.sqrt(var) / mean
                stats_all_dict[cnt] = {
                    "replicate_group_name": replicate_group,
                    "Metabolite": met_name,
                    "Formula": list(
                        list_of_compound_intensities[
                            list_of_compound_intensities["Metabolite"]
                            == met_name
                        ]["Formula"]
                    )[0],
                    "Mean": mean,
                    "Variance": var,
                    "RSD": rsd,
                }
                cnt = cnt + 1
    return pd.DataFrame.from_dict(stats_all_dict, "index")
