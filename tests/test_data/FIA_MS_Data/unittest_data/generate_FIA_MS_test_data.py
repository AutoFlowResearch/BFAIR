# generate test_data
# Last date : 13.12.2020
# By : Matthias Mattanovich (matmat@biosustain.dtu.dk)
# This script is intended to generate sample data and save them into the
# test_data file. The saved objects will then be used to test the
# tecan_od_analyzer using unit testing.
import pickle
import pandas as pd
import pathlib
from AutoFlow_OmicsDataHandling.FIA_MS_tools import (
    extractNamesAndIntensities,
    calculateMeanVarRSD,
)

current_dir = str(pathlib.Path(__file__).parent.absolute())

pd.set_option("mode.chained_assignment", None)
# Use pickle to save python variables
filehandler = open("test_data.obj", "wb")

# load data and assign variables needed for the functions
feature_dir = current_dir + "/../features_AdditionalAdducts"
sequence_triplicates = pd.read_csv(
    current_dir +
    "/../sequence_EColi.csv", sep=";"
    )
print(sequence_triplicates.columns.values)
sample_names_triplicates = sequence_triplicates["sample_group_name"].unique()
database_triplicates = pd.read_csv(
    current_dir + "/../CHEMISTRY/iJO1366_struct.tsv", sep="\t", header=None
)
sequence_single = pd.read_csv(
    current_dir +
    "/../sequence_HumanSerum.csv", sep=";"
    )
sample_names_single = sequence_single["sample_group_name"].unique()
database_single = pd.read_csv(
    current_dir + "/../CHEMISTRY/HMDB_struct.tsv", sep="\t", header=None
)
sequence_standard = pd.read_csv(
    current_dir +
    "/../sequence_Standards.csv", sep=";"
    )
sample_names_standard = sequence_standard["sample_group_name"].unique()
database_standard = pd.read_csv(
    current_dir +
    "/../CHEMISTRY/standard_mix_struct.tsv", sep="\t", header=None
)


# extractNamesAndIntensities, triplicates
intensities_triplicates = extractNamesAndIntensities(
    feature_dir, sample_names_triplicates, database_triplicates
)
intensities_single = extractNamesAndIntensities(
    feature_dir, sample_names_single, database_single
)
intensities_standard = extractNamesAndIntensities(
    feature_dir, sample_names_standard, database_standard
)

# calculateMeanVarRSD
stats_triplicates = calculateMeanVarRSD(
    intensities_triplicates,
    sequence_triplicates.drop_duplicates(
        ["sample_group_name", "replicate_group_name"]
    ),
)
stats_single = calculateMeanVarRSD(
    intensities_single,
    sequence_single.drop_duplicates(
        ["sample_group_name", "replicate_group_name"]
    ),
)
stats_standard = calculateMeanVarRSD(
    intensities_standard,
    sequence_standard.drop_duplicates(
        ["sample_group_name", "replicate_group_name"]),
)

pickle.dump(
    [
        intensities_triplicates,
        stats_triplicates,
        intensities_single,
        stats_single,
        intensities_standard,
        stats_standard,
    ],
    filehandler,
)

filehandler.close()
