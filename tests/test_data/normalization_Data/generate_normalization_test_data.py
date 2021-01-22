# generate test_data
# Last date : ##.##.2021
# By : Matthias Mattanovich (matmat@biosustain.dtu.dk)
# This script is intended to generate sample data and save them into the
# test_data file. The saved objects will then be used to test the
# tecan_od_analyzer using unit testing.
import pickle
import pandas as pd
import pathlib
import BFAIR.normlization as normalization

current_dir = str(pathlib.Path(__file__).parent.absolute())

pd.set_option("mode.chained_assignment", None)
# Use pickle to save python variables
filehandler = open("test_data.obj", "wb")

min_max = normalization.min_max_norm(
    df, columnname="Intensity", groupname_colname="sample_group_name"
)
tsi = normalization.tsi_norm(
    df, columnname="Intensity", groupname_colname="sample_group_name"
)
biomass_tsi = normalization.biomass_tmi_norm(
    biomass_substrate_df,
    df,
    columnname="Intensity",
    groupname_colname="sample_group_name",
)
biomass_formula_tsi = normalization.biomass_formula_tmi_norm(
    biomass_substrate_df,
    biomass_product_df,
    biomass_value,
    df,
    columnname="Intensity",
    groupname_colname="sample_group_name",
)
amino_acid_tsi = normalization.amino_acid_tmi_norm(
    amino_acids,
    df,
    columnname="Intensity",
    groupname_colname="sample_group_name",
)
pqn = normalization.pqn_norm(
    df,
    groupname_colname="sample_group_name",
    value_colname="Intensity",
    corr_type="median",
    qc_vector=None,
)

pickle.dump(
    [

    ],
    filehandler,
)

filehandler.close()
