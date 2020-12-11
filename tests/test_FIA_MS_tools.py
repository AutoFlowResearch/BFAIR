import unittest
import pandas as pd
import pickle
import sys
import pathlib

sys.path.insert(1, "../")
from AutoFlow_OmicsDataHandling.FIA_MS_tools import (  # noqa E402
    extractNamesAndIntensities,
    calculateMeanVarRSD,
)

current_dir = str(pathlib.Path(__file__).parent.absolute())


class test_methods(unittest.TestCase):
    def setUp(self):
        file_obj = open(
            current_dir + "/test_data/FIA_MS_Data/unittest_data/test_data.obj",
            "rb",
        )
        (
            intensities_triplicates,
            stats_triplicates,
            intensities_single,
            stats_single,
            intensities_standard,
            stats_standard,
        ) = pickle.load(file_obj)
        file_obj.close()

        self.intensities_triplicates = intensities_triplicates
        self.stats_triplicates = stats_triplicates
        self.intensities_single = intensities_single
        self.stats_single = stats_single
        self.intensities_standard = intensities_standard
        self.stats_standard = stats_standard

    def test_extractNamesAndIntensities_triplicates(self):
        intensities_triplicates = self.intensities_triplicates
        feature_dir = current_dir + "/test_data/FIA_MS_Data/features_AdditionalAdducts"
        sequence_triplicates = pd.read_csv(
            current_dir + "/test_data/FIA_MS_Data/sequence_EColi.csv"
        )
        sample_names_triplicates = sequence_triplicates[
            "sample_group_name"
        ].unique()
        database_triplicates = pd.read_csv(
            current_dir + "/test_data/FIA_MS_Data/CHEMISTRY/iJO1366_struct.tsv",
            sep="\t",
            header=None,
        )
        intensities_triplicates_ = extractNamesAndIntensities(
            feature_dir, sample_names_triplicates, database_triplicates
        )
        self.assertEqual(
            len(intensities_triplicates), len(intensities_triplicates_)
        )
        self.assertEqual(
            len(intensities_triplicates[0]),
            len(intensities_triplicates_.loc[0]),
        )

    def test_calculateMeanVarRSD_triplicates(self):
        stats_triplicates = self.stats_triplicates
        intensities_triplicates = self.intensities_triplicates
        sequence_triplicates = pd.read_csv(
            current_dir + "/test_data/FIA_MS_Data/sequence_EColi.csv"
        )
        stats_triplicates_ = calculateMeanVarRSD(
            intensities_triplicates,
            sequence_triplicates.drop_duplicates(
                ["sample_group_name", "replicate_group_name"]
            ),
            min_reps=3,
        )
        self.assertEqual(
            len(stats_triplicates), len(stats_triplicates_)
        )
        self.assertEqual(
            len(stats_triplicates[0]),
            len(stats_triplicates_.loc[0]),
        )

    def test_extractNamesAndIntensities_single(self):
        intensities_single = self.intensities_single
        feature_dir = current_dir + "/test_data/FIA_MS_Data/features_AdditionalAdducts"
        sequence_single = pd.read_csv(
            current_dir + "/test_data/FIA_MS_Data/sequence_HumanSerum.csv"
        )
        sample_names_single = sequence_single["sample_group_name"].unique()
        database_single = pd.read_csv(
            current_dir + "/test_data/FIA_MS_Data/CHEMISTRY/HMDB_struct.tsv",
            sep="\t",
            header=None,
        )
        intensities_single_ = extractNamesAndIntensities(
            feature_dir, sample_names_single, database_single
        )
        self.assertEqual(len(intensities_single), len(intensities_single_))
        self.assertEqual(
            len(intensities_single[0]), len(intensities_single_.loc[0])
        )

    def test_calculateMeanVarRSD_single(self):
        stats_single = self.stats_single
        intensities_single = self.intensities_single
        sequence_single = pd.read_csv(
            current_dir + "/test_data/FIA_MS_Data/sequence_HumanSerum.csv"
        )
        stats_single_ = calculateMeanVarRSD(
            intensities_single,
            sequence_single.drop_duplicates(
                ["sample_group_name", "replicate_group_name"]
            ),
            min_reps=1,
        )
        self.assertEqual(len(stats_single), len(stats_single_))
        self.assertEqual(
            len(stats_single[0]), len(stats_single_.loc[0])
        )

    def test_extractNamesAndIntensities_standard(self):
        intensities_standard = self.intensities_standard
        feature_dir = current_dir + "/test_data/FIA_MS_Data/features_AdditionalAdducts"
        sequence_standard = pd.read_csv(
            current_dir + "/test_data/FIA_MS_Data//sequence_Standards.csv"
        )
        sample_names_standard = sequence_standard["sample_group_name"].unique()
        database_standard = pd.read_csv(
            "/test_data/FIA_MS_Data//CHEMISTRY/standard_mix_struct.tsv",
            sep="\t",
            header=None,
        )
        intensities_standard_ = extractNamesAndIntensities(
            feature_dir, sample_names_standard, database_standard
        )
        self.assertEqual(len(intensities_standard), len(intensities_standard_))
        self.assertEqual(
            len(intensities_standard[0]), len(intensities_standard_.loc[0])
        )

    def test_calculateMeanVarRSD_standard(self):
        stats_standard = self.stats_standard
        intensities_standard = self.intensities_standard
        sequence_standard = pd.read_csv(
            current_dir + "/test_data/FIA_MS_Data//sequence_Standards.csv"
        )
        stats_standard_ = calculateMeanVarRSD(
            intensities_standard,
            sequence_standard.drop_duplicates(
                ["sample_group_name", "replicate_group_name"]
            ),
            min_reps=1,
        )
        self.assertEqual(len(stats_standard), len(stats_standard_))
        self.assertEqual(
            len(stats_standard[0]), len(stats_standard_.loc[0])
        )
