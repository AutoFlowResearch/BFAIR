import unittest
import pickle
import pathlib
import os

# import numpy as np
# from freezegun import freeze_time
# import datetime
import pandas as pd
from BFAIR.INCA import INCA_reimport

current_dir = str(pathlib.Path(__file__).parent.absolute())

os.chdir(current_dir + "/test_data/MFA_modelInputsData")


class test_methods(unittest.TestCase):

    maxDiff = None

    # Create method to compare dataframes
    def assertDataframeEqual(self, a, b, msg):
        try:
            pd.testing.assert_frame_equal(a, b)
        except AssertionError as e:
            raise self.failureException(msg) from e

    def setUp(self):
        file_obj = open(
            current_dir
            + "/test_data/MFA_modelInputsData/sampling_test_data.obj",
            "rb",
        )
        (
            name
        ) = pickle.load(file_obj)
        file_obj.close()

        self.name = name

        # Add the method to compare dataframes in the class
        self.addTypeEqualityFunc(pd.DataFrame, self.assertDataframeEqual)







if __name__ == "__main__":
    unittest.main()