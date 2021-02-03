import unittest
import pickle
import pathlib
import os
# import numpy as np
from freezegun import freeze_time
from BFAIR.INCA import INCA_reimport

current_dir = str(pathlib.Path(__file__).parent.absolute())

os.chdir(current_dir + "/test_data/MFA_modelInputsData")


@freeze_time("2021-01-15 19:24:19")
class test_methods(unittest.TestCase):
    def setUp(self):
        file_obj = open(
            current_dir
            + "/test_data/MFA_modelInputsData/reimport_test_data.obj",
            "rb",
        )
        (
            filename,
            simulation_info,
            simulation_id,
            info,
            parallel,
            non_stationary,
            m,
            f,
            model_info,
            f_mnt_info,
            f_mnt_res_info,
            f_par_info,
            fittedData,
            fittedFluxes,
            fittedFragments,
            fittedMeasuredFluxes,
            fittedMeasuredFragments,
            fittedMeasuredFluxResiduals,
            fittedMeasuredFragmentResiduals,
            simulationParameters,
        ) = pickle.load(file_obj)
        file_obj.close()

        self.filename = filename
        self.simulation_info = simulation_info
        self.simulation_id = simulation_id
        self.info = info
        self.parallel = parallel
        self.non_stationary = non_stationary
        self.m = m
        self.f = f
        self.model_info = model_info
        self.f_mnt_info = f_mnt_info
        self.f_mnt_res_info = f_mnt_res_info
        self.f_par_info = f_par_info
        self.fittedData = fittedData
        self.fittedFluxes = fittedFluxes
        self.fittedFragments = fittedFragments
        self.fittedMeasuredFluxes = fittedMeasuredFluxes
        self.fittedMeasuredFragments = fittedMeasuredFragments
        self.fittedMeasuredFluxResiduals = fittedMeasuredFluxResiduals
        self.fittedMeasuredFragmentResiduals = fittedMeasuredFragmentResiduals
        self.simulationParameters = simulationParameters
        self.INCA_reimport = INCA_reimport()

    """
    The general instance tests will be added/activated once the problem in
    the other INCA module is solved.
    """

    def test_extract_file_info(self):
        info = self.INCA_reimport.extract_file_info(
            self.filename
        )
        info_ = self.info
        self.assertEqual(info, info_)

    def test_det_simulation_type(self):
        parallel, non_stationary = self.INCA_reimport.det_simulation_type(
            self.simulation_info
        )
        parallel_ = self.parallel
        non_stationary_ = self.non_stationary
        self.assertEqual(parallel, parallel_)
        self.assertEqual(non_stationary, non_stationary_)

    def test_data_extraction(self):
        """
        comparing m fails
        """
        m, f = self.INCA_reimport.data_extraction(
            self.filename
        )
        m_ = self.m
        f_ = self.f
        self.assertEqual(m.tolist(), m_.tolist())
        # self.np.testing.assert_array_equal(m, m_)
        self.assertEqual(f, f_)

    """
    # datetime issue
    def test_extract_model_info(self):
        model_info = self.INCA_reimport.extract_model_info(self.m)
        model_info_ = self.model_info
        self.assertEqual(model_info, model_info_)
    """

    def test_extract_sim_params(self):
        simulationParameters = self.INCA_reimport.extract_sim_params(
            self.simulation_id, self.info, self.m,
            self.filename
        )
        simulationParameters_ = self.simulationParameters
        self.assertEqual(simulationParameters, simulationParameters_)

    # def test_simulationParameters_testdata(self):
    #     simulationParameters_ = self.simulationParameters
    #     self.assertTrue(len(simulationParameters_) > 0)

    def test_extract_base_stats(self):
        fittedData = self.INCA_reimport.extract_base_stats(
            self.f, self.simulation_id, self.info
        )
        fittedData_ = self.fittedData
        self.assertEqual(fittedData, fittedData_)

    # def test_fittedData_testdata(self):
    #     fittedData_ = self.fittedData
    #     self.assertTrue(len(fittedData_) > 0)

    def test_get_fit_info(self):
        f_mnt_info = self.INCA_reimport.get_fit_info(self.f)
        f_mnt_info_ = self.f_mnt_info
        self.assertEqual(f_mnt_info, f_mnt_info_)

    def test_sort_fit_info(self):
        (
            fittedMeasuredFluxes,
            fittedMeasuredFragments,
        ) = self.INCA_reimport.sort_fit_info(
            self.f_mnt_info, self.simulation_info, self.fittedData
        )
        fittedMeasuredFluxes_ = self.fittedMeasuredFluxes
        fittedMeasuredFragments_ = self.fittedMeasuredFragments
        self.assertEqual(fittedMeasuredFluxes, fittedMeasuredFluxes_)
        self.assertEqual(fittedMeasuredFragments, fittedMeasuredFragments_)

    # def test_fittedMeasuredFluxes_testdata(self):
    #     fittedMeasuredFluxes_ = self.fittedMeasuredFluxes
    #     self.assertTrue(len(fittedMeasuredFluxes_) > 0)

    # def test_fittedMeasuredFragments_testdata(self):
    #     fittedMeasuredFragments_ = self.fittedMeasuredFragments
    #     self.assertTrue(len(fittedMeasuredFragments_) > 0)

    def test_get_residuals_info(self):
        f_mnt_res_info = self.INCA_reimport.get_residuals_info(
            self.f, self.simulation_info
        )
        f_mnt_res_info_ = self.f_mnt_res_info
        self.assertEqual(f_mnt_res_info, f_mnt_res_info_)

    def test_sort_residual_info(self):
        (
            fittedMeasuredFluxResiduals,
            fittedMeasuredFragmentResiduals,
        ) = self.INCA_reimport.sort_residual_info(
            self.f_mnt_res_info, self.simulation_info, self.fittedData
        )
        fittedMeasuredFluxResiduals_ = self.fittedMeasuredFluxResiduals
        fittedMeasuredFragmentResiduals_ = self.fittedMeasuredFragmentResiduals
        self.assertEqual(
            fittedMeasuredFluxResiduals, fittedMeasuredFluxResiduals_
        )
        self.assertEqual(
            fittedMeasuredFragmentResiduals, fittedMeasuredFragmentResiduals_
        )

    # def test_fittedMeasuredFluxResiduals_testdata(self):
    #     fittedMeasuredFluxResiduals_ = self.fittedMeasuredFluxResiduals
    #     self.assertTrue(len(fittedMeasuredFluxResiduals_) > 0)

    # def test_fittedMeasuredFragmentResiduals_testdata(self):
    #     fittedMeasuredFragmentResiduals_ =
    # self.fittedMeasuredFragmentResiduals
    #     self.assertTrue(len(fittedMeasuredFragmentResiduals_) > 0)

    def test_get_fitted_parameters(self):
        f_par_info = self.INCA_reimport.get_fitted_parameters(
            self.f, self.simulation_info
        )
        f_par_info_ = self.f_par_info
        self.assertEqual(f_par_info, f_par_info_)

    def test_sort_parameter_info(self):
        fittedFluxes, fittedFragments = self.INCA_reimport.sort_parameter_info(
            self.f_par_info, self.simulation_info, self.fittedData
        )
        fittedFluxes_ = self.fittedFluxes
        fittedFragments_ = self.fittedFragments
        self.assertEqual(fittedFluxes, fittedFluxes_)
        self.assertEqual(fittedFragments, fittedFragments_)

    # def test_fittedFluxes_testdata(self):
    #     fittedFluxes_ = self.fittedFluxes
    #     self.assertTrue(len(fittedFluxes_) > 0)

    # def test_fittedFragments_testdata(self):
    #     fittedFragments_ = self.fittedFragments
    #     self.assertTrue(len(fittedFragments_) > 0)

    """
    # datetime issue
    def test_reimport(self):
        (
            fittedData,
            fittedFluxes,
            fittedFragments,
            fittedMeasuredFluxes,
            fittedMeasuredFragments,
            fittedMeasuredFluxResiduals,
            fittedMeasuredFragmentResiduals,
            simulationParameters,
        ) = self.INCA_reimport.reimport(
            self.filename,
            self.simulation_info,
            self.simulation_id
        )
        fittedData_ = self.fittedData
        fittedFluxes_ = self.fittedFluxes
        fittedFragments_ = self.fittedFragments
        fittedMeasuredFluxes_ = self.fittedMeasuredFluxes
        fittedMeasuredFragments_ = self.fittedMeasuredFragments
        fittedMeasuredFluxResiduals_ = self.fittedMeasuredFluxResiduals
        fittedMeasuredFragmentResiduals_ = self.fittedMeasuredFragmentResiduals
        simulationParameters_ = self.simulationParameters
        self.assertEqual(fittedData, fittedData_)
        self.assertEqual(fittedFluxes, fittedFluxes_)
        self.assertEqual(fittedFragments, fittedFragments_)
        self.assertEqual(fittedMeasuredFluxes, fittedMeasuredFluxes_)
        self.assertEqual(fittedMeasuredFragments, fittedMeasuredFragments_)
        self.assertEqual(
            fittedMeasuredFluxResiduals, fittedMeasuredFluxResiduals_
        )
        self.assertEqual(
            fittedMeasuredFragmentResiduals, fittedMeasuredFragmentResiduals_
        )
        self.assertEqual(simulationParameters, simulationParameters_)
    """
