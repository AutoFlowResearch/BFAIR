import os.path

import unittest
import pandas as pd

import BFAIR.thermo as thermo


def _check_table(test: unittest.TestCase, table: pd.DataFrame, colnames):
    test.assertIsInstance(table, pd.DataFrame)
    test.assertGreater(len(table), 0)
    test.assertTrue(table.shape[1], 2)
    test.assertTrue(table.notnull().values.any())
    test.assertTrue((table.iloc[:, 0] != 0).any())
    test.assertEqual(tuple(table.columns), colnames)


class TestThermo(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # load data once
        cls.tdata = thermo.load_data("small_ecoli")
        bounds = pd.read_csv(
            os.path.join(os.path.dirname(__file__), "test_data", "test_bounds.csv")
        )
        cls.rxn_bounds = bounds[bounds["bound_type"] == "flux"]
        cls.lc_bounds = bounds[bounds["bound_type"] == "log_concentration"]


class TestIO(unittest.TestCase):
    def test_load_cbm(self):
        actual = thermo.load_cbm("small_ecoli")
        self.assertTrue(abs(actual.slim_optimize() - 0.8109621653343296) < 1e-5)

    def test_load_data(self):
        actual = thermo.load_data("small_ecoli")
        self.assertEqual(len(actual), 3)

        self.assertIsInstance(actual[0], dict)
        self.assertEqual(len(actual[0]), 4)

        self.assertIsInstance(actual[1], pd.DataFrame)
        self.assertEqual(actual[1].shape, (304, 1))

        self.assertIsInstance(actual[2], dict)
        self.assertEqual(len(actual[2]), 3)

    def test_create_model(self):
        tdata = thermo.load_data("small_ecoli")
        actual = thermo.create_model("small_ecoli", *tdata)
        self.assertTrue(abs(actual.slim_optimize() - 0.8109972502600706) < 1e-5)
        self.assertRaises(ValueError, thermo.create_model, "small_ecoli", thermo_data=tdata[0])
        self.assertRaises(ValueError, thermo.create_model, "small_ecoli", lexicon=tdata[1])
        self.assertRaises(ValueError, thermo.create_model, "small_ecoli", compartment_data=tdata[2])

    def test_adjust_model(self):
        from math import log

        model = thermo.create_model("small_ecoli")

        rxn_bounds = pd.DataFrame([{"id": "PGI", "lb": 0, "ub": 0}])
        lc_bounds = pd.DataFrame([{"id": "atp", "lb": log(5e-3), "ub": log(3e-2)}])

        thermo.adjust_model(model, rxn_bounds, lc_bounds)

        self.assertEqual(model.reactions.PGI.bounds, (0, 0))
        self.assertEqual(model.log_concentration.atp_c.variable.lb, log(5e-3))
        self.assertEqual(model.log_concentration.atp_c.variable.ub, log(3e-2))

    def test_factory(self):
        actual = thermo.models.small_ecoli
        expected = thermo.create_model("small_ecoli")
        self.assertTrue(abs(actual.slim_optimize() - expected.slim_optimize()) < 1e-5)


class TestRelaxation(TestThermo):
    def setUp(self):
        self.test_model = thermo.create_model("small_ecoli", *self.tdata)
        # add constraints to require relaxation
        thermo.adjust_model(self.test_model, self.rxn_bounds, self.lc_bounds)

    def test_relax_dgo(self):
        try:
            actual = thermo.relax_lc(self.test_model, [], True)
            _check_table(self, actual, ("bound_change", "subsystem"))
        except Exception as exception:
            self.assertIsInstance(exception, ModuleNotFoundError)
            self.assertEqual(str(exception), "Requires Gurobi.")

    def test_relax_lc(self):
        try:
            actual = thermo.relax_lc(self.test_model, [], True)
            _check_table(self, actual, ("bound_change", "compartment"))
        except Exception as exception:
            self.assertIsInstance(exception, ModuleNotFoundError)
            self.assertEqual(str(exception), "Requires Gurobi.")


class TestUtils(TestThermo):
    @classmethod
    def setUpClass(cls):
        # create one model only
        super(TestUtils, cls).setUpClass()
        cls.test_model = thermo.create_model("small_ecoli", *cls.tdata)
        cls.test_model.slim_optimize()
        cls.flux_table = thermo.get_flux(cls.test_model)
        cls.dg_table = thermo.get_delta_g(cls.test_model)
        cls.lc_table = thermo.get_log_concentration(cls.test_model)
        cls.ratio_table = thermo.get_mass_action_ratio(cls.test_model)

    def test_get_flux(self):
        _check_table(self, self.flux_table, ("flux", "subsystem"))

    def test_get_delta_g(self):
        _check_table(self, self.dg_table, ("delta_g", "subsystem"))

    def test_get_log_concentration(self):
        _check_table(self, self.lc_table, ("log_concentration", "compartment"))

    def test_get_mass_action_ratio(self):
        _check_table(self, self.ratio_table, ("mass_action_ratio", "subsystem"))

    # def test_get_dgo_bound_change(self):
    #     colnames = ["lb_in", "ub_in", "lb_change", "ub_change", "lb_out", "ub_out"]
    #     dummy_table = pd.DataFrame.from_dict(
    #         {"ACONTa": dict(zip(colnames, range(0, 6)))}, orient="index"
    #     )
    #     actual = thermo.get_dgo_bound_change(self.test_model, dummy_table)
    #     _check_table(self, actual, ("bound_change", "subsystem"))
