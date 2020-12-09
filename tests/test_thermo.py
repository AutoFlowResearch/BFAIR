import os.path

import unittest
import pandas as pd

import AutoFlow_OmicsDataHandling.thermo as thermo


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
        cls.tdata = thermo.load_data("iJO1366")
        bounds = pd.read_csv(
            os.path.join(os.path.dirname(__file__), "data", "test_bounds.csv")
        )
        cls.rxn_bounds = bounds[bounds["bound_type"] == "flux"]
        cls.lc_bounds = bounds[bounds["bound_type"] == "log_concentration"]


class TestIO(unittest.TestCase):
    def test_load_cbm(self):
        actual = thermo.load_cbm("iJO1366")
        self.assertTrue(abs(actual.slim_optimize() - 0.8139991066760914) < 1e-5)

    def test_load_data(self):
        actual = thermo.load_data("iJO1366")
        self.assertEqual(len(actual), 3)

        self.assertIsInstance(actual[0], dict)
        self.assertEqual(len(actual[0]), 4)

        self.assertIsInstance(actual[1], pd.DataFrame)
        self.assertEqual(actual[1].shape, (1807, 1))

        self.assertIsInstance(actual[2], dict)
        self.assertEqual(len(actual[2]), 3)

    def test_create_model(self):
        tdata = thermo.load_data("iJO1366")
        actual = thermo.create_model("iJO1366", *tdata)

        self.assertEqual(actual.slim_optimize(), 0.0)

    def test_adjust_model(self):
        from math import log

        tdata = thermo.load_data("iJO1366")
        model = thermo.create_model("iJO1366", *tdata)

        rxn_bounds = pd.DataFrame([{"id": "PGI", "lb": 0, "ub": 0}])
        lc_bounds = pd.DataFrame([{"id": "atp", "lb": log(5e-3), "ub": log(3e-2)}])

        thermo.adjust_model(model, rxn_bounds, lc_bounds)

        self.assertEqual(model.reactions.PGI.bounds, (0, 0))
        self.assertEqual(model.log_concentration.atp_c.variable.lb, log(5e-3))
        self.assertEqual(model.log_concentration.atp_c.variable.ub, log(3e-2))


class TestRelaxation(TestThermo):
    def setUp(self):
        self.test_model = thermo.create_model("iJO1366", *self.tdata)
        thermo.adjust_model(self.test_model, self.rxn_bounds, self.lc_bounds)

    def test_relax_dgo(self):
        actual = thermo.relax_dgo(self.test_model, [])
        self.assertEqual(len(actual), 2)
        _check_table(self, actual[1], ("bound_change", "subsystem"))

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
        cls.test_model = thermo.create_model("iJO1366", *cls.tdata)
        thermo.adjust_model(cls.test_model, cls.rxn_bounds, cls.lc_bounds)
        cls.test_model.slim_optimize()

    def test_get_flux(self):
        actual = thermo.get_flux(self.test_model)
        _check_table(self, actual, ("flux", "subsystem"))

    def test_get_delta_g(self):
        actual = thermo.get_delta_g(self.test_model)
        _check_table(self, actual, ("delta_g", "subsystem"))

    def test_get_log_concentration(self):
        actual = thermo.get_log_concentration(self.test_model)
        _check_table(self, actual, ("log_concentration", "compartment"))

    def test_get_dgo_bound_change(self):
        thermo.adjust_model(self.test_model, self.rxn_bounds, self.lc_bounds)
        actual = thermo.get_dgo_bound_change(*thermo.relax_dgo(self.test_model))
        _check_table(self, actual, ("bound_change", "subsystem"))
