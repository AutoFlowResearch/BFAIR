import unittest
import pickle
import pathlib
import os

# import numpy as np
# from freezegun import freeze_time
# import datetime
import pandas as pd
from BFAIR.INCA.sampling import (
    add_constraints,
    find_biomass_reaction,
    get_min_solution_val,
    replace_biomass_rxn_name,
    add_feasible_constraints,
    reshape_fluxes_escher,
    bound_relaxation,
)

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
            + "/sampling_test_data.obj",
            "rb",
        )
        (
            fittedFluxes,
            unconstraint_bounds,
            constrained_bounds,
            min_val,
            adj_fittedFluxes,
            problems,
            feasible_constrained_bounds,
            sampled_fluxes,
            fluxes_sampling,
            solution,
            fluxes_solution,
            cons_table,
            relaxed_bounds,
        ) = pickle.load(file_obj)
        file_obj.close()

        self.simulation_info = pd.read_csv(
            current_dir + '/experimentalMS_data_I.csv'
            )
        self.simulation_id = 'WTEColi_113C80_U13C20_01'
        self.model = cobra.io.load_json_model(
            current_dir + '/iJO1366.json'
            )
        self.constrained_model = self.add_constraints(
            model.copy(),
            adj_fittedFluxes
        )
        self.fittedFluxes = fittedFluxes
        self.unconstraint_bounds = unconstraint_bounds
        self.constrained_bounds = constrained_bounds
        self.min_val = min_val
        self.adj_fittedFluxes = adj_fittedFluxes
        self.problems = problems
        self.feasible_constrained_bounds = feasible_constrained_bounds
        self.sampled_fluxes = sampled_fluxes
        self.fluxes_sampling = fluxes_sampling
        self.solution = solution
        self.fluxes_solution = fluxes_solution
        self.cons_table = cons_table
        self.relaxed_bounds = relaxed_bounds

        # Add the method to compare dataframes in the class
        self.addTypeEqualityFunc(pd.DataFrame, self.assertDataframeEqual)

    @staticmethod
    def get_bounds_df(model):
        # Helper function to have a way to compare the bounds
        bounds_temp = {}
        for cnt, rxn in enumerate(model.reactions):
            bounds_temp[cnt] = {
                "rxn_id": rxn.id,
                "lb": rxn.lower_bound,
                "ub": rxn.upper_bound
            }
        return pd.DataFrame.from_dict(bounds_temp, "index")

    def test_add_constraints(self):
        unconstraint_bounds = self.unconstraint_bounds
        constrained_bounds = self.constrained_bounds
        constrained_model = self.add_constraints(
            self.model.copy(),
            self.adj_fittedFluxes,
        )
        constrained_bounds_ = self.get_bounds_df(constrained_model)
        self.assertNotEqual(unconstraint_bounds, constrained_bounds_)
        self.assertEqual(constrained_bounds, constrained_bounds_)

    def test_find_biomass_reaction(self):
        biomass_reaction_ids = self.find_biomass_reaction(
            self.constrained_model,
            biomass_string=['Biomass', 'BIOMASS', 'biomass']
        )
        # This one basically just checks if there is an output
        self.assertIsInstance(biomass_reaction_ids, list)
        # This one makes sure that the output list contains a reaction name
        self.assertIsInstance(biomass_reaction_ids[0], str)

    def test_get_min_solution_val(self):
        # Find in fittedFluxes
        # Do not find in fake example

    def test_replace_biomass_rxn_name(self):
        # Replace in fittedFluxes
        # Replace in a one liner

    def test_add_feasible_constraints(self):
        # Check if model is the same?
        # If that doesn't work, check if optimized value is almost Equal?
        # That kind of sucks
        # Maybe create a dataframe of the model bounds and compare those

    def test_reshape_fluxes_escher_sampling(self):
        # Test if a samples result has more than one value per reaction coming in and the right shape coming out
        # Compare to reference?
        # If that doesn't work, compare shape to reference

    def test_reshape_fluxes_escher_solution(self):
        # Test if a samples result has only one value per reaction coming in and the right shape coming out
        # Compare to reference?
        # If that doesn't work, compare shape to reference

    def test_bound_relaxation(self):
        try:
            # regular test
            cons_table = self.cons_table
            cons_table_ = self.bound_relaxation(
                infeasible_model,
                fittedFluxes,
                destructive=True,
                fluxes_to_ignore=[],
            )
            self.assertEqual(cons_table, cons_table_)
        except Exception as exception:
            # Exception if can't find Gurobi
            self.assertIsInstance(exception, ModuleNotFoundError)
            self.assertEqual(str(exception), "Requires Gurobi.")


if __name__ == "__main__":
    unittest.main()
