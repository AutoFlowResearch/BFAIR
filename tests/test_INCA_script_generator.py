import unittest
import pickle
import sys
import pathlib

sys.path.insert(1, "../")
from AutoFlow_OmicsDataHandling.INCA_script_generator import (  # noqa E402
    limit_to_one_model,
    limit_to_one_experiment,
    # prepare_input,
    initiate_MATLAB_script,
    # reaction_mapping,
    add_reactions_to_script,
    initialize_model,
    symmetrical_metabolites,
    unbalanced_reactions,
    add_reaction_parameters,
    verify_and_estimate,
    add_experimental_parameters,
    mapping,
    script_generator,
    # save_INCA_script,
    runner_script_generator,
    # save_runner_script,
    # run_INCA_in_MATLAB,
)

current_dir = str(pathlib.Path(__file__).parent.absolute())


class test_methods(unittest.TestCase):
    def setUp(self):
        file_obj = open(
            current_dir + "/test_data/MFA_modelInputsData/test_data.obj", "rb"
        )
        (
            modelReaction_data_I,
            atomMappingReactions_data_I,
            atomMappingMetabolite_data_I,
            measuredFluxes_data_I,
            experimentalMS_data_I,
            tracer_I,
            initiated_MATLAB_script,
            model_reactions,
            model_rxn_ids,
            initialized_model,
            symmetrical_metabolites_script,
            unbalanced_reactions_script,
            reaction_parameters,
            verify_and_estimate_script,
            experimental_parameters,
            fragments_used,
            mapping_script,
            script,
            runner,
        ) = pickle.load(file_obj)
        file_obj.close()

        self.modelReaction_data_I = modelReaction_data_I
        self.atomMappingReactions_data_I = atomMappingReactions_data_I
        self.atomMappingMetabolite_data_I = atomMappingMetabolite_data_I
        self.measuredFluxes_data_I = measuredFluxes_data_I
        self.experimentalMS_data_I = experimentalMS_data_I
        self.tracer_I = tracer_I
        self.initiated_MATLAB_script = initiated_MATLAB_script
        self.model_reactions = model_reactions
        self.model_rxn_ids = model_rxn_ids
        self.initialized_model = initialized_model
        self.symmetrical_metabolites_script = symmetrical_metabolites_script
        self.unbalanced_reactions_script = unbalanced_reactions_script
        self.reaction_parameters = reaction_parameters
        self.verify_and_estimate_script = verify_and_estimate_script
        self.experimental_parameters = experimental_parameters
        self.fragments_used = fragments_used
        self.mapping_script = mapping_script
        self.script = script
        self.runner = runner

    """
    For now, all the tests are based on comparing the output of the custom
    functions to manually checked previously generated instances
    """

    def test_initiate_MATLAB_script(self):
        """
        Tests initiate_MATLAB_script() function.
        Checks if the output is the same as in a previously generated instance
        """
        initiated_MATLAB_script_ = initiate_MATLAB_script()
        initiated_MATLAB_script = self.initiated_MATLAB_script
        self.assertEqual(initiated_MATLAB_script_, initiated_MATLAB_script)

    def test_add_reactions_to_script(self):
        """
        Tests add_reactions_to_script() function.
        This function depends on prepare_input() and reaction_mapping().
        This test also serves as a test for them.
        Checks if the output is the same as in a previously generated instance
        """
        model_reactions_, model_rxn_ids_ = add_reactions_to_script(
            self.modelReaction_data_I,
            self.atomMappingReactions_data_I,
        )
        model_reactions = self.model_reactions
        model_rxn_ids = self.model_rxn_ids
        self.assertEqual(model_reactions_, model_reactions)
        self.assertEqual(model_rxn_ids_, model_rxn_ids)

    def test_initialize_model(self):
        """
        Tests  function.
        Checks if the output is the same as in a previously generated instance
        """
        initialized_model_ = initialize_model()
        initialized_model = self.initialized_model
        self.assertEqual(initialized_model_, initialized_model)

    def test_symmetrical_metabolites(self):
        """
        Tests symmetrical_metabolites() function.
        Checks if the output is the same as in a previously generated instance
        """
        symmetrical_metabolites_script_ = symmetrical_metabolites(
            self.atomMappingMetabolite_data_I
        )
        symmetrical_metabolites_script = self.symmetrical_metabolites_script
        self.assertEqual(
            symmetrical_metabolites_script_, symmetrical_metabolites_script
        )

    # this one does not have an output in this set up, change for future
    # version
    def test_unbalanced_reactions(self):
        """
        Tests unbalanced_reactions() function.
        Checks if the output is the same as in a previously generated instance
        """
        unbalanced_reactions_script_ = unbalanced_reactions(
            self.atomMappingMetabolite_data_I,
        )
        unbalanced_reactions_script = self.unbalanced_reactions_script
        self.assertEqual(
            unbalanced_reactions_script_, unbalanced_reactions_script
        )

    def test_add_reaction_parameters(self):
        """
        Tests add_reaction_parameters() function.
        Checks if the output is the same as in a previously generated instance
        """
        reaction_parameters_ = add_reaction_parameters(
            self.modelReaction_data_I,
            self.measuredFluxes_data_I,
            self.model_rxn_ids,
            fluxes_present=True,
        )
        reaction_parameters = self.reaction_parameters
        self.assertEqual(reaction_parameters_, reaction_parameters)

    def test_verify_and_estimate(self):
        """
        Tests verify_and_estimate() function.
        Checks if the output is the same as in a previously generated instance
        """
        verify_and_estimate_script_ = verify_and_estimate()
        verify_and_estimate_script = self.verify_and_estimate_script
        self.assertEqual(
            verify_and_estimate_script_, verify_and_estimate_script
        )

    def test_add_experimental_parameters(self):
        """
        Tests add_experimental_parameters() function.
        Checks if the output is the same as in a previously generated instance
        """
        (
            experimental_parameters_,
            fragments_used_,
        ) = add_experimental_parameters(
            self.experimentalMS_data_I,
            self.tracer_I,
            self.measuredFluxes_data_I,
            self.atomMappingMetabolite_data_I,
        )
        experimental_parameters = self.experimental_parameters
        fragments_used = self.fragments_used
        self.assertEqual(experimental_parameters_, experimental_parameters)
        self.assertEqual(fragments_used_, fragments_used)

    def test_mapping(self):
        """
        Tests mapping() function.
        Checks if the output is the same as in a previously generated instance
        """
        mapping_script_ = mapping(
            self.experimentalMS_data_I, self.fragments_used
        )
        mapping_script = self.mapping_script

        self.assertEqual(mapping_script_, mapping_script)

    def test_script_generator(self):
        """
        Compares the output of the generated script with the manually curated
        one.
        """
        script_ = script_generator(
            self.modelReaction_data_I,
            self.atomMappingReactions_data_I,
            self.atomMappingMetabolite_data_I,
            self.measuredFluxes_data_I,
            self.experimentalMS_data_I,
            self.tracer_I,
        )
        script = self.script
        self.assertEqual(script_, script)

    def test_runner_generator(self):
        """
        Compares the output of the generated runner script with the manually
        curated one.
        """
        runner_ = runner_script_generator("TestFile", n_estimates=10)
        runner = self.runner
        self.assertEqual(runner_, runner)
