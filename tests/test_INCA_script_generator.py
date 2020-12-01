import unittest
import pickle
import sys
import pathlib
sys.path.insert(1, '../')
from INCA_script_generator.INCA_script_generator import *

current_dir = str(pathlib.Path(__file__).parent.absolute())


class test_methods(unittest.TestCase):

    def setUp(self):
        file_obj = open(current_dir + '/test_data/MFA_modelInputsData/test_data.obj', 'rb')
        modelReaction_data_I, atomMappingReactions_data_I, \
            atomMappingMetabolite_data_I, measuredFluxes_data_I, \
            experimentalMS_data_I, tracer_I, script, \
            runner = pickle.load(file_obj)
        file_obj.close()

        self.modelReaction_data_I = modelReaction_data_I
        self.atomMappingReactions_data_I = atomMappingReactions_data_I
        self.atomMappingMetabolite_data_I = atomMappingMetabolite_data_I
        self.measuredFluxes_data_I = measuredFluxes_data_I
        self.experimentalMS_data_I = experimentalMS_data_I
        self.tracer_I = tracer_I
        self.script = script
        self.runner = runner

    def test_script_generator(self):
        script_ = script_generator(self.modelReaction_data_I,
                                   self.atomMappingReactions_data_I,
                                   self.atomMappingMetabolite_data_I,
                                   self.measuredFluxes_data_I,
                                   self.experimentalMS_data_I, self.tracer_I)
        script = self.script
        self.assertEqual(script, script_)

    def test_runner_generator(self):
        runner_ = runner_script_generator('TestFile', n_estimates=10)
        runner = self.runner
        self.assertEqual(runner, runner_)
        