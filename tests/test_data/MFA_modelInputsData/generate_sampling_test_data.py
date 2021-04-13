# generate test_data
# Last date : 13.04.2021
# By : Matthias Mattanovich (matmat@biosustain.dtu.dk)
# This script is intended to generate sample data and save them into the
# test_data file. The saved objects will then be used to test the
# MFA_sampling functions using unit testing.
import pickle
import pandas as pd
from BFAIR.INCA.sampling import (

)

pd.set_option("mode.chained_assignment", None)
# Use pickle to save python variables
filehandler = open("sampling_test_data.obj", "wb")



pickle.dump(
    [
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
    ],
    filehandler,
)

filehandler.close()
