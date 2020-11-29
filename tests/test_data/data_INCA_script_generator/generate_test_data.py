# generate test_data
# Last date : 27.11.2020
# By : Matthias Mattanovich (matmat@biosustain.dtu.dk)
# This script is intended to generate sample data and save them into the
# test_data file. The saved objects will then be used to test the
# INCA_script_generator using unit testing.
import pickle
import pandas as pd
from INCA_script_generator.INCA_script_generator import *

pd.set_option('mode.chained_assignment', None)
# Use pickle to save python variables
filehandler = open("test_data.obj", 'wb')

# measured fragments/MS data, tracers and measured fluxes should be limited to
# one experiment

atomMappingReactions_data_I = pd.read_csv(
    'data_stage02_isotopomer_atomMappingReactions2.csv')
modelReaction_data_I = pd.read_csv(
    'data_stage02_isotopomer_modelReactions.csv')
atomMappingMetabolite_data_I = pd.read_csv(
    'data_stage02_isotopomer_atomMappingMetabolites.csv')
measuredFluxes_data_I = pd.read_csv(
    'data_stage02_isotopomer_measuredFluxes.csv')
experimentalMS_data_I = pd.read_csv(
    'data-1604345289079.csv')
tracer_I = pd.read_csv(
    'data_stage02_isotopomer_tracers.csv')

# The files need to be limited by model id and mapping id, I picked
# "ecoli_RL2013_02" here
atomMappingReactions_data_I = limit_to_one_model(
    atomMappingReactions_data_I, 'mapping_id', 'ecoli_RL2013_02')
modelReaction_data_I = limit_to_one_model(
    modelReaction_data_I, 'model_id', 'ecoli_RL2013_02')
atomMappingMetabolite_data_I = limit_to_one_model(
    atomMappingMetabolite_data_I, 'mapping_id', 'ecoli_RL2013_02')
measuredFluxes_data_I = limit_to_one_model(
    measuredFluxes_data_I, 'model_id', 'ecoli_RL2013_02')

# Limiting fluxes, fragments and tracers to one experiment
measuredFluxes_data_I = limit_to_one_experiment(
    measuredFluxes_data_I, 'experiment_id', 'WTEColi_113C80_U13C20_01')
experimentalMS_data_I = limit_to_one_experiment(
    experimentalMS_data_I, 'experiment_id', 'WTEColi_113C80_U13C20_01')
tracer_I = limit_to_one_experiment(
    tracer_I, 'experiment_id', 'WTEColi_113C80_U13C20_01')

# Generate variables to save
script = script_generator(modelReaction_data_I, atomMappingReactions_data_I,
                          atomMappingMetabolite_data_I,
                          measuredFluxes_data_I, experimentalMS_data_I,
                          tracer_I)
runner = runner_script_generator('TestFile', n_estimates=10)

pickle.dump([modelReaction_data_I, atomMappingReactions_data_I,
             atomMappingMetabolite_data_I, measuredFluxes_data_I,
             experimentalMS_data_I, tracer_I, script, runner], filehandler)

filehandler.close()
