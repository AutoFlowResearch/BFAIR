# generate test_data
# Last date : 13.04.2021
# By : Matthias Mattanovich (matmat@biosustain.dtu.dk)
# This script is intended to generate sample data and save them into the
# test_data file. The saved objects will then be used to test the
# MFA_sampling functions using unit testing.
import pickle
import pandas as pd
import cobra
from BFAIR.INCA import INCA_reimport
from BFAIR.INCA.sampling import (
    add_constraints,
    find_biomass_reaction,
    get_min_solution_val,
    replace_biomass_rxn_name,
    add_feasible_constraints,
    reshape_fluxes_escher,
    bound_relaxation,
)

pd.set_option("mode.chained_assignment", None)
# Use pickle to save python variables
filehandler = open("sampling_test_data.obj", "wb")

filename = 'data/MFA_modelInputsData/TestFile.mat'
simulation_info = pd.read_csv(
    'data/MFA_modelInputsData/Re-import/experimentalMS_data_I.csv'
    )
simulation_id = 'WTEColi_113C80_U13C20_01'
model = cobra.io.load_json_model(
    'data/FIA_MS_example/database_files/iJO1366.json'
    )
reimport_data = INCA_reimport()
(fittedData,
 fittedFluxes,
 fittedFragments,
 fittedMeasuredFluxes,
 fittedMeasuredFragments,
 fittedMeasuredFluxResiduals,
 fittedMeasuredFragmentResiduals,
 simulationParameters) = reimport_data.reimport(
    filename,
    simulation_info,
    simulation_id
)

model_input = model
constrained_model = add_constraints(model_input, fittedFluxes)
find_biomass_reaction(
    model,
    biomass_string=['Biomass', 'BIOMASS', 'biomass']
)
min_val = get_min_solution_val(fittedFluxes, biomass_string='Biomass')
fittedFluxes = replace_biomass_rxn_name(
    fittedFluxes,
    biomass_rxn_name='BIOMASS_Ec_iJO1366_core_53p95M',
    biomass_string='Biomass',
)
model_input = model
feasible_constrained_model, problems = add_feasible_constraints(
    model_input, fittedFluxes, min_val=0)
### add one sampling solution and one solution type solution
fluxes_sampling = reshape_fluxes_escher(sampled_fluxes)
fluxes_solution = reshape_fluxes_escher(solution)
infeasible_model = model
cons_table = bound_relaxation(
    infeasible_model,
    fittedFluxes,
    destructive=True,
    fluxes_to_ignore=[],
)

pickle.dump(
    [
        constrained_model,
        min_val,
        fittedFluxes,
        feasible_constrained_model,
        problems,
        fluxes_sampling,
        cons_table,
    ],
    filehandler,
)

filehandler.close()
