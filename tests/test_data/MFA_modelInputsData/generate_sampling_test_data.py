# generate test_data
# Last date : 13.04.2021
# By : Matthias Mattanovich (matmat@biosustain.dtu.dk)
# This script is intended to generate sample data and save them into the
# test_data file. The saved objects will then be used to test the
# MFA_sampling functions using unit testing.
import pickle
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

pd.set_option("mode.chained_assignment", None)
# Use pickle to save python variables
filehandler = open("sampling_test_data.obj", "wb")

add_constraints(model_input, fittedFluxes)
find_biomass_reaction(
    model,
    biomass_string=['Biomass', 'BIOMASS', 'biomass']
)
get_min_solution_val(fittedFluxes, biomass_string='Biomass')
replace_biomass_rxn_name(
    fittedFluxes,
    biomass_rxn_name,
    biomass_string='Biomass',
)
add_feasible_constraints(model_input, fittedFluxes, min_val=0)
reshape_fluxes_escher(sampled_fluxes)
bound_relaxation(
    infeasible_model,
    fittedFluxes,
    destructive=True,
    fluxes_to_ignore=[],
)

pickle.dump(
    [
        name,
    ],
    filehandler,
)

filehandler.close()
