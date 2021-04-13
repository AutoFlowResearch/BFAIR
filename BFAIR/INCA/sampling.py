import pandas as pd
import numpy as np
import time
import ast
import sys
import escher
import cobra
from gurobipy import Model as GRBModel
import re
from BFAIR.INCA import INCA_reimport
import functools

def timer(func):
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        print("--- start ---")
        start_time = time.time()
        output = func(*args, **kwargs)
        elapsed_time = (time.time() - start_time)
        min_, sec = divmod(round(elapsed_time), 60)
        hour, min_ = divmod(min_, 60)
        print(f'{func.__name__} takes {hour}h: {min_}min: {sec}sec to run')
        print("--- end ---")
        return output
    return wrapper


def _reshape_fluxes(fittedFluxes):
    MFA_flux_bounds = {}
    for i, rxn in fittedFluxes.iterrows():
        MFA_flux_bounds[rxn['rxn_id']] = [rxn['flux_ub'], rxn['flux_lb']]
    return MFA_flux_bounds

def _adjust_bounds(model, rxn, bounds):
    skip = False
    if bounds[0] < bounds[1]: # to fix the issue with negaive values above
        try:
            model.reactions.get_by_id(rxn).lower_bound = round(bounds[0], 1)
            model.reactions.get_by_id(rxn).upper_bound = round(bounds[1], 1)
        except KeyError:
            #print(f'Did not work for {rxn}')
            skip = True
    else:
        try:
            model.reactions.get_by_id(rxn).upper_bound = round(bounds[0], 1)
            model.reactions.get_by_id(rxn).lower_bound = round(bounds[1], 1)
        except KeyError:
            #print(f'Did not work for', rxn)
            skip = True
    return model, skip

def _restart_message(restart_counter):
    extra_minus = '-' * len(str(restart_counter))
    print(f'--------------------------{extra_minus}')
    print(f'Total number of restarts: {restart_counter}')

@timer
def add_constraints(model_input, fittedFluxes):
    MFA_fluxes = _reshape_fluxes(fittedFluxes)
    model = model_input.copy()
    for rxn, bounds in MFA_fluxes.items():
        model, skip = _adjust_bounds(model, rxn, bounds)
    return model

def find_biomass_reaction(
    model,
    biomass_string=['Biomass', 'BIOMASS', 'biomass']
):
    if type(biomass_string) == 'list':
        biomass = biomass_string
    else:
        biomass = list(biomass_string)
    biomass_reaction_ids = []
    for reaction in model.reactions:
        for biomass_spelling in biomass:
            if biomass_spelling in reaction.id:
                biomass_reaction_ids.append(reaction.id)
    return biomass_reaction_ids

def get_min_solution_val(fittedFluxes, biomass_string='Biomass'):
    min_val = 0
    for cnt, name in enumerate(fittedFluxes['rxn_id']):
        if biomass_string in name:
            min_val = fittedFluxes.at[cnt, 'flux']
    return min_val

def replace_biomass_rxn_name(
    fittedFluxes,
    biomass_rxn_name,
    biomass_string='Biomass',
):
    for cnt, name in enumerate(fittedFluxes['rxn_id']):
        if biomass_string in name:
            fittedFluxes.at[cnt, 'rxn_id'] = biomass_rxn_name

@timer
def add_feasible_constraints(model_input, fittedFluxes, min_val=0):
    no_restart = False
    problems = []
    restart_counter = 0
    MFA_fluxes = _reshape_fluxes(fittedFluxes)
    while no_restart == False:
        model = model_input.copy()
        for rxn, bounds in MFA_fluxes.items():
            if rxn in problems:
                continue
            model, skip = _adjust_bounds(model, rxn, bounds)
            if skip: continue
            solution_after_adj = model.optimize()
            if solution_after_adj.objective_value is not None and\
                solution_after_adj.objective_value >= min_val:
                no_restart = True
            else:
                print(f'Solution infeasible if adding {rxn}')
                problems.append(rxn)
                no_restart = False
                restart_counter += 1
                break
    _restart_message(restart_counter)
    return model, problems

def reshape_fluxes_escher(sampled_fluxes):
    if type(sampled_fluxes) is pd.core.frame.DataFrame:
        fluxes_escher = {}
        for col in sampled_fluxes.columns:
            fluxes_escher[col] = sampled_fluxes[col].mean()
    elif type(sampled_fluxes) is cobra.core.solution.Solution:
        fluxes_escher = {}
        reactions = list(sampled_fluxes.fluxes.index)
        for cnt, row in enumerate(sampled_fluxes.fluxes):
            fluxes_escher[reactions[cnt]] = row
    else:
        raise TypeError(f"The input is a '{type(sampled_fluxes)}', this type of object cannot be used here")
    return fluxes_escher

def _extract_bound_violations(grb_model, variables, negative_slack_prefix, positive_slack_prefix):
    # Extracts the magnitude of the bound violations on the relaxed Gurobi model
    epsilon = grb_model.Params.FeasibilityTol
    violations = []
    for var in variables:
        lb_change = -grb_model.getVarByName(negative_slack_prefix + var).X
        ub_change = grb_model.getVarByName(positive_slack_prefix + var).X
        if lb_change < -epsilon:
            violations.append((var, lb_change - epsilon, 0))
        elif ub_change > epsilon:
            violations.append((var, 0, ub_change + epsilon))
    return violations

@timer
def bound_relaxation(infeasible_model, fittedFluxes, destructive=True):
    model = infeasible_model
    model.solver = 'optlang-gurobi'
    if model.solver.interface.__name__ != "optlang.gurobi_interface":
        raise ModuleNotFoundError("Requires Gurobi solver.")

    # copy Gurobi model
    grb_model: GRBModel = model.solver.problem.copy()
    MFA_fluxes = _reshape_fluxes(fittedFluxes)
    reactions_to_relax = list(MFA_fluxes.keys())

    # There are no metabolites to relax to we leave relax_cons empty
    relax_cons = []

    # Get forward and reverse reaction separately
    relax_vars = [
        grb_model.getVarByName(con.id) for con in model.reactions if con.id in reactions_to_relax
    ]
    relax_vars_reverse = [
        grb_model.getVarByName(con.reverse_id) for con in model.reactions if con.id in reactions_to_relax
    ]

    # Combine the lists
    relax_vars = relax_vars + relax_vars_reverse

    cons_penalties = relax_cons
    vars_penalties = [1] * len(relax_vars)

    # perform relaxation of variable bounds
    relax_obj = grb_model.feasRelax(
        relaxobjtype=0,
        minrelax=True,
        vars=relax_vars,
        lbpen=vars_penalties,
        ubpen=vars_penalties,
        constrs=relax_cons,
        rhspen=cons_penalties,
    )

    # check if relaxation was successful
    if relax_obj <= 0:
        raise Infeasible("Failed to create the feasibility relaxation!")

    grb_model.optimize()

    # transfer/record the lower/upper bound changes to the cobra model
    if relax_vars:
        variables = [var.VarName for var in relax_vars]
        rows = []
        violations = _extract_bound_violations(grb_model, variables, "ArtL_", "ArtU_")
        for rxn_id, lb_change, ub_change in violations:
            if destructive:
                if '_reverse_' in rxn_id:
                    rxn_id = re.match('.+?(?=_reverse_)', rxn_id)[0]
                    lb_change = -ub_change
                    ub_change = 0
                lb = model.reactions.get_by_id(rxn_id).lower_bound
                ub = model.reactions.get_by_id(rxn_id).upper_bound
                model.reactions.get_by_id(rxn_id).lower_bound = lb + lb_change
                model.reactions.get_by_id(rxn_id).upper_bound = ub + ub_change
            rows.append(
                {
                    "reaction": rxn_id,
                    "bound_change": lb_change + ub_change,
                    "subsystem": model.reactions.get_by_id(rxn_id).subsystem,
                }
            )
        cons_table = pd.DataFrame.from_records(rows, index="reaction")
    else:
        cons_table = None

    return cons_table
