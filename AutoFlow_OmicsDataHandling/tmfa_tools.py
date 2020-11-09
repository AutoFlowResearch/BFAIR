import pandas as pd
from pytfa.optim import DeltaG, LogConcentration


def get_flux(tmodel, reaction_or_id):
    """
    Returns the flux of the specified reaction.
    
    Parameters:
        tmodel (pyTFA model): Model of interest
        reaction_or_id (cobra.Reaction or String): Reaction of interest

    Returns:
        flux: Resulting flux

    """
    rxn = (
        reaction_or_id
        if not isinstance(reaction_or_id, str)
        else tmodel.reactions.get_by_id(reaction_or_id)
    )
    return tmodel.variables[rxn.id].primal - tmodel.variables[rxn.reverse_id].primal


def get_fluxes(tmodel):
    """
    Returns all the calculated fluxes from the pyTFA model
    """
    return pd.Series(
        {
            rxn.id: tmodel.variables[rxn.id].primal
            - tmodel.variables[rxn.reverse_id].primal
            for rxn in tmodel.reactions
        }
    )


def get_free_energy(tmodel, reaction_or_id):
    """
    Returns the delta G of the specified reaction.

    Parameters:
        tmodel (pyTFA model): Model of interest
        reaction_or_id (cobra.Reaction or String): Reaction of interest

    Returns:
        energy: Resulting Gibbs energy of the reaction
    """
    rxn_id = reaction_or_id if isinstance(reaction_or_id, str) else reaction_or_id.id
    var_id = DeltaG.prefix + rxn_id
    if not tmodel.variables.has_key(var_id):
        raise KeyError(rxn_id)
    return tmodel.variables[var_id].primal


def get_free_energies(tmodel):
    """
    Returns all the calculated delta G.
    """
    sr = tmodel.get_primal(DeltaG)
    sr.index = [name[len(DeltaG.prefix) :] for name in sr.index]
    return sr


def get_log_concentration(tmodel, metabolite_or_id):
    """
    Returns the log concentration of the specified metabolite.

    Parameters:
        tmodel (pyTFA model): Model of interest
        metabolite_or_id (cobra.Metabolite or String): Metabolite of interest

    Returns:
        log concentration
    """
    met_id = (
        metabolite_or_id if isinstance(metabolite_or_id, str) else metabolite_or_id.id
    )
    var_id = LogConcentration.prefix + met_id
    if not tmodel.variables.has_key(var_id):
        raise KeyError(met_id)
    return tmodel.variables[var_id].primal


def get_log_concentrations(tmodel):
    """
    Returns all the calculated log concentrations.
    """
    sr = tmodel.get_primal(LogConcentration)
    sr.index = [name[len(LogConcentration.prefix) :] for name in sr.index]
    return sr
