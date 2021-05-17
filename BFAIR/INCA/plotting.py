import numpy as np
import math
import pandas as pd
import scipy.stats as stats
import matplotlib.pyplot as plt
import seaborn as sns


def sampled_fluxes_minrange(sampled_fluxes, min_val, max_val):
    """
    Sometimes the sampled fluxes include very small predicted fluxes
    for a lot of different reactions. They might be relevant but can
    easily overwhealm plots. This function reduces the subset of
    fluxes regarded for plots.

    Parameters
    ----------
    sampled_fluxes : pandas.Dataframe
        The calculated fluxes, output of sampling. For each reaction,
        n fluxes will be calculated (n = number of samples taken).
    min_val : float, int
        Min value for excluded subset.
    max_val : float, int
        Max value for excluded subset.

    Returns
    -------
    sampled_fluxes : pandas.Dataframe
        Reduced subset of the input.
    """

    def in_between(column):
        return not all(column.between(min_val, max_val, inclusive=False))
    criteria = list(sampled_fluxes.apply(in_between, axis=0, result_type='reduce'))
    return sampled_fluxes[sampled_fluxes.columns[criteria]]


def show_reactions(reactions):
    """
    A function to inspect the reactions in a given subset,
    e.g. a subsytem or pathway.

    Parameters
    ----------
    reactions : list
        List of reactions, e.g. output of `get_subsytem_reactions()`.
    """
    for cnt, rxn in enumerate(reactions):
        print(f'{cnt}: {rxn}')


def _sampled_reaction_fit(sampled_values, alpha=1e-3):
    """
    Helper function to prepare the input for the sampled reaction
    fluxes plots.

    Parameters
    ----------
    sampled_values : pandas.Series
        Sampled values for one reaction.

    sampled_value = sampled_fluxes[reactions[reaction_id]]
    if all(v == 0 for v in sampled_value):
        print(f'All sample values are 0 for {reactions[reaction_id]}')
    Returns
    -------
    sampled_reaction : list
        List of the values sampled for a given reaction.
        Sorted ascendingly.
    fit : list
        Normal distribution fit (y-Values) over the sampled
        data.
    """
    else:
        sampled_reaction, fit = _sampled_reaction_fit(
            sampled_value,
        )
        fig = plt.plot(sampled_reaction, fit,'-o', color='darksalmon')
        plt.hist(sampled_reaction, bins=bins, 
                 rwidth=.95, color='darkcyan', density=True)
        plt.ylabel('Distribution [%]')
        plt.xlabel('Flux [mmol * gDCW$^{-1}$ * h$^{-1}$]')
        plt.title(reactions[reaction_id], size = 20)
        return fig        
        
        
    sampled_value = sampled_fluxes[reactions[reaction_id]]
    if all(v == 0 for v in sampled_value):
def plot_sampled_reaction_fluxes(
    sampled_fluxes, reactions, reaction_id, bins=10, alpha=1e-3
):
    """
    Plots the distribution of the sampled values for one selected
    reaction as a histogram and a normal distribution.

    Parameters
    ----------
    sampled_fluxes : pandas.Dataframe
        The calculated fluxes, output of sampling. For each reaction,
        n fluxes will be calculated (n = number of samples taken).
    reactions : list
        List of reactions, e.g. output of `get_subsytem_reactions()`.
    reaction_id : int
        ID of the reaction in the `reactions` list.
    bins : int
        Number of bins used for the histogram. Defaults to 10.

    Returns
    -------
    fig : matplotlib.Figure
        Figure of sampled value distribution for one reaction.
    """
def _subplot_sampled_reaction_fluxes(
    sampled_fluxes, reactions, reaction_id, bins=10, alpha=1e-3
):
    """
    Produces a subplot of the distribution of the sampled values
    for one selected reaction as a histogram and a normal distribution.

    Parameters
    ----------
    sampled_fluxes : pandas.Dataframe
        The calculated fluxes, output of sampling. For each reaction,
        n fluxes will be calculated (n = number of samples taken).
    reactions : list
        List of reactions, e.g. output of `get_subsytem_reactions()`.
    reaction_id : int
        ID of the reaction in the `reactions` list.
    bins : int
        Number of bins used for the histogram. Defaults to 10.

    Returns
    -------
    fig : matplotlib.Figure
        Figure of sampled value distribution for one reaction.
    """
        fig = plt.plot(figsize=(15, 20))
        plt.axis('off')
        plt.text(0, -0.02, f'All sampled\nvalues are 0\nfor {reactions[reaction_id]}', ha='center', fontsize=20)
    else:
        sampled_reaction, fit = _sampled_reaction_fit(
            sampled_value,
        )
        fig = plt.plot(sampled_reaction, fit,'-o', color='darksalmon')
        plt.hist(sampled_reaction, bins=bins, 
                 rwidth=.95, color='darkcyan', density=True)
        plt.title(reactions[reaction_id], size = 20)
        return fig        
        
        
def plot_all_subsystem_fluxes(sampled_fluxes, reactions, bins=10):
    rows = math.ceil(len(reactions)/3)
    fig, ax = plt.subplots(figsize=(15, 5 * (rows/2)))
    """
    Plots the distribution of the sampled values for all
    reactions in a selected subsytem as a histogram and a
    normal distribution.

    Parameters
    ----------
    sampled_fluxes : pandas.Dataframe
        The calculated fluxes, output of sampling. For each reaction,
        n fluxes will be calculated (n = number of samples taken).
    reactions : list
        List of reactions, e.g. output of `get_subsytem_reactions()`.
    bins : int
        Number of bins used for the histogram. Defaults to 10.

    Returns
    -------
    fig : matplotlib.Figure
        Figure of sampled value distribution for one reaction.
    """
    for i in range(0, len(reactions)):
        plt.subplot(rows,3,i+1)
        _subplot_sampled_reaction_fluxes(sampled_fluxes, reactions, reaction_id=i, bins=bins)
    plt.tight_layout()
    fig.text(0.5, -0.02, 'Flux [mmol * gDCW$^{-1}$ * h$^{-1}$]', ha='center', fontsize=20)
    fig.text(-0.02, 0.5, 'Distribution [%]', va='center', rotation='vertical', fontsize=20)
    return fig          
        
        
def show_subsystems(model):
    """
    Prints the names of all the subsystems in the model with
    their IDs.

    Parameters
    ----------
    model : cobra.Model
        Metabolic model.
    """
    subsystems = []
    for rxn in model.reactions:
        subsystem = rxn.subsystem
        if subsystem not in subsystems:
            subsystems.append(subsystem)
    for cnt, sub in enumerate(subsystems):
        print(f'{cnt}: {sub}')        
        
        
def get_subsytem_reactions(model, subsystem_id):
    """
    Extracts the names of the reactions in a selected subsytem. Also
    outputs a list of all the subsytems in the model.

    Parameters
    ----------
    model : cobra.Model
        Metabolic model.
    subsystem_id : int
        ID of the subsystem the reactions should be extracted from.
        Can be selected after inspection using `show_subsystems()`.

    Returns
    -------
    reactions : list
        List of reactions in a selected subsytem.
    subsystems : list
        List of subsystems in a model. Input for `plot_subsystem_fluxes()`.
    """
    subsystems = []
    subsystems_reactions = {}
    for rxn in model.reactions:
        subsystem = rxn.subsystem
        if subsystem not in subsystems:
            subsystems.append(subsystem)
        subsystems_reactions[rxn.id] = subsystem
    reactions = [
        value[0]
        for value in subsystems_reactions.items()
        if value[1] == subsystems[subsystem_id]
    ]
    return reactions, subsystems


def _reduce_sampled_fluxes(sampled_fluxes, reactions):
    reactions_mask = [True if col in reactions else False for col in sampled_fluxes.columns]
    """
    Reduces the input dataframe of sampled fluxes to a subset of
    selected reactions.

    Parameters
    ----------
    sampled_fluxes : pandas.Dataframe
        The calculated fluxes, output of sampling. For each reaction,
        n fluxes will be calculated (n = number of samples taken).
    reactions : list
        List of reactions, e.g. output of `get_subsytem_reactions()`.

    Returns
    -------
    re_arranged_df : pandas.Dataframe
        Sampled fluxes reduced to a specified set of reactions.
    """
    subsystem_df = sampled_fluxes[sampled_fluxes.columns[reactions_mask]]
    re_arranged_df = subsystem_df.stack().reset_index(level=[1], name='Sampled Fluxes')
    re_arranged_df = re_arranged_df.rename(columns = {'level_1':'Reaction'}).reset_index(drop=True)
    return re_arranged_df


def plot_subsystem_fluxes(
    model, sampled_fluxes, subsystem_id, no_zero_cols=False
):
    """
    Visualizes the distribution of sampled fluxes in a selected
    subsytem.

    Parameters
    ----------
    model : cobra.Model
        Metabolic model.
    sampled_fluxes : pandas.Dataframe
        The calculated fluxes, output of sampling. For each reaction,
        n fluxes will be calculated (n = number of samples taken).
    subsystem_id : int
        ID of the subsystem the reactions should be extracted from.
        Can be selected after inspection using `show_subsystems()`.
    no_zero_cols : boolean
        Indicates if columns containing only 0s should be excluded.

    Returns
    -------
    fig : matplotlib.Figure
        Figure of sampled value distribution for one reaction.
    """
    fig = plt.plot(figsize=(15, 20))
    if no_zero_cols:
        sampled_fluxes = sampled_fluxes.loc[:, (sampled_fluxes != 0).any(axis=0)]
    reactions, subsystems = get_subsytem_reactions(model, sampled_fluxes, subsystem_id)
    re_arranged_df = _reduce_sampled_fluxes(sampled_fluxes, reactions)
    sns.boxplot(x="Reaction", y="Sampled Fluxes", data=re_arranged_df, orient = 'v')
    plt.xticks(rotation=70)
    plt.title(subsystems[subsystem_id], size = 20)

# Flux split


def _prepare_input_fluxes(sampled_fluxes, flux):
    """
    Makes sure that the provided fluxes are in the right format
    for other functions. Input can either be a flux id or a
    cobbination of fluxes provided as a Series.

    Parameters
    ----------
    sampled_fluxes : pandas.Dataframe
        The calculated fluxes, output of sampling. For each reaction,
        n fluxes will be calculated (n = number of samples taken).
    flux : pandas.Series or str
        Either a selected column (or combination of fluxes) or the
        name of a flux.

    Returns
    -------
    flux : pandas.Series
        Selected flux in the relevant format.
    """
    if isinstance(flux, pd.Series):
        return flux
    elif isinstance(flux, str):
        flux = sampled_fluxes[flux]
        return flux
    else:
        raise TypeError(
            f"The input is a '{type(flux)}', this type of object"
            " cannot be used here, please provide a pd.Series or"
            "a column identifier (string)")
def calculate_split_ratio(
    sampled_fluxes,
    influx,
    split_flux1,
    split_flux2=None,
    branch_point_name="Branch point",
):
    """
    Calculates the ratios of the in- and effluxes at a split point.

    Parameters
    ----------
    sampled_fluxes : pandas.Dataframe
        The calculated fluxes, output of sampling. For each reaction,
        n fluxes will be calculated (n = number of samples taken).
    influx : string or pandas.Series
        Name of the reaction flowing to the branch point or pre-processed
        combined flux coming to the branch point.
    split_flux1 : string or pandas.Series
        Name of flux #1 coming out of the branch point or pre-processed
        combined flux coming out of the branch point.
    split_flux2 : string or pandas.Series
        Name of flux #2 coming out of the branch point or pre-processed
        combined flux coming out of the branch point.
        Defaults to `None`.
    branch_point_name : string
        name of the branch point to be added to the output.
        Defaults to "Branch point".

    Returns
    -------
    output_df : pandas.Dataframe
        Presents the split flux ratio.
    """
    influx = _prepare_input_fluxes(sampled_fluxes, influx)
    split_flux1 = _prepare_input_fluxes(sampled_fluxes, split_flux1)
    mean_list = [np.abs(np.mean(split_flux1 / influx))]
    stdev_list = [np.std(split_flux1 / influx)]
    index_list = [influx.name + '/' + split_flux1.name]
    if split_flux2 is not None:
        split_flux2 = _prepare_input_fluxes(sampled_fluxes, split_flux2)
        mean_list.append(np.abs(np.mean(split_flux2 / influx)))
        stdev_list.append(np.std(split_flux2 / influx))
        index_list.append(influx.name + '/' + split_flux2.name)
    fluxes = {'Mean': mean_list,
        'Stdev': stdev_list,
        }    
    output_df = pd.DataFrame(fluxes, columns = ['Mean','Stdev'], index = [index_list])
    return output_df


def plot_split_ratio(
    sampled_fluxes,
    influx,
    split_flux1,
    split_flux2=None,
    branch_point_name="Branch point",
):
    """
    Visualization of the in- and effluxes at a split point.

    Parameters
    ----------
    sampled_fluxes : pandas.Dataframe
        The calculated fluxes, output of sampling. For each reaction,
        n fluxes will be calculated (n = number of samples taken).
    influx : string or pandas.Series
        Name of the reaction flowing to the branch point or pre-processed
        combined flux coming to the branch point.
    split_flux1 : string or pandas.Series
        Name of flux #1 coming out of the branch point or pre-processed
        combined flux coming out of the branch point.
    split_flux2 : string or pandas.Series
        Name of flux #2 coming out of the branch point or pre-processed
        combined flux coming out of the branch point.
        Defaults to `None`.
    branch_point_name : string
        name of the branch point to be added to the output.
        Defaults to "Branch point".

    Returns
    -------
    fig : matplotlib.Figure
        Figure of the sampled value distribution of the in-
        and effluxes at a split point.
    """
    influx = _prepare_input_fluxes(sampled_fluxes, influx)
    split_flux1 = _prepare_input_fluxes(sampled_fluxes, split_flux1)
    if split_flux2 is not None:
        split_flux2 = _prepare_input_fluxes(sampled_fluxes, split_flux2)
        sum_ratios = abs(split_flux1 / influx) + abs(split_flux2 / influx)
        split_df = pd.DataFrame({
            'Sum': sum_ratios,
            split_flux1.name: abs(split_flux1 / influx),
            split_flux2.name: abs(split_flux2 / influx),
        })
        plot_title = f'{influx.name} / {split_flux1.name} vs.\n{influx.name} / {split_flux2.name}'
    else:
        split_df = pd.DataFrame({
            split_flux1.name: abs(split_flux1 / influx),
        })
        plot_title = f'{influx.name} / {split_flux1.name}'
    re_arranged_df = _reduce_sampled_fluxes(split_df, split_df.columns)
    fig = sns.boxplot(x="Reaction", y="Sampled Fluxes", data=re_arranged_df, orient = 'v')
    plt.xticks(rotation=70)
    plt.title(f'{branch_point_name}:\n{plot_title}', size = 15)    """
    Checks which calculated fluxes qualify as being 'observable'
    (flux value at least 4 times the confidence interval and 0
    not included in the confidence interval).

    Parameters
    ----------
    fittedFluxes : pandas.DataFrame
        Dataframe (reimported output of an INCA simulation)
        that contains the confidence intervals predicted for
        the model.

    Returns
    -------
    observable_fluxes : pandas.DataFrame
        Input dataframe reduced to the fluxes that qualify as
        observable fluxes.
    """
    big_enough_fluxes = fittedFluxes[
    """
    Calculates the percentage of how many of the fluxes in
    the input dataframe qualify as observable fluxes.

    Parameters
    ----------
    fittedFluxes : pandas.DataFrame
        Dataframe (reimported output of an INCA simulation)
        that contains the confidence intervals predicted for
        the model.

    Returns
    -------
    percentage of observable fluxes : string
        relative to total number of fluxes.
    """
    percent_float = (
    """
    Calculates the mean standard deviation of the observable
    fluxes.

    Parameters
    ----------
    fittedFluxes : pandas.DataFrame
        Dataframe (reimported output of an INCA simulation)
        that contains the confidence intervals predicted for
        the model.

    Returns
    -------
    flux precision : float
        Mean of stdev of observable fluxes.
    """
    return get_observable_fluxes(fittedFluxes)["flux_stdev"].mean()
