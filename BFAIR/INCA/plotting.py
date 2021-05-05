import numpy as np
import math
import pandas as pd
import scipy.stats as stats
import matplotlib.pyplot as plt
import seaborn as sns


def sampled_fluxes_minrange(sampled_fluxes, min_val, max_val):
    def in_between(column):
        return not all(column.between(min_val, max_val, inclusive=False))
    criteria = list(sampled_fluxes.apply(in_between, axis=0, result_type='reduce'))
    return sampled_fluxes[sampled_fluxes.columns[criteria]]


def show_reactions(reactions):
    for cnt, rxn in enumerate(reactions):
        print(f'{cnt}: {rxn}')


def _sampled_reaction_fit(sampled_value):
    sampled_reaction = sorted(sampled_value)
    fit = stats.norm.pdf(
        sampled_reaction,
        np.mean(sampled_reaction),
        np.std(sampled_reaction),
    )
    return sampled_reaction, fit

        
def plot_sampled_reaction_fluxes(sampled_fluxes, reactions, reaction_id, bins=10):
    sampled_value = sampled_fluxes[reactions[reaction_id]]
    if all(v == 0 for v in sampled_value):
        print(f'All sample values are 0 for {reactions[reaction_id]}')
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
        
        
def _subplot_sampled_reaction_fluxes(sampled_fluxes, reactions, reaction_id, bins=10):
    sampled_value = sampled_fluxes[reactions[reaction_id]]
    if all(v == 0 for v in sampled_value):
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
    for i in range(0, len(reactions)):
        plt.subplot(rows,3,i+1)
        _subplot_sampled_reaction_fluxes(sampled_fluxes, reactions, reaction_id=i, bins=bins)
    plt.tight_layout()
    fig.text(0.5, -0.02, 'Flux [mmol * gDCW$^{-1}$ * h$^{-1}$]', ha='center', fontsize=20)
    fig.text(-0.02, 0.5, 'Distribution [%]', va='center', rotation='vertical', fontsize=20)
    return fig          
        
        
def show_subsystems(model):
    subsystems = []
    for rxn in model.reactions:
        subsystem = rxn.subsystem
        if subsystem not in subsystems:
            subsystems.append(subsystem)
    for cnt, sub in enumerate(subsystems):
        print(f'{cnt}: {sub}')        
        
        
def get_subsytem_reactions(model, sampled_fluxes, subsystem_id):
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
    subsystem_df = sampled_fluxes[sampled_fluxes.columns[reactions_mask]]
    re_arranged_df = subsystem_df.stack().reset_index(level=[1], name='Sampled Fluxes')
    re_arranged_df = re_arranged_df.rename(columns = {'level_1':'Reaction'}).reset_index(drop=True)
    return re_arranged_df


def plot_subsystem_fluxes(model, sampled_fluxes, subsystem_id, no_zero_cols=False):
    if no_zero_cols:
        sampled_fluxes = sampled_fluxes.loc[:, (sampled_fluxes != 0).any(axis=0)]
    reactions, subsystems = get_subsytem_reactions(model, sampled_fluxes, subsystem_id)
    re_arranged_df = _reduce_sampled_fluxes(sampled_fluxes, reactions)
    sns.boxplot(x="Reaction", y="Sampled Fluxes", data=re_arranged_df, orient = 'v')
    plt.xticks(rotation=70)
    plt.title(subsystems[subsystem_id], size = 20)

    
# Flux split

def _prepare_input_fluxes(sampled_fluxes, flux):
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

        
def calculate_split_ratio(sampled_fluxes, influx, split_flux1, split_flux2=None, branch_point_name="Branch point"):
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


def plot_split_ratio(sampled_fluxes, influx, split_flux1, split_flux2=None, branch_point_name="Branch point"):
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
    plt.title(f'{branch_point_name}:\n{plot_title}', size = 15)