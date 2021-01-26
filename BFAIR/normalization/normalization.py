import pandas as pd
import numpy as np
import copy


def min_max_norm(
    df, columnname="Intensity", groupname_colname="sample_group_name"
):
    """
    Function that applies min/max scaling to the input dataframe; the data
    will range from 0 for the min- to 1 for the max value

    Parameters
    ----------
    df: pandas.DataFrame
        input dataframe, output of either the extractNamesAndIntensities()
        or the calculateMeanVarRSD() function
    columnname: string
        the name of the column with the data that needs to be
        normalized, defaults to 'Intensity'
    groupname_colname: string
        the name of the column with the sample group names,
        defaults to 'sample_group_name'

    Returns
    -------
    output_df: pandas.DataFrame
        the output dataframe. It follows the same architecture as
        the input dataframe, just with normalized values
    """
    output_df = pd.DataFrame()
    sample_group_names = df[groupname_colname].unique()
    for i, sample_group_name in enumerate(sample_group_names):
        new_df = copy.deepcopy(df[df[groupname_colname] == sample_group_name])
        min_val = min(new_df[columnname])
        max_val = max(new_df[columnname])
        new_df[columnname] = new_df[columnname].sub(min_val)
        new_df[columnname] = new_df[columnname].div((max_val - min_val))
        if i == 0:
            output_df = new_df
        else:
            output_df = output_df.append(new_df)
    return output_df


def tsi_norm(
    df, columnname="Intensity", groupname_colname="sample_group_name"
):
    """
    Applies Total Sum Intensity normalization; all values will be divided by
    the sum of all values. Like that, summing all the values up would
    amout to '1'

    Parameters
    ----------
    df: pandas.DataFrame
        input dataframe, output of either the extractNamesAndIntensities()
        or the calculateMeanVarRSD() function
    columnname: string
        the name of the column with the data that needs to be
        normalized, defaults to 'Intensity'
    groupname_colname: string
        the name of the column with the sample group names,
        defaults to 'sample_group_name'

    Returns
    -------
    output_df: pandas.DataFrame
        the output dataframe. It follows the same architecture as
        the input dataframe, just with normalized values
    """
    output_df = pd.DataFrame()
    sample_group_names = df[groupname_colname].unique()
    for i, sample_group_name in enumerate(sample_group_names):
        new_df = copy.deepcopy(df[df[groupname_colname] == sample_group_name])
        tsi = new_df[columnname].sum(skipna=True)
        new_df[columnname] = new_df[columnname].div(tsi)
        if i == 0:
            output_df = new_df
        else:
            output_df = output_df.append(new_df)
    return output_df


def lim_tsi_norm(
    metabolite_input,
    df,
    lim_type,
    product_df=None,
    biomass_value=None,
    columnname="Intensity",
    groupname_colname="sample_group_name",
):
    """
    Applies a modified version of Total Sum Intensity normalization; all
    values will be divided by the sum of all all the values of metabolites
    that are a) amino acids, b) part of the biomass function of a model of the
    organism of interest or c) are part of the biomass function of a model of
    the organism of interest weighted by their multipliers in said model.
    Like that, summing all the values for the metabolites
    contributing to the biomass function up would amout to '1'
    
    Parameters
    ----------
    metabolite_input: pandas.DataFrame or list
        either a dataframe with one column denoting the
        metabolites that are part of the biomass function and one
        column with their corresponding multipliers or a list
        with the metabolites of interest (e.g. amino acids)
    df: pandas.DataFrame
        input dataframe, output of either the extractNamesAndIntensities()
        or the calculateMeanVarRSD() function
    lim_type: string 
        which type of notmalization should be run, based on a) amino acids,
        b) metabolites of the biomass function or c) metabolites of the
        biomass function weighted by their multipliers
    product_df: pandas.DataFrame
        a dataframe with one column denoting the
        metabolites that are on the product side of the biomass
        function and one column with their corresponding multipliers.
        Defaults to 'None'
    biomass_value: float
        the multiplier of biomass on the product side of the
        biomass function in the corresponding metabolic model.
        Defaults to 'None'
    columnname: string
        the name of the column with the data that needs to be
        normalized, defaults to 'Intensity'.
        Defaults to 'Intensity'
    groupname_colname: string
        the name of the column with the sample group names,
        defaults to 'sample_group_name'.
        Defaults to "sample_group_name"

    Returns
    -------
    output_df: pandas.DataFrame
        the output dataframe. It follows the same architecture as
        the input dataframe, just with normalized values

    Raises
    ------
    ValueError
        Either if the wrong 'lim_type' is used or 'product_df' or
        'biomass_value' is missing
    """
    lim_tsi = 0
    output_df = pd.DataFrame()
    sample_group_names = df[groupname_colname].unique()
    for i, sample_group_name in enumerate(sample_group_names):
        new_df = copy.deepcopy(df[df[groupname_colname] == sample_group_name])
        if lim_type == 'amino_acid':
            for amino_acid in metabolite_input:
                met_tsi = sum(df[df["Metabolite"] == amino_acid][columnname])
                lim_tsi += met_tsi
        elif lim_type == 'biomass':
            for biomass_met in metabolite_input["Metabolite"]:
                met_tsi = sum(
                    new_df[new_df["Metabolite"] == biomass_met][columnname]
                )
                lim_tsi += met_tsi
        elif lim_type == 'bm_function':
            if product_df is None:
                raise ValueError("'product_df' is missing!")
            if biomass_value is None:
                raise ValueError("'biomass_value' is missing!")
            for cnt, biomass_met in enumerate(metabolite_input["Metabolite"]):
                met_tsi = sum(
                    new_df[new_df["Metabolite"] == biomass_met][columnname]
                )
                norm_met_tsi = met_tsi * (
                    metabolite_input["Value"][cnt] / biomass_value
                )
                lim_tsi += norm_met_tsi
            for cnt_p, biomass_mets_product in enumerate(
                product_df["Metabolite"]
            ):
                met_tsi = sum(
                    new_df[new_df["Metabolite"] == biomass_mets_product][
                        columnname
                    ]
                )
                norm_met_tsi = met_tsi * (
                    product_df["Value"][cnt_p] / biomass_value
                )
                lim_tsi -= norm_met_tsi
        else:
            raise ValueError(
                "'lim_type' must be either 'amino_acid', "
                "'biomass' or 'bm_function', not '" + lim_type
                + "'"
            )
        norm_vals = new_df[columnname].div(lim_tsi)
        new_df[columnname] = norm_vals
        if i == 0:
            output_df = new_df
        else:
            output_df = output_df.append(new_df)
    return output_df


def pqn_norm(
    df,
    groupname_colname="sample_group_name",
    value_colname="Intensity",
    corr_type="median",
    qc_vector=None,
):
    """
    Probabilistic Quotient Normalization: This method adjusts for dilutions.
    This is a modified version of Total Sum Intensity normalization

    Parameters
    ----------
    df: pandas.DataFrame
        input dataframe, output of either the extractNamesAndIntensities()
        or the calculateMeanVarRSD() function
    groupname_colname: string
        the name of the column with the sample group names,
        defaults to 'sample_group_name'
    value_colname: string
        the name of the column with the data that needs to be
        normalized, defaults to 'Intensity'
    corr_type: string
        type of midpoint finding, 'median' or 'mean'
    qc_vector: list
        an optional QC vector that can be provided

    Returns
    -------
    output_df: pandas.DataFrame
        the output dataframe. It follows the same architecture as
        the input dataframe, just with normalized values
    """
    # 0)
    grouped_df = pd.DataFrame()
    grouped_df["Metabolite"] = df["Metabolite"].unique()
    sample_group_names = df[groupname_colname].unique()
    for sample_group_name in sample_group_names:
        group_df = df[df[groupname_colname] == sample_group_name]
        new_col = []
        for metabolite in list(grouped_df["Metabolite"]):
            if metabolite in list(group_df["Metabolite"]):
                new_col.append(
                    float(
                        group_df[group_df["Metabolite"] == metabolite][
                            value_colname
                        ]
                    )
                )
            else:
                new_col.append(None)
        grouped_df[sample_group_name] = new_col
    # 1)
    tsi_df = copy.deepcopy(grouped_df)
    for colname in grouped_df.columns[1:]:
        new_df = copy.deepcopy(grouped_df)
        tsi = grouped_df[colname].sum(skipna=True)
        new_col = new_df[colname].div(tsi)
        tsi_df[colname] = new_col
    # 2)
    if qc_vector is not None:
        qc_vector = qc_vector
    else:
        qc_vector = []
        for x, row in tsi_df.iterrows():
            if corr_type == "median":
                qc_vector.append(row[1:].median(skipna=True))
            elif corr_type == "mean":
                qc_vector.append(row[1:].mean(skipna=True))
    # 3)
    dilution_factor_list = []
    for colname in tsi_df.columns[1:]:
        new_col = np.array(tsi_df[colname])
        new_col_div = list(new_col / np.array(qc_vector))
        # 4)
        if corr_type == "median":
            dilution_factor = pd.Series(new_col_div).median(skipna=True)
        elif corr_type == "mean":
            dilution_factor = pd.Series(new_col_div).mean(skipna=True)
        dilution_factor_list.append(dilution_factor)
    # print(dilution_factor_list)
    # 5)
    norm_df = copy.deepcopy(tsi_df)
    for position, colname in enumerate(tsi_df.columns[1:]):
        norm_col = tsi_df[colname] * dilution_factor_list[position]
        norm_df[colname] = norm_col
    # 6)
    output_df = copy.deepcopy(df)
    for x, row in norm_df.iterrows():
        for sample_num, sample in enumerate(row.index[1:]):
            if not pd.isna(row[sample_num + 1]):
                orig_df_reduced = df[df[groupname_colname] == sample]
                index_pos = orig_df_reduced[
                    orig_df_reduced["Metabolite"] == row[0]
                ].index[0]
                output_df.at[index_pos, value_colname] = row[sample_num + 1]
    return output_df
