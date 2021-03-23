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
    the sum of all values. Like that, summing up all the values would
    amount to '1'

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
    tsi = df[[groupname_colname, columnname]].groupby(groupname_colname).transform(np.sum)
    output_df = df.copy()
    output_df[columnname] /= tsi[columnname]
    return output_df


def lim_tsi_norm(
    metabolite_input,
    df,
    biomass_value=None,
    columnname="Intensity",
    groupname_colname="sample_group_name",
):
    """
    Applies a modified version of Total Sum Intensity normalization; all
    values will be divided by the sum of all the values of metabolites
    that are a) amino acids, b) part of the biomass function of a model of the
    organism of interest or c) are part of the biomass function of a model of
    the organism of interest weighted by their multipliers in said model.
    Like that, summing up all the values for the metabolites
    contributing to the biomass function would amount to '1'

    Parameters
    ----------
    metabolite_input: list, pandas.Series, np.ndarray or
        pandas.DataFrame
        either an object listing the metabolites of interest, e.g
        metabolites that are part of the biomass function, or a
        dataframe with the metabolties in one columns and another
        column with their corresponding multipliers
    df: pandas.DataFrame
        input dataframe, output of either the extractNamesAndIntensities()
        or the calculateMeanVarRSD() function
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
        Alerts user if 'biomass_value' is needed but missing
    ValueError
        Wrong type of input
    """
    lim_tsi = 0
    output_df = pd.DataFrame()
    sample_group_names = df[groupname_colname].unique()
    for i, sample_group_name in enumerate(sample_group_names):
        new_df = copy.deepcopy(df[df[groupname_colname] == sample_group_name])
        if isinstance(metabolite_input,
                      (list, pd.core.series.Series, np.ndarray)):
            for metabolite in metabolite_input:
                met_tsi = sum(new_df[new_df["Metabolite"] == metabolite][columnname])
                lim_tsi += met_tsi
        elif type(metabolite_input) == pd.core.frame.DataFrame:

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
        else:
            raise ValueError(
                "Wrong type of input! Input must be either "
                "'list', 'pd.core.series.Series', 'np.ndarray' or "
                "'pd.core.frame.DataFrame', not '"
                + type(metabolite_input) + "'"
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
        defaults to ‘sample_group_name’
    value_colname: string
        the name of the column with the data that needs to be
        normalized, defaults to ‘Intensity’
    corr_type: string
        type of midpoint determination, ‘median’ or ‘mean’
    qc_vector: list
        an optional QC vector that can be provided
    Returns
    -------
    output_df: pandas.DataFrame
        the output dataframe. It follows the same architecture as
        the input dataframe, just with normalized values
    """
    # TODO: handle missing metabolites in provided qc vector
    # 0, pivot table
    pivot_df = df.pivot(index="Metabolite", columns=groupname_colname, values=value_colname)
    # 1, divide each column by its total intensity
    pivot_df = pivot_df.div(pivot_df.sum(axis=0), axis=1)
    # 2, calculate QC vector as mean/median per row (i.e., metabolite)
    qc_df = qc_vector if qc_vector is not None else getattr(pivot_df, corr_type)(axis=1, skipna=True)
    # 3, normalize TSI by diving each row by corresponding value in QC vector
    # 4, calculate dilution factor as mean/median per column (i.e., sample group name)
    dilution = getattr(pivot_df.div(qc_df, axis=0), corr_type)(axis=0, skipna=True)
    # 5, multiply each column by its dilution factor
    pivot_df = pivot_df.multiply(dilution, axis=1)
    # 6, output
    pqn = (
        pivot_df.stack()  # back to original shape
        .reset_index()  # reset index, making a "0" column with intensities
        .merge(df.drop(columns=value_colname), how="right")  # add missing data
        .rename(columns={0: value_colname})  # renaming the "0" column
    )[df.columns]  # match original column order
    return pqn
