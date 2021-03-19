"""Struct.

Contains basic functions to read, merge, and write FIA MS struct files.
"""

__all__ = ["merge", "read", "write"]

from collections import namedtuple

import pandas as pd

FILE_COLUMNS = ["id", "formula", "unused_smiles", "unused_inchi"]

Row = namedtuple("StructRow", FILE_COLUMNS)


def merge(struct_a, struct_b):
    """
    Merges two mapping files.

    Parameters
    ----------
    struct_a, struct_b : pandas.DataFrame
        Dataframes containing 4 columns of metabolite ID, formula, SMILES, and InChI.

    Returns
    -------
    pandas.DataFrame
    """
    return struct_a.append(struct_b, ignore_index=True)


def read(pathname):
    """
    Reads a FIA MS struct file.

    Parameters
    ----------
    pathname : Path or str
        Pathname corresponding to a FIA MS struct file in TSV format.

    Returns
    -------
    pandas.DataFrame
        Dataframe containing 4 columns of metabolite ID, formula, SMILES, and InChI.
    """
    return pd.read_csv(pathname, names=FILE_COLUMNS, delimiter="\t")


def write(struct, pathname):
    """
    Writes a struct file to disk.

    Parameters
    ----------
    struct : pandas.DataFrame
        Dataframe containing 4 columns of metabolite ID, formula, SMILES, and InChI.
    """
    struct.to_csv(pathname, sep="\t", index=None, header=None)
