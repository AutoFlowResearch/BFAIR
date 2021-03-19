"""Database.

This module hosts functions to manipulate FIA MS mapping and struct files.
"""

__all__ = ["merge_mappings", "merge_structs", "read_mapping", "read_struct", "write_mapping", "write_struct"]

import csv
from collections import Iterable

import pandas as pd

MAPPING_FILE_COLUMNS = ["unused", "formula", "ids"]
STRUCT_FILE_COLUMNS = ["id", "formula", "smiles", "inchi"]


def _get_metadata(df, attr_name):
    # Extracts metadata from a dataframe and raises exception if not found.
    try:
        return df.attrs[attr_name]
    except KeyError:
        raise ValueError(f"Attribute '{attr_name}' missing.")


def merge_mappings(mapping_a, mapping_b):
    """
    Merges two mapping files.

    Parameters
    ----------
    mapping_a, mapping_b : pandas.DataFrame
        3-column dataframes representing mapping files. First column is unused. Chemical formulas and lists of
        metabolite IDs are stored in the remaining columns.
    """
    mapping_merged = mapping_a.merge(mapping_b.iloc[:, 1:], how="outer", on="formula")
    rows = []
    for left_list, right_list in zip(mapping_merged["ids_x"], mapping_merged["ids_y"]):
        row = []
        if isinstance(left_list, Iterable):
            row.extend(left_list)
        if isinstance(right_list, Iterable):
            row.extend([item for item in right_list if item not in row])
        rows.append(row)
    mapping_merged.drop(columns=["ids_x", "ids_y"], inplace=True)
    mapping_merged["unused"] = 0
    mapping_merged["ids"] = rows
    return mapping_merged


def merge_structs(struct_a, struct_b):
    """
    Merges two mapping files.

    Parameters
    ----------
    struct_a, struct_b : pandas.DataFrame
        4-column dataframes representing struct files. First and second columns correspond to metabolite ID and chemical
        formula, respectively.
    """
    return struct_a.append(struct_b, ignore_index=True)


def read_mapping(pathname):
    """
    Reads a FIA MS mapping file.

    Parameters
    ----------
    pathname : Path or str
        Pathname corresponding to a FIA MS mapping file in TSV format.

    Returns
    -------
    pandas.DataFrame
        A dataframe with 3-columns. First column is unused. Second column contains chemical formulas and last column
        stores a list of metabolite IDs that share a chemical formula.
    """
    database_name = ""
    database_version = ""
    rows = []
    with open(pathname, "r") as file:
        reader = csv.reader(file, delimiter="\t")
        database_name = next(reader)[1]
        database_version = next(reader)[1]
        for line in reader:
            # elements 2 onwards are metabolite IDs
            rows.append([*line[:2], line[2:]])
    mapping = pd.DataFrame(rows, columns=MAPPING_FILE_COLUMNS)
    mapping.attrs["database_name"] = database_name
    mapping.attrs["database_version"] = database_version
    return mapping


def read_struct(pathname):
    """
    Reads a FIA MS struct file.

    Parameters
    ----------
    pathname : Path or str
        Pathname corresponding to a FIA MS struct file in TSV format.

    Returns
    -------
    pandas.DataFrame
        A dataframe with 4-columns. First and second columns correspond to metabolite ID and chemical formula,
        respectively.
    """
    return pd.read_csv(pathname, names=STRUCT_FILE_COLUMNS, delimiter="\t")


def write_struct(struct, pathname):
    """
    Writes a struct file to disk.

    Parameters
    ----------
    struct : pandas.DataFrame
        A dataframe with 4-columns. First and second columns correspond to metabolite ID and chemical formula,
        respectively.
    """
    struct.to_csv(pathname, sep="\t", index=None, header=None)


def write_mapping(mapping, pathname, database_name=None, database_version=None):
    """
    Writes a mapping file to disk.

    Parameters
    ----------
    mapping : pandas.DataFrame
        A 3-column dataframe representing a mapping file. First column is unused. Chemical formulas and lists of
        metabolite IDs that are stored in the remaining columns.
    pathname : Path or str
        Destination pathname.
    database_name : str
        Name of the mapping file.
    database_version : str
        Version of the mapping file.

    Raises
    ------
    ValueError
        If `database_name` or `database_version` are not provided and cannot be found as metadata attributes of the
        input dataframe.
    """
    if database_name is None:
        database_name = _get_metadata(mapping, "database_name")
    if database_version is None:
        database_version = _get_metadata(mapping, "database_version")
    with open(pathname, "w", newline="", encoding="utf-8") as file:
        writer = csv.writer(file, delimiter="\t")
        writer.writerow(["database_name", database_name])
        writer.writerow(["database_version", database_version])
        for row in mapping.itertuples(index=False, name=None):
            writer.writerow([*row[:2], *row[2]])
