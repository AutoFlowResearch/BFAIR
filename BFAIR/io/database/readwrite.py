"""Read/Write.

Contains basic functions to read, merge, and write FIA MS database files.
"""

__all__ = ["merge_mappings", "merge_structs", "read_mapping", "read_struct", "write_mapping", "write_struct"]

import csv
from collections import Iterable, namedtuple

import pandas as pd

MAPPING_FILE_COLUMNS = ["unused_mass", "formula", "ids"]
STRUCT_FILE_COLUMNS = ["id", "formula", "unused_smiles", "unused_inchi"]

StructRow = namedtuple("StructRow", STRUCT_FILE_COLUMNS)


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
        Dataframes containing 3 columns of mass, formula, and metabolite IDs.

    Returns
    -------
    pandas.DataFrame
    """
    mapping_merged = mapping_a.merge(mapping_b.iloc[:, 1:], how="outer", on=MAPPING_FILE_COLUMNS[1])
    rows = []
    dup_colname = MAPPING_FILE_COLUMNS[2]
    for left_list, right_list in zip(
        mapping_merged[dup_colname + "_x"], mapping_merged[dup_colname + "_y"]
    ):
        row = []
        if isinstance(left_list, Iterable):
            row.extend(left_list)
        if isinstance(right_list, Iterable):
            row.extend([item for item in right_list if item not in row])
        rows.append(row)
    mapping_merged[MAPPING_FILE_COLUMNS[0]] = 0
    mapping_merged[dup_colname] = rows
    return mapping_merged[MAPPING_FILE_COLUMNS]


def merge_structs(struct_a, struct_b):
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
        Dataframe containing 3 columns of mass, formula, and metabolite IDs.
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
        Dataframe containing 4 columns of metabolite ID, formula, SMILES, and InChI.
    """
    return pd.read_csv(pathname, names=STRUCT_FILE_COLUMNS, delimiter="\t")


def write_struct(struct, pathname):
    """
    Writes a struct file to disk.

    Parameters
    ----------
    struct : pandas.DataFrame
        Dataframe containing 4 columns of metabolite ID, formula, SMILES, and InChI.
    """
    struct.to_csv(pathname, sep="\t", index=None, header=None)


def write_mapping(mapping, pathname, database_name=None, database_version=None):
    """
    Writes a mapping file to disk.

    Parameters
    ----------
    mapping : pandas.DataFrame
        Dataframe containing 3 columns of mass, formula, and metabolite IDs.
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
