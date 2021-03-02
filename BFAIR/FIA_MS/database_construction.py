import numpy as np
import pandas as pd
from cobra.core.formula import Formula

__version__ = "0.0.1"


INVALID_FORMULA_STR = ["(", "Generic", "R", "X"]


def print_formula(elements):
    """
    Parameters
    ----------
    elements: dict
        The elements that form the metabolite and their corresponding amount

    Returns
    -------
    formula: str
        The formula of the metabolite
    """
    formula = "".join([f"{k}{int(v)}" for k, v in elements.items()])
    return formula


def zero_charge(metabolite):
    """
    Parameters
    ----------
    metabolite: cobra.Metabolite
        A metabolite in the cobra format including additional information
        (formula, charge, elements)

    Returns
    -------
    formula: str
        The formula of the metabolite
    """
    formula = Formula(metabolite.formula)
    if metabolite.charge != 0:
        if "desulfurated" in metabolite.name:
            formula.elements["S"] += metabolite.charge / 2
            formula = Formula(print_formula(formula.elements))
        else:
            if "H" in formula.elements and not np.isnan(metabolite.charge):
                formula.elements["H"] += -metabolite.charge
                formula = Formula(print_formula(formula.elements))
    return formula


def is_valid(metabolite):
    """
    Parameters
    ----------
    metabolite: cobra.Metabolite
        A metabolite in the cobra format including additional information
        (formula, charge, elements)

    Returns
    -------
    boolean
        Valid or invalid metabolite
    """
    if not metabolite.formula:
        return False
    for string in INVALID_FORMULA_STR:
        if string in metabolite.formula:
            return False
    return True


def make_struct(formulas, ids):
    """
    Parameters
    ----------
    formulas: list
        List of the formulas of the metabolites in the model
    ids: list
        List if the ids of the metabolites in the model

    Returns
    -------
    df_formulas: pandas.DataFrame
        Dataframe with information about the metabolites in the model.
        The masses are all set to 0
    """
    masses = [0] * len(formulas)
    df_formulas = pd.DataFrame(
        {"mass": masses, "formula": formulas, "id": ids}
    )
    return df_formulas


def make_mapping(df_formulas):
    """
    Parameters
    ----------
    df_formulas: pandas.DataFrame
        Dataframe with information about the metabolites in the model.
        The masses are all set to 0

    Returns
    -------
    df_mapping: pandas.DataFrame
        Dataframe with information about overlapping metabolite formulas in
        the model. These include the same metabolites in different
        compartments and metabolites with the same composition but different
        conformations
    """
    map_id = {}
    for formula, df in df_formulas.groupby("formula"):
        map_id[formula] = "\t".join(df["id"])
    df_mapping = df_formulas.drop_duplicates(subset=["formula"])
    df_mapping.loc[:, "id"] = df_mapping["formula"].map(map_id)
    return df_mapping


def store_struct(df_formulas, name, dirpath):
    """
    Parameters
    ----------
    df_formulas: pandas.DataFrame
        Dataframe with information about the metabolites in the model.
        The masses are all set to 0
    name: str
        Name of the database
    dirpath: path
        Path to where the database should be saved
    """
    df_struct_mapping = pd.DataFrame(
        {
            "id": df_formulas["id"],
            "formula": df_formulas["formula"],
            "smiles": ["smiles"] * len(df_formulas["id"]),
            "inchi": ["inchi"] * len(df_formulas["id"]),
        }
    )
    df_struct_mapping.to_csv(
        f"{dirpath}/{name}_struct.tsv", sep="\t", header=None, index=None
    )


def store_mapping(df_mapping, name, dirpath):
    """
    Parameters
    ----------
    df_mapping: pandas.DataFrame
        Dataframe with information about overlapping metabolite formulas
        in the model. These include the same metabolites in different
        compartments and metabolites with the same composition but
        different conformations
    name: str
        Name of the database
    dirpath: path
        Path to where the database should be saved
    """
    filename_mapping_csv = f"{dirpath}/{name}_mapping_csv.tsv"
    filename_mapping = f"{dirpath}/{name}_mapping.tsv"
    df_mapping.to_csv(
        filename_mapping_csv, sep="\t", header=None, index=None, quoting=0
    )
    with open(filename_mapping_csv) as f:
        with open(filename_mapping, "w") as fm:
            fm.write(f"database_name	{name}\ndatabase_version	{name}\n")
            for l in f:
                fm.write(l.replace('"', ""))


def create_database(metabolites, name, dirpath):
    """
    Parameters
    ----------
    metabolites: cobra.DictList or list
        List of the metabolites in the model
    name: str
        Name of the database
    dirpath: path
        Path to where the database should be saved
    """
    formulas = []
    ids = []
    for m in metabolites:
        if not is_valid(m):
            continue
        formula = zero_charge(m)
        formulas.append(str(formula))
        ids.append(m.id)

    df_formulas = make_struct(formulas, ids)
    df_mapping = make_mapping(df_formulas)

    store_struct(df_formulas, name, dirpath)
    store_mapping(df_mapping, name, dirpath)
