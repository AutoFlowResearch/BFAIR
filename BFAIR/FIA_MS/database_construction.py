import numpy as np
import pandas as pd
from cobra.core.formula import Formula

__version__ = "0.0.1"


INVALID_FORMULA_STR = ["(", "Generic", "R", "X"]


def print_formula(elements):
    return "".join([f"{k}{int(v)}" for k, v in elements.items()])


def zero_charge(metabolite):
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
    if not metabolite.formula:
        return False
    for string in INVALID_FORMULA_STR:
        if string in metabolite.formula:
            return False


def make_struct(formulas, ids):
    masses = [0] * len(formulas)
    df_formulas = pd.DataFrame(columns=["mass", "formula", "id"])
    df_formulas["mass"] = masses
    df_formulas["formula"] = formulas
    df_formulas["id"] = ids
    return df_formulas


def make_mapping(df_formulas):
    map_id = {}
    for formula, df in df_formulas.groupby("formula"):
        map_id[formula] = "\t".join(df["id"])
    df_mapping = df_formulas.drop_duplicates(subset=["formula"])
    df_mapping["id"] = df_mapping["formula"].map(map_id)
    return df_mapping


def store_struct(df_formulas, name, dirpath):
    df_struct_mapping = pd.DataFrame(
        columns=["id", "formula", "smiles", "inchi"]
    )
    df_struct_mapping["formula"] = df_formulas["formula"]
    df_struct_mapping["id"] = df_formulas["id"]
    df_struct_mapping["smiles"] = "smiles"
    df_struct_mapping["inchi"] = "inchi"
    df_struct_mapping.to_csv(
        f"{dirpath}/{name}_struct.tsv", sep="\t", header=None, index=None
    )


def store_mapping(df_mapping, name, dirpath):
    filename_mapping_csv = f"{dirpath}/{name}_mapping_csv.tsv"
    filename_mapping = f"{dirpath}/{name}_mapping.tsv"
    df_mapping.to_csv(
        filename_mapping_csv, sep="\t", header=None, index=None, quoting=0
    )
    with open(filename_mapping_csv) as f:
        with open(filename_mapping, "w") as fm:
            fm.write(f"database_name	{name}\ndatabase_version	{name}\n")
            for line in f:
                fm.write(line.replace('"', ""))


def create_database(metabolites, name, dirpath):
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
