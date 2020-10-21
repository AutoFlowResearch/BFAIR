import pytfa.io
from pathlib import Path

from pytfa.io.enrichment import (
    read_lexicon,
    annotate_from_lexicon,
    read_compartment_data,
    apply_compartment_data,
)
from cobra.test import create_test_model

import pytest


@pytest.fixture
def data_directory():
    tests_folder = Path(__file__).parents[0]
    return tests_folder / "test_data/"


@pytest.fixture
def create_small_model_for_tesing(data_directory):

    import logging

    # Small model for simpler tests

    small_model = create_test_model("textbook")

    # Make your computations on it
    # tmodel = pytfa.ThermoModel(thermo_data, cobra_model)

    lexicon = read_lexicon(str(data_directory / "iJO1366" / "lexicon.csv"))

    compartment_data = read_compartment_data(
        str(data_directory / "iJO1366" / "compartment_data.json")
    )

    thermo_data = pytfa.io.load_thermoDB(str(data_directory / "thermo_data.thermodb"))
    # Initialize the cobra_model
    small_tmodel = pytfa.ThermoModel(thermo_data, small_model)
    logging.getLogger("thermomodel_").setLevel(logging.ERROR)

    # Annotate the cobra_model
    annotate_from_lexicon(small_tmodel, lexicon)
    apply_compartment_data(small_tmodel, compartment_data)

    small_tmodel.solver = "optlang-glpk"
    small_tmodel.prepare()
    small_tmodel.convert()
    small_tmodel.optimize()

    return small_tmodel
