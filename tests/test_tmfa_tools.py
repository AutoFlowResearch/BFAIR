# these test is to ensure that the object is not None, has non-zero lengths and that not all the values are 0
def test_get_fluxes(create_small_model_for_tesing):
    from AutoFlow_OmicsDataHandling.tmfa_tools import get_fluxes

    calc = get_fluxes(create_small_model_for_tesing)
    assert (
        (len(calc) > 0)
        & (calc is not None)
        & (calc.notnull().any())
        & ((calc != 0).any())
    )


# these test is to ensure that the object is not None, has non-zero lengths and that not all the values are 0
def test_get_log_concentrations(create_small_model_for_tesing):
    from AutoFlow_OmicsDataHandling.tmfa_tools import get_log_concentrations

    calc = get_log_concentrations(create_small_model_for_tesing)
    assert (
        (len(calc) > 0)
        & (calc is not None)
        & (calc.notnull().any())
        & ((calc != 0).any())
    )


# these test is to ensure that the object is not None, has non-zero lengths and that not all the values are 0
def test_get_free_energies(create_small_model_for_tesing):
    from AutoFlow_OmicsDataHandling.tmfa_tools import get_free_energies

    calc = get_free_energies(create_small_model_for_tesing)
    assert (
        (len(calc) > 0)
        & (calc is not None)
        & (calc.notnull().any())
        & ((calc != 0).any())
    )
