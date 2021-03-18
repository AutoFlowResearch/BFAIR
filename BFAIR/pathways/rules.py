"""Rules.

This module contains functions to apply reaction rules to chemical compounds.
"""

import json
from functools import reduce
from concurrent.futures import ThreadPoolExecutor, TimeoutError

import numpy as np
import pandas as pd

from cobra.core.singleton import Singleton
from rdkit import Chem
from rdkit.Chem import AllChem

from BFAIR.logger import get_logger
from BFAIR.pathways.standardization import standardize
from BFAIR.pathways.utils import get_compound, get_molecular_fingerprint, calculate_similarity


class _RuleSimulator():
    def __init__(self, compound: AllChem.Mol, reaction: AllChem.ChemicalReaction):
        self._compound = compound
        self._reaction = reaction
        self._interrupted = False

    def __call__(self):
        # can_run = any([compound.HasSubstructMatch(substrate) for substrate in reaction.GetReactants()])
        results = {}
        for products in self._reaction.RunReactants((self._compound,)):
            if self._interrupted:
                break
            products = [
                standardize(frag)
                for product in products for frag in Chem.GetMolFrags(product, asMols=True, sanitizeFrags=False)
            ]
            inchis = [Chem.MolToInchi(product) for product in products]
            inchikeys = '.'.join([Chem.InchiToInchiKey(inchi) for inchi in inchis])
            if inchikeys not in results:
                results[inchikeys] = inchis
        return [*results.values()]

    def interrupt(self):
        self._interrupted = True


class RuleLibrary(metaclass=Singleton):
    def __init__(self, data: pd.DataFrame):
        self._data = data
        for col in ["substrate_fps", "product_fps"]:
            self._data[col] = self._data[col].apply(json.loads)
        self._logger = get_logger(__name__)
        self.reset()

    def __iter__(self):
        return self.available.itertuples(index=False, name=None)

    def __len__(self):
        return len(self.available)

    def __repr__(self):
        return repr(self.available)

    def _repr_html_(self):
        return self.available._repr_html_()

    @property
    def available(self):
        if self._filters:
            return self._data[reduce(np.logical_and, self._filters)]
        else:
            return self._data

    def apply_to(self, input_compound, input_type="inchi", timeout=60, **kwargs):
        """
        Applies the available rules to the input compound.

        Parameters
        ----------
        input_compound : str
            String representation of a chemical compound.
        input_type : {'inchi', 'smiles'}
            Type of notation describing the input compound.
        timeout : int
            Number of seconds after which a rule simulation will be stopped.
        **kwargs
            Standardization parameters applied to the input compound, see `BFAIR.pathways.standardization.standardize`.
            Thorough standardization is enforced.

        Yields
        ------
        tuple
            A tuple containing the applied rule ID and a list of InChI and SMILES depictions of the obtained products.
        """
        # Force thorough standardization
        kwargs["thorough"] = True
        compound = get_compound(input_compound, input_type, **kwargs)
        with ThreadPoolExecutor(max_workers=1) as executor:
            tasks = [
                _RuleSimulator(compound, reaction) for reaction in self.available["rule_smarts"]
            ]
            for i, (task, future) in enumerate([(task, executor.submit(task)) for task in tasks]):
                try:
                    yield (self.available.index[i], future.result(timeout))
                except TimeoutError:
                    self._logger.warn(f"Timed out processing rule '{self.available.index[i]}'.")
                    task.interrupt()
            del tasks

    def filter_by_diameter(self, cutoff):
        """
        Applies a filter to exclude rules that do not meet the diameter criteria.

        Parameters
        ----------
        cutoff : {2, 4, 6, 8, 10, 12, 14, 16}
            Diameter below which rules are filtered out. The diameter is a proxy for enzyme promiscuity. Formally, it is
            defined as the distance in molecular bonds around the reaction center. Thus, small diameters capture general
            reactions; whereas, bigg diameters correspond to specific reactions.
        """
        if cutoff not in tuple(range(2, 18, 2)):
            raise ValueError(f"Unsupported diameter: {cutoff}. See allowed values in documentation.")
        self._filters.append(self._data["diameter"] >= cutoff)

    def filter_by_metanetx_id(self, *args):
        """
        Applies a filter to the rule library that excludes reactions that do not match the given MNXref identifiers.
        """
        raise NotImplementedError()

    def filter_by_organism(self, model_name, additional_identifiers=None):
        """
        Applies a filter to the rule library that excludes reactions that are not available in the input genome-scale
        model.
        """
        raise NotImplementedError()

    def filter_by_compound(
        self, input_compound, input_type="inchi", bio_score_cutoff=0.1, chem_score_cutoff=0.7, **kwargs
    ):
        """
        Applies a filter to excludes rules that cannot be applied to an input compound. Applicability is based on
        biochemical uncertainty and chemical similarity.

        Parameters
        ----------
        input_compound : str or rdkit.Chem.rdchem.Mol
            String representation of a compound or a RDKit molecule object.
        input_type : {'inchi', 'smiles', 'rdkit'}
            Type of notation describing the input compound.
        bio_score_cutoff : float
            Biochemical uncertainty score below which reaction rules are filtered out. Each reaction rule has an
            associated score that functions as a proxy for uncertainty and specificity. Scores closer to 1 indicate less
            uncertainty and higher specificity with respect to the annotated enzyme sequences. Note that, due to how it
            is calculated, the score is biased against rarely studied enzymes.
        chem_score_cutoff : float
            Chemical similarity score (between the input compound and the native substrate of a reaction) below which a
            reaction rule is considered inapplicable. The chemical similarity score between two compounds ranges from 0
            to 1, the latter meaning both molecules are identical.
        **kwargs
            Standardization parameters applied to the input compound, see `BFAIR.pathways.standardization.standardize`.
        """
        input_fp = get_molecular_fingerprint(input_compound, input_type, **kwargs)
        # temporarily add filter to limit the number of calculated similarity scores
        self._filters.append(self._data["bio_score"] >= bio_score_cutoff)
        chem_score_filter = self._data.index.isin(
            [
                rule_id
                for rule_id, substrate_fps in zip(self.available.index, self.available["substrate_fps"])
                if max([calculate_similarity(input_fp, native_fp) for native_fp in substrate_fps]) >= chem_score_cutoff
            ]
        )
        bio_score_filter = self.pop_filter()
        self._filters.append(np.logical_and(bio_score_filter, chem_score_filter))
        if len(self.available) == 0:
            self._logger.warn(f"No reaction rules match the given compound:\n{input_compound}")

    def pop_filter(self, index=-1):
        """
        Removes the filter (a boolean array) at the given index from the list of filters.

        Parameters
        ----------
        index : int
            Index of filter to be popped.

        Returns
        -------
        array-like
            Popped filter.
        """
        return self._filters.pop()

    def reset(self):
        """
        Resets all the applied filters.
        """
        self._filters = []
