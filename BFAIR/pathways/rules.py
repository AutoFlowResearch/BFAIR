"""Rules.

This module contains functions to apply reaction rules to chemical compounds.
"""

import json
import sqlite3 as sql
from concurrent.futures import ThreadPoolExecutor, TimeoutError
from itertools import chain
from typing import Generator

import pandas as pd
from cobra.core.singleton import Singleton
from rdkit import Chem
from rdkit.Chem import AllChem
from BFAIR import models
from BFAIR.io._path import static_path
from BFAIR.logger import get_logger
from BFAIR.pathways.standardization import standardize
from BFAIR.pathways.utils import get_compound, get_molecular_fingerprint
from BFAIR.pathways import _queries as q


class _RuleSimulator():
    # simulates a reaction

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
    """Rule Library.

    The Rule Library object connects to the reaction rules database and exposes functions to filter the reaction rules
    by different criteria.

    Parameters
    ----------
    **kwargs
        Standardization parameters applied to input compounds, see `BFAIR.pathways.standardization.standardize`.
        Thorough standardization is enforced.

    Attributes
    ----------
    available : pandas.DataFrame
        A dataframe with the available rules (after filters applied). It has columns for rule ID, associated MetaNetX
        ID, and SMARTS expression.

    Notes
    -----
    The Rule Library class is a singleton, meaning only one instance of it exists at any time.
    """

    def __init__(self, **kwargs):
        kwargs.setdefault('add_hs', True)
        kwargs["remove_stereo"] = True
        kwargs["thorough"] = True
        if not kwargs['add_hs']:
            raise NotImplementedError()
        self._logger = get_logger(__name__)
        self._conn = sql.connect(static_path(f"rules_{'' if kwargs['add_hs'] else 'no'}hs.db"))
        self._conn.create_aggregate(q.ChemicalSimilarity.__name__, 2, q.ChemicalSimilarity)
        self._cursor = self._conn.cursor()
        self._filters = []
        self._std_args = kwargs.copy()

    def __iter__(self):
        return self.available.itertuples(index=False, name=None)

    def __len__(self):
        return len(self.available)

    def __repr__(self):
        return repr(self.available)

    def _repr_html_(self):
        return self.available._repr_html_()

    @property
    def available(self) -> pd.DataFrame:
        return pd.read_sql(q.select_by(self._filters), self._conn, index_col="rule_id")

    def _fetch_rules(self, query):
        # Returns a set of unique rules obtained by executing the specified SQL statement
        return set(chain.from_iterable(self._cursor.execute(query)))

    def apply_to(self, input_compound, input_type="inchi", timeout=60.) -> Generator:
        """
        Applies the available rules to the input compound.

        Parameters
        ----------
        input_compound : str or rdkit.Chem.rdchem.Mol
            String representation of a compound or a RDKit molecule object.
        input_type : {'inchi', 'smiles', 'rdkit'}
            Type of notation describing the input compound.
        timeout : float
            Number of seconds after which a rule simulation will be stopped.

        Yields
        ------
        tuple
            A tuple containing the applied rule ID and a list of InChI and SMILES depictions of the obtained products.
        """
        compound = get_compound(input_compound, input_type, **self._std_args)
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

        Returns
        -------
        same type as caller

        Raises
        ------
        ValueError
            If an unsupported diameter type is specified.
        """
        if cutoff not in tuple(range(2, 18, 2)):
            raise ValueError(f"Unsupported diameter: {cutoff}. See allowed values in documentation.")
        self._filters.append(q.rules.diameter >= cutoff)
        return self

    def filter_by_organism(self, model_name, additional_identifiers=None, additional_ec_numbers=None):
        """
        Applies a filter to the rule library that excludes reactions that are not available in an input genome-scale
        model.

        Parameters
        ----------
        model_name : str
            The name of a genome-scale model, from which MetaNetX reaction identifiers will be extracted. Currently,
            only pre-curated models are supported. See `BFAIR.models` for a list of available models.
        additional_identifiers, additional_ec_numbers : list of str
            List of additional MetaNetX reaction identifiers or E.C. numbers that are considered as part of the
            organism. These can be, for instance, used to include heterologous enzymatic steps not described in the
            genome-scale model.

        Returns
        -------
        same type as caller
        """
        model = getattr(models, model_name)
        reaction_ids = set(additional_identifiers)
        ec_numbers = set(additional_ec_numbers)
        for rxn in model.reactions:
            reaction_id = rxn.annotation.get('metanetx.reaction', None)
            if reaction_id:
                reaction_ids.add(reaction_id)
            ec_num = rxn.annotation.get('ec-code', [])
            if isinstance(ec_num, str):
                ec_numbers.add(ec_num)
            else:
                ec_numbers.update(ec_num)

        fetched_rules = set()
        if reaction_ids:
            fetched_rules.update(self._fetch_rules(q.select_by_reaction_id(reaction_ids, self._filters)))
            # because internal or external reaction IDs can be deprecated, we should look for synonyms in the thesaurus
            fetched_rules.update(self._fetch_rules(q.select_by_synonymous_id(reaction_ids, self._filters)))

        if ec_numbers:
            fetched_rules.update(self._fetch_rules(q.select_by_ec_number(ec_numbers, self._filters)))

        if not fetched_rules:
            self._logger.error("The model has no identifiable reactions.")

        self._filters.append(q.rules.rule_id.isin(fetched_rules))

        return self

    def filter_by_compound(self, input_compound, input_type="inchi", cutoff=0.7):
        """
        Applies a filter to excludes rules that cannot be applied to an input compound. Applicability is based on
        chemical similarity.

        Parameters
        ----------
        input_compound : str or rdkit.Chem.rdchem.Mol
            String representation of a compound or a RDKit molecule object.
        input_type : {'inchi', 'smiles', 'rdkit'}
            Type of notation describing the input compound.
        cutoff : float
            Chemical similarity score (between the input compound and the native substrate of a reaction) below which a
            reaction rule is considered inapplicable. The chemical similarity score between two compounds ranges from 0
            to 1, the latter meaning both molecules are identical.

        Returns
        -------
        same type as caller
        """
        input_fp = json.dumps(get_molecular_fingerprint(input_compound, input_type, **self._std_args))

        # execute query and create new criterion for rules matching the resulting IDs
        # we use this new criterion instead of a sub-query to avoid overhead costs of fingerprint calculation
        fetched_rules = self._fetch_rules(q.select_by_similarity(input_fp, cutoff, self._filters))

        if not fetched_rules:
            self._logger.warn(f"No reaction rules match the given compound:\n{input_compound}")

        self._filters.append(q.rules.rule_id.isin(fetched_rules))

        return self

    def filter_by_uncertainty(self, cutoff=0.1):
        """
        Applies a filter to exclude rules that do not meet the biochemical uncertainty criteria.

        Parameters
        ----------
        cutoff : float
            Biochemical uncertainty score below which reaction rules are filtered out. Each reaction rule has an
            associated score that functions as a proxy for uncertainty and specificity. Scores closer to 1 indicate less
            uncertainty and higher specificity with respect to the annotated enzyme sequences. Note that, due to how it
            is calculated, the score is biased against rarely studied enzymes.

        Returns
        -------
        same type as caller
        """
        self._filters.append(q.rules.score >= cutoff)
        return self

    def pop_filter(self, index=-1):
        """
        Removes a filter criterion from the rule library.

        Parameters
        ----------
        index : int, default -1
            Index of filter to pop (the default is -1, which points to the last filter applied).

        Returns
        -------
        same type as caller
        """
        self._filters.pop()
        return self

    def reset(self):
        """
        Resets all the applied filters.

        Returns
        -------
        same type as caller
        """
        self._filters = []
        return self
