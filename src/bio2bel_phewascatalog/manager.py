# -*- coding: utf-8 -*-

"""Manager for Bio2BEL PheWAS Catalog."""

import logging
from typing import Dict, Mapping, Optional

import pandas as pd
from bio2bel.manager.bel_manager import BELManagerMixin
from protmapper.api import hgnc_name_to_id
from tqdm import tqdm

import pybel.dsl
from pybel import BELGraph
from .constants import MODULE_NAME
from .parser import get_df, make_dict

logger = logging.getLogger(__name__)

columns = [
    "chromosome",
    "snp",
    "phewas phenotype",
    "cases",
    "p-value",
    "odds-ratio",
    "gene_name",
    "phewas code",
    "gwas-associations",
]


def make_graph(
    df: Optional[pd.DataFrame] = None,
    use_tqdm: bool = True,
) -> BELGraph:
    """Convert the data to a BEL graph."""
    if df is None:
        df = get_df()

    graph = BELGraph(
        name="PheWAS gene-phenotype relationships",
        version="1.0.0",
    )

    it = df[["snp", 'gene_name', 'phewas phenotype', 'odds-ratio']].iterrows()
    if use_tqdm:
        it = tqdm(it, total=len(df.index), desc='PheWAS Catalog - generating BEL')
    for i, (snp, gene_symbol, phenotype, odds_ratio) in it:
        if not snp or not gene_symbol or not phenotype or pd.isna(phenotype):
            logger.debug('Skipping', i, snp, gene_symbol, phenotype, odds_ratio)
            continue

        graph.add_association(
            pybel.dsl.Gene("dbsnp", identifier=snp),
            pybel.dsl.Pathology("mesh", name=phenotype),
            citation="24270849",
            evidence="from PheWAS database",
            annotations={
                'bio2bel': MODULE_NAME,
                'OR': odds_ratio,
            }
        )

        if pd.notna(gene_symbol):
            hgnc_id = hgnc_name_to_id.get(gene_symbol)
            if hgnc_id is None:
                it.write(f'Missing identifier for {gene_symbol}')
            else:
                graph.add_association(
                    pybel.dsl.Gene(
                        namespace="hgnc",
                        name=gene_symbol,
                        identifier=hgnc_id,
                    ),
                    pybel.dsl.Pathology(
                        namespace="mesh",
                        name=phenotype,
                    ),
                    citation="24270849",
                    evidence="from PheWAS database",
                    annotations={
                        'bio2bel': MODULE_NAME,
                        'OR': odds_ratio,
                    }
                )

    return graph


class Manager(BELManagerMixin):
    """Gene-disease relationships."""

    module_name = MODULE_NAME

    def __init__(self, *args, **kwargs):
        self._graph = None
        self._dict = None

    @property
    def graph(self) -> BELGraph:
        if self._graph is None:
            self._graph = make_graph()
        return self._graph

    @property
    def dict(self) -> Dict:
        if self._dict is None:
            self._dict = make_dict()
        return self._dict

    @classmethod
    def _get_connection(cls):
        pass

    def summarize(self) -> Mapping[str, int]:
        """Summarize the database."""
        return dict(
            associations=self.count_relations(),
            diseases=self.count_diseases(),
        )

    @staticmethod
    def is_populated() -> bool:
        """Return if the database is populated."""
        return True

    def count_diseases(self) -> int:
        """Count the number of diseases in the database."""
        return sum(isinstance(node, pybel.dsl.Pathology) for node in self.graph)

    def count_genes(self) -> int:
        """Count the number of diseases in the database."""
        return sum(isinstance(node, pybel.dsl.Gene) for node in self.graph)

    def count_relations(self) -> int:
        """Count the number of gene-phenotypes relationships."""
        return self.graph.number_of_edges()

    def to_bel(self) -> BELGraph:
        """Export as a BEL graph."""
        return self.graph

    def to_dict(self) -> Dict:
        """Export as a BEL dictionary with gene as keys."""
        return self.dict
