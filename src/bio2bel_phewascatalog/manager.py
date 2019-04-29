# -*- coding: utf-8 -*-

"""Manager for Bio2BEL PheWAS Catalog."""

import os
from typing import Mapping, Optional
from zipfile import ZipFile

import pandas as pd
from tqdm import tqdm

from bio2bel.downloading import make_downloader, make_zipped_df_getter
from bio2bel.manager.bel_manager import BELManagerMixin
from pybel import BELGraph
from pybel.dsl import Gene, Pathology
from .constants import DATA_PATH, DATA_URL, MODULE_NAME

columns = [
    "chromosome",
    "snp",
    "phewas phenotype",
    "cases",
    "p-value",
    "odds-ratio",
    "gene_name",
    "phewas code",
    "gwas-associations"
]

download_data = make_downloader(DATA_URL, DATA_PATH)
"""Download the data from Denny JC, *et al.* 2013."""

ZIP_INTERNAL_PATH = 'phewas-catalog.csv'

get_df = make_zipped_df_getter(DATA_URL, DATA_PATH, ZIP_INTERNAL_PATH, sep=',', header=0)


def extract_data(path: Optional[str] = None) -> pd.DataFrame:
    """Load the data."""
    with ZipFile(path or DATA_PATH) as zip_file:
        with zip_file.open(ZIP_INTERNAL_PATH) as file:
            return pd.read_csv(file, sep=",", header=0)


def make_graph() -> BELGraph:
    """Convert the data to a BEL graph."""
    if not os.path.exists(DATA_PATH):
        download_data()
    df = extract_data()
    return _make_graph(df)


def _make_graph(
        df: pd.DataFrame,
        use_tqdm: bool = True,
) -> BELGraph:
    graph = BELGraph(
        name="PheWAS gene-phenotype relationships",
        version="1.0.0",
    )
    it = df[["snp", 'gene_name', 'phewas phenotype', 'odds-ratio']].iterrows()
    if use_tqdm:
        it = tqdm(it, total=len(df.index), desc='PheWAS Catalog - generating BEL')
    for i, (snp, gene_symbol, phenotype, odds_ratio) in it:

        if not snp or not gene_symbol or not phenotype or pd.isna(phenotype):
            print('Skipping', i, snp, gene_symbol, phenotype, odds_ratio)
            continue

        graph.add_association(
            Gene("dbsnp", snp),
            Pathology("mesh", phenotype),
            citation="24270849",
            evidence="from PheWAS database",
            annotations={
                'bio2bel': MODULE_NAME,
                'OR': odds_ratio,
            }
        )

        if pd.notna(gene_symbol):
            graph.add_association(
                Gene("hgnc", gene_symbol),
                Pathology("mesh", phenotype),
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

    def __init__(self, *args, **kwargs):
        self._graph = None

    @property
    def graph(self) -> BELGraph:
        if self._graph is None:
            self._graph = make_graph()
        return self._graph

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
        return sum(isinstance(node, Pathology) for node in self.graph)

    def count_genes(self) -> int:
        """Count the number of diseases in the database."""
        return sum(isinstance(node, Gene) for node in self.graph)

    def count_relations(self) -> int:
        """Count the number of gene-phenotypes relationships."""
        return self.graph.number_of_edges()

    def to_bel(self) -> BELGraph:
        """Export as a BEL graph."""
        return self.graph
