# -*- coding: utf-8 -*-

"""Download PheWAS phenotype-gene relationship data and convert to BELGraph."""

import os
from typing import Mapping, Optional
from zipfile import ZipFile

import pandas as pd
from tqdm import tqdm

from bio2bel.downloading import make_downloader
from bio2bel.manager.bel_manager import BELManagerMixin
from pybel import BELGraph
from pybel.dsl import gene, Pathology
from .constants import DATA_DIR, MODULE_NAME

DATA_URL = 'https://phewascatalog.org/files/phewas-catalog.csv.zip'
DATA_PATH = os.path.join(DATA_DIR, 'data.zip')

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


def extract_data(path: Optional[str] = None) -> pd.DataFrame:
    """Load the data."""
    with ZipFile(path or DATA_PATH) as zip_file:
        with zip_file.open('data/DataS4_disease_pairs.csv') as file:
            return pd.read_csv(file, sep=",", header=0)


def make_graph():
    """Convert the data to a BEL graph."""
    if not os.path.exists(DATA_PATH):
        download_data()
    df = extract_data()
    return _make_graph(df)


def _make_graph(df: pd.DataFrame,
                use_tqdm: bool = False,
                min_network_separation: float = 0,
                ) -> BELGraph:
    graph = BELGraph(
        name="PheWAS gene-phenotype relationships",
        version="1.0.0",
    )
    it = df[['gene_name', 'odds-ratio', 'phewas phenotype']].iterrows()
    if use_tqdm:
        it = tqdm(it, total=len(df.index), desc='generating BEL')
    for _, (gene_name, odds_ratio, phenotype) in it:
        if not gene_name or \
                not phenotype or \
                phenotype > min_network_separation:
            continue
        graph.add_association(
            gene("HGNC", gene_name),
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
    """Disease-disease relationships."""

    def __init__(self, *args, **kwargs):
        self.graph = make_graph()

    @classmethod
    def _get_connection(self):
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
        return len([1 for x in self.graph.nodes if isinstance(x, Pathology)])

    def count_genes(self) -> int:
        """Count the number of diseases in the database."""
        return len([1 for x in self.graph.nodes if isinstance(x, gene)])

    def count_relations(self) -> int:
        """Count the number of gene-phenotypes relationships."""
        return self.graph.number_of_edges()

    def to_bel(self) -> BELGraph:
        """Export as a BEL graph"""
        return self.graph


main = Manager.get_cli()

if __name__ == '__main__':
    main()
