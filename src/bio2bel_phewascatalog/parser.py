# -*- coding: utf-8 -*-

"""Parsers and downloaders for Bio2BEL PheWAS Catalog."""

import logging
from typing import Dict

import pandas as pd
from tqdm import tqdm

from bio2bel.downloading import make_zipped_df_getter
from .constants import DATA_PATH, DATA_URL

__all__ = [
    'get_df',
    'make_dict',
]

ZIP_INTERNAL_PATH = 'phewas-catalog.csv'

get_df = make_zipped_df_getter(DATA_URL, DATA_PATH, ZIP_INTERNAL_PATH, sep=',', header=0)
"""Download the data from Denny JC, *et al.* 2013."""


def make_dict() -> Dict:
    """Convert the data to a dictionary."""
    return _make_dict(get_df())


def _make_dict(
        df: pd.DataFrame,
        use_tqdm: bool = True,
) -> Dict:
    _dict = dict()
    it = df[["snp", 'gene_name', 'phewas phenotype', 'odds-ratio', 'phewas code']].iterrows()
    if use_tqdm:
        it = tqdm(it, total=len(df.index), desc='PheWAS Catalog - generating Dict')
    for i, (snp, gene_symbol, phenotype, odds_ratio, icd_code) in it:

        if not snp or not gene_symbol or not phenotype or pd.isna(phenotype):
            logging.debug('Skipping', i, snp, gene_symbol, phenotype, odds_ratio)
            continue

        if pd.notna(gene_symbol):
            if gene_symbol in _dict:
                _dict[gene_symbol] += [(odds_ratio, icd_code)]
            else:
                _dict[gene_symbol] = [(odds_ratio, icd_code)]

    return _dict
