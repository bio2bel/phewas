# -*- coding: utf-8 -*-

"""Constants for Bio2BEL PheWAS Catalog."""

import os

from bio2bel import get_data_dir

__all__ = [
    'VERSION',
    'MODULE_NAME',
    'DATA_DIR',
    'DATA_URL',
    'DATA_PATH',
    'get_version',
]

VERSION = '0.0.1-dev'
MODULE_NAME = 'phewascatalog'
DATA_DIR = get_data_dir(MODULE_NAME)

DATA_URL = 'https://phewascatalog.org/files/phewas-catalog.csv.zip'
DATA_PATH = os.path.join(DATA_DIR, 'phewas-catalog.csv.zip')


def get_version() -> str:
    """Get the software version."""
    return VERSION
