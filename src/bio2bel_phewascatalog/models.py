# -*- coding: utf-8 -*-

"""SQLAlchemy models for Bio2BEL PheWAS Catalog."""

import logging

import pybel.dsl
from sqlalchemy.ext.declarative import DeclarativeMeta, declarative_base

from .constants import MODULE_NAME

__all__ = [
    'Base',
    'Gene',
    'Variant',
    'Phenotype'
]

logger = logging.getLogger(__name__)

GENE_TABLE_NAME = f'{MODULE_NAME}_gene'
VARIANT_TABLE_NAME = f'{MODULE_NAME}_variant'
PHENOTYPE_TABLE_NAME = f'{MODULE_NAME}_phenotype'
GENE_PHENOTYPE_TABLE_NAME = f'{MODULE_NAME}_gene_phenotype'
VARIANT_PHENOTYPE_TABLE_NAME = f'{MODULE_NAME}_variant_phenotype'
GENE_VARIANT_TABLE_NAME = f'{MODULE_NAME}_variant_gene'

Base: DeclarativeMeta = declarative_base()


class Gene(Base):
    """Represents a gene."""

    def to_bel(self) -> pybel.dsl.Gene:
        """Convert to BEL."""
        return pybel.dsl.Gene()


class Variant(Base):
    """Represents a genetic variant."""

    def to_bel(self) -> pybel.dsl.Gene:
        """Convert to BEL."""
        return pybel.dsl.Gene()


class Phenotype(Base):
    """Represents a phenotype."""

    def to_bel(self) -> pybel.dsl.Pathology:
        """Convert to BEL."""
        return pybel.dsl.Pathology()
