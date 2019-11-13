# -*- coding: utf-8 -*-

"""Disease associations to variants and genes."""

from .constants import get_version
from .manager import Manager
from .parser import get_df

__all__ = [
    'Manager',
    'get_version',
]
