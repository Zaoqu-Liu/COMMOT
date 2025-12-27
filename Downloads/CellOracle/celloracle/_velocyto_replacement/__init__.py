"""
Velocyto Replacement Module for CellOracle

This module provides drop-in replacements for velocyto functions that CellOracle depends on,
eliminating the need to install the unmaintained velocyto package.

The original velocyto code used Cython for performance. This implementation uses numba instead,
which provides similar performance without requiring compilation during installation.

Components:
- diffusion.py: Markov chain diffusion (pure Python from velocyto)
- estimation.py: Correlation calculations (reimplemented with numba)  
- serialization.py: HDF5 object serialization (pure Python from velocyto)

Author: Adapted for CellOracle by removing velocyto dependency
License: BSD-2-Clause (same as original velocyto)
"""

from .diffusion import Diffusion
from .estimation import colDeltaCor, colDeltaCorpartial, colDeltaCorLog10, colDeltaCorLog10partial, colDeltaCorSqrt, colDeltaCorSqrtpartial
from .serialization import dump_hdf5, load_hdf5

__all__ = [
    'Diffusion',
    'colDeltaCor', 
    'colDeltaCorpartial',
    'colDeltaCorLog10',
    'colDeltaCorLog10partial', 
    'colDeltaCorSqrt',
    'colDeltaCorSqrtpartial',
    'dump_hdf5',
    'load_hdf5'
]

