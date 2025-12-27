"""
Estimation module reimplemented with numba (replacing velocyto's Cython implementation)

Original velocyto.estimation used Cython+OpenMP for performance.
This implementation uses numba.jit with parallel=True for similar performance
without requiring compilation during installation.

Author: Adapted for CellOracle 
License: BSD-2-Clause (same as original velocyto)
"""

import numpy as np
import multiprocessing
from numba import jit, prange
from typing import *


@jit(nopython=True, parallel=True, cache=True)
def _colDeltaCor_numba(e: np.ndarray, d: np.ndarray, out: np.ndarray) -> None:
    """
    Calculate correlation between displacement and cell-to-cell differences
    
    Numba implementation of velocyto's Cython colDeltaCor
    
    Arguments:
        e: gene expression matrix (ngenes, ncells)
        d: gene velocity/displacement matrix (ngenes, ncells)
        out: output matrix (ncells, ncells) - modified in place
    """
    rows, cols = e.shape
    
    for c in prange(cols):
        # Compute A = e - e[:, c]
        A = np.zeros((rows, cols))
        for j in range(rows):
            for i in range(cols):
                A[j, i] = e[j, i] - e[j, c]
        
        # Compute mean of A along axis 0
        muA = np.zeros(cols)
        for j in range(rows):
            for i in range(cols):
                muA[i] += A[j, i]
        for i in range(cols):
            muA[i] /= rows
        
        # A_mA = A - muA
        A_mA = np.zeros((rows, cols))
        for j in range(rows):
            for i in range(cols):
                A_mA[j, i] = A[j, i] - muA[i]
        
        # Compute mean of b (d[:, c])
        mub = 0.0
        for j in range(rows):
            mub += d[j, c]
        mub /= rows
        
        # b_mb = b - mub
        b_mb = np.zeros(rows)
        for j in range(rows):
            b_mb[j] = d[j, c] - mub
        
        # Compute ssA = sum of squares of A_mA
        ssA = np.zeros(cols)
        for j in range(rows):
            for i in range(cols):
                ssA[i] += A_mA[j, i] * A_mA[j, i]
        for i in range(cols):
            if ssA[i] > 0:
                ssA[i] = 1.0 / np.sqrt(ssA[i])
            else:
                ssA[i] = 0.0
        
        # Compute ssb = sum of squares of b_mb
        ssb = 0.0
        for j in range(rows):
            ssb += b_mb[j] * b_mb[j]
        if ssb > 0:
            ssb = 1.0 / np.sqrt(ssb)
        else:
            ssb = 0.0
        
        # Compute correlation: dot(b_mb, A_mA) / (sqrt(ssA) * sqrt(ssb))
        for j in range(rows):
            tmp = b_mb[j] * ssb
            for i in range(cols):
                out[c, i] += (A_mA[j, i] * ssA[i]) * tmp


@jit(nopython=True, parallel=True, cache=True)
def _colDeltaCorpartial_numba(e: np.ndarray, d: np.ndarray, out: np.ndarray, ixs: np.ndarray) -> None:
    """
    Calculate correlation only for specified neighbors (partial version)
    
    Numba implementation of velocyto's Cython colDeltaCorpartial
    
    Arguments:
        e: gene expression matrix (ngenes, ncells)
        d: gene velocity/displacement matrix (ngenes, ncells)
        out: output matrix (ncells, ncells) - modified in place
        ixs: neighborhood matrix (ncells, nneighbors)
    """
    rows, cols = e.shape
    n_neighbors = ixs.shape[1]
    
    for c in prange(cols):
        neighbors = ixs[c, :]
        
        # Compute A for neighbors only
        A = np.zeros((rows, n_neighbors))
        for j in range(rows):
            for k in range(n_neighbors):
                i = neighbors[k]
                A[j, k] = e[j, i] - e[j, c]
        
        # Compute mean of A
        muA = np.zeros(n_neighbors)
        for j in range(rows):
            for k in range(n_neighbors):
                muA[k] += A[j, k]
        for k in range(n_neighbors):
            muA[k] /= rows
        
        # A_mA = A - muA
        A_mA = np.zeros((rows, n_neighbors))
        for j in range(rows):
            for k in range(n_neighbors):
                A_mA[j, k] = A[j, k] - muA[k]
        
        # Compute mean of b
        mub = 0.0
        for j in range(rows):
            mub += d[j, c]
        mub /= rows
        
        # b_mb = b - mub
        b_mb = np.zeros(rows)
        for j in range(rows):
            b_mb[j] = d[j, c] - mub
        
        # Compute ssA
        ssA = np.zeros(n_neighbors)
        for j in range(rows):
            for k in range(n_neighbors):
                ssA[k] += A_mA[j, k] * A_mA[j, k]
        for k in range(n_neighbors):
            if ssA[k] > 0:
                ssA[k] = 1.0 / np.sqrt(ssA[k])
            else:
                ssA[k] = 0.0
        
        # Compute ssb
        ssb = 0.0
        for j in range(rows):
            ssb += b_mb[j] * b_mb[j]
        if ssb > 0:
            ssb = 1.0 / np.sqrt(ssb)
        else:
            ssb = 0.0
        
        # Compute correlation for neighbors
        for j in range(rows):
            tmp = b_mb[j] * ssb
            for k in range(n_neighbors):
                i = neighbors[k]
                out[c, i] += (A_mA[j, k] * ssA[k]) * tmp


def colDeltaCor(emat: np.ndarray, dmat: np.ndarray, threads: int = None) -> np.ndarray:
    """Calculate the correlation between the displacement (d[:,i])
    and the difference between a cell and every other (e - e[:, i])
    
    Numba parallel implementation (replacement for velocyto's Cython version)

    Arguments
    ---------
    emat: np.ndarray (ngenes, ncells)
        gene expression matrix
    dmat: np.ndarray (ngenes, ncells)
        gene velocity/displacement matrix
    threads: int
        number of parallel threads to use (Note: numba manages this automatically)
    """
    out = np.zeros((emat.shape[1], emat.shape[1]))
    _colDeltaCor_numba(emat, dmat, out)
    return out


def colDeltaCorpartial(emat: np.ndarray, dmat: np.ndarray, ixs: np.ndarray, threads: int = None) -> np.ndarray:
    """Calculate the correlation between the displacement (d[:,i])
    and the difference between a cell and every other (e - e[:, i])
    
    Numba parallel implementation (replacement for velocyto's Cython version)

    Arguments
    ---------
    emat: np.ndarray (ngenes, ncells)
        gene expression matrix
    dmat: np.ndarray (ngenes, ncells)
        gene velocity/displacement matrix
    ixs: np.ndarray (ncells, nneighbours)
        ixs[i, k] is the kth neighbour to the cell i
    threads: int
        number of parallel threads to use (Note: numba manages this automatically)
    """
    out = np.zeros((emat.shape[1], emat.shape[1]))
    emat = np.require(emat, requirements="C")
    dmat = np.require(dmat, requirements="C") 
    ixs = np.require(ixs, requirements="C").astype(np.intp)
    _colDeltaCorpartial_numba(emat, dmat, out, ixs)
    return out


# Stub implementations for other functions (not used by CellOracle but part of velocyto API)

def colDeltaCorLog10(emat: np.ndarray, dmat: np.ndarray, threads: int = None, psc: float = 1.0) -> np.ndarray:
    """Stub: Not implemented (not used by CellOracle)"""
    raise NotImplementedError("colDeltaCorLog10 is not used by CellOracle and not implemented in this replacement")


def colDeltaCorLog10partial(emat: np.ndarray, dmat: np.ndarray, ixs: np.ndarray, threads: int = None, psc: float = 1.0) -> np.ndarray:
    """Stub: Not implemented (not used by CellOracle)"""
    raise NotImplementedError("colDeltaCorLog10partial is not used by CellOracle and not implemented in this replacement")


def colDeltaCorSqrt(emat: np.ndarray, dmat: np.ndarray, threads: int = None, psc: float = 0.0) -> np.ndarray:
    """Stub: Not implemented (not used by CellOracle)"""
    raise NotImplementedError("colDeltaCorSqrt is not used by CellOracle and not implemented in this replacement")


def colDeltaCorSqrtpartial(emat: np.ndarray, dmat: np.ndarray, ixs: np.ndarray, threads: int = None, psc: float = 0.0) -> np.ndarray:
    """Stub: Not implemented (not used by CellOracle)"""
    raise NotImplementedError("colDeltaCorSqrtpartial is not used by CellOracle and not implemented in this replacement")

