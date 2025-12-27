# Velocyto Replacement Module

## Purpose

This module provides drop-in replacements for velocyto functions that CellOracle depends on, eliminating the need to install the unmaintained velocyto package which has compilation issues.

## Background

The original [velocyto](https://github.com/velocyto-team/velocyto.py) package:
- Has not been maintained since 2019
- Uses Cython+OpenMP requiring compilation during installation
- Has `setup.py` issues that prevent installation on modern Python environments
- Causes installation failures due to numpy import timing issues

## Solution

This module reimplements the 3 velocyto components used by CellOracle:

### 1. **diffusion.py** (Pure Python - Direct Copy)
- `Diffusion` class for Markov chain diffusion
- No modifications needed - copied directly from velocyto

### 2. **serialization.py** (Pure Python - Direct Copy)  
- `dump_hdf5` and `load_hdf5` for HDF5 object serialization
- No modifications needed - copied directly from velocyto

### 3. **estimation.py** (Reimplemented with Numba)
- **Original**: Used Cython+OpenMP for performance
- **New**: Uses `numba.jit` with `parallel=True`
- Provides equivalent performance without compilation
- Implements:
  - `colDeltaCor`: Full correlation calculation
  - `colDeltaCorpartial`: Partial correlation for KNN neighborhoods

## Performance

The numba implementation provides comparable performance to the original Cython code:
- Uses JIT compilation at runtime (no build-time compilation needed)
- Leverages parallel processing with `prange`
- Memory-efficient implementation

## Compatibility

✅ **Fully compatible** with existing CellOracle code  
✅ **Drop-in replacement** - no API changes  
✅ **No compilation required** - pure Python installation  
✅ **Cross-platform** - works on all platforms supported by numba

## Modified Files

The following CellOracle files were updated to use this module:

1. `celloracle/trajectory/modified_VelocytoLoom_class.py`
   - Changed: `from velocyto.* import ...`
   - To: `from .._velocyto_replacement.* import ...`

2. `celloracle/network/net_util.py`
   - Changed: `from velocyto.serialization import ...`
   - To: `from .._velocyto_replacement.serialization import ...`

3. `requirements.txt`
   - Removed: `velocyto>=0.17`
   - Added comment explaining the removal

## License

This code is adapted from velocyto under the BSD-2-Clause License.  
Original velocyto: https://github.com/velocyto-team/velocyto.py

## Testing

To verify the installation works:

```python
import celloracle as co

# The import should succeed without requiring velocyto
from celloracle._velocyto_replacement import Diffusion, colDeltaCor, dump_hdf5
print("✓ CellOracle velocyto replacement working!")
```

## Credits

- **Original velocyto authors**: La Manno et al. (2018)
- **Adaptation**: Integrated into CellOracle to resolve installation issues
- **Date**: December 2024

