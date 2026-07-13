# Changelog

All notable changes to PA3Py will be documented here.

---

## [1.0.0] — 2026-07-13

First stable release.

### Physics
- Pebble accretion engine following Ormel (2017) and Drążkowska et al. (2023)
- Headwind/Shear regime transition via $M_{\rm hw/sh} = v_{\rm hw}^3 / (8G\Omega_K t_{\rm stop})$
- Smooth 2D–3D turbulence transition (Ormel 2017 Eq. 7.24)
- Safronov ballistic onset below $M_{\rm PA,onset}$
- Isolation mass cap: $M_{\rm iso} = 25\,M_\oplus\,(h/0.05)^3\,(M_\star/M_\odot)$ — embryo history continues after reaching $M_{\rm iso}$ for consistency
- Dynamic snowline from Oka et al. (2011) + Hartmann accretion rate decay

### Architecture
- `PA3Py` facade class: single entry point for load → run → plot
- `PebbleAccretionModule3`: physics engine, agnostic to HDF5 format
- `CompositionModel` ABC with three implementations: `SimpleWaterComposition`, `MultiSnowlineComposition`, `FunctionComposition`
- `DiskData` dataclass for structured disk snapshot storage
- All physical constants centralized in `constants.py` (CGS)

### Code quality
- Removed all hardcoded constants from module scope; everything references `constants.py`
- Eliminated `self.M_EARTH / M_SUN / AU / G_CGS` instance aliases
- Dynamic HDF5 array loading via dict maps (no repetitive `.append()` chains)
- Array broadcasting for `calculate_isolation_mass_map` (no Python loop over time)

### Tests
- 22 tests across 5 test files
- Pure unit tests for all composition models (no HDF5 required)
- Physical invariant tests: history length, $M_{\rm core} \leq M_{\rm iso}$, isolation mass map shape
- Constant sanity checks including Keplerian dimensional consistency
- Hovmöller smoke test for all three field modes

### Known limitations
- Only compatible with TripodPy HDF5 output format
- No N-body interactions between embryos (pebble flux shadowing only)
- No gap-opening feedback on disk structure

---

## Ideas / Future directions
- Implicit gap-opening feedback on $\Sigma_{\rm gas}$ near embryo location
- Support for DustPy HDF5 output
- Population synthesis mode (batch embryo grid)
- PyPI release (`pip install pa3py`)
