# PA3Py — Pebble Accretion Post-Processing for TripodPy

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT) [![Tests](https://github.com/mmmaxii/PA3Py/actions/workflows/tests.yml/badge.svg)](https://github.com/mmmaxii/PA3Py/actions) [![Documentation Status](https://readthedocs.org/projects/pa3py/badge/?version=latest)](https://pa3py.readthedocs.io/en/latest/?badge=latest)

**PA3Py** (Pebble Accretion 3 in Python) is a post-processing framework for computing the growth and composition of planetary embryos via pebble accretion on top of [TripodPy](https://github.com/tripod-code/tripodpy) hydrodynamic disk simulations. It is **not** a standalone disk evolution code — it reads TripodPy HDF5 snapshots and integrates embryo growth over the fully evolved gas and dust structure, including **structured disks** with planet-carved gaps and pressure bumps.

PA3Py was developed as part of the undergraduate thesis *"Regulación del Flujo de Pebbles por Subestructuras en Discos Protoplanetarios: Implicancias para el Crecimiento Planetario y la Formación de Waterworlds"* at the Instituto de Astrofísica, Pontificia Universidad Católica de Chile, where it was used to compute over 5,000 planetary growth trajectories across 1,300+ TripodPy disk simulations spanning different turbulence levels, substructure configurations, and dust fragmentation properties.

Please read the [documentation](https://pa3py.readthedocs.io) for a full description.

## Features of the Model

Users can easily set up and run pebble accretion models on top of realistic disk simulations. Major features include:

- **Realistic pebble flux from full dust evolution.** Instead of prescribing an analytic pebble flux, PA3Py reads the gas and dust surface densities, temperature, and three-population dust size distribution evolved by TripodPy — so the pebble supply is regulated self-consistently by disk substructures (gaps, pressure traps, sinusoidal perturbations) and by the fragmentation physics of the dust.

- **A dynamic water snowline based on Oka et al. (2011) and Hartmann et al. (1998).** The snowline position is computed from the radiative disk models of Oka et al. (2011), coupled to the decay of the stellar accretion rate ($\dot{M}_\star \propto t^{-\eta}$, Hartmann et al. 1998). The interpolation table ships with the package, and the snowline can also be set explicitly or replaced by any user prescription.

- **Composition-agnostic chemistry with evolving zones.** By default PA3Py tracks dry silicates inside the snowline and a 50/50 silicate–ice mixture beyond it, but any number of species and radial zones can be injected as a plain Python function (`FunctionComposition`) or as multiple migrating snowlines (`MultiSnowlineComposition`). Species are auto-discovered and tracked through the full growth history.

- **Pebble accretion physics following Ormel (2017) and Drążkowska et al. (2023):** onset (Safronov) mass check, headwind/shear regime switching, smooth 2D–3D transition set by the pebble scale height, and growth capped by the pebble isolation mass (Lambrechts & Johansen 2014), which acts as the physical ceiling of every growth track.

- **Any number of embryos, anywhere in the disk,** with configurable seed mass — from a single tracked embryo to a synthetic population in one call.

- **Built-in, publication-ready plotting:** Hovmöller diagrams of the dust evolution, mass growth curves (log-scale, with the isolation-mass ceiling overlaid), single- and multi-species composition evolution, and synthetic population diagrams (final mass vs. initial position, colored by water fraction).

- **HDF5 I/O for growth tracks,** preserving species metadata for later analysis.

Because PA3Py runs in post-processing, accreted mass is not removed from the disk (no pebble filtering between embryos); growth is instead bounded by the isolation mass. A dynamic TripodPy–PA3Py coupling is planned future work.

## Installation

### Option 1: Install from PyPI (Recommended)
You can easily install the latest stable version of PA3Py using pip:

```bash
pip install pa3py
```

### Option 2: Install from Source
If you want to modify the source code or get the latest development version, you can clone the repository:

```bash
git clone https://github.com/mmmaxii/PA3Py.git
cd PA3Py
pip install .
```

For an editable installation (recommended during development):

```bash
pip install -e ".[dev]"
```

## Requirements

- Python ≥ 3.9
- [TripodPy](https://github.com/tripod-code/tripodpy) simulation output (HDF5 snapshots)
- `numpy`, `scipy`, `h5py`, `matplotlib`

## Quick Start

```python
from pa3py import PA3Py

sim = PA3Py("path/to/tripodpy/output/")

# Run pebble accretion for embryos at 2, 5, 10, 20 AU
results = sim.run_growth([2.0, 5.0, 10.0, 20.0])

# Hovmöller diagram (dust-to-gas ratio vs time)
sim.plot_hovmoller(field='epsilon')

# Mass growth curves (log scale, isolation-mass ceiling; works for 1 or many embryos)
sim.plot_growth_curves(results)

# Species mass-fraction evolution: one species, a list, or None for all tracked species
sim.plot_species_fraction(results, species='H2O')

# Synthetic population: final mass vs initial position, with isolation-mass limit
sim.plot_population(results)
```

For custom chemistry (multi-snowline, 4-species):

```python
from pa3py import PA3Py, FunctionComposition
from pa3py import constants as c

def my_chemistry(r_cm, t_sec):
    if r_cm < 3.0 * c.AU:
        return {'silicates': 1.0}
    elif r_cm < 5.0 * c.AU:
        return {'silicates': 0.5, 'H2O': 0.5}
    else:
        return {'silicates': 0.3, 'H2O': 0.3, 'CO2': 0.4}

sim = PA3Py("path/to/output/", comp_model=FunctionComposition(my_chemistry))
results = sim.run_growth([5.0, 10.0, 20.0])

# Full composition history — colors per species for one embryo,
# line styles per species when comparing several embryos
sim.plot_species_fraction(results, species=None)
```

## Related Codes

- **[PPOLS](https://spmccloat.github.io/thePPOLSmodel/)** (McCloat et al. 2025) — a standalone pebble accretion model with an evolving dust reservoir, a temperature-based or self-consistent snowline, and pebble filtering between seeds. PA3Py takes a complementary approach: it post-processes full TripodPy dust-evolution simulations (capturing substructure-regulated pebble fluxes that analytic disks cannot) and parameterizes the snowline migration through the Oka et al. (2011) + Hartmann et al. (1998) prescription.
- **[Pebble Predictor](https://github.com/astrojoanna/pebble-predictor)** (Drążkowska et al. 2021) — analytic estimates of the pebble flux and pebble sizes in smooth disks.
- **[TripodPy](https://github.com/tripod-code/tripodpy)** (Pfeil et al. 2024; Kaufmann et al. 2025) — the gas + three-population dust evolution framework that produces the input snapshots for PA3Py.

## Physics & References

- **Ormel (2017)** — *The Emerging Paradigm of Pebble Accretion*: accretion regimes and rates
- **Drążkowska et al. (2023)** — *Planet Formation Theory in the Era of ALMA and Kepler* (Protostars & Planets VII)
- **Lambrechts & Johansen (2014)** — pebble isolation mass
- **Oka, Nakamoto & Ida (2011)** — radiative disk models for the snowline position
- **Hartmann et al. (1998)** — stellar accretion-rate decay driving the snowline migration
- **Drążkowska, Stammler & Birnstiel (2021)** — the *pebble snow* scenario
- **Pfeil, Birnstiel & Klahr (2024)** — TripodPy three-population dust evolution

## Citing & Acknowledgements

PA3Py was developed by **Maximiliano Valderrama Vargas** at the Instituto de Astrofísica, Pontificia Universidad Católica de Chile, under the supervision of **Prof. Gijs Mulders**, as part of the Licenciatura en Astronomía thesis *"Regulación del Flujo de Pebbles por Subestructuras en Discos Protoplanetarios"* (2026).

If you use PA3Py in your research, please cite the thesis and link to this repository.

Contact: [mvalderrav@uc.cl](mailto:mvalderrav@uc.cl)
