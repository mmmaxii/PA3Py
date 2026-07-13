# PA3Py — Pebble Accretion Post-Processing for TripodPy

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Tests](https://github.com/mmmaxii/PA3Py/actions/workflows/tests.yml/badge.svg)](https://github.com/mmmaxii/PA3Py/actions)

**PA3Py** is a Python post-processing package for computing pebble accretion growth of planetary embryos from [TripodPy](https://github.com/tripod-code/tripodpy) hydrodynamic disk simulations. It is **not** a standalone disk evolution code — it reads TripodPy HDF5 output snapshots and computes embryo growth on top of them.

The physics follows Ormel (2017) and Drążkowska et al. (2023), including headwind/shear regime switching, 2D–3D turbulence transition, Safronov ballistic onset, and a dynamic isolation mass cap.

Please read the [documentation](https://pa3py.readthedocs.io) for a full description.

## Installation

Clone the repository and install via pip:

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
```

## Physics & References

- **Ormel (2017)**: *The Emerging Paradigm of Pebble Accretion*
- **Drążkowska et al. (2023)**: *Planet Formation Theory in the Era of ALMA and Kepler*

## Acknowledgements

PA3Py was developed at the Pontificia Universidad Católica de Chile. 
Contact: [mvalderrav@uc.cl](mailto:mvalderrav@uc.cl)
