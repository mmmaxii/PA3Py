Welcome to PA3Py's documentation!
===================================

**PA3Py** (Pebble Accretion 3 in Python) is a post-processing framework for computing the growth and composition of planetary embryos via pebble accretion on top of `TripodPy <https://github.com/tripod-code/tripodpy>`_ hydrodynamic disk simulations. It is **not** a standalone disk evolution code — it reads TripodPy HDF5 snapshots and integrates embryo growth over the fully evolved gas and dust structure, including structured disks with planet-carved gaps and pressure bumps.

The pebble accretion physics follows Ormel (2017) and Drążkowska et al. (2023) — headwind/shear regime switching, 2D–3D turbulence transition, Safronov ballistic onset, and a dynamic isolation mass cap (Lambrechts & Johansen 2014). The water snowline migrates following the radiative disk models of Oka et al. (2011) coupled to the Hartmann et al. (1998) stellar accretion-rate decay, and the chemistry module is fully composition-agnostic: any number of species and radial zones can be injected as a plain Python function.

PA3Py was developed as part of the undergraduate thesis *"Regulación del Flujo de Pebbles por Subestructuras en Discos Protoplanetarios: Implicancias para el Crecimiento Planetario y la Formación de Waterworlds"* at the Instituto de Astrofísica, Pontificia Universidad Católica de Chile.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   notebooks/00_Implemented_Physics
   notebooks/01_Getting_Started
   notebooks/02_Advanced_Chemistry
   notebooks/03_Synthetic_Populations_and_IO
   api

Installation
------------

Install the latest stable release from PyPI:

.. code-block:: bash

    pip install pa3py

Or clone the repository and install from source:

.. code-block:: bash

    git clone https://github.com/mmmaxii/PA3Py.git
    cd PA3Py
    pip install .

Requirements
------------

* Python ≥ 3.9
* `TripodPy <https://github.com/tripod-code/tripodpy>`_ simulation output (HDF5 snapshots)
* ``numpy``, ``scipy``, ``h5py``, ``matplotlib``
