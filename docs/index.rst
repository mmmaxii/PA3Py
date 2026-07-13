Welcome to PA3Py's documentation!
===================================

**PA3Py** is a Python post-processing package for computing pebble accretion growth of planetary embryos from `TripodPy <https://github.com/tripod-code/tripodpy>`_ hydrodynamic disk simulations. It is **not** a standalone disk evolution code — it reads TripodPy HDF5 output snapshots and computes embryo growth on top of them.

The physics follows Ormel (2017) and Drążkowska et al. (2023), including headwind/shear regime switching, 2D–3D turbulence transition, Safronov ballistic onset, and a dynamic isolation mass cap.

.. toctree::
   :maxdepth: 2
   :caption: Contents:
   
   notebooks/00_Fisica_Implementada
   notebooks/01_Primeros_Pasos
   notebooks/02_Quimica_Avanzada
   notebooks/03_Poblaciones_Sinteticas_e_IO
   api

Installation
------------

Clone the repository and install via pip:

.. code-block:: bash

    git clone https://github.com/mmmaxii/PA3Py.git
    cd PA3Py
    pip install .

Requirements
------------

* Python ≥ 3.9
* `TripodPy <https://github.com/tripod-code/tripodpy>`_ simulation output (HDF5 snapshots)
* ``numpy``, ``scipy``, ``h5py``, ``matplotlib``
