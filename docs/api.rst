API Reference
=============

.. module:: pa3py

Core Simulation Interface
-------------------------

.. autoclass:: PA3Py
   :members:
   :undoc-members:
   :show-inheritance:

Physics Engine
--------------

.. autoclass:: pa3py.pebble_accretion.PebbleAccretionModule3
   :members:
   :undoc-members:
   :show-inheritance:

Data Handling
-------------

.. autofunction:: load_tripodpy_hdf5

.. autoclass:: pa3py.data.DiskData
   :members:
   :undoc-members:
   :show-inheritance:

Plotting
--------

.. autofunction:: pa3py.plot_hovmoller

.. autofunction:: pa3py.plotting.plot_population

.. autofunction:: pa3py.plot_growth_curves

.. autofunction:: pa3py.plot_species_fraction

.. autofunction:: pa3py.plot_growth_summary

Composition Models
------------------

.. autoclass:: CompositionModel
   :members:
   :undoc-members:
   :show-inheritance:

.. autoclass:: SimpleWaterComposition
   :members:
   :undoc-members:
   :show-inheritance:

.. autoclass:: MultiSnowlineComposition
   :members:
   :undoc-members:
   :show-inheritance:

.. autoclass:: FunctionComposition
   :members:
   :undoc-members:
   :show-inheritance:
