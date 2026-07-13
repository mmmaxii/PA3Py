"""
PA3Py - Pebble Accretion 3 in Python
"""

from .data import load_tripodpy_hdf5, DiskData
from .composition import CompositionModel, SimpleWaterComposition, MultiSnowlineComposition, FunctionComposition
from .pebble_accretion import PebbleAccretionModule3
from .snowline import generate_rsnow_array, get_rsnow_from_mdot_au, mdot_time, r_snow_time_cgs
from .plotting import plot_hovmoller
from .core import PA3Py
from . import constants

def easy_run(data_dir: str, embryos_au: list, m_seed_me: float = 1e-4) -> dict:
    """
    Ruta rápida para correr el simulador PA3Py con el modelo por defecto.
    
    Asume:
      - Snowline dinámica basada en Oka et al. (2011) y Hartmann.
      - Composición: 100% Silicatos (dentro de snowline) y 50% H2O / 50% Silicatos (fuera).
      
    Parámetros:
    -----------
    data_dir : str
        Ruta al directorio con los HDF5 de TripodPy.
    embryos_au : list
        Lista con las posiciones iniciales de los embriones en AU (ej: [3.0, 5.0, 10.0]).
    m_seed_me : float
        Masa semilla inicial de los embriones en Masas Terrestres.
        
    Retorna:
    --------
    dict
        Diccionario con la historia de crecimiento de cada embrión.
    """
    # 1. Cargar datos
    disk = load_tripodpy_hdf5(data_dir)
    
    # 2. Configurar snowline y composición clásica
    rsnow = generate_rsnow_array(disk.times)
    comp = SimpleWaterComposition(rsnow)
    
    # 3. Inicializar simulador y correr
    sim = PebbleAccretionModule3(disk, comp_model=comp)
    results = sim.run_growth(embryos_au, M0_g=m_seed_me * sim.M_EARTH)
    
    # 4. Mostrar resumen automáticamente
    sim.summary(results)
    return results

__all__ = [
    "load_tripodpy_hdf5",
    "DiskData",
    "CompositionModel",
    "SimpleWaterComposition",
    "MultiSnowlineComposition",
    "FunctionComposition",
    "PebbleAccretionModule3",
    "generate_rsnow_array",
    "easy_run",
    "plot_hovmoller",
    "PA3Py",
    "constants"
]
