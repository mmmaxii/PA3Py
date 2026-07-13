"""
core.py - Interfaz principal unificada (Facade) para PA3Py.
"""

from typing import List, Callable, Optional
from .data import load_tripodpy_hdf5
from .composition import CompositionModel, SimpleWaterComposition, FunctionComposition
from .pebble_accretion import PebbleAccretionModule3
from .snowline import generate_rsnow_array
from .plotting import plot_hovmoller

class PA3Py:
    """
    Clase principal que centraliza el flujo de trabajo de PA3Py:
    Carga de datos HDF5, configuración química, simulación de acreción y gráficas.
    """
    def __init__(self, data_dir: str, comp_model: Optional[CompositionModel] = None):
        """
        Inicializa la simulación cargando los datos en memoria.
        
        Parámetros:
        -----------
        data_dir : str
            Ruta al directorio que contiene los archivos HDF5 de la simulación del disco.
        comp_model : CompositionModel, opcional
            Modelo de composición. Si no se provee, asume el modelo clásico de 
            agua con migración de snowline (Oka + Hartmann).
        """
        self.data_dir = data_dir
        self.disk = load_tripodpy_hdf5(data_dir)
        
        if comp_model is None:
            rsnow = generate_rsnow_array(self.disk.times)
            self.comp = SimpleWaterComposition(rsnow)
        else:
            self.comp = comp_model
            
        self._init_engine()

    def _init_engine(self):
        """Inicializa el motor de física subyacente."""
        self.engine = PebbleAccretionModule3(self.disk, comp_model=self.comp)
        
    def set_custom_chemistry(self, user_func: Callable, species: List[str]):
        """
        Atajo para redefinir la química de la simulación usando una función de Python.
        
        Ejemplo:
            sim.set_custom_chemistry(mi_funcion, ["silicatos", "H2O", "CO2"])
        """
        self.comp = FunctionComposition(user_func, species)
        self._init_engine()
        
    def run_growth(self, embryos_au: list, m_seed_me: float = 1e-3) -> dict:
        """
        Corre la simulación de acreción para los embriones dados.
        
        Parámetros:
        -----------
        embryos_au : list
            Lista de radios iniciales en AU (ej: [1.0, 5.0, 15.0]).
        m_seed_me : float
            Masa semilla en Masas Terrestres. Default: 1e-3.
            
        Retorna:
        --------
        dict
            Resultados de la evolución de masa en el tiempo.
        """
        results = self.engine.run_growth(embryos_au, M0_g=m_seed_me * self.engine.M_EARTH)
        self.engine.summary(results)
        return results

    def plot_hovmoller(self, field: str = 'dust_Sigma', show_snowlines: bool = True, **kwargs):
        """
        Genera el diagrama de Hovmöller (Radio vs Tiempo).
        field puede ser 'dust_Sigma', 'gas_Sigma', o 'eps'.
        """
        return plot_hovmoller(self.disk, field=field, show_snowlines=show_snowlines, **kwargs)
        
    def calculate_isolation_mass_map(self):
        """Calcula el mapa teórico de masa de aislamiento en todo el disco."""
        return self.engine.calculate_isolation_mass_map()
