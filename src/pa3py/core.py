"""
core.py - Interfaz principal (Facade).
"""

from typing import List, Callable, Optional
import h5py

from .data import load_tripodpy_hdf5
from .composition import CompositionModel, SimpleWaterComposition, FunctionComposition
from .pebble_accretion import PebbleAccretionModule3
from .snowline import generate_rsnow_array
from .plotting import plot_hovmoller
from . import constants as c

class PA3Py:
    """
    Centraliza flujo de PA3Py: Datos HDF5, química, acreción y gráficas.
    """
    def __init__(self, data_dir: str, comp_model: Optional[CompositionModel] = None):
        """
        Inicializa la simulación y carga los datos.
        
        Parámetros:
        -----------
        data_dir : str
            Ruta a archivos HDF5.
        comp_model : CompositionModel, opcional
            Modelo de composición. Por defecto: Agua simple + migración snowline.
        """
        self.data_dir = data_dir
        self.disk = load_tripodpy_hdf5(data_dir)
        self.comp = comp_model or SimpleWaterComposition(generate_rsnow_array(self.disk.times))
        self.engine = PebbleAccretionModule3(self.disk, comp_model=self.comp)
        
    def set_custom_chemistry(self, user_func: Callable, species: Optional[List[str]] = None):
        """
        Atajo para redefinir la química con una función Python.
        
        Ejemplo:
            sim.set_custom_chemistry(mi_funcion, ["silicatos", "H2O"])
        """
        self.comp = FunctionComposition(user_func, species)
        self.engine = PebbleAccretionModule3(self.disk, comp_model=self.comp)
        
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
        results = self.engine.run_growth(embryos_au, M0_g=m_seed_me * c.M_EARTH)
        self.engine.summary(results)
        return results

    def plot_hovmoller(self, field: str = 'dust_Sigma', show_snowlines: bool = True, **kwargs):
        """
        Genera el diagrama de Hovmöller (Radio vs Tiempo).
        field puede ser 'dust_Sigma', 'gas_Sigma', o 'epsilon'.
        """
        return plot_hovmoller(self.disk, field=field, show_snowlines=show_snowlines, **kwargs)
        
    def calculate_isolation_mass_map(self):
        """Calcula el mapa teórico de masa de aislamiento en todo el disco."""
        return self.engine.calculate_isolation_mass_map()

    def save_results(self, results: dict, filename: str):
        """Guarda tracks de crecimiento en HDF5."""
        with h5py.File(filename, 'w') as f:
            # Encode species como bytes para compatibilidad h5py
            f.attrs['tracked_species'] = [s.encode('ascii') for s in self.engine.tracked_species]
            f.attrs['M_EARTH'] = c.M_EARTH
            
            tracks = f.create_group('tracks')
            for r_au, track_data in results.items():
                dset = tracks.create_dataset(f"embryo_{r_au:.4f}_AU", data=track_data)
                dset.attrs['r_au'] = float(r_au)

    @staticmethod
    def load_results(filename: str) -> tuple[dict, list]:
        """Carga tracks de crecimiento desde HDF5. Retorna (results_dict, species_list)."""
        results = {}
        with h5py.File(filename, 'r') as f:
            raw_species = list(f.attrs['tracked_species'])
            tracked_species = [s.decode('ascii') if isinstance(s, bytes) else s for s in raw_species]
            
            tracks = f['tracks']
            for key in tracks.keys():
                r_au = float(tracks[key].attrs['r_au'])
                results[r_au] = tracks[key][:]
                
        return results, tracked_species
