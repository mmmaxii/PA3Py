"""
composition.py - Modelos para determinar la composición de la materia acretada.
"""

from typing import Dict, List, Optional
import numpy as np
from abc import ABC, abstractmethod

class CompositionModel(ABC):
    """Clase base para todos los modelos de composición química."""
    
    @abstractmethod
    def get_fractions(self, r: float, t_idx: int) -> Dict[str, float]:
        """
        Dado un radio orbital `r` (en cm) y un índice temporal `t_idx`, devuelve un diccionario 
        con las fracciones de masa de cada elemento químico que se está acretando. 
        Las fracciones deben sumar 1.0.
        """
        pass

    @abstractmethod
    def get_species(self) -> List[str]:
        """Devuelve la lista de nombres de las especies químicas rastreadas."""
        pass


class SimpleWaterComposition(CompositionModel):
    """
    Modelo por defecto clásico de PA3Py:
    - 100% silicatos dentro de la snowline de agua.
    - 50% silicatos, 50% H2O fuera de la snowline de agua.
    
    La snowline se provee como un arreglo dependiente del tiempo.
    """
    def __init__(self, rsnow_h2o_array: np.ndarray):
        """
        Parámetros:
        -----------
        rsnow_h2o_array : np.ndarray
            Arreglo 1D de tamaño (Nt,) con la posición radial de la snowline de agua en cm.
        """
        self.rsnow = rsnow_h2o_array
        
    def get_fractions(self, r: float, t_idx: int) -> Dict[str, float]:
        r_snow = float(self.rsnow[t_idx]) if not np.isnan(self.rsnow[t_idx]) else 0.0
        
        if r >= r_snow:
            return {"silicates": 0.5, "H2O": 0.5}
        else:
            return {"silicates": 1.0, "H2O": 0.0}
            
    def get_species(self) -> List[str]:
        return ["silicates", "H2O"]


class MultiSnowlineComposition(CompositionModel):
    """
    Modelo de composición personalizable por el usuario.
    Permite definir un número arbitrario de snowlines dependientes del tiempo,
    y especificar las abundancias en cada zona separada por ellas.
    """
    def __init__(self, snowlines: Dict[str, np.ndarray], zone_abundances: List[Dict[str, float]]):
        """
        Parámetros:
        -----------
        snowlines : Dict[str, np.ndarray]
            Diccionario de snowlines, ordenadas de MENOR radio a MAYOR radio.
            Ej: {"H2O": rsnow_h2o_array, "CO2": rsnow_co2_array}
        zone_abundances : List[Dict[str, float]]
            Lista de abundancias para las zonas delimitadas por las snowlines.
            Debe tener exactamente len(snowlines) + 1 elementos.
            - zone_abundances[0] = región interna a la primera snowline.
            - zone_abundances[-1] = región externa a la última snowline.
        """
        self.snowline_names = list(snowlines.keys())
        self.snowlines = snowlines
        self.zone_abundances = zone_abundances
        
        if len(zone_abundances) != len(snowlines) + 1:
            raise ValueError("Debes proveer exactamente N+1 diccionarios de abundancias para N snowlines.")
            
        # Extraer todas las especies mencionadas
        all_sps = set()
        for abund in self.zone_abundances:
            all_sps.update(abund.keys())
        self.species = list(all_sps)

    def get_fractions(self, r: float, t_idx: int) -> Dict[str, float]:
        # Encontrar en qué zona estamos
        zone_idx = 0
        for name in self.snowline_names:
            r_snow = float(self.snowlines[name][t_idx]) if not np.isnan(self.snowlines[name][t_idx]) else 0.0
            if r >= r_snow:
                zone_idx += 1
            else:
                break
                
        # Normalizar fracciones
        raw_fracs = self.zone_abundances[zone_idx]
        total = sum(raw_fracs.values())
        if total > 0:
            return {k: v / total for k, v in raw_fracs.items()}
        return {sp: 0.0 for sp in self.species}
        
    def get_species(self) -> List[str]:
        return self.species
