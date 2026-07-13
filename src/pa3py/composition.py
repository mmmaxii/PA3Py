"""
composition.py - Modelos de composición química.
"""

from typing import Dict, List, Optional, Callable
import numpy as np
from abc import ABC, abstractmethod

class CompositionModel(ABC):
    """Clase base para modelos de composición."""
    
    @abstractmethod
    def get_fractions(self, r: float, t_sec: float, t_idx: int) -> Dict[str, float]:
        """
        Dado un radio orbital `r` (en cm), el tiempo físico `t_sec` (s), y el índice temporal `t_idx`, 
        devuelve un diccionario con las fracciones de masa de cada elemento químico que se está acretando. 
        Las fracciones deben sumar 1.0.
        """
        pass

    @abstractmethod
    def get_species(self) -> List[str]:
        """Devuelve la lista de nombres de las especies químicas rastreadas."""
        pass


class SimpleWaterComposition(CompositionModel):
    """
    Modelo clásico PA3Py:
    - 100% silicatos interior a la snowline.
    - 50% silicatos, 50% H2O exterior a la snowline.
    """
    def __init__(self, rsnow_h2o_array: np.ndarray):
        self.rsnow = rsnow_h2o_array
        
    def get_fractions(self, r: float, t_sec: float, t_idx: int) -> Dict[str, float]:
        r_snow = float(self.rsnow[t_idx]) if not np.isnan(self.rsnow[t_idx]) else 0.0
        
        if r >= r_snow:
            return {"silicates": 0.5, "H2O": 0.5}
        else:
            return {"silicates": 1.0, "H2O": 0.0}
            
    def get_species(self) -> List[str]:
        return ["silicates", "H2O"]


class MultiSnowlineComposition(CompositionModel):
    """Zonas estáticas separadas por snowlines pre-calculadas."""
    def __init__(self, snowlines: Dict[str, np.ndarray], zone_abundances: List[Dict[str, float]]):
        self.snowline_names = list(snowlines.keys())
        self.snowlines = snowlines
        self.zone_abundances = zone_abundances
        
        if len(zone_abundances) != len(snowlines) + 1:
            raise ValueError("Debes proveer exactamente N+1 diccionarios de abundancias para N snowlines.")
            
        self.species = list(set().union(*(abund.keys() for abund in self.zone_abundances)))

    def get_fractions(self, r: float, t_sec: float, t_idx: int) -> Dict[str, float]:
        zone_idx = 0
        for name in self.snowline_names:
            r_snow = float(self.snowlines[name][t_idx]) if not np.isnan(self.snowlines[name][t_idx]) else 0.0
            if r >= r_snow:
                zone_idx += 1
            else:
                break
                
        raw_fracs = self.zone_abundances[zone_idx]
        total = sum(raw_fracs.values())
        if total > 0:
            return {k: v / total for k, v in raw_fracs.items()}
        return {sp: 0.0 for sp in self.species}
        
    def get_species(self) -> List[str]:
        return self.species


class FunctionComposition(CompositionModel):
    """Química dinámica definida por el usuario mediante una función Python."""
    def __init__(self, user_func: Callable[[float, float], Dict[str, float]], species: Optional[List[str]] = None):
        """
        Parámetros:
        -----------
        user_func: callable
            Función del usuario. Debe recibir `r` (en cm) y `t` (en segundos)
            y retornar un diccionario con las fracciones, ej: {"silicates": 0.4, "H2O": 0.6}.
        species: List[str], opcional
            Lista explícita de las especies químicas. Si no se provee, PA3Py la deducirá
            automáticamente ejecutando la función con un valor de prueba (r=1.0, t=0.0).
        """
        self.user_func = user_func
        
        if species is None:
            try:
                res_sets = [set(self.user_func(r, t).keys()) 
                            for r in [1e11, 1e12, 1e13, 1e14, 1e15, 1e16] 
                            for t in [0.0, 1e10, 1e12, 1e14]]
                self.species = list(set.union(*res_sets)) if res_sets else []
                if not self.species:
                    raise ValueError("La función devolvió diccionarios vacíos en todos los tests.")
            except AttributeError:
                raise ValueError("La función debe retornar un dict (con el método .keys()).")
            except Exception as e:
                raise RuntimeError(f"Fallo en la auto-detección de especies: {e}\nProvee 'species' explícitamente.")
        else:
            self.species = species
        
    def get_fractions(self, r: float, t_sec: float, t_idx: int) -> Dict[str, float]:
        raw_fracs = self.user_func(r, t_sec)
        total = sum(raw_fracs.values())
        if total > 0:
            # Normalizar y devolver, rellenando con ceros las especies que no aparezcan esta vez
            norm_fracs = {sp: 0.0 for sp in self.species}
            for k, v in raw_fracs.items():
                if k not in norm_fracs:
                    # Si la función retornó una especie nueva dinámicamente, la ignoramos o fallamos.
                    # Aquí la ignoramos pero lo ideal es detectarlas todas al inicio.
                    pass
                else:
                    norm_fracs[k] = v / total
            return norm_fracs
        return {sp: 0.0 for sp in self.species}

    def get_species(self) -> List[str]:
        return self.species
