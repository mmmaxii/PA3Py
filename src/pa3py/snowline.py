"""
snowline.py - Cálculo analítico de la posición de la snowline en base a Oka et al. (2011) y Hartmann.
"""

import os
import numpy as np
from scipy.interpolate import interp1d
from . import constants as c

# ============================================================
# INICIALIZACIÓN DE LA INTERPOLACIÓN (Ejecutada al importar)
# ============================================================

# Ruta robusta al CSV (empaquetado junto al módulo)
_current_dir = os.path.dirname(__file__)
_csv_path = os.path.join(_current_dir, "data_files", "oka_et_al_2011_data_corrected.csv")

# Leer datos velozmente con NumPy
data = np.loadtxt(_csv_path, delimiter=';', skiprows=1)
_mdot_raw, _rsnow_raw = data[:, 0], data[:, 1]

# Asegurar orden monotónico para la interpolación log-log
_sort_idx = np.argsort(_mdot_raw)
_mdot_sorted = _mdot_raw[_sort_idx]
_rsnow_sorted = _rsnow_raw[_sort_idx]

_mdot_unique, _unique_indices = np.unique(_mdot_sorted, return_index=True)
_rsnow_unique = _rsnow_sorted[_unique_indices]

_log_mdot = np.log10(_mdot_unique)
_log_rsnow = np.log10(_rsnow_unique)

# Interpolador principal
_interp_log_rsnow = interp1d(
    _log_mdot, 
    _log_rsnow, 
    kind='linear', 
    bounds_error=False, 
    fill_value="extrapolate"
)

# ============================================================
# FUNCIONES FÍSICAS PÚBLICAS
# ============================================================

def get_rsnow_from_mdot_au(mdot_val):
    """Retorna la posición de snowline en AU dada una tasa de acreción."""
    log_m = np.log10(mdot_val)
    log_r = _interp_log_rsnow(log_m)
    return 10**log_r

def mdot_time(t_myr, eta=1.5):
    """Calcula la tasa de acreción de gas en el tiempo (Hartmann)."""
    mdot_1myr = 1e-8
    # Proteger contra tiempo 0 para evitar singularidades matemáticas
    t_myr = np.maximum(t_myr, 1e-6) 
    return mdot_1myr * (t_myr)**(-eta)

def r_snow_time_cgs(t_sec, eta=1.5, r_min_au=0.5):
    """
    Dada el tiempo en segundos (sim.t), calcula la posición de la snowline
    en centímetros (CGS), aplicando un corte inferior estricto dictado por la grilla.
    """
    # 1. Convertir t_sec a Myr
    t_myr = t_sec / (1e6 * c.YEAR)
    
    # 2. Obtener Mdot
    m = mdot_time(t_myr, eta)
    
    # 3. Obtener r_snow en AU
    # Si Mdot > 1e-7, capamos Mdot para que r_snow no exceda ~2.73 AU
    m_eff = min(m, 1e-7)
    r_au = get_rsnow_from_mdot_au(m_eff)
    
    # 4. Forzar límite inferior dinámico (r_min de la grilla) y límite superior de 2.73 AU
    r_au_limited = np.clip(r_au, r_min_au, 2.73)
    
    # 5. Retornar en CGS
    return float(r_au_limited * c.AU)

def generate_rsnow_array(times_sec: np.ndarray, eta: float = 1.5, r_min_au: float = 0.5) -> np.ndarray:
    """
    Genera un arreglo de posiciones de la snowline en función del tiempo.
    Ideal para inyectar en los CompositionModels.
    """
    return np.array([r_snow_time_cgs(t, eta=eta, r_min_au=r_min_au) for t in times_sec])
