"""
Pruebas para el cálculo de la snowline (Oka et al. 2011 & Hartmann).
"""

import os
import sys
import numpy as np

# Permitir correr el test sin instalar el paquete (pip install -e .)
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'src')))

from pa3py.snowline import get_rsnow_from_mdot_au, mdot_time, r_snow_time_cgs, generate_rsnow_array
from pa3py import constants as c

def test_mdot_time():
    """Prueba que la tasa de acreción de Hartmann decaiga en el tiempo."""
    print("Testing mdot_time()...")
    m_early = mdot_time(1.0, eta=1.5) # 1 Myr
    m_late = mdot_time(3.0, eta=1.5)  # 3 Myr
    assert m_early > m_late
    assert np.isclose(m_early, 1e-8)

def test_r_snow_time_cgs():
    """Prueba la evolución de la snowline a través del tiempo."""
    print("Testing r_snow_time_cgs()...")
    # A t=0 o temprano, debe limitarse al máximo de 2.73 AU por el Mdot_cap
    r_cgs_early = r_snow_time_cgs(1e5 * c.YEAR)
    assert r_cgs_early / c.AU <= 2.73
    
    # A t muy tardío (ej. 10 Myr), debe haber migrado hacia adentro, 
    # pero no bajar del r_min_au (ej 0.5 AU)
    r_cgs_late = r_snow_time_cgs(10e6 * c.YEAR, r_min_au=0.5)
    assert r_cgs_late / c.AU >= 0.5
    assert r_cgs_late / c.AU < 2.73

def test_generate_rsnow_array():
    """Prueba la generación de arreglos dinámicos para el módulo de acreción."""
    print("Testing generate_rsnow_array()...")
    times_sec = np.array([1e5 * c.YEAR, 2e6 * c.YEAR, 10e6 * c.YEAR])
    rsnow_arr = generate_rsnow_array(times_sec, r_min_au=0.7)
    
    assert len(rsnow_arr) == 3
    # Debe migrar hacia adentro
    assert rsnow_arr[0] > rsnow_arr[1]
    assert rsnow_arr[1] >= rsnow_arr[2]
    # Límite interior
    assert rsnow_arr[2] >= 0.7 * c.AU

if __name__ == "__main__":
    print("="*60)
    print("RUNNING SNOWLINE TESTS")
    print("="*60)
    test_mdot_time()
    test_r_snow_time_cgs()
    test_generate_rsnow_array()
    print("[OK] Todos los tests de la Snowline pasaron correctamente.")
