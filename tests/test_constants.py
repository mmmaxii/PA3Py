"""
Tests de sanidad para las constantes físicas de PA3Py.
No requieren datos HDF5.
"""

import sys
import os
import math

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'src')))

from pa3py import constants as c


def test_constants_physical_values():
    """Las constantes tienen valores dentro del 0.1% del valor NIST/IAU."""
    assert abs(c.AU     - 1.496e13) / 1.496e13 < 1e-3, "AU incorrecto"
    assert abs(c.G      - 6.674e-8) / 6.674e-8 < 1e-3, "G incorrecto"
    assert abs(c.M_SUN  - 1.989e33) / 1.989e33 < 1e-3, "M_SUN incorrecto"
    assert abs(c.M_EARTH- 5.972e27) / 5.972e27 < 1e-3, "M_EARTH incorrecto"
    assert abs(c.YEAR   - 3.156e7)  / 3.156e7  < 1e-3, "YEAR incorrecto"


def test_constants_dimensional_consistency():
    """G*M_sun/AU^3 debe reproducir T_orbital ≈ 1 año en r=1 AU (Kepler)."""
    omega_k = math.sqrt(c.G * c.M_SUN / c.AU**3)   # rad/s a 1 AU
    T_orbital = 2 * math.pi / omega_k                # segundos
    # Debe ser ≈ 1 año con ≤ 1% de error
    assert abs(T_orbital - c.YEAR) / c.YEAR < 0.01, \
        f"Inconsistencia dimensional: T={T_orbital:.3e} s ≠ {c.YEAR:.3e} s"
