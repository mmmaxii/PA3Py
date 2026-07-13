"""
Tests unitarios puros para los modelos de composición química.
No requieren datos HDF5.
"""

import sys
import os
import pytest
import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'src')))

from pa3py.composition import SimpleWaterComposition, FunctionComposition, MultiSnowlineComposition
from pa3py import constants as c


def test_simple_water_inside_outside():
    """SimpleWaterComposition: 100% silicatos dentro, 50/50 fuera de la snowline."""
    rsnow = np.full(5, 3.0 * c.AU)
    comp  = SimpleWaterComposition(rsnow)

    fracs_in = comp.get_fractions(1.0 * c.AU, 0.0, 0)
    assert fracs_in['silicates'] == 1.0
    assert fracs_in['H2O']       == 0.0

    fracs_out = comp.get_fractions(5.0 * c.AU, 0.0, 0)
    assert fracs_out['silicates'] == 0.5
    assert fracs_out['H2O']       == 0.5


def test_fractions_always_sum_to_one():
    """Las fracciones de los tres modelos suman exactamente 1.0 en todas las zonas."""
    rsnow = np.full(5, 3.0 * c.AU)

    # SimpleWaterComposition
    comp_s = SimpleWaterComposition(rsnow)
    for r_au in [1.0, 2.5, 3.5, 8.0]:
        fracs = comp_s.get_fractions(r_au * c.AU, 0.0, 0)
        assert abs(sum(fracs.values()) - 1.0) < 1e-12, f"SimpleWater: suma ≠ 1 en r={r_au} AU"

    # FunctionComposition (3 zonas)
    def quimica(r, t):
        if r < 3.0 * c.AU:
            return {'silicates': 1.0}
        elif r < 6.0 * c.AU:
            return {'silicates': 0.5, 'H2O': 0.5}
        else:
            return {'silicates': 0.3, 'H2O': 0.3, 'CO2': 0.4}

    comp_f = FunctionComposition(quimica)
    for r_au in [1.0, 4.0, 10.0]:
        fracs = comp_f.get_fractions(r_au * c.AU, 0.0, 0)
        assert abs(sum(fracs.values()) - 1.0) < 1e-12, f"FunctionComp: suma ≠ 1 en r={r_au} AU"

    # MultiSnowlineComposition
    snowlines = {
        'H2O': np.full(5, 3.0 * c.AU),
        'CO2': np.full(5, 6.0 * c.AU),
    }
    zones = [
        {'silicates': 1.0},
        {'silicates': 0.5, 'H2O': 0.5},
        {'silicates': 0.3, 'H2O': 0.3, 'CO2': 0.4},
    ]
    comp_m = MultiSnowlineComposition(snowlines, zones)
    for r_au in [1.0, 4.0, 10.0]:
        fracs = comp_m.get_fractions(r_au * c.AU, 0.0, 0)
        assert abs(sum(fracs.values()) - 1.0) < 1e-12, f"MultiSnowline: suma ≠ 1 en r={r_au} AU"


def test_function_composition_bad_return():
    """FunctionComposition lanza ValueError si la función no retorna un dict."""
    with pytest.raises(ValueError):
        FunctionComposition(lambda r, t: 42.0)


def test_function_composition_empty_dict():
    """FunctionComposition lanza RuntimeError si la función retorna siempre dict vacío."""
    with pytest.raises(RuntimeError):
        FunctionComposition(lambda r, t: {})


def test_multisnowline_wrong_zone_count():
    """MultiSnowlineComposition lanza ValueError si el número de zonas es incorrecto."""
    snowlines = {
        'H2O': np.full(5, 3.0 * c.AU),
        'CO2': np.full(5, 6.0 * c.AU),
    }
    too_few = [{'silicates': 1.0}, {'silicates': 0.5, 'H2O': 0.5}]  # necesita 3, tiene 2
    with pytest.raises(ValueError):
        MultiSnowlineComposition(snowlines, too_few)


def test_multisnowline_species_discovered():
    """MultiSnowlineComposition descubre todas las especies de todas las zonas."""
    snowlines = {'H2O': np.full(5, 3.0 * c.AU)}
    zones = [
        {'silicates': 1.0},
        {'silicates': 0.5, 'H2O': 0.5},
    ]
    comp = MultiSnowlineComposition(snowlines, zones)
    assert set(comp.get_species()) == {'silicates', 'H2O'}
