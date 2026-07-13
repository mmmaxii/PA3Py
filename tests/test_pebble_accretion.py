"""
Tests principales para validar la física modular de PA3Py con datos reales.
"""

import os
import sys
import numpy as np

# Permitir correr el test sin instalar el paquete (pip install -e .)
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'src')))

from pa3py.data import load_tripodpy_hdf5
from pa3py.composition import SimpleWaterComposition
from pa3py.pebble_accretion import PebbleAccretionModule3
from pa3py.snowline import generate_rsnow_array

# Rutas de datos (asumiendo ejecución local en la máquina de desarrollo)
DATA_DIR_SMOOTH = r"C:\astro\Codigos practica + docs + papers\codigos\data\runs\vf_10ms\smooth\run_smooth_a0.001_v10"
DATA_DIR_SIN = r"C:\astro\Codigos practica + docs + papers\codigos\data\runs\vf_10ms\sinusoidal\run_ngap5_A1.0_a0.001_rmin0.7"

def test_smooth_run():
    """Prueba el módulo sobre un disco suave (alpha 0.001, v_frag 10ms)."""
    if not os.path.exists(DATA_DIR_SMOOTH):
        print(f"Data not found: {DATA_DIR_SMOOTH}")
        return

    print("\n" + "="*80)
    print("TEST: RUN SMOOTH (alpha=0.001, vf=10ms)")
    print("="*80)

    # 1. Cargar datos
    disk = load_tripodpy_hdf5(DATA_DIR_SMOOTH)

    # 2. Generar snowline dinámica (Oka + Hartmann)
    rsnow_cgs = generate_rsnow_array(disk.times)
    comp = SimpleWaterComposition(rsnow_cgs)

    # 3. Correr acreción
    sim = PebbleAccretionModule3(disk, comp_model=comp)
    
    # 4. Probar en un rango de posiciones (1, 3, 5, 10, 15 AU)
    embryos = [1.0, 3.0, 5.0, 10.0, 15.0]
    results = sim.run_growth(embryos, M0_g=1e-3 * sim.M_EARTH)

    sim.summary(results)
    
    # Validaciones básicas
    for r in embryos:
        assert r in results
        hist = results[r]
        if len(hist) > 0:
            assert hist[-1][1] > 0 # Masa final mayor a cero


def test_sinusoidal_run():
    """Prueba el módulo sobre un disco con 5 gaps sinusoidales de amplitud 1.0 (alpha 0.001, v_frag 10ms)."""
    if not os.path.exists(DATA_DIR_SIN):
        print(f"Data not found: {DATA_DIR_SIN}")
        return

    print("\n" + "="*80)
    print("TEST: RUN SINUSOIDAL (5 gaps, A=1.0, alpha=0.001, vf=10ms)")
    print("="*80)

    # 1. Cargar datos
    disk = load_tripodpy_hdf5(DATA_DIR_SIN)

    # 2. Generar snowline dinámica (Oka + Hartmann)
    rsnow_cgs = generate_rsnow_array(disk.times)
    comp = SimpleWaterComposition(rsnow_cgs)

    # 3. Correr acreción
    sim = PebbleAccretionModule3(disk, comp_model=comp)
    
    # 4. Probar embriones en los gaps (aprox 3, 5, 7, 10 AU)
    embryos = [3.0, 5.0, 7.0, 10.0]
    results = sim.run_growth(embryos, M0_g=1e-3 * sim.M_EARTH)

    sim.summary(results)

    # Validaciones básicas
    for r in embryos:
        assert r in results
        
if __name__ == "__main__":
    test_smooth_run()
    test_sinusoidal_run()
