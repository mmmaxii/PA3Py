"""
Tests principales para validar la física modular de PA3Py con datos reales (localizados en test_data).
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

# Rutas locales relativas a la carpeta del test
_current_dir = os.path.dirname(__file__)
DATA_DIR_SMOOTH = os.path.join(_current_dir, "test_data", "run_smooth_a0.001_v10")
DATA_DIR_SIN = os.path.join(_current_dir, "test_data", "run_ngap5_A1.0_a0.001_rmin0.7")
DATA_DIR_GAP = os.path.join(_current_dir, "test_data", "run_r10.0_m0.01_a0.001")

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

def test_sinusoidal_run():
    """Prueba el módulo sobre un disco con 5 gaps sinusoidales de amplitud 1.0."""
    if not os.path.exists(DATA_DIR_SIN):
        print(f"Data not found: {DATA_DIR_SIN}")
        return

    print("\n" + "="*80)
    print("TEST: RUN SINUSOIDAL (5 gaps, A=1.0, alpha=0.001, vf=10ms)")
    print("="*80)

    disk = load_tripodpy_hdf5(DATA_DIR_SIN)
    comp = SimpleWaterComposition(generate_rsnow_array(disk.times))
    sim = PebbleAccretionModule3(disk, comp_model=comp)
    
    embryos = [3.0, 5.0, 7.0, 10.0]
    results = sim.run_growth(embryos, M0_g=1e-3 * sim.M_EARTH)
    sim.summary(results)

def test_strong_gap_run():
    """Prueba el módulo sobre un disco con un gap fuerte en 10 AU."""
    if not os.path.exists(DATA_DIR_GAP):
        print(f"Data not found: {DATA_DIR_GAP}")
        return

    print("\n" + "="*80)
    print("TEST: RUN STRONG GAP (r=10 AU, masa=0.01)")
    print("="*80)

    disk = load_tripodpy_hdf5(DATA_DIR_GAP)
    comp = SimpleWaterComposition(generate_rsnow_array(disk.times))
    sim = PebbleAccretionModule3(disk, comp_model=comp)
    
    embryos = [8.0, 10.0, 12.0]
    results = sim.run_growth(embryos, M0_g=1e-3 * sim.M_EARTH)
    sim.summary(results)

if __name__ == "__main__":
    test_smooth_run()
    test_sinusoidal_run()
    test_strong_gap_run()
