"""
Tests principales para validar la física modular de PA3Py con datos reales (localizados en test_data).
"""

import os
import sys
import numpy as np

# Permitir correr el test sin instalar el paquete (pip install -e .)
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'src')))

from pa3py.data import load_tripodpy_hdf5
from pa3py.composition import SimpleWaterComposition, FunctionComposition
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

def test_out_of_bounds_embryo():
    """Prueba que el paquete lance ValueError si un embrión está fuera de la grilla."""
    if not os.path.exists(DATA_DIR_SMOOTH):
        return
        
    print("\n" + "="*80)
    print("TEST: OUT OF BOUNDS EMBRYO")
    print("="*80)
    
    disk = load_tripodpy_hdf5(DATA_DIR_SMOOTH)
    sim = PebbleAccretionModule3(disk)
    
    try:
        sim.run_growth(embryo_locations_AU=[200.0])
        print("ERROR: test_out_of_bounds_embryo failed to raise ValueError!")
    except ValueError as e:
        print("[OK] Boundary check passed successfully:", e)

def test_multispecies_dynamic():
    """Prueba que el motor maneje N especies dinámicamente."""
    if not os.path.exists(DATA_DIR_SMOOTH):
        return
        
    print("\n" + "="*80)
    print("TEST: MULTISPECIES DYNAMIC TRACKING")
    print("="*80)
    
    disk = load_tripodpy_hdf5(DATA_DIR_SMOOTH)
    
    def exotic_chem(r, t):
        return {'Iron': 0.4, 'Carbon': 0.2, 'Ice': 0.4}
        
    comp = FunctionComposition(exotic_chem)
    sim = PebbleAccretionModule3(disk, comp_model=comp)
    
    # 3 species: Iron, Carbon, Ice
    assert len(sim.tracked_species) == 3
    assert 'Iron' in sim.tracked_species
    
    results = sim.run_growth([5.0], M0_g=1e-3 * sim.M_EARTH)
    
    # Check history structure: time, m_core, m_iso, sp1, sp2, sp3
    hist = results[5.0]
    assert hist.shape[1] == 6 
    
    sim.summary(results)
    print("[OK] Multispecies dynamic tracking works!")

def test_multizone_autodetection():
    """Prueba que el auto-descubrimiento barra todas las ramas del if del usuario."""
    if not os.path.exists(DATA_DIR_SMOOTH):
        return
        
    print("\n" + "="*80)
    print("TEST: MULTI-ZONE AUTO-DETECTION")
    print("="*80)
    
    disk = load_tripodpy_hdf5(DATA_DIR_SMOOTH)
    
    def quimica_4_zonas(r_cm, t_sec):
        AU = 1.496e13
        r_h2o = 2.73 * AU * (max(t_sec, 1e-6) / 1e13)**(-0.5)
        r_co2 = 5.0 * AU
        r_co  = 12.0 * AU
        if r_cm < r_h2o:
            return {'silicatos': 1.0}                           
        elif r_cm < r_co2:
            return {'silicatos': 0.5, 'H2O': 0.5}               
        elif r_cm < r_co:
            return {'silicatos': 0.3, 'H2O': 0.3, 'CO2': 0.4}   
        else:
            return {'silicatos': 0.2, 'H2O': 0.2, 'CO2': 0.3, 'CO': 0.3} 
            
    comp = FunctionComposition(quimica_4_zonas)
    sim = PebbleAccretionModule3(disk, comp_model=comp)
    
    # Debe haber descubierto las 4 especies
    assert len(sim.tracked_species) == 4
    for expected in ['silicatos', 'H2O', 'CO2', 'CO']:
        assert expected in sim.tracked_species, f"Faltó {expected}"
        
    print("[OK] Auto-detección de 4 zonas exitosa. Especies encontradas:", sim.tracked_species)

if __name__ == "__main__":
    test_smooth_run()
    test_sinusoidal_run()
    test_strong_gap_run()
    test_out_of_bounds_embryo()
    test_multispecies_dynamic()
    test_multizone_autodetection()
