"""
Tests principales para validar la física modular de PA3Py con datos reales (localizados en test_data).
"""

import os
import sys
import numpy as np

# Permitir correr el test sin instalar el paquete (pip install -e .)
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'src')))

from pa3py import PebbleAccretionModule3, load_tripodpy_hdf5, FunctionComposition, SimpleWaterComposition, generate_rsnow_array
from pa3py import constants as c

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
    results = sim.run_growth(embryos, M0_g=1e-3 * c.M_EARTH)

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
    results = sim.run_growth(embryos, M0_g=1e-3 * c.M_EARTH)
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
    results = sim.run_growth(embryos, M0_g=1e-3 * c.M_EARTH)
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
    
    results = sim.run_growth([5.0], M0_g=1e-3 * c.M_EARTH)
    
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
        r_h2o = 2.73 * c.AU * (max(t_sec, 1e-6) / 1e13)**(-0.5)
        r_co2 = 5.0 * c.AU
        r_co  = 12.0 * c.AU
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


def test_history_full_length():
    """Todos los embriones retornan exactamente Nt filas, incluso post M_iso."""
    if not os.path.exists(DATA_DIR_SMOOTH):
        return

    disk = load_tripodpy_hdf5(DATA_DIR_SMOOTH)
    sim  = PebbleAccretionModule3(disk, comp_model=SimpleWaterComposition(generate_rsnow_array(disk.times)))
    results = sim.run_growth([3.0, 10.0], M0_g=1e-3 * c.M_EARTH)

    for r_au, hist in results.items():
        assert hist.shape[0] == disk.Nt, \
            f"Embrión {r_au} AU: {hist.shape[0]} filas, esperaba {disk.Nt}"


def test_isolation_mass_saturation():
    """M_core nunca supera M_iso en ningún timestep."""
    if not os.path.exists(DATA_DIR_SMOOTH):
        return

    disk = load_tripodpy_hdf5(DATA_DIR_SMOOTH)
    sim  = PebbleAccretionModule3(disk, comp_model=SimpleWaterComposition(generate_rsnow_array(disk.times)))
    results = sim.run_growth([3.0, 5.0, 10.0], M0_g=1e-3 * c.M_EARTH)

    for r_au, hist in results.items():
        M_core = hist[:, 1]
        M_iso  = hist[:, 2]
        overshoot = (M_core - M_iso).max()
        assert overshoot <= 1e-20, \
            f"Embrión {r_au} AU: M_core supera M_iso en {overshoot:.2e} g"


def test_isolation_mass_map_shape():
    """calculate_isolation_mass_map retorna shape (Nt, Nr) con todos los valores positivos."""
    if not os.path.exists(DATA_DIR_SMOOTH):
        return

    disk     = load_tripodpy_hdf5(DATA_DIR_SMOOTH)
    sim      = PebbleAccretionModule3(disk)
    M_iso_map = sim.calculate_isolation_mass_map()

    assert M_iso_map.shape == (disk.Nt, disk.Nr), \
        f"Shape incorrecto: {M_iso_map.shape} != ({disk.Nt}, {disk.Nr})"
    assert np.all(M_iso_map > 0), "M_iso_map contiene valores no positivos"


def test_hovmoller_smoke():
    """plot_hovmoller corre sin errores en los tres modos de campo."""
    if not os.path.exists(DATA_DIR_SMOOTH):
        return

    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from pa3py import PA3Py

    sim = PA3Py(DATA_DIR_SMOOTH)
    for field in ['dust_Sigma', 'gas_Sigma', 'epsilon']:
        fig, ax = sim.plot_hovmoller(field=field, show_snowlines=False)
        assert fig is not None
        assert ax  is not None
        plt.close(fig)


if __name__ == "__main__":
    test_smooth_run()
    test_sinusoidal_run()
    test_strong_gap_run()
    test_out_of_bounds_embryo()
    test_multispecies_dynamic()
    test_multizone_autodetection()
    test_history_full_length()
    test_isolation_mass_saturation()
    test_isolation_mass_map_shape()
    test_hovmoller_smoke()
