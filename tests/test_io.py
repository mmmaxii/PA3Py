import os
import sys
import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'src')))

from pa3py import PA3Py

DATA_DIR_SMOOTH = os.path.abspath(os.path.join(os.path.dirname(__file__), 'test_data', 'run_smooth_a0.001_v10'))

def test_hdf5_io():
    """Prueba que el módulo pueda guardar y cargar datos desde HDF5 de forma idéntica."""
    if not os.path.exists(DATA_DIR_SMOOTH):
        print("Faltan datos de prueba, saltando test HDF5 IO...")
        return
        
    print("\n" + "="*80)
    print("TEST: HDF5 INPUT / OUTPUT")
    print("="*80)
    
    sim = PA3Py(DATA_DIR_SMOOTH)
    
    # Usaremos una química custom para asegurar que lea bien la lista de especies
    def test_chem(r, t):
        return {'Iron': 0.3, 'Water': 0.7}
    
    sim.set_custom_chemistry(test_chem)
    
    # Corremos 3 embriones
    results_original = sim.run_growth([5.0, 10.0, 15.0])
    
    # Guardamos a HDF5
    test_file = "test_io_output.h5"
    if os.path.exists(test_file):
        os.remove(test_file)
        
    sim.save_results(results_original, test_file)
    print(f"[OK] Resultados guardados en {test_file}")
    
    # Cargamos desde HDF5
    results_loaded, loaded_species = PA3Py.load_results(test_file)
    print(f"[OK] Resultados cargados desde {test_file}")
    
    # Verificaciones
    assert loaded_species == ['Iron', 'Water'], f"Error: Especies cargadas {loaded_species} != ['Iron', 'Water']"
    assert len(results_loaded) == 3, "Error: Faltan embriones"
    assert 5.0 in results_loaded, "Error: Llave 5.0 no se preservó"
    
    # Verificamos que los numpy arrays sean idénticos bit a bit
    for r_au in results_original.keys():
        arr_orig = results_original[r_au]
        arr_load = results_loaded[r_au]
        
        # Comparación estricta de numpy arrays
        np.testing.assert_array_equal(arr_orig, arr_load)
        
    print("[OK] Las matrices HDF5 cargadas son EXACTAMENTE iguales bit a bit a la memoria RAM original.")
    
    # Limpiamos el archivo temporal
    os.remove(test_file)

if __name__ == "__main__":
    test_hdf5_io()
