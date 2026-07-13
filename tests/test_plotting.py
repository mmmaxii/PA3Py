import pytest
import os
import matplotlib
matplotlib.use('Agg')  # Usar backend no interactivo para evitar fallos en CI
import matplotlib.pyplot as plt

from pa3py import PA3Py

def test_plot_hovmoller(tmp_path):
    data_dir = 'tests/test_data/run_smooth_a0.001_v10'
    if not os.path.exists(data_dir):
        return
        
    sim = PA3Py(data_dir)
    fig, ax = sim.plot_hovmoller(field='dust_Sigma')
    assert fig is not None
    assert ax is not None
    
    out_file = tmp_path / "hovmoller.png"
    fig.savefig(out_file)
    assert out_file.exists()
    plt.close(fig)

def test_plot_population(tmp_path):
    data_dir = 'tests/test_data/run_smooth_a0.001_v10'
    if not os.path.exists(data_dir):
        return
        
    sim = PA3Py(data_dir)
    
    # Hacemos una simulación corta para el test
    results = sim.run_growth([1.0, 3.0])
    
    fig, ax = sim.plot_population(results)
    assert fig is not None
    assert ax is not None
    
    out_file = tmp_path / "population.png"
    fig.savefig(out_file)
    assert out_file.exists()
    plt.close(fig)
