import pytest
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')  # Usar backend no interactivo para evitar fallos en CI
import matplotlib.pyplot as plt

from pa3py import PA3Py
from pa3py.plotting import plot_growth_curves, plot_species_fraction

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

def test_plot_growth_curves(tmp_path):
    data_dir = 'tests/test_data/run_smooth_a0.001_v10'
    if not os.path.exists(data_dir):
        return

    sim = PA3Py(data_dir)
    results = sim.run_growth([1.0, 3.0, 5.0])

    fig, ax = sim.plot_growth_curves(results)
    assert fig is not None
    assert ax is not None
    # La masa siempre en escala log, con M_iso como techo del eje
    assert ax.get_yscale() == 'log'

    out_file = tmp_path / "growth_curves.png"
    fig.savefig(out_file)
    assert out_file.exists()
    plt.close(fig)

def test_plot_species_fraction_default(tmp_path):
    data_dir = 'tests/test_data/run_smooth_a0.001_v10'
    if not os.path.exists(data_dir):
        return

    sim = PA3Py(data_dir)
    results = sim.run_growth([1.0, 3.0, 5.0])

    fig, ax = sim.plot_species_fraction(results)
    assert fig is not None
    assert ax is not None

    out_file = tmp_path / "species_fraction.png"
    fig.savefig(out_file)
    assert out_file.exists()
    plt.close(fig)

def test_plot_species_fraction_invalid_species_raises():
    data_dir = 'tests/test_data/run_smooth_a0.001_v10'
    if not os.path.exists(data_dir):
        return

    sim = PA3Py(data_dir)
    results = sim.run_growth([1.0])

    with pytest.raises(ValueError):
        sim.plot_species_fraction(results, species='not_a_species')

def test_plot_growth_summary(tmp_path):
    data_dir = 'tests/test_data/run_smooth_a0.001_v10'
    if not os.path.exists(data_dir):
        return

    sim = PA3Py(data_dir)
    results = sim.run_growth([1.0, 3.0])

    fig, axes = sim.plot_growth_summary(results)
    assert fig is not None
    assert len(axes) == 2

    out_file = tmp_path / "growth_summary.png"
    fig.savefig(out_file)
    assert out_file.exists()
    plt.close(fig)

def test_plot_growth_curves_many_embryos_colorbar(tmp_path):
    # Diccionario `results` sintético (sin correr la física real) para ejercitar
    # de forma barata la rama de colormap continuo / colorbar (>8 embriones),
    # y para verificar que las funciones funcionan standalone (sin PA3Py).
    tracked_species = ["silicates", "H2O"]
    results = {}
    for r_au in np.linspace(1.0, 10.0, 12):
        t = np.linspace(0, 1e6 * 3.15576e7, 20)
        m_core = np.linspace(1e26, 5e27, 20)
        m_iso = np.full(20, 6e27)
        sil = m_core * 0.7
        h2o = m_core * 0.3
        results[float(r_au)] = np.column_stack([t, m_core, m_iso, sil, h2o])

    fig, ax = plot_growth_curves(results)
    assert fig is not None
    assert ax is not None
    fig.savefig(tmp_path / "growth_curves_many.png")
    assert (tmp_path / "growth_curves_many.png").exists()
    plt.close(fig)

    fig2, ax2 = plot_species_fraction(results, tracked_species, species='H2O')
    assert fig2 is not None
    assert ax2 is not None
    fig2.savefig(tmp_path / "species_fraction_many.png")
    assert (tmp_path / "species_fraction_many.png").exists()
    plt.close(fig2)

def test_plot_species_fraction_multi_species(tmp_path):
    data_dir = 'tests/test_data/run_smooth_a0.001_v10'
    if not os.path.exists(data_dir):
        return

    from pa3py import constants as c

    def chemistry_4_zones(r_cm, t_sec):
        if r_cm < 2.7 * c.AU:
            return {'silicates': 1.0}
        elif r_cm < 5.0 * c.AU:
            return {'silicates': 0.5, 'H2O': 0.5}
        else:
            return {'silicates': 0.3, 'H2O': 0.3, 'CO2': 0.4}

    sim = PA3Py(data_dir)
    sim.set_custom_chemistry(chemistry_4_zones)
    results = sim.run_growth([2.0, 8.0])

    # Lista de especies: color por embrión, estilo de línea por especie
    fig, ax = sim.plot_species_fraction(results, species=['H2O', 'CO2'])
    assert fig is not None
    fig.savefig(tmp_path / "multi_species.png")
    assert (tmp_path / "multi_species.png").exists()
    plt.close(fig)

    # None = todas las especies rastreadas, un solo embrión (color por especie)
    fig2, ax2 = sim.plot_species_fraction(results, species=None, embryos=[8.0])
    assert fig2 is not None
    fig2.savefig(tmp_path / "all_species_single_embryo.png")
    assert (tmp_path / "all_species_single_embryo.png").exists()
    plt.close(fig2)
