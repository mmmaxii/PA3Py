import json
import sys

def create_nb(filename, cells):
    nb = {
        "cells": cells,
        "metadata": {
            "kernelspec": {
                "display_name": "Python 3",
                "language": "python",
                "name": "python3"
            }
        },
        "nbformat": 4,
        "nbformat_minor": 4
    }
    with open(filename, "w", encoding="utf-8") as f:
        json.dump(nb, f, indent=2, ensure_ascii=False)

def markdown_cell(source):
    return {"cell_type": "markdown", "metadata": {}, "source": [s + "\n" for s in source.split("\n")]}

def code_cell(source):
    return {"cell_type": "code", "execution_count": None, "metadata": {}, "outputs": [], "source": [s + "\n" for s in source.split("\n")]}

cells_04 = [
    markdown_cell("# Poblaciones Sintéticas y Diagramas de Hovmöller\nEste notebook enseña cómo simular muchos embriones simultáneamente a lo largo del disco y cómo usar las herramientas avanzadas de visualización de PA3Py."),
    code_cell("import sys\nimport os\nimport numpy as np\nimport matplotlib.pyplot as plt\n\nsys.path.insert(0, os.path.abspath('../../src'))\nfrom pa3py.data import load_tripodpy_hdf5\nfrom pa3py.pebble_accretion import PebbleAccretionModule3\nfrom pa3py.plotting import plot_hovmoller\nfrom pa3py import constants as c"),
    markdown_cell("## 1. El Disco (Diagrama de Hovmöller)\nEl diagrama de Hovmöller nos permite ver cómo evoluciona la superficie del polvo en el espacio (Radio) y el tiempo."),
    code_cell("data_path = '../../tests/test_data/run_smooth_a0.001_v10'\ndisk = load_tripodpy_hdf5(data_path)\n\n# Graficamos la densidad de superficie del polvo\nfig, ax = plot_hovmoller(disk, field='dust_Sigma', show_snowlines=True)\nplt.show()"),
    markdown_cell("## 2. Simulando una Población de 100 Embriones\nVamos a colocar 100 semillas de $10^{-3} M_\oplus$ repartidas entre 1 y 30 AU para ver la huella final de la acreción de pebbles."),
    code_cell("# Instanciamos el simulador con la física estándar\nsim = PebbleAccretionModule3(disk)\n\n# 100 embriones desde 1 AU hasta 30 AU\nembryos_au = np.linspace(1.0, 30.0, 100).tolist()\n\n# Corremos la simulación de toda la población (¡PA3Py es muy rápido!)\nresultados = sim.run_growth(embryos_au, M0_g=1e-3 * sim.M_EARTH)"),
    markdown_cell("## 3. Análisis Poblacional: Masa Final vs Distancia y Isolation Mass\nExtraemos la masa de aislamiento teórica en todo el disco y graficamos cómo se compara con las masas que realmente alcanzaron nuestros planetas."),
    code_cell("# Masa de aislamiento para todo el disco en el tiempo t=0\nm_iso_map = sim.calculate_isolation_mass_map()\nm_iso_t0 = m_iso_map[0, :] / c.M_EARTH\nr_disco = disk.r / c.AU\n\n# Extraer masas finales de la población\nradios = []\nmasas_finales = []\n\nfor r_au, hist in resultados.items():\n    if len(hist) > 0:\n        radios.append(r_au)\n        masas_finales.append(hist[-1][1] / c.M_EARTH)\n\nplt.figure(figsize=(10, 5))\n\n# Línea teórica de masa de aislamiento inicial\nplt.plot(r_disco, m_iso_t0, 'k--', alpha=0.5, label='Isolation Mass teórica (t=0)')\n\n# Resultados de la simulación\nplt.scatter(radios, masas_finales, c=radios, cmap='viridis', s=50, edgecolor='k', label='Masa Final Alcanzada')\n\nplt.xscale('log')\nplt.yscale('log')\nplt.xlabel('Radio Orbital [AU]')\nplt.ylabel('Masa del Planeta [$M_\\\\oplus$]')\nplt.title('Población Sintética: Masa Final vs Distancia')\nplt.xlim(0.5, 50)\nplt.legend()\nplt.grid(True, alpha=0.3)\nplt.show()")
]

create_nb('docs/notebooks/04_Poblaciones_Sinteticas.ipynb', cells_04)
print("Notebook 04 generado.")
