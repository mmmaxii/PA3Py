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

# 01 Quickstart
cells_01 = [
    markdown_cell("# PA3Py Quickstart: Crecimiento de Waterworlds\nEste notebook muestra la forma más sencilla de correr PA3Py usando la interfaz unificada `PA3Py`."),
    code_cell("import sys\nimport os\nimport numpy as np\nimport matplotlib.pyplot as plt\n\n# Añadir PA3Py al path local\nsys.path.insert(0, os.path.abspath('../../src'))\nfrom pa3py import PA3Py"),
    markdown_cell("## 1. Inicializar la Simulación\nSolo necesitamos la ruta a un archivo HDF5 válido de TripodPy/DustPy."),
    code_cell("# NOTA: Asegúrate de apuntar a la ruta correcta donde tienes tus datos\ndata_path = '../../tests/test_data/run_smooth_a0.001_v10'\n\n# Un solo objeto centraliza todo el paquete\nsim = PA3Py(data_path)"),
    markdown_cell("## 2. Correr la Acreción\nCorremos la simulación para embriones en 1, 5 y 15 AU."),
    code_cell("resultados = sim.run_growth(embryos_au=[1.0, 5.0, 15.0], m_seed_me=1e-3)"),
    markdown_cell("## 3. Visualización\nGraficaremos la masa total vs el tiempo para cada planeta."),
    code_cell("plt.figure(figsize=(10, 6))\n\nfor r_au, hist in resultados.items():\n    if len(hist) == 0: continue\n    \n    tiempos_myr = hist[:, 0] / (1e6 * 3.15576e7)  # Convertir seg a Myr\n    masa_total = hist[:, 1] / 5.97e27             # Convertir g a Masas Terrestres\n    \n    plt.plot(tiempos_myr, masa_total, label=f'Embrión en {r_au} AU')\n\nplt.xlabel('Tiempo [Myr]')\nplt.ylabel('Masa del Núcleo [$M_\\\\oplus$]')\nplt.yscale('log')\nplt.legend()\nplt.grid(True, alpha=0.3)\nplt.title('Crecimiento de Embriones Planetarios')\nplt.show()")
]

# 02 Física y Snowline
cells_02 = [
    markdown_cell("# La Física del Disco y la Snowline Dinámica\nEn este notebook profundizaremos en cómo PA3Py modela la migración térmica de la Snowline utilizando los resultados de transferencia radiativa de **Oka et al. (2011)** combinados con la evolución de la tasa de acreción de gas de **Hartmann et al. (1998)**."),
    code_cell("import sys\nimport os\nimport numpy as np\nimport matplotlib.pyplot as plt\nimport matplotlib.ticker as ticker\nfrom scipy.interpolate import interp1d\n\n# Añadir PA3Py al path local\nsys.path.insert(0, os.path.abspath('../../src'))\nfrom pa3py import PA3Py\nfrom pa3py.snowline import generate_rsnow_array, mdot_time, _csv_path\nfrom pa3py import constants as c\n\n# Estilo de Plotting (Paper)\nplt.rcParams.update({\n    'text.usetex': False,\n    'font.family': 'serif',\n    'font.size': 14,\n    'axes.labelsize': 16,\n    'lines.linewidth': 2.0,\n    'figure.dpi': 150,\n})"),
    markdown_cell("## 1. Evolución de la Tasa de Acreción (Hartmann et al.)\nAsumimos que la tasa de acreción del gas decae en el tiempo como $\\dot{M}(t) = \\dot{M}_0 (t / 1 \\text{Myr})^{-\\eta}$."),
    code_cell("t_eval = np.logspace(-1, 1, 500) # 0.1 a 10 Myr\neta_plot = 1.5\nm_t = mdot_time(t_eval, eta=eta_plot)\nlog_t_eval = np.log10(t_eval)\nlog_m_t = np.log10(m_t)"),
    markdown_cell("## 2. Interpolación de la Snowline (Oka et al. 2011)\nCargamos los datos corregidos de Oka et al. 2011 y generamos una interpolación log-log para mapear cualquier $\\dot{M}$ a su respectiva $r_{snow}$."),
    code_cell("mdot_list, rsnow_list = [], []\nwith open(_csv_path, 'r') as f:\n    lines = f.readlines()\n    for line in lines[1:]:\n        line = line.strip()\n        if not line: continue\n        parts = line.split(';')\n        if len(parts) == 2:\n            mdot_list.append(float(parts[0]))\n            rsnow_list.append(float(parts[1]))\n\nmdot_raw, rsnow_raw = np.array(mdot_list), np.array(rsnow_list)\nsort_idx = np.argsort(mdot_raw)\nmdot_unique, unique_indices = np.unique(mdot_raw[sort_idx], return_index=True)\nrsnow_unique = rsnow_raw[sort_idx][unique_indices]\n\ninterp_log_rsnow = interp1d(np.log10(mdot_unique), np.log10(rsnow_unique), kind='linear', fill_value='extrapolate')\n\ndef r_snow_from_mdot(mdot_val):\n    return 10**interp_log_rsnow(np.log10(mdot_val))\n\ndef r_snow_time(t_myr, eta=1.5):\n    return r_snow_from_mdot(mdot_time(t_myr, eta))"),
    markdown_cell("## 3. Diagnóstico Completo (El Gráfico Oficial)\nConstruimos un panel triple mostrando toda la física de la snowline."),
    code_cell("fig, axes = plt.subplots(1, 3, figsize=(16, 5), constrained_layout=True)\nax1, ax2, ax3 = axes\ncolor_main, color_raw = '#d95f02', '#1f78b4'\n\n# --- PANEL 1: Evolución de Mdot(t) ---\nax1.plot(log_t_eval, log_m_t, color=color_main, linewidth=2.5, label=fr'$\\eta = {eta_plot}$')\nax1.set_xlim(-1, 1)\nax1.set_xticks([-1, -0.5, 0, 0.5, 1])\nax1.set_xlabel(r'$\\log_{10}(t\\ [\\mathrm{Myr}])$')\nax1.set_ylabel(r'$\\log_{10}(\\dot{M}\\ [M_\\odot\\ \\mathrm{yr}^{-1}])$')\nax1.legend(loc='upper right', framealpha=0.85)\n\n# --- PANEL 2: Interpolación Oka et al ---\nm_test = np.logspace(-12, -7, 500)\nr_test = r_snow_from_mdot(m_test)\nax2.plot(np.log10(mdot_raw), rsnow_raw, color=color_raw, linestyle='-', linewidth=3, alpha=0.5, label='Oka et al. 2011')\nax2.plot(np.log10(m_test), r_test, color=color_main, linestyle='--', linewidth=2.5, label='Interpolación')\nax2.set_xlim(-7, -12)\nax2.set_yscale('log')\nax2.set_ylim(0.1, 10)\nax2.set_yticks([0.1, 1, 10])\nax2.get_yaxis().set_major_formatter(ticker.ScalarFormatter())\nax2.set_xlabel(r'$\\log_{10}(\\dot{M}\\ [M_\\odot\\ \\mathrm{yr}^{-1}])$')\nax2.set_ylabel(r'$r_{\\mathrm{snow}}\\ [\\mathrm{AU}]$')\n\n# Linea en 10 Myr\nlog_m_10myr = np.log10(mdot_time(10, eta=eta_plot))\nax2.axvline(log_m_10myr, color='k', linestyle=':', linewidth=1.5, alpha=0.6, label=r'$\\dot{M}(10\\,\\mathrm{Myr})$')\nax2.legend(loc='lower right', framealpha=0.85)\n\n# --- PANEL 3: Evolución de r_snow(t) ---\nax3.plot(log_t_eval, r_snow_time(t_eval, eta=eta_plot), color=color_main, linewidth=2.5, label=fr'$\\eta = {eta_plot}$')\nax3.set_xlim(-1, 1)\nax3.set_xticks([-1, -0.5, 0, 0.5, 1])\nax3.set_yscale('log')\nax3.set_ylim(0.1, 10)\nax3.set_yticks([0.1, 1, 10])\nax3.get_yaxis().set_major_formatter(ticker.ScalarFormatter())\nax3.set_xlabel(r'$\\log_{10}(t\\ [\\mathrm{Myr}])$')\nax3.set_ylabel(r'$r_{\\mathrm{snow}}\\ [\\mathrm{AU}]$')\n\nfrom scipy.optimize import brentq\ntry:\n    t_1au = brentq(lambda t: r_snow_time(t, eta=eta_plot) - 1.0, 0.1, 10.0)\n    ax3.plot(np.log10(t_1au), 1.0, 'ko', markersize=6, zorder=5, label=f'{t_1au:.2f} Myr')\n    ax3.axvline(np.log10(t_1au), color='k', linestyle=':', alpha=0.6)\n    ax3.axhline(1.0, color='k', linestyle=':', alpha=0.6)\nexcept ValueError:\n    pass\nax3.legend(loc='lower left', framealpha=0.85)\n\nfor ax in [ax1, ax2, ax3]:\n    ax.grid(True, which='both', linestyle=':', alpha=0.3, color='gray')\nplt.show()")
]

# 03 Advanced Chemistry
cells_03 = [
    markdown_cell("# Química Compleja y Zonas Múltiples\nEste tutorial muestra cómo inyectar funciones de Python para modelar química altamente compleja de manera simple usando la fachada `PA3Py`."),
    code_cell("import sys\nimport os\nimport numpy as np\nimport matplotlib.pyplot as plt\n\nsys.path.insert(0, os.path.abspath('../../src'))\nfrom pa3py import PA3Py\nfrom pa3py import constants as c"),
    markdown_cell("## 1. Definiendo tu propia física química\nPA3Py te permite definir una función simple `f(r, t)` que retorne diccionarios de abundancias."),
    code_cell("def quimica_multizona(r_cm, t_sec):\n    r_h2o = 2.73 * c.AU * (max(t_sec, 1e-6) / 1e13)**(-0.5)\n    r_co2 = 5.0 * c.AU\n    r_co  = 12.0 * c.AU\n    \n    if r_cm < r_h2o:\n        return {'silicatos': 1.0}  # Zona rocosa pura\n    elif r_cm < r_co2:\n        return {'silicatos': 0.5, 'H2O': 0.5}  # Zona agua\n    elif r_cm < r_co:\n        return {'silicatos': 0.3, 'H2O': 0.3, 'CO2': 0.4} # Zona CO2\n    else:\n        return {'silicatos': 0.2, 'H2O': 0.2, 'CO2': 0.3, 'CO': 0.3} # Zona CO"),
    markdown_cell("## 2. Inyectar el modelo en la simulación\nUsamos `sim.set_custom_chemistry` para cambiar la dinámica al vuelo."),
    code_cell("data_path = '../../tests/test_data/run_smooth_a0.001_v10'\nsim = PA3Py(data_path)\n\n# Cambiamos la química\nsim.set_custom_chemistry(quimica_multizona, ['silicatos', 'H2O', 'CO2', 'CO'])\n\n# Corremos en varias regiones para ver los resultados distintos\nresultados = sim.run_growth([2.0, 4.0, 8.0, 15.0], m_seed_me=1e-3)")
]

# 04 Synthetic Populations
cells_04 = [
    markdown_cell("# Poblaciones Sintéticas y Diagramas de Hovmöller\nEste notebook enseña cómo simular muchos embriones simultáneamente y usar los gráficos integrados de `PA3Py`."),
    code_cell("import sys\nimport os\nimport numpy as np\nimport matplotlib.pyplot as plt\n\nsys.path.insert(0, os.path.abspath('../../src'))\nfrom pa3py import PA3Py\nfrom pa3py import constants as c"),
    markdown_cell("## 1. El Disco (Diagrama de Hovmöller)\n`PA3Py` trae gráficos integrados. Solo llama a `plot_hovmoller`."),
    code_cell("data_path = '../../tests/test_data/run_smooth_a0.001_v10'\nsim = PA3Py(data_path)\n\n# Graficamos la densidad de superficie del polvo\nfig, ax = sim.plot_hovmoller(field='dust_Sigma', show_snowlines=True)\nplt.show()"),
    markdown_cell("## 2. Simulando una Población de 100 Embriones"),
    code_cell("# 100 embriones desde 1 AU hasta 30 AU\nembryos_au = np.linspace(1.0, 30.0, 100).tolist()\n\n# Corremos la simulación de toda la población (¡PA3Py es muy rápido!)\nresultados = sim.run_growth(embryos_au, m_seed_me=1e-3)"),
    markdown_cell("## 3. Análisis Poblacional: Masa Final vs Masa de Aislamiento"),
    code_cell("# Extraemos la masa de aislamiento del objeto central\nm_iso_map = sim.calculate_isolation_mass_map()\nm_iso_t0 = m_iso_map[0, :] / c.M_EARTH\nr_disco = sim.disk.r / c.AU\n\n# Extraer masas finales de la población\nradios = []\nmasas_finales = []\n\nfor r_au, hist in resultados.items():\n    if len(hist) > 0:\n        radios.append(r_au)\n        masas_finales.append(hist[-1][1] / c.M_EARTH)\n\nplt.figure(figsize=(10, 5))\nplt.plot(r_disco, m_iso_t0, 'k--', alpha=0.5, label='Isolation Mass teórica (t=0)')\nplt.scatter(radios, masas_finales, c=radios, cmap='viridis', s=50, edgecolor='k', label='Masa Final Alcanzada')\n\nplt.xscale('log')\nplt.yscale('log')\nplt.xlabel('Radio Orbital [AU]')\nplt.ylabel('Masa del Planeta [$M_\\\\oplus$]')\nplt.title('Población Sintética: Masa Final vs Distancia')\nplt.xlim(0.5, 50)\nplt.legend()\nplt.grid(True, alpha=0.3)\nplt.show()")
]

create_nb('../docs/notebooks/01_Quickstart_Waterworlds.ipynb', cells_01)
create_nb('../docs/notebooks/02_Fisica_del_Disco_y_Snowline.ipynb', cells_02)
create_nb('../docs/notebooks/03_Quimica_Compleja_y_Zonas.ipynb', cells_03)
create_nb('../docs/notebooks/04_Poblaciones_Sinteticas.ipynb', cells_04)

print("All notebooks updated with Facade API.")
