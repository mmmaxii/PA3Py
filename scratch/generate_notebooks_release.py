import os
import nbformat as nbf

DOCS_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'docs', 'notebooks'))
os.makedirs(DOCS_DIR, exist_ok=True)

# Delete old notebooks to clean up
import glob
for f in glob.glob(os.path.join(DOCS_DIR, "*.ipynb")):
    os.remove(f)

def create_nb_00():
    nb = nbf.v4.new_notebook()
    nb.cells = [
        nbf.v4.new_markdown_cell("# 0. Física Implementada en PA3Py\n\nEste notebook contiene la documentación científica completa del modelo físico subyacente de `PA3Py`."),
        nbf.v4.new_markdown_cell("## 1. Acreción de Pebbles (Pebble Accretion)\n\nEl núcleo del simulador resuelve el crecimiento del núcleo planetario gobernado por la acreción de *pebbles* (Lambrechts & Johansen 2012, 2014). Dependiendo de si la capa de *pebbles* es más gruesa o más delgada que el radio de Hill del embrión, el régimen puede ser 2D o 3D."),
        nbf.v4.new_markdown_cell("### Régimen 3D\n\nCuando el embrión es pequeño, su radio de acreción $R_{acc}$ es menor que la escala de altura de los *pebbles* ($H_d$). La tasa de acreción volumétrica es:\n\n$$ \dot{M}_{3D} = \\pi R_{acc}^2 \\Delta v \\rho_{d} $$\n\nDonde $\\rho_d = \\frac{\\Sigma_d}{\\sqrt{2\\pi} H_d}$ es la densidad volumétrica de los *pebbles* en el plano medio."),
        nbf.v4.new_markdown_cell("### Régimen 2D\n\nA medida que el embrión crece y su $R_{acc}$ excede $H_d$, la acreción se vuelve puramente bidimensional (disco aplanado desde la perspectiva del planeta):\n\n$$ \dot{M}_{2D} = 2 R_{acc} \\Sigma_d \\Delta v $$\n\nEl simulador interpola suavemente entre ambos regímenes."),
        nbf.v4.new_markdown_cell("## 2. Masa de Transición ($M_t$)\n\nLa eficiencia de acreción depende del parámetro aerodinámico de Stokes ($St$). Cuando el número de Stokes excede $0.1$, el régimen entra en un límite de acreción óptima.\n\nLa masa crítica donde la acreción es máxima para un $St$ dado es la Masa de Transición (Bitsch et al. 2015):\n\n$$ M_t = \\sqrt{\\frac{1}{3}} \\frac{\\Delta v^3}{G \\Omega} $$"),
        nbf.v4.new_markdown_cell("## 3. Masa de Aislamiento de Pebbles ($M_{iso}$)\n\nA medida que el planeta gana masa, perturba el gas a su alrededor creando un gradiente de presión positivo (un *pressure bump*) fuera de su órbita. Este gradiente atrapa los *pebbles* que derivan hacia adentro, cortando el suministro.\n\n`PA3Py` utiliza la prescripción de masa de aislamiento definida por Bitsch et al. (2018):\n\n$$ M_{iso} = 25 M_\\oplus \\left( \\frac{H/r}{0.05} \\right)^3 \\left( \\frac{a}{1 \\times 10^{-3}} \\right) \\times \\Pi $$\n\nDonde $\\Pi$ es una función compleja de los parámetros del disco y el gradiente de presión local."),
        nbf.v4.new_markdown_cell("## 4. Evolución de la Snowline y Química\n\nLas *snowlines* dictan qué especies se condensan sobre los *pebbles*. En el modelo por defecto, la posición de la línea de hielo de agua ($r_{snow}$) evoluciona en el tiempo debido a los cambios de temperatura del disco (radiación pasiva + calentamiento viscoso), de acuerdo a los tratamientos teóricos tipo Oka et al. (2011) o extrayendo la temperatura directamente de simulaciones numéricas.")
    ]
    with open(os.path.join(DOCS_DIR, '00_Fisica_Implementada.ipynb'), 'w', encoding='utf-8') as f:
        nbf.write(nb, f)

def create_nb_01():
    nb = nbf.v4.new_notebook()
    nb.cells = [
        nbf.v4.new_markdown_cell("# 1. Primeros Pasos con PA3Py\n\nEn este notebook aprenderemos a inicializar `PA3Py`, correr una simulación de crecimiento para un embrión y graficar los resultados."),
        nbf.v4.new_code_cell("""import sys
import os
import matplotlib.pyplot as plt

# Agregar la ruta del paquete (si corres desde el repositorio local)
sys.path.insert(0, os.path.abspath('../../src'))
from pa3py import PA3Py"""),
        nbf.v4.new_markdown_cell("## Inicializando la Simulación\n\nCargamos la simulación simplemente pasando la ruta de nuestro disco HDF5 a la clase principal `PA3Py`."),
        nbf.v4.new_code_cell("""# Ruta a tus datos HDF5 de TripodPy
data_path = '../../tests/test_data/run_smooth_a0.001_v10'

# Inicializamos el motor
sim = PA3Py(data_path)"""),
        nbf.v4.new_markdown_cell("## Corriendo el Crecimiento\n\nLe diremos a `PA3Py` que corra el crecimiento para un embrión ubicado a 5 AU, con una masa inicial semilla de $10^{-3} M_\\oplus$."),
        nbf.v4.new_code_cell("""# Embriones a correr
embriones = [5.0]

# Ejecutamos el módulo
resultados = sim.run_growth(embriones, m_seed_me=1e-3)"""),
        nbf.v4.new_markdown_cell("## Graficando la Evolución de Masa\n\nEl diccionario `resultados` contiene matrices donde la primera columna es el tiempo (en segundos) y la segunda es la masa total del embrión."),
        nbf.v4.new_code_cell("""# Extraemos los datos del embrión en 5 AU
historia_5au = resultados[5.0]

tiempo_años = historia_5au[:, 0] / (3.154e7)  # Convertir s a años
masa_total = historia_5au[:, 1] / sim.engine.M_EARTH # Convertir a Masas Terrestres
masa_iso = historia_5au[:, 2] / sim.engine.M_EARTH

plt.figure(figsize=(8,5))
plt.plot(tiempo_años, masa_total, lw=2, label="Masa del Núcleo")
plt.plot(tiempo_años, masa_iso, '--', color='gray', label="Masa de Aislamiento (Límite)")
plt.xlabel("Tiempo [Años]")
plt.ylabel("Masa [M_\\oplus]")
plt.title("Crecimiento de un Embrión a 5 AU")
plt.legend()
plt.grid(True)
plt.show()"""),
        nbf.v4.new_markdown_cell("## Fracción de Agua\n\nPor defecto, PA3Py corre el modelo clásico `silicates` y `H2O` midiendo la posición de la snowline.\nPodemos extraer la masa de agua (columna 4) y ver cómo subió la fracción."),
        nbf.v4.new_code_cell("""masa_silicatos = historia_5au[:, 3]
masa_agua = historia_5au[:, 4]
masa_suma = masa_silicatos + masa_agua

fraccion_agua = 100 * masa_agua / masa_suma

plt.figure(figsize=(8,5))
plt.plot(tiempo_años, fraccion_agua, color='blue', lw=2)
plt.xlabel("Tiempo [Años]")
plt.ylabel("Fracción de Agua [%]")
plt.title("Evolución Química del Planeta a 5 AU")
plt.ylim(0, 60)
plt.grid(True)
plt.show()""")
    ]
    with open(os.path.join(DOCS_DIR, '01_Primeros_Pasos.ipynb'), 'w', encoding='utf-8') as f:
        nbf.write(nb, f)

def create_nb_02():
    nb = nbf.v4.new_notebook()
    nb.cells = [
        nbf.v4.new_markdown_cell("# 2. Química Avanzada (Zonas Dinámicas)\n\n`PA3Py` fue diseñado para ser completamente agnóstico en términos de composición. Tú puedes inyectar tus propias funciones químicas de Python y él automáticamente descubrirá las especies y construirá las tablas."),
        nbf.v4.new_code_cell("""import sys
import os
sys.path.insert(0, os.path.abspath('../../src'))
from pa3py import PA3Py

# Cargar simulación
sim = PA3Py('../../tests/test_data/run_smooth_a0.001_v10')"""),
        nbf.v4.new_markdown_cell("## Función de Composición Personalizada\n\nDefiniremos un disco protoplanetario con 4 especies químicas (`silicatos`, `H2O`, `CO2`, `CO`) particionado en 4 zonas dinámicas."),
        nbf.v4.new_code_cell("""from pa3py import constants as c

def quimica_4_zonas(r_cm, t_sec):
    # 1. Definimos las snowlines (algunas pueden migrar)
    r_h2o = 2.73 * c.AU * (max(t_sec, 1e-6) / 1e13)**(-0.5)
    r_co2 = 5.0 * c.AU
    r_co  = 12.0 * c.AU
    
    # 2. Asignamos abundancias relativas según la distancia al sol
    if r_cm < r_h2o:
        return {'silicatos': 1.0}                           # 100% rocoso
    elif r_cm < r_co2:
        return {'silicatos': 0.5, 'H2O': 0.5}               # Zona de Hielo
    elif r_cm < r_co:
        return {'silicatos': 0.3, 'H2O': 0.3, 'CO2': 0.4}   # Zona CO2
    else:
        return {'silicatos': 0.2, 'H2O': 0.2, 'CO2': 0.3, 'CO': 0.3} # Zona fría

# Inyectamos nuestra función al motor. PA3Py auto-detectará las especies.
sim.set_custom_chemistry(quimica_4_zonas)"""),
        nbf.v4.new_markdown_cell("## Corriendo la Química Personalizada\n\nColocaremos embriones en varias regiones y dejaremos que PA3Py se encargue del resto."),
        nbf.v4.new_code_cell("""# Embriones en la zona rocosa, hielo, CO2 y CO
resultados = sim.run_growth([2.0, 4.0, 8.0, 15.0])""")
    ]
    with open(os.path.join(DOCS_DIR, '02_Quimica_Avanzada.ipynb'), 'w', encoding='utf-8') as f:
        nbf.write(nb, f)

def create_nb_03():
    nb = nbf.v4.new_notebook()
    nb.cells = [
        nbf.v4.new_markdown_cell("# 3. Poblaciones Sintéticas y HDF5 I/O\n\nSi vas a simular la acreción de miles de planetas, usar un arreglo `numpy` es computacionalmente vital. Además, `PA3Py` provee métodos nativos para guardar toda la población en disco (`.h5`)."),
        nbf.v4.new_code_cell("""import sys
import os
import numpy as np
sys.path.insert(0, os.path.abspath('../../src'))
from pa3py import PA3Py

# Inicializamos el motor
sim = PA3Py('../../tests/test_data/run_smooth_a0.001_v10')"""),
        nbf.v4.new_markdown_cell("## Graficando el Diagrama de Hovmöller\n\nAntes de colocar embriones, es útil visualizar las snowlines sobre el disco."),
        nbf.v4.new_code_cell("""import matplotlib.pyplot as plt
sim.plot_hovmoller(field='dust_Sigma', show_snowlines=True)
plt.show()"""),
        nbf.v4.new_markdown_cell("## Definiendo la Población Sintética\n\nUsaremos `np.linspace` para generar 100 embriones a lo largo del disco."),
        nbf.v4.new_code_cell("""# Importante: para evitar decimales infinitos (ej. 1.98989), usamos 99 puntos 
# para un rango de 98 unidades (1 a 99) y así asegurar pasos de 1 entero.
embryos = np.linspace(1.0, 99.0, 99).tolist()

# Si no necesitas imprimir las keys, usar np.linspace(1, 99, 100) es perfectamente válido.

print("Simulando", len(embryos), "embriones...")
resultados = sim.run_growth(embryos)"""),
        nbf.v4.new_markdown_cell("## Guardando Resultados (HDF5)\n\nGuardaremos toda la corrida a disco duro."),
        nbf.v4.new_code_cell("""sim.save_results(resultados, "poblacion_sintetica_100.h5")
print("¡Datos guardados!")"""),
        nbf.v4.new_markdown_cell("## Cargando Resultados\n\nMañana, cuando abras este notebook de nuevo, no tendrás que volver a correr la física."),
        nbf.v4.new_code_cell("""# Carga la matriz y también la lista de química usada en esa corrida
resultados_cargados, quimica = PA3Py.load_results("poblacion_sintetica_100.h5")

print("Especies de la simulación original:", quimica)
print("Planetas cargados:", len(resultados_cargados))""")
    ]
    with open(os.path.join(DOCS_DIR, '03_Poblaciones_Sinteticas_e_IO.ipynb'), 'w', encoding='utf-8') as f:
        nbf.write(nb, f)

if __name__ == '__main__':
    create_nb_00()
    create_nb_01()
    create_nb_02()
    create_nb_03()
    print("Notebooks generados con exito en docs/notebooks/")
