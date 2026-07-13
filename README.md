# PA3Py (Pebble Accretion 3)

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**PA3Py** is a modular Python package for post-processing protoplanetary disk hydrodynamic simulations (like those from `dustpy` or `tripodpy`) to compute the 2D and 3D pebble accretion growth of planetary embryos.

## Installation

```bash
git clone https://github.com/tu-usuario/pa3py.git
cd pa3py
pip install -e .
```

## Physics & References

Este módulo integra estricta y directamente las formulaciones físicas obtenidas en:
- **Ormel (2017)**: *The Emerging Paradigm of Pebble Accretion*.
- **Drążkowska et al. (2023)**: *Planet Formation Theory in the Era of ALMA and Kepler: from Pebbles to Exoplanets*.

### 1. Masa Crítica para iniciar la Acreción de Pebbles (*Onset Mass*)
Aparece una exclusión física para los cuerpos muy livianos. Para poder capturar el gas necesario de arrastre, un cuerpo debe tener una masa mínima. Implementamos:
$$ M_{\rm PA\ onset} = {\rm St} \eta^3 M_{\star} $$

*Implementación:* Si $M_{\rm pl} < M_{\rm PA\ onset}$, el modelo no asume tasa nula, sino que cruza al modo de **Acreción de Planetesimales (Safronov)** (Ec. 7.14 de Ormel 2017), asumiendo régimen gas-free balístico bajo densidad propia rocosa. Esto produce que haya una transición real entre el crecimiento primigenio de embriones pequeños al runaway estrepitoso provisto por pebbles.

### 2. Ecuaciones Dinámicas de Acreción
Se dividen orgánicamente los regímenes del acercamiento orbital según las propiedades fluidas analizadas por **Ormel (2017)**. 
- **Headwind (Bondi-like drift limit)** ($M \lesssim M_{\rm hw/sh}$):
  $$ \dot{M}_{\rm 2D, hw} = \sqrt{8 G M_{\rm pl} t_{\rm stop} v_{\rm hw}} \Sigma_{\rm peb} $$
- **Shear (Hill limit)** ($M \gtrsim M_{\rm hw/sh}$):
  $$ \dot{M}_{\rm 2D, sh} = 2 R_{\rm Hill}^2 \Omega_{\rm K} {\rm St}^{2/3} \Sigma_{\rm peb} $$
  
Con la masa de transición evaludaba matemáticamente con:
$$ M_{\rm hw/sh} = \frac{v_{\rm hw}^3}{8 G \Omega_{\rm K} t_{\rm stop}} $$

### 3. Transición Analítica a disco Extenso 3D
Ya no se evalúa `erf` numéricamente para forzar transiciones. Se utiliza la aproximación turbulenta del factor de decaimiento natural propuesta por **Dubrulle (1995) / Ormel 2017 (Eq 7.24)**:
$$ \dot{M} = \dot{M}_{\rm 2D} \left( \frac{b_{\rm col}}{b_{\rm col} + h_{\rm peb} \sqrt{8/\pi}} \right) $$
Donde $h_{\rm peb} = \sqrt{\frac{\alpha_T}{\alpha_T + {\rm St}}} h_{\rm gas}$. Esta transición converge al límite 3D natural en discos muy agitados verticalmente.

### 4. Masa Reguladora de Aislamiento
Las barreras numéricas y el *overshoot* quedan protegidos por el rediseño más actualizado sobre retroalimentación entre aislamiento del gap y gas proporcionado por Bitsch, modificado bajo la última simplificación analítica de la comparación de **Drążkowska et al. 2023**:
$$ M_{\rm iso, peb} = 25 M_{\oplus} \left(\frac{H/r}{0.05}\right)^3 \left(\frac{M_{\star}}{M_{\odot}}\right) $$

### 5. Consumo y Drift
Los Pebbles viajan hacia la estrella limitados por:
$$ v_{\rm r,solid} \approx \frac{- 2\eta v_{\rm K} {\rm St}}{1 + {\rm St}^2} $$
Para mantener el módulo de acreción puramente parasitario y conservador, se extrae la traza dinámica de los archivos HDF5 asumiendo compatibilidad directa.

### 6. Composición Planetaria (Tracking Físico)
Los embriones se inicializan asumiendo una semilla `100%` rocosa (silicatos). El tracking de composición evalúa la posición del embrión respecto a la *snowline* dinámica (`r_snow`) de agua extraída de la simulación hidrodinámica base.
- Si $r \ge r_{\rm snow}$: El material acretado asume composición de hielo (50% H2O, 50% Silicatos).
- Si $r < r_{\rm snow}$: El material acretado asume composición de polvo seco (100% Silicatos).
