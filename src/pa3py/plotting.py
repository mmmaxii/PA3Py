"""
plotting.py - Herramientas de visualización (ej: Diagrama Hovmöller)
"""

from typing import Optional, Union, List

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import LogNorm, Normalize
from matplotlib.lines import Line2D
from .data import DiskData
from . import constants as c

_DISCRETE_LEGEND_MAX = 8
_SPECIES_LINESTYLES = ['-', '--', '-.', ':']


def _resolve_time_unit(t_unit: str) -> tuple:
    """Devuelve (factor_de_conversion, etiqueta_eje) para un t_unit dado."""
    t_map = {'Myr': (c.YEAR * 1e6, "Time [Myr]"), 'kyr': (c.YEAR * 1e3, "Time [kyr]")}
    return t_map.get(t_unit, (c.YEAR, "Time [Years]"))


def _radius_colors(radii_au, cmap: str = 'viridis'):
    """
    Decide cómo colorear N trayectorias identificadas por su radio inicial.

    Retorna (colors, mode, extra) donde:
      - mode == 'discrete': colors es una lista de colores (uno por radio, colormap
        cualitativo 'tab10'); el caller debe añadir una leyenda por línea.
      - mode == 'continuous': colors es una lista de colores tomados de `cmap`
        normalizados sobre [min(radii), max(radii)]; extra es el ScalarMappable
        que el caller debe usar para dibujar un colorbar.
    """
    radii_au = np.asarray(radii_au, dtype=float)
    if len(radii_au) <= _DISCRETE_LEGEND_MAX:
        tab10 = plt.get_cmap('tab10')
        colors = [tab10(i % 10) for i in range(len(radii_au))]
        return colors, 'discrete', None
    else:
        norm = Normalize(vmin=np.min(radii_au), vmax=np.max(radii_au))
        sm = cm.ScalarMappable(norm=norm, cmap=plt.get_cmap(cmap))
        colors = [sm.to_rgba(r) for r in radii_au]
        return colors, 'continuous', sm

def plot_hovmoller(disk: DiskData, field: str = 'dust_Sigma', 
                   cmap: str = 'magma', vmin: float = None, vmax: float = None,
                   t_unit: str = 'Myr', show_snowlines: bool = True):
    """Genera un diagrama de Hovmöller (Radio vs Tiempo) para una propiedad."""
    # 1. Extraer los datos 2D (Tiempo, Radio)
    if field == 'dust_Sigma':
        # Sumamos sobre todos los tamaños de polvo (último índice)
        Z = np.sum(disk.dust_Sigma, axis=-1)
        title = r"Hovmöller: $\Sigma_{dust}$ [g/cm$^2$]"
    elif field == 'gas_Sigma':
        Z = disk.gas_Sigma
        title = r"Hovmöller: $\Sigma_{gas}$ [g/cm$^2$]"
    elif field in ['eps', 'epsilon']:
        dust = np.sum(disk.dust_Sigma, axis=-1)
        gas = disk.gas_Sigma
        # Le sumamos 1e-300 al gas para evitar división por cero.
        Z = dust / (gas + 1e-300)
        Z = np.clip(Z, 1e-10, None)  # Limitar valores muy pequeños
        title = r"Hovmöller: Dust-to-gas ratio ($\epsilon = \Sigma_d / \Sigma_g$)"
    else:
        raise ValueError("Unsupported field. Options: 'dust_Sigma', 'gas_Sigma', 'epsilon'.")

    # Limpiar ceros para escala log
    Z = np.clip(Z, 1e-20, None)
    
    t_factor, t_label = _resolve_time_unit(t_unit)

    t_array = disk.times / t_factor
    r_array = disk.r / c.AU

    # Ejes extendidos robustos para pcolormesh
    def _log_edges(arr):
        arr = np.maximum(arr, 1e-10)
        log_arr = np.log10(arr)
        mid = (log_arr[:-1] + log_arr[1:]) / 2.0
        return 10**np.concatenate([[log_arr[0] - (mid[0] - log_arr[0])], mid, [log_arr[-1] + (log_arr[-1] - mid[-1])]])

    # Prevenir log10(0)
    safe_t = np.maximum(t_array, 1e-6)
    t_edges = _log_edges(safe_t)
    r_edges = _log_edges(r_array)

    # Escala de colores automática
    logZ = np.log10(Z[Z > 1e-20])
    v_lo = vmin if vmin is not None else 10**np.percentile(logZ, 2)
    v_hi = vmax if vmax is not None else 10**np.percentile(logZ, 98)

    # 2. Dibujar
    fig, ax = plt.subplots(figsize=(10, 5))
    
    # Transponemos Z porque pcolormesh recibe Y, X -> Z debe ser (Ny, Nx)
    pcm = ax.pcolormesh(t_edges, r_edges, Z.T, 
                        norm=LogNorm(vmin=v_lo, vmax=v_hi), 
                        cmap=cmap, shading='flat')
    
    cbar = fig.colorbar(pcm, ax=ax, pad=0.02)
    cbar.set_label("Logarithmic Scale")
    
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim(safe_t[0], safe_t[-1])
    ax.set_ylim(r_array[0], r_array[-1])
    ax.set_xlabel(t_label)
    ax.set_ylabel("Radius [AU]")
    ax.set_title(title, pad=10)
    
    if show_snowlines and disk.hdf5_snowlines:
        colors = ['cyan', 'white', 'lightgreen']
        for idx, (name, rsnow) in enumerate(disk.hdf5_snowlines.items()):
            r_au = rsnow / c.AU
            if np.any(r_au > 0):
                ax.plot(safe_t, r_au, ls='--', color=colors[idx % len(colors)], lw=2, label=f"Snowline {name}")
        ax.legend(loc='lower left')
                
    fig.tight_layout()
    return fig, ax



def plot_population(disk, results: dict, M_iso_map: np.ndarray = None, **kwargs):
    """
    Grafica la población sintética de planetas (Masa Final vs Posición Inicial).
    
    Parámetros:
    -----------
    disk : DiskData
        El objeto que contiene los datos del disco.
    results : dict
        Diccionario con los historiales de crecimiento de cada embrión.
    M_iso_map : np.ndarray, opcional
        Mapa 2D de masa de aislamiento calculada para todo el disco.
    """
    radii = []
    final_masses = []
    water_fractions = []

    for r_au, hist in results.items():
        if len(hist) == 0:
            continue
        
        radii.append(float(r_au))
        final_state = hist[-1]
        total_m = final_state[1] / c.M_EARTH
        
        silicate_m = final_state[3]
        water_m = final_state[4] if len(final_state) > 4 else 0.0
        M_total = silicate_m + water_m
        f_water = water_m / M_total if M_total > 0 else 0.0
        
        final_masses.append(total_m)
        water_fractions.append(f_water * 100)
        
    radii = np.array(radii)
    final_masses = np.array(final_masses)
    water_fractions = np.array(water_fractions)
    
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # 1. Snowline region
    if disk.hdf5_snowlines and 'H2O' in disk.hdf5_snowlines:
        rsnow_array = disk.hdf5_snowlines['H2O'] / c.AU
        rsnow_valid = rsnow_array[rsnow_array > 0]
        if len(rsnow_valid) > 0:
            rsnow_min = np.min(rsnow_valid)
            rsnow_max = np.max(rsnow_array)
            ax.axvspan(rsnow_min, rsnow_max, color='cyan', alpha=0.1, label='Snowline Migration')
            ax.axvline(rsnow_min, color='cyan', linestyle='--', alpha=0.5)
            ax.axvline(rsnow_max, color='cyan', linestyle='--', alpha=0.5)
            
    # 2. Isolation Mass from the full disk grid
    if M_iso_map is not None:
        # Usamos el límite en el último paso de tiempo (estado final del disco)
        iso_mass_final = M_iso_map[-1, :] / c.M_EARTH
        r_grid_au = disk.r / c.AU
        ax.plot(r_grid_au, iso_mass_final, color='gray', linestyle='--', alpha=0.7, zorder=1, label='Isolation Mass Limit')

    # 3. Planetas
    sc = ax.scatter(radii, final_masses, c=water_fractions, cmap='Blues', edgecolor='black', 
                    s=80, alpha=0.9, zorder=2, vmin=0, vmax=50)
                    
    cbar = fig.colorbar(sc, ax=ax)
    cbar.set_label('Water Mass Fraction [%]', rotation=270, labelpad=15)
    
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('Initial Radius [AU]')
    ax.set_ylabel('Final Mass [$M_\\oplus$]')
    ax.set_title('Synthetic Population: Final Mass vs Initial Position')
    ax.grid(True, which='both', linestyle=':', alpha=0.6)
    ax.legend(loc='upper right')

    fig.tight_layout()
    return fig, ax


def plot_growth_curves(results: dict, embryos: Optional[list] = None,
                        time_unit: str = 'Myr', show_isolation_mass: bool = True,
                        cmap: str = 'viridis', ax=None):
    """
    Grafica la masa acretada (M_core) vs tiempo para una o varias trayectorias
    de embriones, con la masa de aislamiento (M_iso) opcionalmente superpuesta.

    Parámetros:
    -----------
    results : dict
        Diccionario {r_au: historial} como el que retorna PA3Py.run_growth /
        PebbleAccretionModule3.run_growth. Columnas del historial:
        [tiempo_s, M_core_g, M_iso_g, ...especies].
    embryos : list, opcional
        Subconjunto de radios iniciales (claves de `results`) a graficar.
        Por defecto se grafican todas las trayectorias.
    time_unit : str
        'Myr', 'kyr', o cualquier otro valor para usar Años.
    show_isolation_mass : bool
        Si es True, superpone M_iso(t) como línea punteada del mismo color
        que la trayectoria de M_core correspondiente. La masa de aislamiento
        actúa como tope físico del crecimiento y define el límite superior
        del eje Y.
    cmap : str
        Colormap continuo a usar cuando hay más de 8 trayectorias.
    ax : matplotlib.axes.Axes, opcional
        Ejes existentes donde dibujar (para componer figuras). Si es None,
        se crea una figura nueva.

    Retorna:
    --------
    (fig, ax)
    """
    keys = sorted(results.keys()) if embryos is None else list(embryos)
    keys = [r_au for r_au in keys if len(results[r_au]) > 0]

    t_factor, t_label = _resolve_time_unit(time_unit)
    colors, mode, sm = _radius_colors(keys, cmap=cmap)

    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 5))
    else:
        fig = ax.figure

    m_floor = np.inf
    m_ceiling = 0.0
    for r_au, color in zip(keys, colors):
        hist = results[r_au]
        t = hist[:, 0] / t_factor
        m_core = hist[:, 1] / c.M_EARTH
        label = f"{r_au:g} AU" if mode == 'discrete' else None
        ax.plot(t, m_core, color=color, lw=2, label=label)
        m_floor = min(m_floor, np.min(m_core[m_core > 0], initial=np.inf))
        m_ceiling = max(m_ceiling, np.max(m_core))
        if show_isolation_mass:
            m_iso = hist[:, 2] / c.M_EARTH
            ax.plot(t, m_iso, color=color, lw=1, ls='--', alpha=0.6)
            m_ceiling = max(m_ceiling, np.max(m_iso))

    ax.set_xscale('log')
    ax.set_yscale('log')
    # La masa de aislamiento (o la masa máxima alcanzada) define el techo del eje
    if np.isfinite(m_floor) and m_ceiling > 0:
        ax.set_ylim(m_floor * 0.5, m_ceiling * 2.0)
    ax.set_xlabel(t_label)
    ax.set_ylabel(r"Mass [$M_\oplus$]")
    ax.set_title("Embryo Growth")
    ax.grid(True, which='both', linestyle=':', alpha=0.6)

    if mode == 'discrete':
        ax.legend(loc='best')
    elif sm is not None:
        fig.colorbar(sm, ax=ax, label='Initial Radius [AU]')

    fig.tight_layout()
    return fig, ax


def plot_species_fraction(results: dict, tracked_species: list,
                           species: Union[str, List[str], None] = 'H2O',
                           embryos: Optional[list] = None, time_unit: str = 'Myr',
                           cmap: str = 'viridis', ylim: Optional[tuple] = None, ax=None):
    """
    Grafica la fracción de masa (%) de una o varias especies rastreadas vs
    tiempo, para una o varias trayectorias de embriones.

    Parámetros:
    -----------
    results : dict
        Diccionario {r_au: historial}, ver `plot_growth_curves`.
    tracked_species : list
        Nombres de las especies rastreadas, en el mismo orden que las columnas
        3+ de cada historial (ej: PA3Py.engine.tracked_species).
    species : str, list o None
        Especie(s) a graficar. Por defecto 'H2O' (fracción de agua). Puede ser
        una lista (ej: ['H2O', 'CO2']) o None para graficar todas las especies
        rastreadas — útil con química multi-snowline (ver notebook 02).
        Con varias especies, el color identifica al embrión y el estilo de
        línea a la especie; con un solo embrión, el color identifica a la
        especie directamente.
    embryos, time_unit, cmap, ax : ver `plot_growth_curves`.
    ylim : tuple, opcional
        Límites del eje Y (%).

    Retorna:
    --------
    (fig, ax)
    """
    if species is None:
        species_list = list(tracked_species)
    elif isinstance(species, str):
        species_list = [species]
    else:
        species_list = list(species)

    for sp in species_list:
        if sp not in tracked_species:
            raise ValueError(f"'{sp}' not in tracked_species={tracked_species}")

    keys = sorted(results.keys()) if embryos is None else list(embryos)
    keys = [r_au for r_au in keys if len(results[r_au]) > 0]

    t_factor, t_label = _resolve_time_unit(time_unit)

    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 5))
    else:
        fig = ax.figure

    def _fraction(hist, sp):
        col_idx = 3 + tracked_species.index(sp)
        total = hist[:, 3:].sum(axis=1)
        return 100 * np.divide(hist[:, col_idx], total, out=np.zeros_like(total), where=total > 0)

    multi_species = len(species_list) > 1

    if multi_species and len(keys) == 1:
        # Un solo embrión, varias especies: el color identifica a la especie.
        r_au = keys[0]
        hist = results[r_au]
        t = hist[:, 0] / t_factor
        tab10 = plt.get_cmap('tab10')
        for i, sp in enumerate(species_list):
            ax.plot(t, _fraction(hist, sp), color=tab10(i % 10), lw=2, label=sp)
        ax.legend(loc='best', title='Species')
        title = f"Composition Evolution at {r_au:g} AU"
    else:
        # Color por embrión; si hay varias especies, estilo de línea por especie.
        colors, mode, sm = _radius_colors(keys, cmap=cmap)
        for r_au, color in zip(keys, colors):
            hist = results[r_au]
            t = hist[:, 0] / t_factor
            for i, sp in enumerate(species_list):
                ls = _SPECIES_LINESTYLES[i % len(_SPECIES_LINESTYLES)]
                label = f"{r_au:g} AU" if (mode == 'discrete' and i == 0) else None
                ax.plot(t, _fraction(hist, sp), color=color, lw=2, ls=ls, label=label)

        if mode == 'discrete':
            emb_legend = ax.legend(loc='upper left', title='Embryo')
            if multi_species:
                ax.add_artist(emb_legend)
        elif sm is not None:
            fig.colorbar(sm, ax=ax, label='Initial Radius [AU]')

        if multi_species:
            proxies = [Line2D([0], [0], color='gray', lw=2,
                              ls=_SPECIES_LINESTYLES[i % len(_SPECIES_LINESTYLES)])
                       for i in range(len(species_list))]
            ax.legend(proxies, species_list, loc='upper right', title='Species')

        title = (f"{species_list[0]} Fraction Evolution" if not multi_species
                 else "Composition Evolution")

    ax.set_xscale('log')
    ax.set_xlabel(t_label)
    ax.set_ylabel(f"{species_list[0]} Mass Fraction [%]" if not multi_species
                  else "Species Mass Fraction [%]")
    ax.set_title(title)
    ax.grid(True, which='both', linestyle=':', alpha=0.6)
    if ylim is not None:
        ax.set_ylim(*ylim)

    fig.tight_layout()
    return fig, ax


def plot_growth_summary(results: dict, tracked_species: list, species: str = 'H2O',
                         embryos: Optional[list] = None, time_unit: str = 'Myr',
                         show_isolation_mass: bool = True, cmap: str = 'viridis',
                         figsize: tuple = (14, 5)):
    """
    Figura de 2 paneles: crecimiento de masa (izquierda) + fracción de una
    especie (derecha), para una o varias trayectorias de embriones.

    Ver `plot_growth_curves` y `plot_species_fraction` para el detalle de
    los parámetros compartidos.

    Retorna:
    --------
    (fig, (ax_mass, ax_frac))
    """
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=figsize)
    plot_growth_curves(results, embryos=embryos, time_unit=time_unit,
                        show_isolation_mass=show_isolation_mass, cmap=cmap, ax=ax1)
    plot_species_fraction(results, tracked_species, species=species, embryos=embryos,
                           time_unit=time_unit, cmap=cmap, ax=ax2)
    fig.tight_layout()
    return fig, (ax1, ax2)
