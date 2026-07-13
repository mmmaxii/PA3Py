"""
plotting.py - Herramientas de visualización (ej: Diagrama Hovmöller)
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from .data import DiskData
from . import constants as c

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
    
    t_map = {'Myr': (c.YEAR * 1e6, "Time [Myr]"), 'kyr': (c.YEAR * 1e3, "Time [kyr]")}
    t_factor, t_label = t_map.get(t_unit, (c.YEAR, "Time [Years]"))
        
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
