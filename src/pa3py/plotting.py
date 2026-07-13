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
        title = r"Hovmöller: Relación polvo/gas ($\epsilon = \Sigma_d / \Sigma_g$)"
    else:
        raise ValueError("Field no soportado. Opciones: 'dust_Sigma', 'gas_Sigma', 'epsilon'.")

    # Limpiar ceros para escala log
    Z = np.clip(Z, 1e-20, None)
    
    t_map = {'Myr': (c.YEAR * 1e6, "Tiempo [Myr]"), 'kyr': (c.YEAR * 1e3, "Tiempo [kyr]")}
    t_factor, t_label = t_map.get(t_unit, (c.YEAR, "Tiempo [años]"))
        
    t_array = disk.times / t_factor
    r_array = disk.r / c.AU

    # Ejes extendidos robustos para pcolormesh
    def _log_edges(arr):
        arr = np.maximum(arr, 1e-10)
        log_arr = np.log10(arr)
        mid = (log_arr[:-1] + log_arr[1:]) / 2.0
        return 10**np.concatenate([[log_arr[0] - (mid[0] - log_arr[0])], mid, [log_arr[-1] + (log_arr[-1] - mid[-1])]])

    # Prevenir log10(0)
    safe_t = np.maximum(t_array, 1.0)
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
    cbar.set_label("Escala Logarítmica")
    
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim(safe_t[0], safe_t[-1])
    ax.set_ylim(r_array[0], r_array[-1])
    ax.set_xlabel(t_label)
    ax.set_ylabel("Radio [AU]")
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


