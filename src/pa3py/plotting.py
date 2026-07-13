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
    """
    Genera un diagrama de Hovmöller (Radio vs Tiempo) para una propiedad del disco.

    Parámetros:
    -----------
    disk : DiskData
        El objeto que contiene los datos temporales del disco.
    field : str
        El campo a graficar. Opciones: 'dust_Sigma', 'gas_Sigma', 'eps' (relación polvo/gas).
    cmap : str
        Mapa de colores para matplotlib.
    vmin, vmax : float, opcional
        Límites de la escala de color logarítmica. Si no se dan, usan el percentil 2 y 98.
    t_unit : str
        Unidad de tiempo ('yr', 'kyr', 'Myr').
    show_snowlines : bool
        Si True, sobrepone las líneas de hielo extraídas de disk.hdf5_snowlines.
    """
    # 1. Extraer los datos 2D (Tiempo, Radio)
    if field == 'dust_Sigma':
        # Sumamos sobre todos los tamaños de polvo (último índice)
        Z = np.sum(disk.dust_Sigma, axis=-1)
        title = r"Hovmöller: $\Sigma_{dust}$ [g/cm$^2$]"
    elif field == 'gas_Sigma':
        Z = disk.gas_Sigma
        title = r"Hovmöller: $\Sigma_{gas}$ [g/cm$^2$]"
    elif field == 'eps':
        dust = np.sum(disk.dust_Sigma, axis=-1)
        gas = disk.gas_Sigma
        Z = dust / (gas + 1e-30)
        title = r"Hovmöller: $\epsilon = \Sigma_{dust}/\Sigma_{gas}$"
    else:
        raise ValueError("Field no soportado. Opciones: 'dust_Sigma', 'gas_Sigma', 'eps'.")

    # Limpiar ceros para escala log
    Z = np.clip(Z, 1e-20, None)
    
    # Eje Temporal
    if t_unit == 'Myr':
        t_factor = c.YEAR * 1e6
        t_label = "Tiempo [Myr]"
    elif t_unit == 'kyr':
        t_factor = c.YEAR * 1e3
        t_label = "Tiempo [kyr]"
    else:
        t_factor = c.YEAR
        t_label = "Tiempo [años]"
        
    t_array = disk.times / t_factor
    r_array = disk.r / c.AU

    # Ejes extendidos para pcolormesh
    def log_edges(x):
        # Asume x escala log. Retorna bordes.
        logx = np.log10(x)
        dlogx = np.diff(logx)
        dlogx = np.append(dlogx, dlogx[-1])
        edges = 10**(logx - dlogx/2)
        edges = np.append(edges, 10**(logx[-1] + dlogx[-1]/2))
        return edges

    # Para el tiempo evitamos tomar log10(0) si time[0]==0
    safe_t = np.maximum(t_array, 1e-6)
    t_edges = log_edges(safe_t)
    r_edges = log_edges(r_array)

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
    
    # 3. Dibujar Snowlines si están presentes
    if show_snowlines and disk.hdf5_snowlines:
        colors = ['cyan', 'white', 'lightgreen']
        for idx, (name, rsnow) in enumerate(disk.hdf5_snowlines.items()):
            color = colors[idx % len(colors)]
            r_au = rsnow / c.AU
            # Solo dibujar si cae dentro del rango y no es todo cero
            if np.any(r_au > 0):
                ax.plot(safe_t, r_au, ls='--', color=color, lw=2, label=f"Snowline {name}")
        ax.legend(loc='lower left')
                
    fig.tight_layout()
    return fig, ax

def plot_hovmoller_epsilon(t_yr, r_au, Sigma_dust_2D, Sigma_gas_2D, fig_path="hovmoller_eps.png"):
    """
    t_yr: arreglo 1D con los tiempos en años (ej. de shape (Nt,))
    r_au: arreglo 1D con los radios en AU (ej. de shape (Nr,))
    Sigma_dust_2D: matriz de polvo de shape (Nt, Nr)
    Sigma_gas_2D: matriz de gas de shape (Nt, Nr)
    """
    
    # 1. Calcular Epsilon (polvo/gas). 
    # Le sumamos 1e-300 al gas para evitar división por cero.
    epsilon = Sigma_dust_2D / (Sigma_gas_2D + 1e-300)
    epsilon = np.clip(epsilon, 1e-10, None)  # Limitar valores muy pequeños
    
    # 2. Función auxiliar para calcular los BORDES de las celdas logarítmicas
    def _log_edges(arr):
        arr = np.maximum(arr, 1e-10) # Evitar log(0)
        log_arr = np.log10(arr)
        # Bordes intermedios
        mid = (log_arr[:-1] + log_arr[1:]) / 2.0
        # Extrapolar bordes iniciales y finales
        first = log_arr[0] - (mid[0] - log_arr[0])
        last  = log_arr[-1] + (log_arr[-1] - mid[-1])
        edges_log = np.concatenate([[first], mid, [last]])
        return 10**edges_log

    # Bordes de Tiempo (Eje X) y Radio (Eje Y)
    # Reemplazamos 0 por 1.0 al inicio si es necesario para el log
    t_yr_safe = np.maximum(t_yr, 1.0) 
    t_edges = _log_edges(t_yr_safe)
    r_edges = _log_edges(r_au)
    
    # 3. Graficar
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # IMPORTANTE: epsilon.T invierte de (Nt, Nr) a (Nr, Nt) para que cuadre
    # con (X=tiempo, Y=radio).
    pcm = ax.pcolormesh(
        t_edges, r_edges, epsilon.T,
        norm=LogNorm(vmin=1e-4, vmax=1e-1), # Ajusta los límites de color a tu gusto
        cmap='magma', shading='flat'
    )
    
    cbar = fig.colorbar(pcm, ax=ax, pad=0.02)
    cbar.set_label(r"Relación polvo/gas ($\epsilon = \Sigma_d / \Sigma_g$)")

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim(t_yr_safe[0], t_yr_safe[-1])
    ax.set_ylim(r_au[0], r_au[-1])
    
    ax.set_xlabel("Tiempo [años]")
    ax.set_ylabel("Distancia Radial [AU]")
    ax.set_title(r"Diagrama Hovmöller: $\epsilon$")
    
    plt.tight_layout()
    fig.savefig(fig_path, dpi=300, bbox_inches='tight')
    plt.show()
    plt.close(fig)
