"""
data.py - Módulo para la carga y manejo de datos del disco protoplanetario.
"""

import os
import glob
import h5py
import numpy as np
from dataclasses import dataclass
from typing import Dict
from . import constants as c

@dataclass
class DiskData:
    """Contenedor de propiedades físicas del disco (gas y polvo)."""
    times: np.ndarray        # (Nt,) Tiempos en segundos
    r: np.ndarray            # (Nr,) Grilla radial en cm
    
    # Gas: shape (Nt, Nr)
    gas_Sigma: np.ndarray    # Densidad superficial del gas
    gas_T: np.ndarray        # Temperatura del gas
    gas_cs: np.ndarray       # Velocidad del sonido
    gas_eta: np.ndarray      # Gradiente de presión
    gas_nu: np.ndarray       # Viscosidad cinemática
    gas_alpha: np.ndarray    # Parámetro alpha de Shakura-Sunyaev
    gas_Hp: np.ndarray       # Escala de altura del gas
    
    # Dust: shape (Nt, Nr, Nd) donde Nd=5 (o Nd=1 si ya se colapsó)
    dust_Sigma: np.ndarray   # Densidad superficial del polvo
    dust_vr: np.ndarray      # Velocidad radial del polvo
    dust_St: np.ndarray      # Número de Stokes
    dust_H: np.ndarray       # Escala de altura del polvo
    
    # Derivados / Auxiliares
    Omega_K: np.ndarray      # Frecuencia kepleriana local (Nt, Nr) o (Nr,)
    M_star: float            # Masa estelar en masas solares
    
    # Snowlines pre-calculadas presentes en el HDF5 (Opcional, útil para fallback)
    # Formato: {"H2O": array(Nt,)}
    hdf5_snowlines: Dict[str, np.ndarray]

    @property
    def Nt(self):
        return len(self.times)

    @property
    def Nr(self):
        return len(self.r)


def load_tripodpy_hdf5(datadir: str, M_star: float = 1.0, t_min_yr: float = 0.0) -> DiskData:
    """Convierte archivos HDF5 de tripodpy a un objeto DiskData."""
    
    files = sorted(glob.glob(os.path.join(datadir, 'data*.hdf5')))
    if not files:
        raise FileNotFoundError(f"No se encontraron archivos HDF5 en {datadir}")

    print(f"[load_tripodpy_hdf5] Leyendo {len(files)} snapshots desde {datadir}...")
    
    times_list, OmegaK_list = [], []
    rsnow = {'H2O': []}
    r_grid = None
    
    gas_keys = ['Sigma', 'T', 'cs', 'eta', 'nu', 'alpha', 'Hp']
    dust_keys = ['Sigma', 'v/rad', 'St', 'H']
    gas_data = {k: [] for k in gas_keys}
    dust_data = {k: [] for k in dust_keys}

    for fpath in files:
        with h5py.File(fpath, 'r') as f:
            t_s = float(f['t'][()])
            if t_s < t_min_yr * c.YEAR:
                continue

            times_list.append(t_s)

            if r_grid is None:
                r_grid = f['grid/r'][:]

            # Cargar arrays
            for k in gas_keys: gas_data[k].append(f[f'gas/{k}'][:])
            for k in dust_keys: dust_data[k].append(f[f'dust/{k}'][:])

            # OmegaK
            if 'grid/OmegaK' in f:
                OmegaK_list.append(f['grid/OmegaK'][:])

            # Snowline dinámica de H2O original de la corrida
            if 'dust/r_snow' in f:
                rsnow['H2O'].append(float(f['dust/r_snow'][()]))
            else:
                rsnow_key = 'grid/rsnow_H2O'
                rsnow['H2O'].append(float(f[rsnow_key][()]) if rsnow_key in f else np.nan)

    times = np.array(times_list)
    
    # OmegaK: 2D del HDF5 o analítico 1D
    if OmegaK_list:
        Omega_K = np.array(OmegaK_list)
    else:
        Omega_K = np.sqrt(c.G * M_star * c.M_SUN / r_grid**3)

    rsnow['H2O'] = np.array(rsnow['H2O'])

    return DiskData(
        times=times,
        r=r_grid,
        gas_Sigma=np.array(gas_data['Sigma']),
        gas_T=np.array(gas_data['T']),
        gas_cs=np.array(gas_data['cs']),
        gas_eta=np.array(gas_data['eta']),
        gas_nu=np.array(gas_data['nu']),
        gas_alpha=np.array(gas_data['alpha']),
        gas_Hp=np.array(gas_data['Hp']),
        dust_Sigma=np.array(dust_data['Sigma']),
        dust_vr=np.array(dust_data['v/rad']),
        dust_St=np.array(dust_data['St']),
        dust_H=np.array(dust_data['H']),
        Omega_K=Omega_K,
        M_star=M_star,
        hdf5_snowlines=rsnow
    )
