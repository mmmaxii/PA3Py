"""
data.py - Módulo para la carga y manejo de datos del disco protoplanetario.
"""

import os
import glob
import h5py
import numpy as np
from dataclasses import dataclass
from typing import Dict, List, Optional
from . import constants as c

@dataclass
class DiskData:
    """
    Contenedor agnóstico para las propiedades físicas del disco (gas y polvo).
    Los arrays 1D tienen dimensión (Nr,) y los temporales/radiales (Nt, Nr).
    """
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
    """
    Lee los archivos HDF5 generados por tripodpy y los convierte
    en un objeto DiskData.
    """
    
    files = sorted(glob.glob(os.path.join(datadir, 'data*.hdf5')))
    if not files:
        raise FileNotFoundError(f"No se encontraron archivos HDF5 en {datadir}")

    print(f"[load_tripodpy_hdf5] Leyendo {len(files)} snapshots desde {datadir}...")
    
    times_list = []
    rsnow = {'H2O': []}
    
    gas_S, gas_T, gas_cs, gas_eta, gas_nu, gas_alpha, gas_Hp = [], [], [], [], [], [], []
    dust_Sigma, dust_vr, dust_St, dust_H = [], [], [], []
    OmegaK_list = []
    r_grid = None

    for fpath in files:
        with h5py.File(fpath, 'r') as f:
            t_s = float(f['t'][()])
            if t_s < t_min_yr * 3.156e7:
                continue

            times_list.append(t_s)

            if r_grid is None:
                r_grid = f['grid/r'][:]

            # Gas
            gas_S.append(f['gas/Sigma'][:])
            gas_T.append(f['gas/T'][:])
            gas_cs.append(f['gas/cs'][:])
            gas_eta.append(f['gas/eta'][:])
            gas_nu.append(f['gas/nu'][:])
            gas_alpha.append(f['gas/alpha'][:])
            gas_Hp.append(f['gas/Hp'][:])

            # Dust
            dust_Sigma.append(f['dust/Sigma'][:])
            dust_vr.append(f['dust/v/rad'][:])
            dust_St.append(f['dust/St'][:])
            dust_H.append(f['dust/H'][:])

            # OmegaK
            if 'grid/OmegaK' in f:
                OmegaK_list.append(f['grid/OmegaK'][:])

            # Snowline dinámica de H2O original de la corrida
            if 'dust/r_snow' in f:
                rsnow['H2O'].append(float(f['dust/r_snow'][()]))
            else:
                key = f'grid/rsnow_H2O'
                rsnow['H2O'].append(float(f[key][()]) if key in f else np.nan)

    times = np.array(times_list)
    
    # OmegaK (2D o 1D analítico)
    if OmegaK_list:
        Omega_K = np.array(OmegaK_list)
    else:
        # G está en CGS en constants, M_star está en masas solares
        G_M = c.G * (M_star * c.M_SUN)
        Omega_K = np.sqrt(G_M / r_grid**3)

    rsnow['H2O'] = np.array(rsnow['H2O'])

    return DiskData(
        times=times,
        r=r_grid,
        gas_Sigma=np.array(gas_S),
        gas_T=np.array(gas_T),
        gas_cs=np.array(gas_cs),
        gas_eta=np.array(gas_eta),
        gas_nu=np.array(gas_nu),
        gas_alpha=np.array(gas_alpha),
        gas_Hp=np.array(gas_Hp),
        dust_Sigma=np.array(dust_Sigma),
        dust_vr=np.array(dust_vr),
        dust_St=np.array(dust_St),
        dust_H=np.array(dust_H),
        Omega_K=Omega_K,
        M_star=M_star,
        hdf5_snowlines=rsnow
    )
