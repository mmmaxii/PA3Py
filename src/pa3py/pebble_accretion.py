"""
PA3Py/pebble_accretion.py — Módulo de acreción de pebbles (Ormel 2017 & Drążkowska et al. 2023)
=====================================================================================

Física implementada y referencias:
  - M_PA_onset = St * η³ * M_star             (Drążkowska et al. 2023 Eq. 3, Ormel 2017 Eq 7.11)
  - M_hw/sh   = v_hw³ / (8 G Ω_K t_stop)      (Ormel 2017 Eq 7.9)
  - Ṁ_2D_hw   = √(8GM t_stop v_hw) Σ_peb      (Ormel 2017 Eq 7.13 Headwind)
  - Ṁ_2D_sh   = 2 R_H² Ω_K St^(2/3) Σ_peb     (Ormel 2017 Eq 7.13 Shear, Drążkowska et al. 2023 Eq 5)
  - Ṁ_3D      = 2π G M t_stop ρ_peb           (Ormel 2017 Eq 7.12)
  - Transición = Ṁ_2D * b_col / (b_col + H_peb √(8/π)) (Ormel 2017 Eq 7.24)
  - M_iso_peb = 25 M⊕ (H_gas/r / 0.05)³ (M_star/M_sun) (Drążkowska et al. 2023 Eq 6)
  - M < M_onset: Acreción de Safronov Balística (Ormel 2017 Eq 7.14)

Unidades internas: CGS (g, cm, s).
"""

import numpy as np

from .data import DiskData
from .composition import CompositionModel, SimpleWaterComposition
from . import constants as c


class PebbleAccretionModule3:
    """Motor físico de acreción agnóstico a HDF5 (Ormel 2017 & Drążkowska 2023)."""

    # Índice único de pebbles. (El último bin de polvo representa los pebbles)
    peb_idx = -1

    def __init__(self, disk_data: DiskData, comp_model: CompositionModel = None):
        """
        Inicializa el simulador con los datos del disco y un modelo de composición.
        
        Parámetros:
        -----------
        disk_data : DiskData
            Contenedor con todas las propiedades radiales y temporales del disco.
        comp_model : CompositionModel, opcional
            Modelo que define las snowlines y abundancias por regiones.
            Si no se provee, se usará el modelo clásico de PA3Py (SimpleWaterComposition)
            usando la snowline de agua que venga del HDF5 (si existe).
        """
        self.data = disk_data
        if comp_model is None:
            rsnow_h2o = self.data.hdf5_snowlines.get('H2O', np.zeros(self.data.Nt))
            self.comp = SimpleWaterComposition(rsnow_h2o)
        else:
            self.comp = comp_model
        self.tracked_species = self.comp.get_species()

    # ══════════════════════════════════════════════════════════════════════
    # Helpers Interp
    # ══════════════════════════════════════════════════════════════════════

    def _interp(self, field_1d: np.ndarray, r_emb: float) -> float:
        """
        Interpola logarítmicamente un campo radial 1D.
        """
        return float(np.interp(np.log(r_emb), np.log(self.data.r), field_1d))

    def _local(self, t: int, r_emb: float) -> dict:
        """Extrae e interpola las propiedades locales del disco."""
        I = lambda arr: self._interp(arr[t], r_emb)
        
        Sigma_peb = self._interp(self.data.dust_Sigma[t, :, self.peb_idx], r_emb)
        eta   = I(self.data.gas_eta)
        St    = self._interp(self.data.dust_St[t, :, self.peb_idx], r_emb)
        H_gas = I(self.data.gas_Hp)
        H_peb = self._interp(self.data.dust_H[t, :, self.peb_idx], r_emb)
        
        if self.data.Omega_K.ndim == 2:
            omega_1d = self.data.Omega_K[t]
        else:
            omega_1d = self.data.Omega_K
            
        Omega = float(np.interp(np.log(r_emb), np.log(self.data.r), omega_1d))
        v_K   = Omega * r_emb
        v_hw  = eta * v_K
        
        return dict(
            Sigma_peb=Sigma_peb, eta=eta, St=St, H_gas=H_gas, H_peb=H_peb, 
            Omega=Omega, v_K=v_K, v_hw=v_hw
        )

    # ══════════════════════════════════════════════════════════════════════
    # Física Ormel 2017 & Drążkowska et al. 2023
    # ══════════════════════════════════════════════════════════════════════

    def _pebble_flux(self, t: int, r_emb: float) -> float:
        """Ṁ_peb = 2π r Σ_peb |v_r|"""
        Sigma_peb = self._interp(self.data.dust_Sigma[t, :, self.peb_idx], r_emb)
        v_r       = self._interp(self.data.dust_vr[t, :, self.peb_idx], r_emb)
        return 2 * np.pi * r_emb * Sigma_peb * abs(v_r)

    def _isolation_mass(self, r_emb: float, t: int) -> float:
        """Drążkowska 2023 Eq.6: M_iso = 25 M⊕ (h/0.05)³ (M★/M☉)"""
        h_gas = self._interp(self.data.gas_Hp[t], r_emb) / r_emb
        M_iso = 25 * (h_gas / 0.05)**3 * self.data.M_star * c.M_EARTH
        return max(M_iso, 0.1 * c.M_EARTH)

    def _accretion_rate(self, M_core: float, r_emb: float, t: int) -> float:
        if M_core <= 0: return 0.0
        p = self._local(t, r_emb)
        if p['Sigma_peb'] < 1e-30: return 0.0

        G, M, Omega = c.G, M_core, p['Omega']
        St, v_hw, Sigma = p['St'], p['v_hw'], p['Sigma_peb']
        
        t_stop = St / Omega
        M_PA_onset = St * (p['eta']**3) * (self.data.M_star * c.M_SUN)
        H_peb = max(p['H_peb'], 1e-10 * p['H_gas'])
        rho_peb = Sigma / (np.sqrt(2 * np.pi) * H_peb)

        # Safronov (Balístico)
        if M < M_PA_onset:
            rho_core = 3.0  # g/cm3
            R_pl = (3 * M / (4 * np.pi * rho_core))**(1/3)
            v_impact = max(v_hw, 1.0) 
            return (2 * np.pi * R_pl * G * M / v_impact) * rho_peb

        # Transición Headwind/Shear
        M_hw_sh = (v_hw**3) / (8 * G * Omega * St)

        # Regímenes 2D
        if M < M_hw_sh:
            Mdot_2D = np.sqrt(8 * G * M * t_stop * v_hw) * Sigma
            b_col   = np.sqrt(2 * G * M * t_stop / max(v_hw, 1e-5))
        else:
            R_H     = r_emb * (M / (3 * self.data.M_star * c.M_SUN))**(1/3)
            Mdot_2D = 2 * R_H**2 * Omega * St**(2/3) * Sigma
            b_col   = (St**(1/3)) * R_H
            
        # Transición 2D-3D turbulencia suave
        denominator = b_col + H_peb * np.sqrt(8.0 / np.pi)
        Mdot_eff = Mdot_2D * (b_col / denominator)
        
        return max(Mdot_eff, 0.0)

    # ══════════════════════════════════════════════════════════════════════
    # API / Loop Principal
    # ══════════════════════════════════════════════════════════════════════

    def run_growth(self, embryo_locations_AU: list, M0_g: float = None) -> dict:
        r_min_au, r_max_au = self.data.r[0] / c.AU, self.data.r[-1] / c.AU
        for r_au in embryo_locations_AU:
            if not (r_min_au <= r_au <= r_max_au):
                raise ValueError(f"El embrión en {r_au:.2f} AU fuera del disco ({r_min_au:.2f}–{r_max_au:.2f} AU).")

        if M0_g is None:
            M0_g = 1e-3 * c.M_EARTH
            
        locs_outer_to_inner = sorted(embryo_locations_AU, reverse=True)
        M_core = {r: float(M0_g) for r in locs_outer_to_inner}
        primary_rock = "silicates" if "silicates" in self.tracked_species else self.tracked_species[0]
        
        M_X = {r: {sp: float(M0_g) if sp == primary_rock else 0.0 for sp in self.tracked_species} for r in locs_outer_to_inner}
            
        histories= {r: [] for r in locs_outer_to_inner}

        for i in range(self.data.Nt):
            dt = float(self.data.times[i] - (self.data.times[i-1] if i > 0 else 0.0))
            if dt <= 0: continue

            flux_consumed = 0.0

            for r_au in locs_outer_to_inner:
                r_emb = r_au * c.AU

                M_iso = self._isolation_mass(r_emb, i)
                
                if M_core[r_au] < M_iso:
                    Mdot_peb_disk  = self._pebble_flux(i, r_emb)
                    Mdot_peb_avail = max(0.0, Mdot_peb_disk - flux_consumed)

                    Mdot_core_r = self._accretion_rate(M_core[r_au], r_emb, i)
                    Mdot_core_r = min(Mdot_core_r, Mdot_peb_avail)

                    dM = Mdot_core_r * dt
                    dM = min(dM, max(0.0, M_iso - M_core[r_au]))
                    flux_consumed += Mdot_core_r

                    fractions = self.comp.get_fractions(r_emb, self.data.times[i], i)
                    for sp in self.tracked_species:
                        M_X[r_au][sp] += fractions.get(sp, 0.0) * dM

                    M_core[r_au] += dM

                histories[r_au].append([self.data.times[i], M_core[r_au], M_iso] + [M_X[r_au].get(sp, 0.0) for sp in self.tracked_species])

        results = {
            r_au: (np.array(histories[r_au]) if histories[r_au] else np.empty((0, 3 + len(self.tracked_species))))
            for r_au in embryo_locations_AU
        }
        return results

    def summary(self, results: dict):
        """Imprime tabla resumen de la composición final."""
        col_widths = [max(10, len(f"f_{sp}[%]") + 2) for sp in self.tracked_species]
        SEP = '-' * (35 + sum(col_widths))
        print(f"\n{SEP}")
        header = f"{'r [AU]':>8} {'M_tot [ME]':>11} {'M_iso [ME]':>11}"
        header += "".join(f"{f'f_{sp}[%]':>{w}}" for sp, w in zip(self.tracked_species, col_widths))
        print(header)
        print(SEP)
        
        for r_au, hist in results.items():
            if len(hist) == 0:
                print(f"{r_au:>8.2f}  -- no accretion")
                continue
            
            _, M, M_iso, *M_species = hist[-1]
            
            M_total = sum(M_species)
            
            line = f"{r_au:>8.2f}  {M/c.M_EARTH:>11.3f}  {M_iso/c.M_EARTH:>11.2f}"
            for i, m_sp in enumerate(M_species):
                f_sp = 0.0 if M_total <= 0 else 100 * m_sp / M_total
                line += f"{f_sp:>{col_widths[i]}.1f}"
                
            print(line)
        print(f"{SEP}\n")

    def calculate_isolation_mass_map(self) -> np.ndarray:
        """Calcula el mapa de masa de aislamiento teórico 2D."""
        h_gas = self.data.gas_Hp / self.data.r
        M_iso = 25 * (h_gas / 0.05)**3 * self.data.M_star * c.M_EARTH
        return np.maximum(M_iso, 0.1 * c.M_EARTH)
