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

# Importar constantes y tipos
from . import constants as c
from .data import DiskData
from .composition import CompositionModel, SimpleWaterComposition


class PebbleAccretionModule3:
    """
    Simula la acreción de pebbles sobre embriones planetarios.
    Version actualizada a física pura Ormel 2017 y revisión PA.
    
    Esta clase ahora es completamente agnóstica de los archivos fuente (HDF5), 
    y depende de un objeto `DiskData` y un modelo de composición `CompositionModel`.
    """

    # ── Constantes físicas (CGS) ──────────────────────────────────────────
    G_CGS   = c.G        # cm³ g⁻¹ s⁻²
    M_SUN   = c.M_SUN    # g
    M_EARTH = c.M_EARTH  # g
    AU      = c.AU       # cm

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
            # Fallback a la clásica water snowline de los datos base
            rsnow_h2o = self.data.hdf5_snowlines.get('H2O', np.zeros(self.data.Nt))
            self.comp = SimpleWaterComposition(rsnow_h2o)
        else:
            self.comp = comp_model
            

    # ══════════════════════════════════════════════════════════════════════
    # Helpers Interp
    # ══════════════════════════════════════════════════════════════════════

    def _interp(self, field_1d: np.ndarray, r_emb: float) -> float:
        """
        Interpola logarítmicamente un campo radial 1D.
        """
        return float(np.interp(np.log(r_emb), np.log(self.data.r), field_1d))

    def _local(self, t: int, r_emb: float) -> dict:
        """
        Extrae e interpola las propiedades locales del disco en t_idx y r_emb.
        """
        I = lambda arr: self._interp(arr[t], r_emb)
        
        Sigma_peb = self._interp(self.data.dust_Sigma[t, :, self.peb_idx], r_emb)
        eta   = I(self.data.gas_eta)
        St    = self._interp(self.data.dust_St[t, :, self.peb_idx], r_emb)
        H_gas = I(self.data.gas_Hp)
        H_peb = self._interp(self.data.dust_H[t, :, self.peb_idx], r_emb)
        
        # OmegaK puede ser 1D o 2D
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
        """
        M_iso_peb = 25 M_earth * (H_gas/r / 0.05)^3 * (M_star/M_sun)
        """
        p = self._local(t, r_emb)
        h_gas = p['H_gas'] / r_emb
        M_iso = 25 * (h_gas / 0.05)**3 * self.data.M_star * self.M_EARTH
        return max(M_iso, 0.1 * self.M_EARTH)

    def _accretion_rate(self, M_core: float, r_emb: float, t: int) -> float:
        if M_core <= 0: return 0.0
        p = self._local(t, r_emb)
        if p['Sigma_peb'] < 1e-30: return 0.0

        G, M, Omega = self.G_CGS, M_core, p['Omega']
        St, v_hw, Sigma = p['St'], p['v_hw'], p['Sigma_peb']
        
        t_stop = St / Omega
        M_PA_onset = St * (p['eta']**3) * (self.data.M_star * self.M_SUN)
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
            R_H     = r_emb * (M / (3 * (self.data.M_star * self.M_SUN)))**(1/3)
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
        if M0_g is None:
            M0_g = 1e-3 * self.M_EARTH
            
        locs_outer_to_inner = sorted(embryo_locations_AU, reverse=True)
        M_core   = {r: float(M0_g) for r in locs_outer_to_inner}
        
        # Inicializar composición (semilla 100% de la especie "rocosa" primaria)
        species = self.comp.get_species()
        primary_rock = "silicates" if "silicates" in species else species[0]
        
        M_X = {r: {sp: 0.0 for sp in species} for r in locs_outer_to_inner}
        for r in locs_outer_to_inner:
            M_X[r][primary_rock] = float(M0_g)
            
        active   = {r: True for r in locs_outer_to_inner}
        histories= {r: [] for r in locs_outer_to_inner}

        for i in range(self.data.Nt):
            dt = float(self.data.times[i] - (self.data.times[i-1] if i > 0 else 0.0))
            if dt <= 0: continue

            flux_consumed = 0.0

            for r_au in locs_outer_to_inner:
                if not active[r_au]: continue
                r_emb = r_au * self.AU

                M_iso = self._isolation_mass(r_emb, i)
                if M_core[r_au] >= M_iso:
                    active[r_au] = False
                    continue

                Mdot_peb_disk  = self._pebble_flux(i, r_emb)
                Mdot_peb_avail = max(0.0, Mdot_peb_disk - flux_consumed)

                Mdot_core_r = self._accretion_rate(M_core[r_au], r_emb, i)
                Mdot_core_r = min(Mdot_core_r, Mdot_peb_avail)

                dM = Mdot_core_r * dt
                dM = min(dM, max(0.0, M_iso - M_core[r_au]))
                flux_consumed += Mdot_core_r

                # Calcular partición de masa acretada usando el modelo de composición
                fractions = self.comp.get_fractions(r_emb, i)
                for sp in species:
                    M_X[r_au][sp] += fractions.get(sp, 0.0) * dM

                M_core[r_au] += dM

                # Guardamos: t, M_core, M_H2O, M_silicates (por legacy) o M_X completo
                histories[r_au].append((
                    self.data.times[i], 
                    M_core[r_au], 
                    M_X[r_au].get('H2O', 0.0),
                    M_X[r_au].get('silicates', 0.0), 
                    M_iso
                ))

        results = {
            r_au: (np.array(histories[r_au]) if histories[r_au] else np.empty((0, 5)))
            for r_au in embryo_locations_AU
        }
        return results

    def summary(self, results: dict):
        """Imprime tabla resumen de composición final con M_iso."""
        SEP = '-' * 80
        print(f"\n{SEP}")
        print(f"{'r [AU]':>8} {'M_tot [ME]':>11} {'M_iso [ME]':>11} {'f_H2O [%]':>10} {'f_Sil [%]':>10}  Tipo")
        print(SEP)
        
        for r_au, hist in results.items():
            if len(hist) == 0:
                print(f"{r_au:>8.2f}  -- sin acreción")
                continue
            
            # hist[-1] format: [time, M_core, M_H2O, M_sil, M_iso]
            _, M, M_H2O, M_sil, M_iso = hist[-1]
            
            M_total = M_H2O + M_sil
            if M_total <= 0:
                f_h2o, f_sil = 0.0, 100.0
            else:
                f_h2o = 100 * M_H2O / M_total
                f_sil = 100 * M_sil / M_total
            
            tipo = "[W] Waterworld" if f_h2o > 10 else "[R] Rocoso"
            print(f"{r_au:>8.2f}  {M/self.M_EARTH:>11.3f}  {M_iso/self.M_EARTH:>11.2f}  {f_h2o:>9.1f}  {f_sil:>9.1f}  {tipo}")
        print(f"{SEP}\n")
