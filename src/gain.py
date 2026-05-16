#!/usr/bin/env python3
"""Material gain estimates built on the compact MQW k.p solver."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np

from BasicMQWDesign import HC_EV_UM
from kp_solver import KpSubbands, MQWProfile

KB_EV_K = 8.617333262145e-5
NM2_TO_CM2 = 1.0e14


@dataclass(frozen=True)
class GainSpectrum:
    energy_eV: np.ndarray
    wavelength_nm: np.ndarray
    gain_te_cm: np.ndarray
    gain_tm_cm: np.ndarray
    quasi_fermi_e_eV: float
    quasi_fermi_h_eV: float
    peak_te_gain_cm: float
    peak_te_wavelength_nm: float
    peak_tm_gain_cm: float
    peak_tm_wavelength_nm: float


def _fermi(energy_eV: np.ndarray | float, fermi_eV: float, kbt_eV: float):
    x = (np.asarray(energy_eV) - fermi_eV) / kbt_eV
    x = np.clip(x, -80.0, 80.0)
    return 1.0 / (1.0 + np.exp(x))


def _integrate_density_cm2(
    kt_nm: np.ndarray,
    energies_by_k: list[list[float]],
    fermi_eV: float,
    kbt_eV: float,
) -> float:
    density_nm2 = 0.0
    for band_index in range(min(len(row) for row in energies_by_k)):
        energies = np.array([row[band_index] for row in energies_by_k], dtype=float)
        occ = _fermi(energies, fermi_eV, kbt_eV)
        density_nm2 += float(np.trapezoid(occ * kt_nm, kt_nm)) / np.pi
    return density_nm2 * NM2_TO_CM2


def _solve_fermi_level(
    kt_nm: np.ndarray,
    energies_by_k: list[list[float]],
    target_density_cm2: float,
    kbt_eV: float,
) -> float:
    flat = np.array([energy for row in energies_by_k for energy in row], dtype=float)
    lo = float(np.min(flat) - 1.0)
    hi = float(np.max(flat) + 1.0)
    for _ in range(120):
        mid = 0.5 * (lo + hi)
        density = _integrate_density_cm2(kt_nm, energies_by_k, mid, kbt_eV)
        if density < target_density_cm2:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def _line_shape(delta_eV: np.ndarray, broadening_eV: float, kind: str) -> np.ndarray:
    if broadening_eV <= 0:
        raise ValueError("broadening_eV must be positive")
    if kind == "gaussian":
        sigma = broadening_eV / (2.0 * np.sqrt(2.0 * np.log(2.0)))
        return np.exp(-0.5 * (delta_eV / sigma) ** 2) / (sigma * np.sqrt(2.0 * np.pi))
    if kind == "lorentzian":
        gamma = 0.5 * broadening_eV
        return (gamma / np.pi) / (delta_eV * delta_eV + gamma * gamma)
    raise ValueError("line_shape must be 'lorentzian' or 'gaussian'")


def _transition_terms(
    profile: MQWProfile,
    subbands: KpSubbands,
) -> list[dict[str, float]]:
    terms: list[dict[str, float]] = []
    dz = profile.dz_nm
    for k_index, kt in enumerate(subbands.kt_nm):
        electrons = subbands.electron[k_index]
        holes = subbands.valence[k_index]
        for e_index, estate in enumerate(electrons):
            e_prob = np.abs(estate.psi) ** 2
            for h_index, hstate in enumerate(holes):
                hh_overlap = float(np.sum(estate.psi * hstate.hh) * dz)
                lh_overlap = float(np.sum(estate.psi * hstate.lh) * dz)
                te_strength = hh_overlap * hh_overlap + (lh_overlap * lh_overlap) / 3.0
                tm_strength = 4.0 * (lh_overlap * lh_overlap) / 3.0
                denom = float(
                    np.sum(e_prob * (np.abs(hstate.hh) ** 2 + np.abs(hstate.lh) ** 2))
                    * dz
                )
                if denom <= 0:
                    continue
                eg_eff = float(
                    np.sum(
                        e_prob
                        * (
                            np.abs(hstate.hh) ** 2 * profile.Eg_hh_eV
                            + np.abs(hstate.lh) ** 2 * profile.Eg_lh_eV
                        )
                    )
                    * dz
                    / denom
                )
                terms.append(
                    {
                        "kt_nm": float(kt),
                        "electron_index": float(e_index),
                        "hole_index": float(h_index),
                        "electron_energy_eV": estate.energy_eV,
                        "hole_energy_eV": hstate.energy_eV,
                        "transition_energy_eV": eg_eff + estate.energy_eV + hstate.energy_eV,
                        "te_strength": max(0.0, te_strength),
                        "tm_strength": max(0.0, tm_strength),
                    }
                )
    return terms


def calculate_gain_spectrum(
    profile: MQWProfile,
    subbands: KpSubbands,
    carrier_density_cm3: float,
    temperature_K: float = 300.0,
    energy_min_eV: float | None = None,
    energy_max_eV: float | None = None,
    energy_points: int = 500,
    broadening_eV: float = 0.030,
    line_shape: str = "lorentzian",
    gain_scale_cm: float = 2400.0,
) -> tuple[GainSpectrum, list[dict[str, float]]]:
    if carrier_density_cm3 <= 0:
        raise ValueError("carrier_density_cm3 must be positive")
    if temperature_K <= 0:
        raise ValueError("temperature_K must be positive")
    if energy_points < 2:
        raise ValueError("energy_points must be at least 2")
    kbt = KB_EV_K * temperature_K
    target_sheet_cm2 = carrier_density_cm3 * profile.active_nm * 1.0e-7

    e_energies = [[state.energy_eV for state in row] for row in subbands.electron]
    h_energies = [[state.energy_eV for state in row] for row in subbands.valence]
    ef_e = _solve_fermi_level(subbands.kt_nm, e_energies, target_sheet_cm2, kbt)
    ef_h = _solve_fermi_level(subbands.kt_nm, h_energies, target_sheet_cm2, kbt)

    terms = _transition_terms(profile, subbands)
    transition_energies = np.array([term["transition_energy_eV"] for term in terms])
    if energy_min_eV is None:
        energy_min_eV = float(max(0.05, np.min(transition_energies) - 0.12))
    if energy_max_eV is None:
        energy_max_eV = float(np.max(transition_energies) + 0.12)
    energy = np.linspace(energy_min_eV, energy_max_eV, energy_points)
    gain_te = np.zeros_like(energy)
    gain_tm = np.zeros_like(energy)

    kt = subbands.kt_nm
    if len(kt) == 1:
        weights = np.array([kt[0] if kt[0] > 0 else 1.0], dtype=float)
    else:
        weights = np.gradient(kt) * kt

    for term in terms:
        k_index = int(np.argmin(np.abs(kt - term["kt_nm"])))
        fe = float(_fermi(term["electron_energy_eV"], ef_e, kbt))
        fh = float(_fermi(term["hole_energy_eV"], ef_h, kbt))
        inversion = fe + fh - 1.0
        if abs(inversion) < 1e-12:
            continue
        lines = _line_shape(energy - term["transition_energy_eV"], broadening_eV, line_shape)
        prefactor = gain_scale_cm * weights[k_index] / max(profile.active_nm, 1e-12)
        gain_te += prefactor * term["te_strength"] * inversion * lines
        gain_tm += prefactor * term["tm_strength"] * inversion * lines

    wavelength_nm = HC_EV_UM * 1000.0 / energy
    peak_te_index = int(np.argmax(gain_te))
    peak_tm_index = int(np.argmax(gain_tm))
    spectrum = GainSpectrum(
        energy_eV=energy,
        wavelength_nm=wavelength_nm,
        gain_te_cm=gain_te,
        gain_tm_cm=gain_tm,
        quasi_fermi_e_eV=ef_e,
        quasi_fermi_h_eV=ef_h,
        peak_te_gain_cm=float(gain_te[peak_te_index]),
        peak_te_wavelength_nm=float(wavelength_nm[peak_te_index]),
        peak_tm_gain_cm=float(gain_tm[peak_tm_index]),
        peak_tm_wavelength_nm=float(wavelength_nm[peak_tm_index]),
    )
    return spectrum, terms


def spectrum_to_rows(spectrum: GainSpectrum) -> list[dict[str, float]]:
    return [
        {
            "energy_eV": float(energy),
            "wavelength_nm": float(wavelength),
            "gain_TE_cm-1": float(gain_te),
            "gain_TM_cm-1": float(gain_tm),
        }
        for energy, wavelength, gain_te, gain_tm in zip(
            spectrum.energy_eV,
            spectrum.wavelength_nm,
            spectrum.gain_te_cm,
            spectrum.gain_tm_cm,
            strict=True,
        )
    ]


def gain_summary_dict(
    spectrum: GainSpectrum,
    terms: list[dict[str, float]],
    carrier_density_cm3: float,
    temperature_K: float,
) -> dict[str, Any]:
    return {
        "carrier_density_cm3": carrier_density_cm3,
        "temperature_K": temperature_K,
        "quasi_fermi_e_eV": spectrum.quasi_fermi_e_eV,
        "quasi_fermi_h_eV": spectrum.quasi_fermi_h_eV,
        "peak_te_gain_cm-1": spectrum.peak_te_gain_cm,
        "peak_te_wavelength_nm": spectrum.peak_te_wavelength_nm,
        "peak_tm_gain_cm-1": spectrum.peak_tm_gain_cm,
        "peak_tm_wavelength_nm": spectrum.peak_tm_wavelength_nm,
        "transition_count": len(terms),
        "strongest_transitions": sorted(
            terms,
            key=lambda item: item["te_strength"] + item["tm_strength"],
            reverse=True,
        )[:12],
    }
