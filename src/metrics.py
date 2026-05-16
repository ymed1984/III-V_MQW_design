#!/usr/bin/env python3
"""Reusable spectrum metrics for MQW gain calibration and comparison."""

from __future__ import annotations

import math
from dataclasses import dataclass, asdict
from typing import Literal

import numpy as np

Polarization = Literal["TE", "TM"]


@dataclass(frozen=True)
class SpectrumMetrics:
    peak_wavelength_nm: float
    peak_energy_eV: float
    peak_gain_cm: float
    fwhm_meV: float

    def as_dict(self) -> dict[str, float]:
        return asdict(self)


def _column_name(polarization: Polarization) -> str:
    if polarization == "TE":
        return "gain_TE_cm-1"
    if polarization == "TM":
        return "gain_TM_cm-1"
    raise ValueError("polarization must be 'TE' or 'TM'")


def _arrays(
    rows: list[dict[str, float]],
    polarization: Polarization,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    if not rows:
        raise ValueError("spectrum rows must not be empty")
    gain_key = _column_name(polarization)
    energy = np.array([float(row["energy_eV"]) for row in rows], dtype=float)
    wavelength = np.array([float(row["wavelength_nm"]) for row in rows], dtype=float)
    gain = np.array([float(row[gain_key]) for row in rows], dtype=float)
    if not (np.all(np.isfinite(energy)) and np.all(np.isfinite(wavelength)) and np.all(np.isfinite(gain))):
        raise ValueError("spectrum rows contain non-finite values")
    order = np.argsort(energy)
    return energy[order], wavelength[order], gain[order]


def peak_metrics(
    rows: list[dict[str, float]],
    polarization: Polarization = "TE",
) -> SpectrumMetrics:
    energy, wavelength, gain = _arrays(rows, polarization)
    peak_index = int(np.argmax(gain))
    return SpectrumMetrics(
        peak_wavelength_nm=float(wavelength[peak_index]),
        peak_energy_eV=float(energy[peak_index]),
        peak_gain_cm=float(gain[peak_index]),
        fwhm_meV=fwhm_from_spectrum(rows, polarization),
    )


def fwhm_from_spectrum(
    rows: list[dict[str, float]],
    polarization: Polarization = "TE",
) -> float:
    energy, _, gain = _arrays(rows, polarization)
    peak = float(np.max(gain))
    if peak <= 0.0:
        return float("nan")
    half = 0.5 * peak
    above = np.where(gain >= half)[0]
    if len(above) == 0:
        return float("nan")
    left_index = int(above[0])
    right_index = int(above[-1])

    left = _half_crossing(energy, gain, left_index, half, side="left")
    right = _half_crossing(energy, gain, right_index, half, side="right")
    width = right - left
    return float(width * 1000.0) if width >= 0.0 else float("nan")


def _half_crossing(
    energy: np.ndarray,
    gain: np.ndarray,
    index: int,
    half: float,
    side: Literal["left", "right"],
) -> float:
    if side == "left":
        if index == 0:
            return float(energy[index])
        lo = index - 1
        hi = index
    else:
        if index == len(energy) - 1:
            return float(energy[index])
        lo = index
        hi = index + 1
    g0 = float(gain[lo])
    g1 = float(gain[hi])
    e0 = float(energy[lo])
    e1 = float(energy[hi])
    if math.isclose(g0, g1):
        return e0
    frac = (half - g0) / (g1 - g0)
    return e0 + frac * (e1 - e0)


def gain_at_wavelength(
    rows: list[dict[str, float]],
    wavelength_nm: float,
    polarization: Polarization = "TE",
) -> float:
    if not math.isfinite(wavelength_nm):
        raise ValueError("wavelength_nm must be finite")
    _, wavelength, gain = _arrays(rows, polarization)
    order = np.argsort(wavelength)
    return float(np.interp(wavelength_nm, wavelength[order], gain[order]))


def spectrum_rmse(
    predicted_rows: list[dict[str, float]],
    reference_rows: list[dict[str, float]],
    polarization: Polarization = "TE",
) -> float:
    _, reference_wavelength, reference_gain = _arrays(reference_rows, polarization)
    predicted_at_reference = np.array(
        [
            gain_at_wavelength(predicted_rows, float(wavelength), polarization)
            for wavelength in reference_wavelength
        ],
        dtype=float,
    )
    return float(np.sqrt(np.mean((predicted_at_reference - reference_gain) ** 2)))
