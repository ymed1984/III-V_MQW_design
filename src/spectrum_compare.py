"""Shared comparison logic for predicted and reference gain spectra."""

from __future__ import annotations

from typing import Any, Literal

import numpy as np

from metrics import Polarization, peak_metrics, spectrum_rmse

PolarizationChoice = Literal["TE", "TM", "both"]


def polarizations_from_choice(choice: PolarizationChoice) -> list[Polarization]:
    if choice == "both":
        return ["TE", "TM"]
    return [choice]


def compare_spectra(
    predicted_rows: list[dict[str, float]],
    reference_rows: list[dict[str, float]],
    polarizations: list[Polarization],
) -> dict[str, Any]:
    comparisons: dict[str, Any] = {}
    for pol in polarizations:
        predicted = peak_metrics(predicted_rows, pol)
        reference = peak_metrics(reference_rows, pol)
        comparisons[pol] = {
            "predicted": predicted.as_dict(),
            "reference": reference.as_dict(),
            "delta": {
                "peak_wavelength_nm": predicted.peak_wavelength_nm
                - reference.peak_wavelength_nm,
                "peak_energy_eV": predicted.peak_energy_eV - reference.peak_energy_eV,
                "peak_gain_cm": predicted.peak_gain_cm - reference.peak_gain_cm,
                "fwhm_meV": predicted.fwhm_meV - reference.fwhm_meV,
            },
            "rmse_gain_cm": spectrum_rmse(predicted_rows, reference_rows, pol),
        }
    return comparisons


def flatten_comparison_metrics(
    predicted_rows: list[dict[str, float]],
    reference_rows: list[dict[str, float]],
    polarizations: list[Polarization],
) -> dict[str, float]:
    comparisons = compare_spectra(predicted_rows, reference_rows, polarizations)
    output: dict[str, float] = {}
    rmses = []
    for pol, comparison in comparisons.items():
        delta = comparison["delta"]
        rmse = float(comparison["rmse_gain_cm"])
        rmses.append(rmse)
        output[f"{pol}_rmse_gain_cm"] = rmse
        output[f"{pol}_peak_wavelength_delta_nm"] = delta["peak_wavelength_nm"]
        output[f"{pol}_peak_gain_delta_cm"] = delta["peak_gain_cm"]
        output[f"{pol}_fwhm_delta_meV"] = delta["fwhm_meV"]
    output["score_rmse_gain_cm"] = float(np.mean(rmses))
    return output
