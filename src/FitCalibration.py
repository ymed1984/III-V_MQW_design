#!/usr/bin/env python3
"""Fit a minimal calibration JSON from target MQW gain metrics."""

from __future__ import annotations

import argparse
import json
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

import numpy as np
from scipy.optimize import minimize_scalar

from BasicMQWDesign import design_default
from calibration import load_calibration, resolve_calibration
from gain import calculate_gain_spectrum, spectrum_to_rows
from kp_solver import solve_kp_subbands
from metrics import SpectrumMetrics, peak_metrics


@dataclass(frozen=True)
class FitState:
    Eg_offset_well_eV: float
    broadening_eV: float
    gain_scale_cm: float
    metrics: SpectrumMetrics

    def as_dict(self) -> dict[str, Any]:
        result = asdict(self)
        result["metrics"] = self.metrics.as_dict()
        return result


def _json_safe(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(key): _json_safe(item) for key, item in value.items()}
    if isinstance(value, list | tuple):
        return [_json_safe(item) for item in value]
    if isinstance(value, np.ndarray):
        return [_json_safe(item) for item in value.tolist()]
    if isinstance(value, np.floating):
        return float(value)
    if isinstance(value, np.integer):
        return int(value)
    return value


def build_arg_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        description="Fit a minimal MQW gain calibration JSON from target metrics"
    )
    ap.add_argument("--calibration-in", type=Path, default=None)
    ap.add_argument("--out", type=Path, default=Path("calibrations/fitted/fit.json"))
    ap.add_argument("--name", type=str, default="mqw_gain_fit")

    ap.add_argument("--family", choices=["algainas", "ingaasp"], default="ingaasp")
    ap.add_argument("--wells", type=int, default=5)
    ap.add_argument("--well-nm", type=float, default=7.0)
    ap.add_argument("--barrier-nm", type=float, default=10.0)
    ap.add_argument("--qc", type=float, default=None)
    ap.add_argument("--well-strain", type=float, default=None)
    ap.add_argument("--barrier-strain", type=float, default=None)
    ap.add_argument("--al-well", type=float, default=None)
    ap.add_argument("--al-barrier", type=float, default=None)
    ap.add_argument("--as-well", type=float, default=None)
    ap.add_argument("--as-barrier", type=float, default=None)
    ap.add_argument("--eg-offset-well-eV", type=float, default=None)
    ap.add_argument("--eg-offset-barrier-eV", type=float, default=None)

    ap.add_argument("--carrier-density-cm3", type=float, default=2.0e18)
    ap.add_argument("--temperature", type=float, default=300.0)
    ap.add_argument("--dz-nm", type=float, default=0.2)
    ap.add_argument("--kt-max-nm", type=float, default=0.35)
    ap.add_argument("--kt-points", type=int, default=9)
    ap.add_argument("--electron-states", type=int, default=1)
    ap.add_argument("--hole-states", type=int, default=2)
    ap.add_argument("--energy-min-eV", type=float, default=None)
    ap.add_argument("--energy-max-eV", type=float, default=None)
    ap.add_argument("--energy-points", type=int, default=160)
    ap.add_argument("--broadening-eV", type=float, default=None)
    ap.add_argument("--line-shape", choices=["lorentzian", "gaussian"], default=None)
    ap.add_argument("--gain-scale-cm", type=float, default=None)

    ap.add_argument("--target-peak-wavelength-nm", type=float, required=True)
    ap.add_argument("--target-te-peak-gain-cm", type=float, required=True)
    ap.add_argument("--target-fwhm-meV", type=float, default=None)
    ap.add_argument("--eg-offset-min-eV", type=float, default=-0.10)
    ap.add_argument("--eg-offset-max-eV", type=float, default=0.10)
    ap.add_argument("--broadening-min-eV", type=float, default=0.01)
    ap.add_argument("--broadening-max-eV", type=float, default=0.08)
    ap.add_argument("--fit-maxiter", type=int, default=16)
    return ap


def _resolved_from_args(args: argparse.Namespace):
    calibration = load_calibration(args.calibration_in)
    return resolve_calibration(args, calibration)


def _run_gain(
    args: argparse.Namespace,
    eg_offset_well_eV: float,
    broadening_eV: float,
    gain_scale_cm: float,
) -> tuple[SpectrumMetrics, list[dict[str, float]], dict[str, Any]]:
    resolved = _resolved_from_args(args)
    design = design_default(
        family=args.family,
        wells=args.wells,
        well_nm=args.well_nm,
        barrier_nm=args.barrier_nm,
        q_c=resolved.q_c,
        well_strain=args.well_strain,
        barrier_strain=args.barrier_strain,
        al_well=args.al_well,
        al_barrier=args.al_barrier,
        as_well=args.as_well,
        as_barrier=args.as_barrier,
        eg_offset_well_eV=eg_offset_well_eV,
        eg_offset_barrier_eV=resolved.Eg_offset_barrier_eV,
    )
    profile, subbands = solve_kp_subbands(
        design,
        dz_nm=args.dz_nm,
        kt_max_nm=args.kt_max_nm,
        kt_points=args.kt_points,
        electron_states=args.electron_states,
        hole_states=args.hole_states,
    )
    spectrum, _ = calculate_gain_spectrum(
        profile,
        subbands,
        carrier_density_cm3=args.carrier_density_cm3,
        temperature_K=args.temperature,
        energy_min_eV=args.energy_min_eV,
        energy_max_eV=args.energy_max_eV,
        energy_points=args.energy_points,
        broadening_eV=broadening_eV,
        line_shape=resolved.line_shape,
        gain_scale_cm=gain_scale_cm,
    )
    rows = spectrum_to_rows(spectrum)
    return peak_metrics(rows, "TE"), rows, design


def _fit_eg_offset(args: argparse.Namespace, initial_broadening: float, initial_gain_scale: float) -> float:
    if args.eg_offset_min_eV >= args.eg_offset_max_eV:
        raise ValueError("--eg-offset-min-eV must be smaller than --eg-offset-max-eV")

    def objective(offset: float) -> float:
        metrics, _, _ = _run_gain(args, offset, initial_broadening, initial_gain_scale)
        return (metrics.peak_wavelength_nm - args.target_peak_wavelength_nm) ** 2

    result = minimize_scalar(
        objective,
        bounds=(args.eg_offset_min_eV, args.eg_offset_max_eV),
        method="bounded",
        options={"xatol": 1e-4, "maxiter": args.fit_maxiter},
    )
    if not result.success and not np.isfinite(result.x):
        raise RuntimeError(f"Eg offset fit failed: {result.message}")
    return float(result.x)


def _fit_broadening(args: argparse.Namespace, eg_offset: float, initial_gain_scale: float) -> float:
    resolved = _resolved_from_args(args)
    if args.target_fwhm_meV is None:
        return resolved.broadening_eV
    if args.broadening_min_eV <= 0 or args.broadening_min_eV >= args.broadening_max_eV:
        raise ValueError("invalid broadening bounds")

    def objective(broadening: float) -> float:
        metrics, _, _ = _run_gain(args, eg_offset, broadening, initial_gain_scale)
        if not np.isfinite(metrics.fwhm_meV):
            return float("inf")
        return (metrics.fwhm_meV - args.target_fwhm_meV) ** 2

    result = minimize_scalar(
        objective,
        bounds=(args.broadening_min_eV, args.broadening_max_eV),
        method="bounded",
        options={"xatol": 5e-4, "maxiter": args.fit_maxiter},
    )
    if not result.success and not np.isfinite(result.x):
        raise RuntimeError(f"broadening fit failed: {result.message}")
    return float(result.x)


def fit_calibration(args: argparse.Namespace) -> tuple[dict[str, Any], FitState]:
    resolved = _resolved_from_args(args)
    eg_offset = _fit_eg_offset(args, resolved.broadening_eV, resolved.gain_scale_cm)
    broadening = _fit_broadening(args, eg_offset, resolved.gain_scale_cm)
    metrics_before_scale, _, design = _run_gain(
        args,
        eg_offset,
        broadening,
        resolved.gain_scale_cm,
    )
    if metrics_before_scale.peak_gain_cm <= 0:
        raise RuntimeError("cannot fit gain scale because current TE peak gain is not positive")
    gain_scale = (
        resolved.gain_scale_cm
        * args.target_te_peak_gain_cm
        / metrics_before_scale.peak_gain_cm
    )
    final_metrics, _, _ = _run_gain(args, eg_offset, broadening, gain_scale)
    state = FitState(
        Eg_offset_well_eV=eg_offset,
        broadening_eV=broadening,
        gain_scale_cm=gain_scale,
        metrics=final_metrics,
    )
    output = {
        "name": args.name,
        "description": "Generated by FitCalibration.py; validate against reference data before use.",
        "reference": {
            "target_peak_wavelength_nm": args.target_peak_wavelength_nm,
            "target_te_peak_gain_cm": args.target_te_peak_gain_cm,
            "target_fwhm_meV": args.target_fwhm_meV,
        },
        "design_filter": {
            "family": args.family,
            "wells": args.wells,
            "well_nm": args.well_nm,
            "barrier_nm": args.barrier_nm,
        },
        "band": {
            "qc": float(design["qc"]),
            "Eg_offset_well_eV": eg_offset,
            "Eg_offset_barrier_eV": resolved.Eg_offset_barrier_eV,
        },
        "gain": {
            "broadening_eV": broadening,
            "line_shape": resolved.line_shape,
            "gain_scale_cm": gain_scale,
        },
        "fit_result": state.as_dict(),
    }
    output["fit_result"]["residuals"] = {
        "peak_wavelength_nm": final_metrics.peak_wavelength_nm
        - args.target_peak_wavelength_nm,
        "peak_gain_cm": final_metrics.peak_gain_cm - args.target_te_peak_gain_cm,
        "fwhm_meV": (
            None
            if args.target_fwhm_meV is None
            else final_metrics.fwhm_meV - args.target_fwhm_meV
        ),
    }
    return output, state


def format_summary(path: Path, state: FitState) -> str:
    metrics = state.metrics
    return "\n".join(
        (
            "=== MQW gain calibration fit ===",
            f"Eg offset well     : {state.Eg_offset_well_eV:+.5f} eV",
            f"broadening         : {state.broadening_eV * 1000:.2f} meV",
            f"gain scale         : {state.gain_scale_cm:.6g} cm^-1 scale",
            f"TE peak wavelength : {metrics.peak_wavelength_nm:.2f} nm",
            f"TE peak gain       : {metrics.peak_gain_cm:.3g} cm^-1",
            f"TE FWHM            : {metrics.fwhm_meV:.2f} meV",
            f"wrote calibration  : {path}",
        )
    )


def main(argv: list[str] | None = None) -> None:
    args = build_arg_parser().parse_args(argv)
    output, state = fit_calibration(args)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(
        json.dumps(_json_safe(output), ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    print(format_summary(args.out, state))


if __name__ == "__main__":
    main()
