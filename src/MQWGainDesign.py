#!/usr/bin/env python3
"""Run compact k.p MQW material-gain estimates from the command line."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np

from BasicMQWDesign import design_default
from gain import calculate_gain_spectrum, gain_summary_dict, spectrum_to_rows
from kp_solver import solve_kp_subbands, subband_summary


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
        description="Compact 4x4 k.p MQW material-gain screening helper"
    )
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

    ap.add_argument("--carrier-density-cm3", type=float, default=2.0e18)
    ap.add_argument("--temperature", type=float, default=300.0)
    ap.add_argument("--dz-nm", type=float, default=0.10)
    ap.add_argument("--kt-max-nm", type=float, default=0.35)
    ap.add_argument("--kt-points", type=int, default=31)
    ap.add_argument("--electron-states", type=int, default=2)
    ap.add_argument("--hole-states", type=int, default=4)
    ap.add_argument("--energy-min-eV", type=float, default=None)
    ap.add_argument("--energy-max-eV", type=float, default=None)
    ap.add_argument("--energy-points", type=int, default=500)
    ap.add_argument("--broadening-eV", type=float, default=0.030)
    ap.add_argument("--line-shape", choices=["lorentzian", "gaussian"], default="lorentzian")
    ap.add_argument(
        "--gain-scale-cm",
        type=float,
        default=2400.0,
        help="Empirical oscillator-strength scale for cm^-1 output calibration",
    )

    ap.add_argument("--out-json", type=Path, default=Path("out/gain_result.json"))
    ap.add_argument("--out-csv", type=Path, default=Path("out/gain_spectrum.csv"))
    ap.add_argument("--plot", type=Path, default=Path("out/gain_spectrum.png"))
    return ap


def write_csv(rows: list[dict[str, float]], path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return path


def write_plot(rows: list[dict[str, float]], path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    wavelength = np.array([row["wavelength_nm"] for row in rows])
    gain_te = np.array([row["gain_TE_cm-1"] for row in rows])
    gain_tm = np.array([row["gain_TM_cm-1"] for row in rows])
    order = np.argsort(wavelength)

    fig, ax = plt.subplots(figsize=(7.0, 4.2), constrained_layout=True)
    ax.plot(wavelength[order], gain_te[order], label="TE")
    ax.plot(wavelength[order], gain_tm[order], label="TM")
    ax.axhline(0.0, color="0.35", linewidth=0.8)
    ax.set_xlabel("Wavelength [nm]")
    ax.set_ylabel("Material gain [cm$^{-1}$]")
    ax.legend()
    ax.grid(True, alpha=0.25)
    fig.savefig(path, dpi=180)
    plt.close(fig)
    return path


def format_summary(result: dict[str, Any], json_path: Path, csv_path: Path, plot_path: Path) -> str:
    gain = result["gain"]
    kp = result["kp"]
    return "\n".join(
        (
            "=== MQW compact k.p gain estimate ===",
            f"family             : {result['design']['family']}",
            f"wells              : {result['design']['wells']}",
            f"well/barrier       : {result['design']['well_nm']} nm / {result['design']['barrier_nm']} nm",
            f"carrier density    : {gain['carrier_density_cm3']:.4g} cm^-3",
            f"temperature        : {gain['temperature_K']:.1f} K",
            f"grid               : dz={kp['dz_nm']:.4g} nm, kt points={kp['kt_points']}",
            f"TE peak gain       : {gain['peak_te_gain_cm-1']:.3g} cm^-1 at {gain['peak_te_wavelength_nm']:.1f} nm",
            f"TM peak gain       : {gain['peak_tm_gain_cm-1']:.3g} cm^-1 at {gain['peak_tm_wavelength_nm']:.1f} nm",
            f"wrote JSON         : {json_path}",
            f"wrote CSV          : {csv_path}",
            f"wrote plot         : {plot_path}",
        )
    )


def main(argv: list[str] | None = None) -> None:
    args = build_arg_parser().parse_args(argv)
    design = design_default(
        family=args.family,
        wells=args.wells,
        well_nm=args.well_nm,
        barrier_nm=args.barrier_nm,
        q_c=args.qc,
        well_strain=args.well_strain,
        barrier_strain=args.barrier_strain,
        al_well=args.al_well,
        al_barrier=args.al_barrier,
        as_well=args.as_well,
        as_barrier=args.as_barrier,
    )
    profile, subbands = solve_kp_subbands(
        design,
        dz_nm=args.dz_nm,
        kt_max_nm=args.kt_max_nm,
        kt_points=args.kt_points,
        electron_states=args.electron_states,
        hole_states=args.hole_states,
    )
    spectrum, terms = calculate_gain_spectrum(
        profile,
        subbands,
        carrier_density_cm3=args.carrier_density_cm3,
        temperature_K=args.temperature,
        energy_min_eV=args.energy_min_eV,
        energy_max_eV=args.energy_max_eV,
        energy_points=args.energy_points,
        broadening_eV=args.broadening_eV,
        line_shape=args.line_shape,
        gain_scale_cm=args.gain_scale_cm,
    )
    rows = spectrum_to_rows(spectrum)
    result = {
        "model_note": (
            "Compact screening model: scalar conduction band plus axial 4x4 valence "
            "k.p represented by one 2x2 HH/LH Kramers block. Absolute gain requires "
            "calibration of offsets, broadening, and gain_scale_cm."
        ),
        "design": design,
        "kp": subband_summary(profile, subbands),
        "gain": gain_summary_dict(
            spectrum,
            terms,
            carrier_density_cm3=args.carrier_density_cm3,
            temperature_K=args.temperature,
        ),
    }

    args.out_json.parent.mkdir(parents=True, exist_ok=True)
    args.out_json.write_text(
        json.dumps(_json_safe(result), ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    write_csv(rows, args.out_csv)
    write_plot(rows, args.plot)
    print(format_summary(result, args.out_json, args.out_csv, args.plot))


if __name__ == "__main__":
    main()
