#!/usr/bin/env python3
"""Sweep compact k.p MQW gain peaks versus one design or calibration knob."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np

from BasicMQWDesign import design_default
from gain import calculate_gain_spectrum, spectrum_to_rows
from kp_solver import solve_kp_subbands

SweepName = str


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


def _parse_values(args: argparse.Namespace) -> list[float]:
    if args.values:
        values = [float(item) for item in args.values.split(",")]
    else:
        if args.points < 2:
            raise ValueError("--points must be at least 2 when --values is omitted")
        if args.log:
            if args.start <= 0 or args.stop <= 0:
                raise ValueError("--log requires positive --start and --stop")
            values = np.geomspace(args.start, args.stop, args.points).tolist()
        else:
            values = np.linspace(args.start, args.stop, args.points).tolist()
    return values


def build_arg_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        description="Sweep compact k.p MQW gain peak wavelength and peak gain"
    )
    ap.add_argument(
        "--sweep",
        choices=[
            "carrier-density",
            "well-nm",
            "well-strain",
            "barrier-strain",
            "qc",
            "broadening-eV",
        ],
        default="carrier-density",
    )
    ap.add_argument("--start", type=float, default=1.0e18)
    ap.add_argument("--stop", type=float, default=3.0e18)
    ap.add_argument("--points", type=int, default=7)
    ap.add_argument(
        "--values",
        type=str,
        default=None,
        help="Comma-separated sweep values. Overrides --start/--stop/--points.",
    )
    ap.add_argument("--log", action="store_true", help="Use logarithmic sweep spacing")

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
    ap.add_argument("--kt-points", type=int, default=21)
    ap.add_argument("--electron-states", type=int, default=2)
    ap.add_argument("--hole-states", type=int, default=4)
    ap.add_argument("--energy-min-eV", type=float, default=None)
    ap.add_argument("--energy-max-eV", type=float, default=None)
    ap.add_argument("--energy-points", type=int, default=320)
    ap.add_argument("--broadening-eV", type=float, default=0.030)
    ap.add_argument("--line-shape", choices=["lorentzian", "gaussian"], default="lorentzian")
    ap.add_argument("--gain-scale-cm", type=float, default=2400.0)

    ap.add_argument("--out-json", type=Path, default=Path("out/gain_sweep.json"))
    ap.add_argument("--out-csv", type=Path, default=Path("out/gain_sweep.csv"))
    ap.add_argument("--plot", type=Path, default=Path("out/gain_sweep.png"))
    ap.add_argument(
        "--spectra-csv",
        type=Path,
        default=Path("out/gain_sweep_spectra.csv"),
        help="CSV output for wavelength-domain spectra at every sweep value.",
    )
    ap.add_argument(
        "--spectra-plot",
        type=Path,
        default=Path("out/gain_sweep_spectra.png"),
        help="Overlay plot of gain versus wavelength at every sweep value.",
    )
    return ap


def _design_kwargs(args: argparse.Namespace, sweep: SweepName, value: float) -> dict[str, Any]:
    kwargs = {
        "family": args.family,
        "wells": args.wells,
        "well_nm": args.well_nm,
        "barrier_nm": args.barrier_nm,
        "q_c": args.qc,
        "well_strain": args.well_strain,
        "barrier_strain": args.barrier_strain,
        "al_well": args.al_well,
        "al_barrier": args.al_barrier,
        "as_well": args.as_well,
        "as_barrier": args.as_barrier,
    }
    if sweep == "well-nm":
        kwargs["well_nm"] = value
    elif sweep == "well-strain":
        kwargs["well_strain"] = value
    elif sweep == "barrier-strain":
        kwargs["barrier_strain"] = value
    elif sweep == "qc":
        kwargs["q_c"] = value
    return kwargs


def _gain_kwargs(args: argparse.Namespace, sweep: SweepName, value: float) -> dict[str, Any]:
    kwargs = {
        "carrier_density_cm3": args.carrier_density_cm3,
        "temperature_K": args.temperature,
        "energy_min_eV": args.energy_min_eV,
        "energy_max_eV": args.energy_max_eV,
        "energy_points": args.energy_points,
        "broadening_eV": args.broadening_eV,
        "line_shape": args.line_shape,
        "gain_scale_cm": args.gain_scale_cm,
    }
    if sweep == "carrier-density":
        kwargs["carrier_density_cm3"] = value
    elif sweep == "broadening-eV":
        kwargs["broadening_eV"] = value
    return kwargs


def run_one(
    args: argparse.Namespace, sweep: SweepName, value: float
) -> tuple[dict[str, float], list[dict[str, float]]]:
    design = design_default(**_design_kwargs(args, sweep, value))
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
        **_gain_kwargs(args, sweep, value),
    )
    summary = {
        "sweep_value": value,
        "carrier_density_cm3": float(_gain_kwargs(args, sweep, value)["carrier_density_cm3"]),
        "well_nm": float(design["well_nm"]),
        "barrier_nm": float(design["barrier_nm"]),
        "qc": float(design["qc"]),
        "well_strain": float(design["transition"]["well_strain"]["eps_parallel"]),
        "barrier_strain": float(design["transition"]["barrier_strain"]["eps_parallel"]),
        "broadening_eV": float(_gain_kwargs(args, sweep, value)["broadening_eV"]),
        "peak_TE_gain_cm-1": spectrum.peak_te_gain_cm,
        "peak_TE_wavelength_nm": spectrum.peak_te_wavelength_nm,
        "peak_TM_gain_cm-1": spectrum.peak_tm_gain_cm,
        "peak_TM_wavelength_nm": spectrum.peak_tm_wavelength_nm,
        "quasi_fermi_e_eV": spectrum.quasi_fermi_e_eV,
        "quasi_fermi_h_eV": spectrum.quasi_fermi_h_eV,
        "transition_count": float(len(terms)),
    }
    spectra_rows = []
    for row in spectrum_to_rows(spectrum):
        spectra_rows.append(
            {
                "sweep_value": value,
                "carrier_density_cm3": summary["carrier_density_cm3"],
                "well_nm": summary["well_nm"],
                "well_strain": summary["well_strain"],
                "qc": summary["qc"],
                **row,
            }
        )
    return summary, spectra_rows


def write_csv(rows: list[dict[str, float]], path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return path


def write_spectra_csv(rows: list[dict[str, float]], path: Path) -> Path:
    return write_csv(rows, path)


def write_plot(rows: list[dict[str, float]], sweep: SweepName, path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    x = np.array([row["sweep_value"] for row in rows], dtype=float)
    te_wavelength = np.array([row["peak_TE_wavelength_nm"] for row in rows], dtype=float)
    tm_wavelength = np.array([row["peak_TM_wavelength_nm"] for row in rows], dtype=float)
    te_gain = np.array([row["peak_TE_gain_cm-1"] for row in rows], dtype=float)
    tm_gain = np.array([row["peak_TM_gain_cm-1"] for row in rows], dtype=float)
    order = np.argsort(x)

    fig, (ax_wave, ax_gain) = plt.subplots(
        2,
        1,
        figsize=(7.0, 6.4),
        sharex=True,
        constrained_layout=True,
    )
    ax_wave.plot(x[order], te_wavelength[order], marker="o", label="TE peak")
    ax_wave.plot(x[order], tm_wavelength[order], marker="s", label="TM peak")
    ax_wave.set_ylabel("Peak wavelength [nm]")
    ax_wave.grid(True, alpha=0.25)
    ax_wave.legend()

    ax_gain.plot(x[order], te_gain[order], marker="o", label="TE peak")
    ax_gain.plot(x[order], tm_gain[order], marker="s", label="TM peak")
    ax_gain.axhline(0.0, color="0.35", linewidth=0.8)
    ax_gain.set_xlabel(_axis_label(sweep))
    ax_gain.set_ylabel("Peak material gain [cm$^{-1}$]")
    ax_gain.grid(True, alpha=0.25)
    ax_gain.legend()

    if sweep == "carrier-density":
        ax_gain.ticklabel_format(axis="x", style="sci", scilimits=(0, 0))
    fig.savefig(path, dpi=180)
    plt.close(fig)
    return path


def write_spectra_plot(
    rows: list[dict[str, float]],
    spectra_rows: list[dict[str, float]],
    sweep: SweepName,
    path: Path,
) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    values = [row["sweep_value"] for row in rows]
    cmap = plt.get_cmap("viridis", max(len(values), 2))

    fig, (ax_te, ax_tm) = plt.subplots(
        2,
        1,
        figsize=(7.4, 6.8),
        sharex=True,
        constrained_layout=True,
    )
    for index, value in enumerate(values):
        subset = [row for row in spectra_rows if row["sweep_value"] == value]
        wavelength = np.array([row["wavelength_nm"] for row in subset], dtype=float)
        gain_te = np.array([row["gain_TE_cm-1"] for row in subset], dtype=float)
        gain_tm = np.array([row["gain_TM_cm-1"] for row in subset], dtype=float)
        order = np.argsort(wavelength)
        label = _sweep_label(sweep, value)
        color = cmap(index)
        ax_te.plot(wavelength[order], gain_te[order], color=color, label=label)
        ax_tm.plot(wavelength[order], gain_tm[order], color=color, label=label)

    ax_te.axhline(0.0, color="0.35", linewidth=0.8)
    ax_tm.axhline(0.0, color="0.35", linewidth=0.8)
    ax_te.set_ylabel("TE gain [cm$^{-1}$]")
    ax_tm.set_ylabel("TM gain [cm$^{-1}$]")
    ax_tm.set_xlabel("Wavelength [nm]")
    ax_te.grid(True, alpha=0.25)
    ax_tm.grid(True, alpha=0.25)
    ax_te.legend(title=_axis_label(sweep), fontsize=8)
    ax_tm.legend(title=_axis_label(sweep), fontsize=8)
    fig.savefig(path, dpi=180)
    plt.close(fig)
    return path


def _axis_label(sweep: SweepName) -> str:
    labels = {
        "carrier-density": "Carrier density [cm$^{-3}$]",
        "well-nm": "Well width [nm]",
        "well-strain": "Well strain",
        "barrier-strain": "Barrier strain",
        "qc": "Conduction band offset fraction",
        "broadening-eV": "Broadening FWHM [eV]",
    }
    return labels[sweep]


def _sweep_label(sweep: SweepName, value: float) -> str:
    if sweep == "carrier-density":
        return f"{value:.2g}"
    if sweep in {"well-strain", "barrier-strain"}:
        return f"{value * 100:+.2f}%"
    if sweep == "broadening-eV":
        return f"{value * 1000:.0f} meV"
    return f"{value:.4g}"


def format_summary(
    rows: list[dict[str, float]],
    sweep: SweepName,
    json_path: Path,
    csv_path: Path,
    plot_path: Path,
    spectra_csv_path: Path,
    spectra_plot_path: Path,
) -> str:
    best_te = max(rows, key=lambda row: row["peak_TE_gain_cm-1"])
    best_tm = max(rows, key=lambda row: row["peak_TM_gain_cm-1"])
    return "\n".join(
        (
            "=== MQW compact k.p gain peak sweep ===",
            f"sweep              : {sweep}",
            f"points             : {len(rows)}",
            f"TE max peak        : {best_te['peak_TE_gain_cm-1']:.3g} cm^-1 at "
            f"{best_te['peak_TE_wavelength_nm']:.1f} nm "
            f"(sweep={best_te['sweep_value']:.6g})",
            f"TM max peak        : {best_tm['peak_TM_gain_cm-1']:.3g} cm^-1 at "
            f"{best_tm['peak_TM_wavelength_nm']:.1f} nm "
            f"(sweep={best_tm['sweep_value']:.6g})",
            f"wrote JSON         : {json_path}",
            f"wrote CSV          : {csv_path}",
            f"wrote plot         : {plot_path}",
            f"wrote spectra CSV  : {spectra_csv_path}",
            f"wrote spectra plot : {spectra_plot_path}",
        )
    )


def main(argv: list[str] | None = None) -> None:
    args = build_arg_parser().parse_args(argv)
    values = _parse_values(args)
    rows = []
    spectra_rows = []
    for value in values:
        summary, spectrum_rows = run_one(args, args.sweep, value)
        rows.append(summary)
        spectra_rows.extend(spectrum_rows)
    result = {
        "model_note": (
            "Peak sweep from compact screening model. Absolute gain requires "
            "calibration of offsets, broadening, and gain_scale_cm."
        ),
        "sweep": args.sweep,
        "values": values,
        "rows": rows,
        "spectra_csv": str(args.spectra_csv),
    }

    args.out_json.parent.mkdir(parents=True, exist_ok=True)
    args.out_json.write_text(
        json.dumps(_json_safe(result), ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    write_csv(rows, args.out_csv)
    write_spectra_csv(spectra_rows, args.spectra_csv)
    write_plot(rows, args.sweep, args.plot)
    write_spectra_plot(rows, spectra_rows, args.sweep, args.spectra_plot)
    print(
        format_summary(
            rows,
            args.sweep,
            args.out_json,
            args.out_csv,
            args.plot,
            args.spectra_csv,
            args.spectra_plot,
        )
    )


if __name__ == "__main__":
    main()
