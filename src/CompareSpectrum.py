#!/usr/bin/env python3
"""Compare predicted MQW gain spectra against reference spectra."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np

from json_utils import json_safe
from metrics import Polarization
from spectrum_compare import compare_spectra, polarizations_from_choice
from spectrum_io import filter_wavelength_range, read_spectrum_csv


def build_arg_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        description="Compare predicted and reference TE/TM gain spectra"
    )
    ap.add_argument("--predicted", type=Path, required=True)
    ap.add_argument("--reference", type=Path, required=True)
    ap.add_argument("--out-json", type=Path, default=Path("out/spectrum_comparison.json"))
    ap.add_argument("--out-plot", type=Path, default=Path("out/spectrum_comparison.png"))
    ap.add_argument(
        "--polarization",
        choices=["TE", "TM", "both"],
        default="both",
        help="Polarization to compare.",
    )
    ap.add_argument("--wavelength-min-nm", type=float, default=None)
    ap.add_argument("--wavelength-max-nm", type=float, default=None)
    return ap


def write_comparison_plot(
    predicted_rows: list[dict[str, float]],
    reference_rows: list[dict[str, float]],
    polarizations: list[Polarization],
    path: Path,
) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(
        len(polarizations),
        1,
        figsize=(7.4, 3.8 * len(polarizations)),
        sharex=True,
        constrained_layout=True,
    )
    if len(polarizations) == 1:
        axes = [axes]
    for ax, pol in zip(axes, polarizations, strict=True):
        key = f"gain_{pol}_cm-1"
        pred_wl = np.array([row["wavelength_nm"] for row in predicted_rows], dtype=float)
        pred_gain = np.array([row[key] for row in predicted_rows], dtype=float)
        ref_wl = np.array([row["wavelength_nm"] for row in reference_rows], dtype=float)
        ref_gain = np.array([row[key] for row in reference_rows], dtype=float)
        pred_order = np.argsort(pred_wl)
        ref_order = np.argsort(ref_wl)
        ax.plot(pred_wl[pred_order], pred_gain[pred_order], label=f"Predicted {pol}")
        ax.plot(
            ref_wl[ref_order],
            ref_gain[ref_order],
            linestyle="--",
            label=f"Reference {pol}",
        )
        ax.axhline(0.0, color="0.35", linewidth=0.8)
        ax.set_ylabel("Gain [cm$^{-1}$]")
        ax.grid(True, alpha=0.25)
        ax.legend()
    axes[-1].set_xlabel("Wavelength [nm]")
    fig.savefig(path, dpi=180)
    plt.close(fig)
    return path


def format_summary(result: dict[str, Any], json_path: Path, plot_path: Path) -> str:
    lines = ["=== MQW gain spectrum comparison ==="]
    for pol, comparison in result["comparisons"].items():
        delta = comparison["delta"]
        lines.append(
            f"{pol} peak delta     : {delta['peak_wavelength_nm']:+.2f} nm, "
            f"{delta['peak_gain_cm']:+.3g} cm^-1"
        )
        lines.append(f"{pol} gain RMSE      : {comparison['rmse_gain_cm']:.3g} cm^-1")
    lines.extend((f"wrote JSON         : {json_path}", f"wrote plot         : {plot_path}"))
    return "\n".join(lines)


def main(argv: list[str] | None = None) -> None:
    args = build_arg_parser().parse_args(argv)
    polarizations = polarizations_from_choice(args.polarization)
    predicted_rows = filter_wavelength_range(
        read_spectrum_csv(args.predicted),
        args.wavelength_min_nm,
        args.wavelength_max_nm,
    )
    reference_rows = filter_wavelength_range(
        read_spectrum_csv(args.reference),
        args.wavelength_min_nm,
        args.wavelength_max_nm,
    )
    result = {
        "predicted": str(args.predicted),
        "reference": str(args.reference),
        "wavelength_filter": {
            "min_nm": args.wavelength_min_nm,
            "max_nm": args.wavelength_max_nm,
        },
        "comparisons": compare_spectra(predicted_rows, reference_rows, polarizations),
    }

    args.out_json.parent.mkdir(parents=True, exist_ok=True)
    args.out_json.write_text(
        json.dumps(json_safe(result), ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    write_comparison_plot(predicted_rows, reference_rows, polarizations, args.out_plot)
    print(format_summary(result, args.out_json, args.out_plot))


if __name__ == "__main__":
    main()
