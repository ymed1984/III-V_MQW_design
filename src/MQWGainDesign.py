#!/usr/bin/env python3
"""Run compact k.p MQW material-gain estimates from the command line."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

from BasicMQWDesign import design_default, load_design_input, load_design_json
from calibration import calibration_summary, load_calibration, resolve_calibration
from gain import calculate_gain_spectrum, gain_summary_dict, spectrum_to_rows
from json_utils import json_safe
from kp_solver import solve_kp_subbands, subband_summary
from spectrum_io import write_rows_csv
from visualization import (
    plot_band_diagram,
    plot_gain_spectrum,
    plot_subband_dispersion,
    plot_wavefunctions,
)


def build_arg_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        description="Compact 4x4 k.p MQW material-gain screening helper"
    )
    ap.add_argument("--design-json", type=Path, default=None,
                    help="Load a pre-computed DesignDict JSON instead of design_default()")
    ap.add_argument("--design-input", type=Path, default=None,
                    help="JSON file specifying MQW compositions and geometry")
    ap.add_argument("--calibration", type=Path, default=None)
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
    ap.add_argument("--dz-nm", type=float, default=0.10)
    ap.add_argument("--kt-max-nm", type=float, default=0.35)
    ap.add_argument("--kt-points", type=int, default=31)
    ap.add_argument("--electron-states", type=int, default=2)
    ap.add_argument("--hole-states", type=int, default=4)
    ap.add_argument("--energy-min-eV", type=float, default=None)
    ap.add_argument("--energy-max-eV", type=float, default=None)
    ap.add_argument("--energy-points", type=int, default=500)
    ap.add_argument("--broadening-eV", type=float, default=None)
    ap.add_argument("--line-shape", choices=["lorentzian", "gaussian"], default=None)
    ap.add_argument(
        "--gain-scale-cm",
        type=float,
        default=None,
        help="Empirical oscillator-strength scale for cm^-1 output calibration",
    )

    ap.add_argument("--out-json", type=Path, default=Path("out/gain_result.json"))
    ap.add_argument("--out-csv", type=Path, default=Path("out/gain_spectrum.csv"))
    ap.add_argument("--plot", type=Path, default=Path("out/gain_spectrum.png"))
    ap.add_argument("--band-plot", type=Path, default=None)
    ap.add_argument("--wavefunction-plot", type=Path, default=None)
    ap.add_argument("--dispersion-plot", type=Path, default=None)
    return ap


def write_csv(rows: list[dict[str, float]], path: Path) -> Path:
    return write_rows_csv(rows, path)


def write_plot(rows: list[dict[str, float]], path: Path) -> Path:
    return plot_gain_spectrum(rows, path)


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
    calibration = load_calibration(args.calibration)
    resolved_calibration = resolve_calibration(args, calibration)
    if args.design_json is not None:
        design = load_design_json(args.design_json)
    elif args.design_input is not None:
        input_kwargs = load_design_input(args.design_input)
        # Calibration overrides can still apply
        if resolved_calibration.q_c is not None and "q_c" not in input_kwargs:
            input_kwargs["q_c"] = resolved_calibration.q_c
        input_kwargs.setdefault("eg_offset_well_eV", resolved_calibration.Eg_offset_well_eV)
        input_kwargs.setdefault("eg_offset_barrier_eV", resolved_calibration.Eg_offset_barrier_eV)
        design = design_default(**input_kwargs)
    else:
        design = design_default(
            family=args.family,
            wells=args.wells,
            well_nm=args.well_nm,
            barrier_nm=args.barrier_nm,
            q_c=resolved_calibration.q_c,
            well_strain=args.well_strain,
            barrier_strain=args.barrier_strain,
            al_well=args.al_well,
            al_barrier=args.al_barrier,
            as_well=args.as_well,
            as_barrier=args.as_barrier,
            eg_offset_well_eV=resolved_calibration.Eg_offset_well_eV,
            eg_offset_barrier_eV=resolved_calibration.Eg_offset_barrier_eV,
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
        broadening_eV=resolved_calibration.broadening_eV,
        line_shape=resolved_calibration.line_shape,
        gain_scale_cm=resolved_calibration.gain_scale_cm,
    )
    rows = spectrum_to_rows(spectrum)
    result = {
        "model_note": (
            "Compact screening model: scalar conduction band plus axial 4x4 valence "
            "k.p represented by one 2x2 HH/LH Kramers block. Absolute gain requires "
            "calibration of offsets, broadening, and gain_scale_cm."
        ),
        "design": design,
        "calibration": calibration_summary(
            calibration,
            resolved_calibration,
            applied_qc=float(design["qc"]),
        ),
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
        json.dumps(json_safe(result), ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    write_csv(rows, args.out_csv)
    write_plot(rows, args.plot)
    if args.band_plot is not None:
        plot_band_diagram(profile, subbands, args.band_plot)
    if args.wavefunction_plot is not None:
        plot_wavefunctions(profile, subbands, args.wavefunction_plot)
    if args.dispersion_plot is not None:
        plot_subband_dispersion(subbands, args.dispersion_plot)
    print(format_summary(result, args.out_json, args.out_csv, args.plot))


if __name__ == "__main__":
    main()
