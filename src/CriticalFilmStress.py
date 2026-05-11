#!/usr/bin/env python3
"""Estimate critical film stress and simple critical-thickness metrics on InP.

This helper script reuses the material models in BasicMQWDesign.py and reports:
- in-plane strain on InP,
- biaxial film stress estimate,
- Matthews-Blakeslee critical thickness estimate.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from BasicMQWDesign import (
    BIN,
    make_algainas,
    make_ingaasp,
    matthews_blakeslee_hc_nm,
    strain_parallel,
)


def biaxial_modulus_gpa(c11_gpa: float, c12_gpa: float) -> float:
    """Return biaxial modulus M = C11 + C12 - 2*C12^2/C11 for (001) cubic film."""
    if c11_gpa <= 0:
        raise ValueError("C11 must be positive")
    return c11_gpa + c12_gpa - 2.0 * (c12_gpa * c12_gpa) / c11_gpa


def film_stress_gpa(eps_parallel: float, c11_gpa: float, c12_gpa: float) -> float:
    """In-plane biaxial stress estimate in GPa."""
    return biaxial_modulus_gpa(c11_gpa, c12_gpa) * eps_parallel


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Calculate film stress and critical thickness on InP substrate"
    )
    parser.add_argument("--family", choices=["ingaasp", "algainas"], default="ingaasp")
    parser.add_argument(
        "--strain",
        type=float,
        required=True,
        help="Target strain eps=(a_sub-a_layer)/a_layer (negative=compressive)",
    )
    parser.add_argument("--thickness-nm", type=float, default=None)

    # InGaAsP knobs
    parser.add_argument("--as-frac", type=float, default=0.567)
    parser.add_argument("--ga-frac", type=float, default=None)

    # AlGaInAs knobs
    parser.add_argument("--al-frac", type=float, default=0.14)
    parser.add_argument("--ga-frac-iii", type=float, default=None)

    parser.add_argument("--json", type=Path, default=None)
    args = parser.parse_args()

    if args.family == "ingaasp":
        mat = make_ingaasp(y_As=args.as_frac, strain_target=args.strain, x_Ga=args.ga_frac)
    else:
        mat = make_algainas(
            x_Al=args.al_frac, strain_target=args.strain, y_Ga=args.ga_frac_iii
        )

    eps = strain_parallel(mat, BIN["InP"])
    stress = film_stress_gpa(eps, mat.C11_GPa, mat.C12_GPa)
    hc_nm = matthews_blakeslee_hc_nm(abs(eps), substrate=BIN["InP"])

    result = {
        "family": args.family,
        "material": mat.name,
        "target_strain": args.strain,
        "resolved_strain": eps,
        "biaxial_modulus_GPa": biaxial_modulus_gpa(mat.C11_GPa, mat.C12_GPa),
        "film_stress_GPa": stress,
        "critical_thickness_nm_est": hc_nm,
    }

    if args.thickness_nm is not None:
        if args.thickness_nm <= 0:
            raise ValueError("thickness-nm must be positive")
        result["thickness_nm"] = args.thickness_nm
        result["stress_thickness_GPa_nm"] = stress * args.thickness_nm
        result["over_critical_thickness"] = args.thickness_nm > hc_nm

    print("=== Critical film stress estimate ===")
    for key, value in result.items():
        if isinstance(value, float):
            print(f"{key:28s}: {value:.8g}")
        else:
            print(f"{key:28s}: {value}")

    if args.json is not None:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8")
        print(f"saved json -> {args.json}")


if __name__ == "__main__":
    main()
