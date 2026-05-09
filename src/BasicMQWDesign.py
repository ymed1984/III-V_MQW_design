#!/usr/bin/env python3
"""MQW design helper for O-band InP SOA active regions.

This script is intentionally a compact, inspectable first-pass design tool, not a
replacement for a calibrated 6x6/8x8 k.p simulator. It estimates:

  * alloy composition on InP for InGaAsP or AlGaInAs,
  * pseudomorphic (001) biaxial strain and HH/LH strained bandgaps,
  * finite-well e1/hh1/lh1 confinement levels by a 1D effective-mass solver,
  * approximate e1-hh1/e1-lh1 transition wavelengths,
  * strain balance and a simple Matthews-Blakeslee critical-thickness estimate,
  * Lumerical MQW input snippets for later 4x4/6x6/8x8 k.p gain simulation.

Conventions
-----------
Strain is eps_parallel = (a_substrate - a_layer) / a_layer.
Negative eps_parallel means compressive strain in the layer, consistent with
Ansys Lumerical's MQW strain convention.

Composition notation
--------------------
InGaAsP:     In_{1-x}Ga_x As_y P_{1-y}; user-facing x is x_Ga, y is y_As.
AlGaInAs:    Al_x Ga_y In_{1-x-y} As; user-facing x is x_Al, y is y_Ga.

Limits
------
The material database below uses common 300 K III-V parameter values and simple
Vegard/bowing interpolation. For tape-out or epitaxy release, replace the defaults
with your lab's calibrated PL/XRD/Hall data and validate against an MQW solver
(e.g., Lumerical MQW object or Nextnano/Silvaco/SimuLase style tools).
"""

from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Iterable, Literal, Optional

import numpy as np
from scipy.sparse import diags
from scipy.sparse.linalg import eigsh

# Physical constants
C_LIGHT = 299792458  # m/s
M0 = 9.10938 * 10 ** (-31)  # kg
HC_EV_UM = 1.239841984  # eV*um
H_PRANCK = 4.135667696 * 10 ** (-15)  # eV*s
HBAR = H_PRANCK / (2 * math.pi)  # eV*s
HBAR2_OVER_2M0_EV_NM2 = 0.0380998212  # hbar^2/(2*m0) in eV*nm^2

DesignFamily = Literal["algainas", "ingaasp"]
MaterialDict = dict[str, Any]
DesignDict = dict[str, Any]
Layer = tuple["Material", float]


@dataclass(frozen=True)
class Material:
    name: str
    family: str
    # Composition fields are family-specific.
    x_Ga: Optional[float] = None
    y_As: Optional[float] = None
    x_Al: Optional[float] = None
    y_Ga: Optional[float] = None
    z_In: Optional[float] = None
    # Material parameters at 300 K unless otherwise stated.
    a_A: float = 0.0  # lattice constant in Angstrom
    Eg_eV: float = 0.0  # unstrained Gamma bandgap
    me: float = 0.0  # electron effective mass / m0
    gamma1: float = 0.0
    gamma2: float = 0.0
    gamma3: float = 0.0
    ac_eV: float = 0.0
    av_eV: float = 0.0
    b_eV: float = 0.0
    C11_GPa: float = 0.0
    C12_GPa: float = 0.0
    eps_static: float = 13.0
    source_note: str = "interpolated"

    @property
    def mhh_z(self) -> float:
        # Heavy-hole mass along [001]. Guard avoids division by zero for rough data.
        return 1.0 / max(1e-6, self.gamma1 - 2.0 * self.gamma2)

    @property
    def mlh_z(self) -> float:
        return 1.0 / max(1e-6, self.gamma1 + 2.0 * self.gamma2)


@dataclass(frozen=True)
class FamilyDefaults:
    well_composition: float
    barrier_composition: float
    well_strain: float
    barrier_strain: float
    qc: float


# Common binary parameters. These are a pragmatic 300 K starter set compiled from
# common III-V parameter tables. Replace with your calibrated parameter set if needed.
BIN: dict[str, Material] = {
    "InP": Material(
        "InP",
        "binary",
        a_A=5.8687,
        Eg_eV=1.344,
        me=0.077,
        gamma1=5.08,
        gamma2=1.60,
        gamma3=2.10,
        ac_eV=-6.0,
        av_eV=0.6,
        b_eV=-1.7,
        C11_GPa=101.1,
        C12_GPa=56.1,
        eps_static=12.5,
        source_note="starter binary parameter",
    ),
    "InAs": Material(
        "InAs",
        "binary",
        a_A=6.0583,
        Eg_eV=0.354,
        me=0.023,
        gamma1=20.0,
        gamma2=8.5,
        gamma3=9.2,
        ac_eV=-5.08,
        av_eV=1.0,
        b_eV=-1.8,
        C11_GPa=83.29,
        C12_GPa=45.26,
        eps_static=15.15,
        source_note="starter binary parameter",
    ),
    "GaAs": Material(
        "GaAs",
        "binary",
        a_A=5.65325,
        Eg_eV=1.424,
        me=0.067,
        gamma1=6.98,
        gamma2=2.06,
        gamma3=2.93,
        ac_eV=-7.17,
        av_eV=1.16,
        b_eV=-2.0,
        C11_GPa=122.1,
        C12_GPa=56.6,
        eps_static=12.9,
        source_note="starter binary parameter",
    ),
    "GaP": Material(
        "GaP",
        "binary",
        a_A=5.4505,
        Eg_eV=2.78,
        me=0.15,
        gamma1=4.05,
        gamma2=0.49,
        gamma3=1.25,
        ac_eV=-8.2,
        av_eV=1.7,
        b_eV=-1.6,
        C11_GPa=140.5,
        C12_GPa=62.0,
        eps_static=11.1,
        source_note="starter binary parameter",
    ),
    "AlAs": Material(
        "AlAs",
        "binary",
        a_A=5.6611,
        Eg_eV=3.13,
        me=0.15,
        gamma1=3.76,
        gamma2=0.82,
        gamma3=1.42,
        ac_eV=-5.64,
        av_eV=2.47,
        b_eV=-2.3,
        C11_GPa=125.0,
        C12_GPa=53.4,
        eps_static=10.1,
        source_note="starter binary parameter",
    ),
}

PARAM_KEYS = [
    "a_A",
    "Eg_eV",
    "me",
    "gamma1",
    "gamma2",
    "gamma3",
    "ac_eV",
    "av_eV",
    "b_eV",
    "C11_GPa",
    "C12_GPa",
    "eps_static",
]

FAMILY_DEFAULTS: dict[DesignFamily, FamilyDefaults] = {
    # Strain-compensated TE-gain candidate around 1.31 um.
    "algainas": FamilyDefaults(
        well_composition=0.14,
        barrier_composition=0.30,
        well_strain=-0.007,
        barrier_strain=+0.004,
        qc=0.65,
    ),
    # InGaAsP needs a longer-wavelength well baseline because confinement blueshifts e1-hh1.
    "ingaasp": FamilyDefaults(
        well_composition=0.567,
        barrier_composition=0.30,
        well_strain=-0.006,
        barrier_strain=+0.003,
        qc=0.40,
    ),
}


def _weighted_mix(
    name: str, family: str, weights: dict[str, float], **kwargs
) -> Material:
    vals = {
        k: sum(w * getattr(BIN[b], k) for b, w in weights.items()) for k in PARAM_KEYS
    }
    return Material(name=name, family=family, **kwargs, **vals)


def _ingaasp_x_ga_for_lattice(y_As: float, a_target_A: float) -> float:
    a_no_ga = y_As * BIN["InAs"].a_A + (1.0 - y_As) * BIN["InP"].a_A
    a_all_ga = y_As * BIN["GaAs"].a_A + (1.0 - y_As) * BIN["GaP"].a_A
    return (a_target_A - a_no_ga) / (a_all_ga - a_no_ga)


def _ingaasp_linear_gap_eV(y_As: float, x_Ga: float) -> float:
    weights = {
        "InAs": (1.0 - x_Ga) * y_As,
        "GaAs": x_Ga * y_As,
        "InP": (1.0 - x_Ga) * (1.0 - y_As),
        "GaP": x_Ga * (1.0 - y_As),
    }
    return sum(weight * BIN[binary].Eg_eV for binary, weight in weights.items())


def strain_parallel(material: Material, substrate: Material = BIN["InP"]) -> float:
    """Pseudomorphic in-plane strain on substrate. Negative = compressive."""
    return (substrate.a_A - material.a_A) / material.a_A


def strain_shifts_001(
    material: Material, substrate: Material = BIN["InP"]
) -> dict[str, float]:
    """Return biaxial strain state and strained bandgaps for a (001) layer."""
    epar = strain_parallel(material, substrate)
    ezz = -2.0 * material.C12_GPa / material.C11_GPa * epar
    hydro = 2.0 * epar + ezz
    dEc = material.ac_eV * hydro
    dEv_hydro = material.av_eV * hydro
    dEv_hh = dEv_hydro - material.b_eV * (ezz - epar)
    dEv_lh = dEv_hydro + material.b_eV * (ezz - epar)
    return {
        "eps_parallel": epar,
        "eps_zz": ezz,
        "hydrostatic_strain": hydro,
        "dEc_eV": dEc,
        "dEv_hh_eV": dEv_hh,
        "dEv_lh_eV": dEv_lh,
        "Eg_hh_eV": material.Eg_eV + dEc - dEv_hh,
        "Eg_lh_eV": material.Eg_eV + dEc - dEv_lh,
        "hh_lh_split_eV": dEv_hh - dEv_lh,
    }


def make_ingaasp(
    y_As: float, strain_target: Optional[float] = 0.0, x_Ga: Optional[float] = None
) -> Material:
    """Create In_{1-x}Ga_x As_y P_{1-y}; solve x_Ga from target strain if omitted."""
    if not (0.0 <= y_As <= 1.0):
        raise ValueError("y_As must be in [0, 1]")
    if x_Ga is None:
        if strain_target is None:
            strain_target = 0.0
        a_target = BIN["InP"].a_A / (1.0 + strain_target)
        x_Ga = _ingaasp_x_ga_for_lattice(y_As, a_target)
    if not (0.0 <= x_Ga <= 1.0):
        raise ValueError(
            f"x_Ga={x_Ga:.4f} is outside [0, 1]; choose another y_As/strain"
        )
    weights = {
        "InAs": (1.0 - x_Ga) * y_As,
        "GaAs": x_Ga * y_As,
        "InP": (1.0 - x_Ga) * (1.0 - y_As),
        "GaP": x_Ga * (1.0 - y_As),
    }
    mat = _weighted_mix(
        name=f"In{1 - x_Ga:.4f}Ga{x_Ga:.4f}As{y_As:.4f}P{1 - y_As:.4f}",
        family="InGaAsP",
        weights=weights,
        x_Ga=x_Ga,
        y_As=y_As,
    )
    # For the InP-lattice-matched family, the empirical 300 K relation
    # Eg ≈ 1.35 - 0.72*y + 0.12*y^2 is often closer than plain binary mixing.
    # For strained off-lattice cases, keep that lattice-matched baseline and add
    # the linear-mixing gap change caused by moving x_Ga away from the
    # lattice-matched composition at the same y_As.
    x_Ga_lm = _ingaasp_x_ga_for_lattice(y_As, BIN["InP"].a_A)
    eg_lm = 1.35 - 0.72 * y_As + 0.12 * y_As * y_As
    eg_linear_delta = _ingaasp_linear_gap_eV(y_As, x_Ga) - _ingaasp_linear_gap_eV(
        y_As, x_Ga_lm
    )
    mat = Material(
        **{
            **asdict(mat),
            "Eg_eV": eg_lm + eg_linear_delta,
            "source_note": (
                "InGaAsP interpolated; Eg uses common InP-lattice-matched formula "
                "plus off-lattice correction"
            ),
        }
    )
    return mat


def make_algainas(
    x_Al: float, strain_target: Optional[float] = 0.0, y_Ga: Optional[float] = None
) -> Material:
    """Create Al_x Ga_y In_{1-x-y} As; solve y_Ga from target strain if omitted."""
    if not (0.0 <= x_Al <= 1.0):
        raise ValueError("x_Al must be in [0, 1]")
    if y_Ga is None:
        if strain_target is None:
            strain_target = 0.0
        a_target = BIN["InP"].a_A / (1.0 + strain_target)
        y_Ga = (a_target - (1.0 - x_Al) * BIN["InAs"].a_A - x_Al * BIN["AlAs"].a_A) / (
            BIN["GaAs"].a_A - BIN["InAs"].a_A
        )
    z_In = 1.0 - x_Al - y_Ga
    if min(x_Al, y_Ga, z_In) < -1e-9:
        raise ValueError(
            f"Invalid AlGaInAs composition: Al={x_Al:.4f}, Ga={y_Ga:.4f}, In={z_In:.4f}"
        )
    weights = {"AlAs": x_Al, "GaAs": y_Ga, "InAs": z_In}
    mat = _weighted_mix(
        name=f"Al{x_Al:.4f}Ga{y_Ga:.4f}In{z_In:.4f}As",
        family="AlGaInAs",
        weights=weights,
        x_Al=x_Al,
        y_Ga=y_Ga,
        z_In=z_In,
    )
    # Approximate direct-gap bowing terms (eV). These are deliberately exposed in code.
    b_InGaAs = 0.477
    b_InAlAs = 0.70
    b_AlGaAs = 0.127
    eg_bowed = (
        mat.Eg_eV
        - b_InGaAs * z_In * y_Ga
        - b_InAlAs * z_In * x_Al
        - b_AlGaAs * x_Al * y_Ga
    )
    mat = Material(
        **{
            **asdict(mat),
            "Eg_eV": eg_bowed,
            "source_note": "AlGaInAs interpolated with simple pairwise Gamma-gap bowing",
        }
    )
    return mat


def finite_well_levels(
    well_width_nm: float,
    barrier_width_nm: float,
    m_well: float,
    m_barrier: float,
    barrier_height_eV: float,
    dz_nm: float = 0.025,
    n_eigs: int = 4,
) -> list[float]:
    """Solve bound states of a symmetric finite quantum well.

    Energies are positive confinement energies above the well band edge.
    """
    if well_width_nm <= 0 or barrier_width_nm <= 0:
        raise ValueError("well_width_nm and barrier_width_nm must be positive")
    if m_well <= 0 or m_barrier <= 0:
        raise ValueError("m_well and m_barrier must be positive")
    if dz_nm <= 0:
        raise ValueError("dz_nm must be positive")
    if n_eigs < 1:
        raise ValueError("n_eigs must be at least 1")
    if barrier_height_eV <= 0:
        return []
    total_nm = well_width_nm + 2.0 * barrier_width_nm
    # Use interior grid for Dirichlet boundary conditions at both ends.
    n = int(math.ceil(total_nm / dz_nm)) + 1
    z = np.linspace(-total_nm / 2.0, total_nm / 2.0, n)
    dz = z[1] - z[0]
    interior = z[1:-1]
    in_well = np.abs(interior) <= well_width_nm / 2.0
    V = np.where(in_well, 0.0, barrier_height_eV)
    m = np.where(in_well, m_well, m_barrier)
    invm = 1.0 / m
    if len(invm) < 3:
        raise ValueError("dz_nm is too coarse for the requested well/barrier thickness")
    invm_half_right = 0.5 * (invm[:-1] + invm[1:])
    # For N interior points, off-diagonal arrays have N-1 elements.
    main = np.empty_like(invm)
    main[0] = invm[0] + invm_half_right[0]
    main[-1] = invm_half_right[-1] + invm[-1]
    if len(invm) > 2:
        main[1:-1] = invm_half_right[:-1] + invm_half_right[1:]
    pref = HBAR2_OVER_2M0_EV_NM2 / (dz * dz)
    H = diags(
        diagonals=[-pref * invm_half_right, pref * main + V, -pref * invm_half_right],
        offsets=[-1, 0, 1],
        format="csr",
    )
    k = min(n_eigs, max(1, len(interior) - 2))
    vals = eigsh(H, k=k, which="SA", return_eigenvectors=False, tol=1e-10)
    vals = sorted(float(v) for v in vals if v < barrier_height_eV * 0.999)
    return vals


def transition_estimate(
    well: Material,
    barrier: Material,
    well_nm: float,
    barrier_nm: float,
    qc: float,
    dz_nm: float = 0.025,
) -> dict[str, float | list[float] | dict[str, float]]:
    """Estimate e1-hh1 and e1-lh1 transitions using finite QW confinement."""
    if not (0.0 <= qc <= 1.0):
        raise ValueError("qc must be in [0, 1]")
    sw = strain_shifts_001(well)
    sb = strain_shifts_001(barrier)
    # Offset model: use strained HH bandgap difference and supplied conduction-band offset ratio.
    # This is a first-pass approximation. For Lumerical/nextnano, prefer absolute VBO/model-solid offsets.
    dEg_hh = max(0.0, sb["Eg_hh_eV"] - sw["Eg_hh_eV"])
    dEg_lh = max(0.0, sb["Eg_lh_eV"] - sw["Eg_lh_eV"])
    dEc = qc * dEg_hh
    dEv_hh = (1.0 - qc) * dEg_hh
    dEv_lh = (1.0 - qc) * dEg_lh
    e_levels = finite_well_levels(
        well_nm, barrier_nm, well.me, barrier.me, dEc, dz_nm=dz_nm
    )
    hh_levels = finite_well_levels(
        well_nm, barrier_nm, well.mhh_z, barrier.mhh_z, dEv_hh, dz_nm=dz_nm
    )
    lh_levels = finite_well_levels(
        well_nm, barrier_nm, well.mlh_z, barrier.mlh_z, dEv_lh, dz_nm=dz_nm
    )
    Ee1 = e_levels[0] if e_levels else float("nan")
    Ehh1 = hh_levels[0] if hh_levels else float("nan")
    Elh1 = lh_levels[0] if lh_levels else float("nan")
    Etr_hh = sw["Eg_hh_eV"] + Ee1 + Ehh1
    Etr_lh = sw["Eg_lh_eV"] + Ee1 + Elh1
    return {
        "well_strain": sw,
        "barrier_strain": sb,
        "dEg_hh_eV": dEg_hh,
        "dEc_eV": dEc,
        "dEv_hh_eV": dEv_hh,
        "dEv_lh_eV": dEv_lh,
        "electron_levels_eV": e_levels,
        "hh_levels_eV": hh_levels,
        "lh_levels_eV": lh_levels,
        "e1_eV": Ee1,
        "hh1_eV": Ehh1,
        "lh1_eV": Elh1,
        "E_transition_hh_eV": Etr_hh,
        "lambda_transition_hh_um": HC_EV_UM / Etr_hh,
        "E_transition_lh_eV": Etr_lh,
        "lambda_transition_lh_um": HC_EV_UM / Etr_lh,
    }


def strain_balance(layers: Iterable[Layer]) -> dict[str, float]:
    """Return thickness-weighted strain metrics for an MQW stack."""
    total_nm = 0.0
    sum_eps_t = 0.0
    sum_abs_eps_t = 0.0
    for mat, t_nm in layers:
        if t_nm <= 0:
            raise ValueError("layer thicknesses must be positive")
        eps = strain_parallel(mat)
        total_nm += t_nm
        sum_eps_t += eps * t_nm
        sum_abs_eps_t += abs(eps) * t_nm
    return {
        "total_nm": total_nm,
        "average_strain": sum_eps_t / total_nm,
        "absolute_strain_thickness_nm": sum_abs_eps_t,
        "signed_strain_thickness_nm": sum_eps_t,
    }


def matthews_blakeslee_hc_nm(
    abs_mismatch: float,
    substrate: Material = BIN["InP"],
    poisson: float = 0.35,
    alpha_deg: float = 60.0,
) -> float:
    """Approximate single-layer critical thickness in nm for a 60-degree dislocation.

    This is a screening number only; MQW strain balance and growth conditions matter.
    """
    f = abs(abs_mismatch)
    if f <= 0:
        return float("inf")
    a_nm = substrate.a_A * 0.1
    # Approximate Burgers vector for a 60-degree dislocation.
    b_nm = a_nm / math.sqrt(2.0)
    alpha = math.radians(alpha_deg)
    pref = (
        b_nm
        / (4.0 * math.pi * f)
        * (1.0 - poisson * math.cos(alpha) ** 2)
        / (1.0 + poisson)
    )
    h = max(b_nm * 2.0, pref * 10.0)
    for _ in range(100):
        h_new = max(b_nm, pref * max(1.0, math.log(h / b_nm)))
        if abs(h_new - h) < 1e-6:
            break
        h = 0.5 * h + 0.5 * h_new
    return h


def build_stack(
    well: Material,
    barrier: Material,
    wells: int,
    well_nm: float,
    barrier_nm: float,
    outer_barrier_nm: Optional[float] = None,
) -> list[Layer]:
    if wells < 1:
        raise ValueError("wells must be at least 1")
    if well_nm <= 0 or barrier_nm <= 0:
        raise ValueError("well_nm and barrier_nm must be positive")
    if outer_barrier_nm is None:
        outer_barrier_nm = barrier_nm
    if outer_barrier_nm <= 0:
        raise ValueError("outer_barrier_nm must be positive")
    layers: list[Layer] = [(barrier, outer_barrier_nm)]
    for i in range(wells):
        layers.append((well, well_nm))
        layers.append((barrier, outer_barrier_nm if i == wells - 1 else barrier_nm))
    return layers


def _family_materials(
    family: DesignFamily,
    q_c: Optional[float],
    well_strain: Optional[float],
    barrier_strain: Optional[float],
    al_well: Optional[float],
    al_barrier: Optional[float],
    as_well: Optional[float],
    as_barrier: Optional[float],
) -> tuple[Material, Material, float]:
    defaults = FAMILY_DEFAULTS[family]
    resolved_well_strain = defaults.well_strain if well_strain is None else well_strain
    resolved_barrier_strain = (
        defaults.barrier_strain if barrier_strain is None else barrier_strain
    )
    qc = defaults.qc if q_c is None else q_c

    if family == "algainas":
        well = make_algainas(
            x_Al=defaults.well_composition if al_well is None else al_well,
            strain_target=resolved_well_strain,
        )
        barrier = make_algainas(
            x_Al=defaults.barrier_composition if al_barrier is None else al_barrier,
            strain_target=resolved_barrier_strain,
        )
        return well, barrier, qc

    well = make_ingaasp(
        y_As=defaults.well_composition if as_well is None else as_well,
        strain_target=resolved_well_strain,
    )
    barrier = make_ingaasp(
        y_As=defaults.barrier_composition if as_barrier is None else as_barrier,
        strain_target=resolved_barrier_strain,
    )
    return well, barrier, qc


def design_default(
    family: DesignFamily = "ingaasp",
    wells: int = 5,
    well_nm: float = 7.0,
    barrier_nm: float = 10.0,
    q_c: Optional[float] = None,
    well_strain: Optional[float] = None,
    barrier_strain: Optional[float] = None,
    al_well: Optional[float] = None,
    al_barrier: Optional[float] = None,
    as_well: Optional[float] = None,
    as_barrier: Optional[float] = None,
) -> DesignDict:
    """Return a default O-band SOA MQW design candidate."""
    if family not in FAMILY_DEFAULTS:
        raise ValueError("family must be 'algainas' or 'ingaasp'")

    well, barrier, qc = _family_materials(
        family=family,
        q_c=q_c,
        well_strain=well_strain,
        barrier_strain=barrier_strain,
        al_well=al_well,
        al_barrier=al_barrier,
        as_well=as_well,
        as_barrier=as_barrier,
    )
    stack = build_stack(well, barrier, wells, well_nm, barrier_nm)
    tr = transition_estimate(well, barrier, well_nm, barrier_nm, qc)
    return {
        "family": family,
        "wells": wells,
        "well_nm": well_nm,
        "barrier_nm": barrier_nm,
        "qc": qc,
        "well": asdict(well),
        "barrier": asdict(barrier),
        "transition": tr,
        "strain_balance": strain_balance(stack),
        "critical_thickness_well_nm_est": matthews_blakeslee_hc_nm(
            abs(strain_parallel(well))
        ),
        "critical_thickness_barrier_nm_est": matthews_blakeslee_hc_nm(
            abs(strain_parallel(barrier))
        ),
    }


def _lum_material_struct(mat: MaterialDict) -> str:
    if mat["family"] == "AlGaInAs":
        return (
            "struct; "
            f"database_material='Al_{{x}}Ga_{{y}}In_{{1-x-y}}As'; "
            f"x={mat['x_Al']:.8f}; y={mat['y_Ga']:.8f};"
        )
    if mat["family"] == "InGaAsP":
        # Lumerical docs use In_xGa_{1-x}As_yP_{1-y}; our x_Ga -> x_In = 1-x_Ga.
        x_in = 1.0 - float(mat["x_Ga"])
        return (
            "struct; "
            f"database_material='In_{{x}}Ga_{{1-x}}As_{{y}}P_{{1-y}}'; "
            f"x={x_in:.8f}; y={mat['y_As']:.8f};"
        )
    raise ValueError(f"Unsupported family for Lumerical export: {mat['family']}")


def _lumerical_stack_arrays(
    design: DesignDict,
) -> tuple[list[float], list[float], list[str]]:
    well_nm = float(design["well_nm"])
    barrier_nm = float(design["barrier_nm"])
    transition = design["transition"]
    layers = []
    strains = []
    materials = []

    # Symmetric stack: barrier, (well, barrier)xN
    for i in range(int(design["wells"]) * 2 + 1):
        is_well = i % 2 == 1
        strain_key = "well_strain" if is_well else "barrier_strain"
        layers.append(well_nm if is_well else barrier_nm)
        strains.append(float(transition[strain_key]["eps_parallel"]))
        materials.append("well_mat" if is_well else "barrier_mat")

    return layers, strains, materials


def write_lumerical_lsf(design: DesignDict, path: str | Path) -> Path:
    """Write a starter Lumerical script-command MQW snippet.

    Current Ansys versions prefer the MQW solver object for newest features, but this
    script-command version is convenient for automation and benchmarking.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    well = design["well"]
    barrier = design["barrier"]
    wells = int(design["wells"])
    well_nm = float(design["well_nm"])
    barrier_nm = float(design["barrier_nm"])
    layers, strains, materials = _lumerical_stack_arrays(design)
    lengths_lsf = "; ".join(f"{t:.9g}e-9" for t in layers)
    strains_lsf = "; ".join(f"{e:.9g}" for e in strains)
    materials_lsf = "; ".join(materials)
    text = f"""# Generated by BasicMQWDesign.py
# Family: {design["family"]}, wells: {wells}, well={well_nm} nm, barrier={barrier_nm} nm
# Estimated e1-hh1 transition: {design["transition"]["lambda_transition_hh_um"]:.4f} um
# Sign convention: negative strain = compressive.

well_mat = {_lum_material_struct(well)};
barrier_mat = {_lum_material_struct(barrier)};

stack = struct;
stack.length = [{lengths_lsf}];
stack.material = {{{materials_lsf}}};
stack.strain = [{strains_lsf}];
stack.gamma = 0.030; # eV Lorentzian FWHM; adjust to measured linewidth
stack.vb = struct;
stack.vb.method = 'palankovski'; # replace with 'direct' offsets if using calibrated VBOs

sim = struct;
sim.T = 300;
sim.kt = linspace(0, 2*pi/{BIN["InP"].a_A * 1e-10}*0.10, 61);
sim.cden = [1e24, 2e24, 3e24, 4e24]; # average carrier density [m^-3] over full MQW span

cfg = struct;
cfg.pml = 1;

# The script command is useful for automation; for newest features use the MQW solver object.
out = mqwgain(stack, sim, cfg);
visualize(out.emission);
"""
    path.write_text(text, encoding="utf-8")
    return path


def write_design_json(design: DesignDict, path: str | Path) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(design, indent=2), encoding="utf-8")
    return path


def build_arg_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        description="First-pass O-band InP MQW SOA design helper"
    )
    ap.add_argument("--family", choices=["algainas", "ingaasp"], default="ingaasp")
    ap.add_argument("--wells", type=int, default=5)
    ap.add_argument("--well-nm", type=float, default=7.0)
    ap.add_argument("--barrier-nm", type=float, default=10.0)
    ap.add_argument(
        "--qc", type=float, default=None, help="Conduction band offset fraction"
    )
    ap.add_argument(
        "--well-strain",
        type=float,
        default=None,
        help="Target well strain eps=(a_sub-a_layer)/a_layer; negative is compressive",
    )
    ap.add_argument(
        "--barrier-strain",
        type=float,
        default=None,
        help="Target barrier strain eps=(a_sub-a_layer)/a_layer; positive is tensile",
    )
    ap.add_argument(
        "--al-well",
        type=float,
        default=None,
        help="Al fraction in Al_xGa_yIn_1-x-yAs well",
    )
    ap.add_argument(
        "--al-barrier",
        type=float,
        default=None,
        help="Al fraction in Al_xGa_yIn_1-x-yAs barrier",
    )
    ap.add_argument(
        "--as-well",
        type=float,
        default=None,
        help="As fraction y in In_1-xGa_xAs_yP_1-y well",
    )
    ap.add_argument(
        "--as-barrier",
        type=float,
        default=None,
        help="As fraction y in In_1-xGa_xAs_yP_1-y barrier",
    )
    ap.add_argument("--json", type=Path, default=Path("out/ingaasp_design_result.json"))
    ap.add_argument("--lsf", type=Path, default=Path("out/ingaasp_lumerical_input.lsf"))
    return ap


def format_summary(design: DesignDict, json_path: Path, lsf_path: Path) -> str:
    tr = design["transition"]
    sb = design["strain_balance"]
    return "\n".join(
        (
            "=== MQW SOA first-pass design ===",
            f"family             : {design['family']}",
            f"wells              : {design['wells']}",
            f"well/barrier       : {design['well_nm']} nm / {design['barrier_nm']} nm",
            f"well material      : {design['well']['name']}",
            f"barrier material   : {design['barrier']['name']}",
            f"well strain        : {tr['well_strain']['eps_parallel'] * 100:+.3f} %",
            f"barrier strain     : {tr['barrier_strain']['eps_parallel'] * 100:+.3f} %",
            f"avg stack strain   : {sb['average_strain'] * 100:+.3f} %",
            f"DeltaEc/DeltaEv_hh : {tr['dEc_eV']:.3f} / {tr['dEv_hh_eV']:.3f} eV",
            f"e1/hh1/lh1         : {tr['e1_eV']:.4f} / {tr['hh1_eV']:.4f} / {tr['lh1_eV']:.4f} eV",
            f"e1-hh1 wavelength  : {tr['lambda_transition_hh_um'] * 1000:.1f} nm",
            f"e1-lh1 wavelength  : {tr['lambda_transition_lh_um'] * 1000:.1f} nm",
            f"wrote JSON         : {json_path}",
            f"wrote LSF          : {lsf_path}",
        )
    )


def main(argv: Optional[list[str]] = None) -> None:
    args = build_arg_parser().parse_args(argv)
    design = design_default(
        args.family,
        args.wells,
        args.well_nm,
        args.barrier_nm,
        args.qc,
        well_strain=args.well_strain,
        barrier_strain=args.barrier_strain,
        al_well=args.al_well,
        al_barrier=args.al_barrier,
        as_well=args.as_well,
        as_barrier=args.as_barrier,
    )
    json_path = write_design_json(design, args.json)
    lsf_path = write_lumerical_lsf(design, args.lsf)
    print(format_summary(design, json_path, lsf_path))


if __name__ == "__main__":
    main()
