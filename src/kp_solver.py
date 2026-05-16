#!/usr/bin/env python3
"""Compact k.p helpers for strained InP MQW gain screening.

The implementation is deliberately scoped as a first-pass solver:

* conduction states use a scalar effective-mass Hamiltonian,
* valence states use one axial 2x2 HH/LH Kramers block of the 4x4
  Luttinger-Kohn Hamiltonian,
* material profiles come from the same strain and offset model used by
  BasicMQWDesign.py.

This is not a replacement for a calibrated 6x6/8x8 k.p package. It gives a
local, inspectable gain-screening path that can be tightened against Lumerical
or nextnano data as calibrated offsets and momentum matrix elements become
available.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Any

import numpy as np
from scipy.sparse import bmat, diags
from scipy.sparse.linalg import ArpackNoConvergence, eigsh

from BasicMQWDesign import (
    HBAR2_OVER_2M0_EV_NM2,
    DesignDict,
    Material,
    build_stack,
    strain_shifts_001,
)


@dataclass(frozen=True)
class MQWProfile:
    z_nm: np.ndarray
    dz_nm: float
    material_names: list[str]
    layer_kind: np.ndarray
    Eg_hh_eV: np.ndarray
    Eg_lh_eV: np.ndarray
    V_e_eV: np.ndarray
    V_hh_eV: np.ndarray
    V_lh_eV: np.ndarray
    me_z: np.ndarray
    me_xy: np.ndarray
    mhh_z: np.ndarray
    mlh_z: np.ndarray
    mhh_xy: np.ndarray
    mlh_xy: np.ndarray
    gamma2: np.ndarray
    total_nm: float
    active_nm: float


@dataclass(frozen=True)
class ScalarState:
    energy_eV: float
    psi: np.ndarray


@dataclass(frozen=True)
class ValenceState:
    energy_eV: float
    hh: np.ndarray
    lh: np.ndarray
    hh_weight: float
    lh_weight: float


@dataclass(frozen=True)
class KpSubbands:
    kt_nm: np.ndarray
    electron: list[list[ScalarState]]
    valence: list[list[ValenceState]]


def material_from_dict(data: dict[str, Any]) -> Material:
    fields = set(Material.__dataclass_fields__)
    return Material(**{key: value for key, value in data.items() if key in fields})


def _positive_mass(value: np.ndarray | float, floor: float = 0.02) -> np.ndarray:
    return np.maximum(np.asarray(value, dtype=float), floor)


def _scalar_hamiltonian(
    potential_eV: np.ndarray,
    mass_z: np.ndarray,
    dz_nm: float,
    kt_nm: float = 0.0,
    mass_xy: np.ndarray | None = None,
):
    mass_z = _positive_mass(mass_z)
    invm = 1.0 / mass_z
    invm_half = 0.5 * (invm[:-1] + invm[1:])
    main = np.empty_like(invm)
    main[0] = invm[0] + invm_half[0]
    main[-1] = invm_half[-1] + invm[-1]
    if len(invm) > 2:
        main[1:-1] = invm_half[:-1] + invm_half[1:]
    in_plane = 0.0
    if mass_xy is not None and kt_nm > 0.0:
        in_plane = HBAR2_OVER_2M0_EV_NM2 * kt_nm * kt_nm / _positive_mass(mass_xy)
    pref = HBAR2_OVER_2M0_EV_NM2 / (dz_nm * dz_nm)
    return diags(
        diagonals=[
            -pref * invm_half,
            pref * main + potential_eV + in_plane,
            -pref * invm_half,
        ],
        offsets=[-1, 0, 1],
        format="csr",
    )


def _normalize(psi: np.ndarray, dz_nm: float) -> np.ndarray:
    norm = float(np.sqrt(np.sum(np.abs(psi) ** 2) * dz_nm))
    if norm <= 0:
        return psi
    return psi / norm


def _lowest_eigenpairs(hamiltonian, n_states: int) -> tuple[np.ndarray, np.ndarray]:
    k = min(n_states, max(1, hamiltonian.shape[0] - 2))
    try:
        return eigsh(hamiltonian, k=k, sigma=0.0, which="LM", tol=1e-9, maxiter=20000)
    except ArpackNoConvergence as exc:
        if (
            exc.eigenvalues is not None
            and exc.eigenvectors is not None
            and len(exc.eigenvalues) >= k
        ):
            return exc.eigenvalues[:k], exc.eigenvectors[:, :k]
        return eigsh(hamiltonian, k=k, which="SA", tol=1e-8, maxiter=40000)


def build_profile_from_design(design: DesignDict, dz_nm: float = 0.10) -> MQWProfile:
    if dz_nm <= 0:
        raise ValueError("dz_nm must be positive")
    well = material_from_dict(design["well"])
    barrier = material_from_dict(design["barrier"])
    stack = build_stack(
        well,
        barrier,
        int(design["wells"]),
        float(design["well_nm"]),
        float(design["barrier_nm"]),
    )
    total_nm = sum(thickness for _, thickness in stack)
    n = int(np.ceil(total_nm / dz_nm)) + 1
    edges = np.linspace(-total_nm / 2.0, total_nm / 2.0, n)
    z_nm = 0.5 * (edges[:-1] + edges[1:])
    dz_eff = float(edges[1] - edges[0])

    materials: list[Material] = []
    layer_kind = np.empty(len(z_nm), dtype="<U7")
    cursor = -total_nm / 2.0
    for index, (mat, thickness_nm) in enumerate(stack):
        right = cursor + thickness_nm
        if index == len(stack) - 1:
            mask = (z_nm >= cursor) & (z_nm <= right)
        else:
            mask = (z_nm >= cursor) & (z_nm < right)
        materials.extend([mat] * int(np.count_nonzero(mask)))
        layer_kind[mask] = "well" if mat.name == well.name else "barrier"
        cursor = right

    if len(materials) != len(z_nm):
        raise RuntimeError("internal profile construction error")

    shifts = [strain_shifts_001(mat) for mat in materials]
    Eg_hh = np.array([item["Eg_hh_eV"] for item in shifts], dtype=float)
    Eg_lh = np.array([item["Eg_lh_eV"] for item in shifts], dtype=float)
    qc = float(design["qc"])
    if not (0.0 <= qc <= 1.0):
        raise ValueError("qc must be in [0, 1]")

    V_e = qc * (Eg_hh - float(np.min(Eg_hh)))
    V_hh = (1.0 - qc) * (Eg_hh - float(np.min(Eg_hh)))
    V_lh = (1.0 - qc) * (Eg_lh - float(np.min(Eg_lh)))

    gamma1 = np.array([mat.gamma1 for mat in materials], dtype=float)
    gamma2 = np.array([mat.gamma2 for mat in materials], dtype=float)
    me = np.array([mat.me for mat in materials], dtype=float)
    mhh_z = _positive_mass(1.0 / np.maximum(gamma1 - 2.0 * gamma2, 1e-6))
    mlh_z = _positive_mass(1.0 / np.maximum(gamma1 + 2.0 * gamma2, 1e-6))
    mhh_xy = _positive_mass(1.0 / np.maximum(gamma1 + gamma2, 1e-6))
    mlh_xy = _positive_mass(1.0 / np.maximum(gamma1 - gamma2, 1e-6))

    return MQWProfile(
        z_nm=z_nm,
        dz_nm=dz_eff,
        material_names=[mat.name for mat in materials],
        layer_kind=layer_kind,
        Eg_hh_eV=Eg_hh,
        Eg_lh_eV=Eg_lh,
        V_e_eV=V_e,
        V_hh_eV=V_hh,
        V_lh_eV=V_lh,
        me_z=me,
        me_xy=me,
        mhh_z=mhh_z,
        mlh_z=mlh_z,
        mhh_xy=mhh_xy,
        mlh_xy=mlh_xy,
        gamma2=gamma2,
        total_nm=total_nm,
        active_nm=float(design["wells"]) * float(design["well_nm"]),
    )


def solve_scalar_states(
    potential_eV: np.ndarray,
    mass_z: np.ndarray,
    dz_nm: float,
    n_states: int,
    kt_nm: float = 0.0,
    mass_xy: np.ndarray | None = None,
) -> list[ScalarState]:
    if n_states < 1:
        return []
    hamiltonian = _scalar_hamiltonian(
        potential_eV, mass_z, dz_nm, kt_nm=kt_nm, mass_xy=mass_xy
    )
    vals, vecs = _lowest_eigenpairs(hamiltonian, n_states)
    order = np.argsort(vals)
    states = []
    for idx in order:
        states.append(
            ScalarState(
                energy_eV=float(vals[idx]),
                psi=_normalize(np.asarray(vecs[:, idx]), dz_nm),
            )
        )
    return states


def solve_valence_kp_states(
    profile: MQWProfile,
    kt_nm: float,
    n_states: int,
) -> list[ValenceState]:
    if n_states < 1:
        return []
    hh = _scalar_hamiltonian(
        profile.V_hh_eV,
        profile.mhh_z,
        profile.dz_nm,
        kt_nm=kt_nm,
        mass_xy=profile.mhh_xy,
    )
    lh = _scalar_hamiltonian(
        profile.V_lh_eV,
        profile.mlh_z,
        profile.dz_nm,
        kt_nm=kt_nm,
        mass_xy=profile.mlh_xy,
    )
    # Axial 4x4 LK block reduced to one HH/LH Kramers block. The R coupling gives
    # in-plane HH/LH mixing; S terms are deferred until calibrated validation data
    # are available for this compact solver.
    mix = np.sqrt(3.0) * HBAR2_OVER_2M0_EV_NM2 * profile.gamma2 * kt_nm * kt_nm
    coupling = diags(mix, offsets=0, format="csr")
    hamiltonian = bmat([[hh, coupling], [coupling, lh]], format="csr")
    vals, vecs = _lowest_eigenpairs(hamiltonian, n_states)
    order = np.argsort(vals)
    n_grid = len(profile.z_nm)
    states: list[ValenceState] = []
    for idx in order:
        hh_part = np.asarray(vecs[:n_grid, idx])
        lh_part = np.asarray(vecs[n_grid:, idx])
        norm = float(
            np.sqrt((np.sum(np.abs(hh_part) ** 2) + np.sum(np.abs(lh_part) ** 2)) * profile.dz_nm)
        )
        if norm > 0:
            hh_part = hh_part / norm
            lh_part = lh_part / norm
        hh_weight = float(np.sum(np.abs(hh_part) ** 2) * profile.dz_nm)
        lh_weight = float(np.sum(np.abs(lh_part) ** 2) * profile.dz_nm)
        states.append(
            ValenceState(
                energy_eV=float(vals[idx]),
                hh=hh_part,
                lh=lh_part,
                hh_weight=hh_weight,
                lh_weight=lh_weight,
            )
        )
    return states


def solve_kp_subbands(
    design: DesignDict,
    dz_nm: float = 0.10,
    kt_max_nm: float = 0.12,
    kt_points: int = 31,
    electron_states: int = 2,
    hole_states: int = 4,
) -> tuple[MQWProfile, KpSubbands]:
    if kt_points < 1:
        raise ValueError("kt_points must be at least 1")
    profile = build_profile_from_design(design, dz_nm=dz_nm)
    kt_grid = np.linspace(0.0, kt_max_nm, kt_points)
    e_states: list[list[ScalarState]] = []
    h_states: list[list[ValenceState]] = []
    for kt in kt_grid:
        e_states.append(
            solve_scalar_states(
                profile.V_e_eV,
                profile.me_z,
                profile.dz_nm,
                electron_states,
                kt_nm=float(kt),
                mass_xy=profile.me_xy,
            )
        )
        h_states.append(solve_valence_kp_states(profile, float(kt), hole_states))
    return profile, KpSubbands(kt_nm=kt_grid, electron=e_states, valence=h_states)


def subband_summary(profile: MQWProfile, subbands: KpSubbands) -> dict[str, Any]:
    e0 = subbands.electron[0]
    h0 = subbands.valence[0]
    return {
        "dz_nm": profile.dz_nm,
        "total_nm": profile.total_nm,
        "active_nm": profile.active_nm,
        "kt_max_nm": float(subbands.kt_nm[-1]),
        "kt_points": int(len(subbands.kt_nm)),
        "electron_levels_eV": [asdict(state) | {"psi": None} for state in e0],
        "valence_levels_eV": [
            {
                "energy_eV": state.energy_eV,
                "hh_weight": state.hh_weight,
                "lh_weight": state.lh_weight,
            }
            for state in h0
        ],
    }
