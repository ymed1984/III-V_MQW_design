import numpy as np
import pytest

from BasicMQWDesign import design_default
from kp_solver import build_profile_from_design, solve_kp_subbands


def test_profile_shapes_and_thickness() -> None:
    design = design_default(family="ingaasp")
    profile = build_profile_from_design(design, dz_nm=0.2)

    expected_total_nm = (
        design["wells"] * design["well_nm"]
        + (design["wells"] + 1) * design["barrier_nm"]
    )

    assert len(profile.z_nm) > 0
    assert profile.total_nm == pytest.approx(expected_total_nm)
    assert len(profile.V_e_eV) == len(profile.z_nm)
    assert len(profile.V_hh_eV) == len(profile.z_nm)
    assert len(profile.V_lh_eV) == len(profile.z_nm)
    assert np.all(np.isfinite(profile.V_e_eV))
    assert np.all(np.isfinite(profile.V_hh_eV))
    assert np.all(np.isfinite(profile.V_lh_eV))


def test_kp_kt0_matches_finite_well_levels() -> None:
    design = design_default(family="ingaasp")
    _, subbands = solve_kp_subbands(
        design,
        dz_nm=0.2,
        kt_max_nm=0.0,
        kt_points=1,
        electron_states=1,
        hole_states=1,
    )

    e1_reference = design["transition"]["e1_eV"]
    hh1_reference = design["transition"]["hh1_eV"]
    e1_kp = subbands.electron[0][0].energy_eV
    hh1_kp = subbands.valence[0][0].energy_eV

    assert abs(e1_kp - e1_reference) < 0.02
    assert abs(hh1_kp - hh1_reference) < 0.02


def test_electron_dispersion_increases_with_kt() -> None:
    design = design_default(family="ingaasp")
    _, subbands = solve_kp_subbands(
        design,
        dz_nm=0.2,
        kt_max_nm=0.2,
        kt_points=3,
        electron_states=1,
        hole_states=1,
    )

    assert subbands.electron[0][0].energy_eV < subbands.electron[-1][0].energy_eV


def test_wavefunctions_are_normalized() -> None:
    design = design_default(family="ingaasp")
    profile, subbands = solve_kp_subbands(
        design,
        dz_nm=0.2,
        kt_max_nm=0.0,
        kt_points=1,
        electron_states=1,
        hole_states=1,
    )

    electron_norm = np.sum(np.abs(subbands.electron[0][0].psi) ** 2) * profile.dz_nm
    hole = subbands.valence[0][0]
    hole_norm = (
        np.sum(np.abs(hole.hh) ** 2) + np.sum(np.abs(hole.lh) ** 2)
    ) * profile.dz_nm

    assert electron_norm == pytest.approx(1.0, abs=1e-3)
    assert hole_norm == pytest.approx(1.0, abs=1e-3)
