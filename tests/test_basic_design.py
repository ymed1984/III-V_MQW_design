import math

import pytest

from BasicMQWDesign import design_default, make_algainas


def test_default_ingaasp_transition_near_oband() -> None:
    design = design_default(family="ingaasp")

    wavelength_nm = design["transition"]["lambda_transition_hh_um"] * 1000.0

    assert 1250.0 < wavelength_nm < 1370.0
    assert math.isfinite(wavelength_nm)


def test_default_strain_signs_and_composition_ranges() -> None:
    design = design_default(family="ingaasp")

    assert design["transition"]["well_strain"]["eps_parallel"] < 0.0
    assert design["transition"]["barrier_strain"]["eps_parallel"] > 0.0
    assert 0.0 <= design["well"]["x_Ga"] <= 1.0
    assert 0.0 <= design["well"]["y_As"] <= 1.0
    assert 0.0 <= design["barrier"]["x_Ga"] <= 1.0
    assert 0.0 <= design["barrier"]["y_As"] <= 1.0


def test_wider_well_redshifts_hh_transition() -> None:
    narrow = design_default(family="ingaasp", well_nm=6.0)
    nominal = design_default(family="ingaasp", well_nm=7.0)
    wide = design_default(family="ingaasp", well_nm=8.0)

    lambda_narrow = narrow["transition"]["lambda_transition_hh_um"]
    lambda_nominal = nominal["transition"]["lambda_transition_hh_um"]
    lambda_wide = wide["transition"]["lambda_transition_hh_um"]

    assert lambda_narrow < lambda_nominal < lambda_wide


def test_algainas_composition_is_valid() -> None:
    material = make_algainas(x_Al=0.14, strain_target=-0.007)

    assert 0.0 <= material.x_Al <= 1.0
    assert material.y_Ga is not None
    assert material.z_In is not None
    assert material.y_Ga >= 0.0
    assert material.z_In >= 0.0
    assert material.x_Al + material.y_Ga + material.z_In == pytest.approx(1.0)
