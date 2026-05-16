import numpy as np

from BasicMQWDesign import design_default
from gain import calculate_gain_spectrum, spectrum_to_rows
from kp_solver import solve_kp_subbands


def _small_subbands():
    design = design_default(family="ingaasp")
    return solve_kp_subbands(
        design,
        dz_nm=0.2,
        kt_max_nm=0.2,
        kt_points=5,
        electron_states=1,
        hole_states=2,
    )


def test_gain_spectrum_schema_and_finite_values() -> None:
    profile, subbands = _small_subbands()
    spectrum, terms = calculate_gain_spectrum(
        profile,
        subbands,
        carrier_density_cm3=2.0e18,
        energy_points=80,
    )
    rows = spectrum_to_rows(spectrum)

    assert len(terms) > 0
    assert len(spectrum.energy_eV) == 80
    assert len(rows) == 80
    assert set(rows[0]) == {
        "energy_eV",
        "wavelength_nm",
        "gain_TE_cm-1",
        "gain_TM_cm-1",
    }
    assert np.all(np.isfinite(spectrum.energy_eV))
    assert np.all(np.isfinite(spectrum.wavelength_nm))
    assert np.all(np.isfinite(spectrum.gain_te_cm))
    assert np.all(np.isfinite(spectrum.gain_tm_cm))
    assert np.isfinite(spectrum.peak_te_gain_cm)
    assert np.isfinite(spectrum.peak_tm_gain_cm)


def test_higher_density_increases_peak_te_gain() -> None:
    profile, subbands = _small_subbands()
    low_density, _ = calculate_gain_spectrum(
        profile,
        subbands,
        carrier_density_cm3=1.0e16,
        energy_points=80,
    )
    high_density, _ = calculate_gain_spectrum(
        profile,
        subbands,
        carrier_density_cm3=3.0e17,
        energy_points=80,
    )

    assert high_density.peak_te_gain_cm > low_density.peak_te_gain_cm
