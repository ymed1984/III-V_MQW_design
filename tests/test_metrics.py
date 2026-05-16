import pytest

from metrics import fwhm_from_spectrum, gain_at_wavelength, peak_metrics, spectrum_rmse


def _rows():
    return [
        {"energy_eV": 0.9, "wavelength_nm": 1377.6, "gain_TE_cm-1": 0.0, "gain_TM_cm-1": 0.0},
        {"energy_eV": 1.0, "wavelength_nm": 1239.8, "gain_TE_cm-1": 10.0, "gain_TM_cm-1": 4.0},
        {"energy_eV": 1.1, "wavelength_nm": 1127.1, "gain_TE_cm-1": 20.0, "gain_TM_cm-1": 8.0},
        {"energy_eV": 1.2, "wavelength_nm": 1033.2, "gain_TE_cm-1": 10.0, "gain_TM_cm-1": 4.0},
        {"energy_eV": 1.3, "wavelength_nm": 953.7, "gain_TE_cm-1": 0.0, "gain_TM_cm-1": 0.0},
    ]


def test_peak_metrics_and_fwhm() -> None:
    metrics = peak_metrics(_rows(), "TE")

    assert metrics.peak_energy_eV == pytest.approx(1.1)
    assert metrics.peak_gain_cm == pytest.approx(20.0)
    assert metrics.fwhm_meV == pytest.approx(200.0)


def test_gain_at_wavelength_interpolates() -> None:
    gain = gain_at_wavelength(_rows(), 1127.1, "TE")

    assert gain == pytest.approx(20.0)


def test_spectrum_rmse_zero_for_identical_rows() -> None:
    assert spectrum_rmse(_rows(), _rows(), "TE") == pytest.approx(0.0)


def test_fwhm_returns_nan_for_non_positive_peak() -> None:
    rows = [
        {"energy_eV": 1.0, "wavelength_nm": 1239.8, "gain_TE_cm-1": -1.0, "gain_TM_cm-1": 0.0}
    ]

    assert fwhm_from_spectrum(rows, "TE") != fwhm_from_spectrum(rows, "TE")
