from pathlib import Path

from BasicMQWDesign import design_default
from gain import calculate_gain_spectrum, spectrum_to_rows
from kp_solver import solve_kp_subbands
from visualization import (
    plot_band_diagram,
    plot_gain_spectrum,
    plot_subband_dispersion,
    plot_sweep_summary,
    plot_wavefunctions,
)


def _small_result():
    design = design_default(family="ingaasp")
    profile, subbands = solve_kp_subbands(
        design,
        dz_nm=0.2,
        kt_max_nm=0.2,
        kt_points=3,
        electron_states=1,
        hole_states=1,
    )
    spectrum, _ = calculate_gain_spectrum(
        profile,
        subbands,
        carrier_density_cm3=2.0e18,
        energy_points=40,
    )
    return profile, subbands, spectrum_to_rows(spectrum)


def _assert_png(path: Path) -> None:
    assert path.exists()
    assert path.stat().st_size > 0


def test_plot_gain_spectrum_writes_png(tmp_path: Path) -> None:
    _, _, rows = _small_result()
    path = tmp_path / "gain.png"

    plot_gain_spectrum(rows, path)

    _assert_png(path)


def test_plot_band_diagram_writes_png(tmp_path: Path) -> None:
    profile, subbands, _ = _small_result()
    path = tmp_path / "band.png"

    plot_band_diagram(profile, subbands, path)

    _assert_png(path)


def test_plot_wavefunctions_writes_png(tmp_path: Path) -> None:
    profile, subbands, _ = _small_result()
    path = tmp_path / "wavefunctions.png"

    plot_wavefunctions(profile, subbands, path)

    _assert_png(path)


def test_plot_subband_dispersion_writes_png(tmp_path: Path) -> None:
    _, subbands, _ = _small_result()
    path = tmp_path / "dispersion.png"

    plot_subband_dispersion(subbands, path)

    _assert_png(path)


def test_plot_sweep_summary_writes_png(tmp_path: Path) -> None:
    path = tmp_path / "sweep.png"
    rows = [
        {
            "sweep_value": 1.0,
            "peak_TE_wavelength_nm": 1300.0,
            "peak_TM_wavelength_nm": 1280.0,
            "peak_TE_gain_cm-1": 10.0,
            "peak_TM_gain_cm-1": 5.0,
        },
        {
            "sweep_value": 2.0,
            "peak_TE_wavelength_nm": 1310.0,
            "peak_TM_wavelength_nm": 1290.0,
            "peak_TE_gain_cm-1": 20.0,
            "peak_TM_gain_cm-1": 7.0,
        },
    ]

    plot_sweep_summary(rows, "test sweep", path)

    _assert_png(path)
