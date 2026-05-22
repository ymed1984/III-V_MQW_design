"""Tests for --design-input (composition-based MQW input JSON)."""

import json
import subprocess
import sys
from pathlib import Path

import pytest

from BasicMQWDesign import design_default, load_design_input

REPO_ROOT = Path(__file__).resolve().parents[1]
FIXTURES = Path(__file__).resolve().parent / "fixtures"


def _run_script(script: str, *args: str) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [sys.executable, "-B", str(REPO_ROOT / "src" / script), *args],
        check=True,
        cwd=REPO_ROOT,
        text=True,
        capture_output=True,
    )


# --- Unit tests for load_design_input ---


def test_load_design_input_oband() -> None:
    kwargs = load_design_input(FIXTURES / "ingaasp_oband_input.json")
    assert kwargs["family"] == "ingaasp"
    assert kwargs["x_Ga_well"] == 0.1755
    assert kwargs["as_well"] == 0.567
    design = design_default(**kwargs)
    wl_nm = design["transition"]["lambda_transition_hh_um"] * 1000.0
    assert 1250.0 < wl_nm < 1400.0


def test_load_design_input_cband() -> None:
    kwargs = load_design_input(FIXTURES / "ingaasp_cband_input.json")
    design = design_default(**kwargs)
    wl_nm = design["transition"]["lambda_transition_hh_um"] * 1000.0
    assert 1450.0 < wl_nm < 1600.0


def test_strain_auto_calculated_from_composition() -> None:
    kwargs = load_design_input(FIXTURES / "ingaasp_oband_input.json")
    design = design_default(**kwargs)
    eps = design["transition"]["well_strain"]["eps_parallel"]
    # Strain should be auto-calculated, not zero
    assert abs(eps) > 1e-4


def test_load_rejects_unknown_keys(tmp_path: Path) -> None:
    bad = {"family": "ingaasp", "wells": 5, "typo_key": 123}
    path = tmp_path / "bad.json"
    path.write_text(json.dumps(bad), encoding="utf-8")
    with pytest.raises(ValueError, match="unknown top-level keys"):
        load_design_input(path)


def test_load_rejects_unknown_layer_keys(tmp_path: Path) -> None:
    bad = {"family": "ingaasp", "well": {"y_As": 0.5, "bad_key": 1}}
    path = tmp_path / "bad.json"
    path.write_text(json.dumps(bad), encoding="utf-8")
    with pytest.raises(ValueError, match="unknown keys"):
        load_design_input(path)


# --- CLI: BasicMQWDesign.py --design-input ---


def test_basic_design_cli_with_design_input(tmp_path: Path) -> None:
    json_path = tmp_path / "design.json"
    lsf_path = tmp_path / "design.lsf"

    _run_script(
        "BasicMQWDesign.py",
        "--design-input", str(FIXTURES / "ingaasp_oband_input.json"),
        "--json", str(json_path),
        "--lsf", str(lsf_path),
    )

    data = json.loads(json_path.read_text(encoding="utf-8"))
    assert data["family"] == "ingaasp"
    assert json_path.stat().st_size > 0


def test_basic_design_cli_cband_input(tmp_path: Path) -> None:
    json_path = tmp_path / "design.json"
    lsf_path = tmp_path / "design.lsf"

    _run_script(
        "BasicMQWDesign.py",
        "--design-input", str(FIXTURES / "ingaasp_cband_input.json"),
        "--json", str(json_path),
        "--lsf", str(lsf_path),
    )

    data = json.loads(json_path.read_text(encoding="utf-8"))
    wl_nm = data["transition"]["lambda_transition_hh_um"] * 1000.0
    assert 1450.0 < wl_nm < 1600.0


# --- CLI: MQWGainDesign.py --design-input ---


def test_gain_cli_with_design_input(tmp_path: Path) -> None:
    json_path = tmp_path / "gain.json"
    csv_path = tmp_path / "gain.csv"
    png_path = tmp_path / "gain.png"

    _run_script(
        "MQWGainDesign.py",
        "--design-input", str(FIXTURES / "ingaasp_oband_input.json"),
        "--dz-nm", "0.2",
        "--kt-points", "5",
        "--energy-points", "60",
        "--electron-states", "1",
        "--hole-states", "2",
        "--out-json", str(json_path),
        "--out-csv", str(csv_path),
        "--plot", str(png_path),
    )

    data = json.loads(json_path.read_text(encoding="utf-8"))
    assert data["design"]["family"] == "ingaasp"
    assert data["gain"]["transition_count"] > 0


# --- CLI: MQWGainSweep.py --design-input ---


def test_sweep_cli_with_design_input(tmp_path: Path) -> None:
    json_path = tmp_path / "sweep.json"
    csv_path = tmp_path / "sweep.csv"
    plot_path = tmp_path / "sweep.png"
    spectra_csv = tmp_path / "spectra.csv"
    spectra_plot = tmp_path / "spectra.png"

    _run_script(
        "MQWGainSweep.py",
        "--design-input", str(FIXTURES / "ingaasp_oband_input.json"),
        "--sweep", "carrier-density",
        "--values", "1.5e18,2.0e18",
        "--dz-nm", "0.2",
        "--kt-points", "3",
        "--energy-points", "40",
        "--electron-states", "1",
        "--hole-states", "2",
        "--out-json", str(json_path),
        "--out-csv", str(csv_path),
        "--plot", str(plot_path),
        "--spectra-csv", str(spectra_csv),
        "--spectra-plot", str(spectra_plot),
    )

    data = json.loads(json_path.read_text(encoding="utf-8"))
    assert len(data["rows"]) == 2


# --- Suzuki 2018 C-band MQW (paper validation) ---


def test_suzuki2018_design_strain_matches_paper() -> None:
    """Validate strain values against Suzuki JJAP 57 094101 (2018) Table I."""
    kwargs = load_design_input(FIXTURES / "suzuki2018_cband_input.json")
    design = design_default(**kwargs)
    well_eps = design["transition"]["well_strain"]["eps_parallel"]
    barrier_eps = design["transition"]["barrier_strain"]["eps_parallel"]
    # Paper: CS 1.07% (compressive → negative), TS 0.15% (tensile → positive)
    assert -0.012 < well_eps < -0.010
    assert 0.001 < barrier_eps < 0.002


def test_suzuki2018_transition_in_cband() -> None:
    kwargs = load_design_input(FIXTURES / "suzuki2018_cband_input.json")
    design = design_default(**kwargs)
    wl_nm = design["transition"]["lambda_transition_hh_um"] * 1000.0
    assert 1500.0 < wl_nm < 1650.0


def test_suzuki2018_gain_cli(tmp_path: Path) -> None:
    json_path = tmp_path / "gain.json"
    csv_path = tmp_path / "gain.csv"
    png_path = tmp_path / "gain.png"

    _run_script(
        "MQWGainDesign.py",
        "--design-input", str(FIXTURES / "suzuki2018_cband_input.json"),
        "--energy-min-eV", "0.70",
        "--energy-max-eV", "0.90",
        "--energy-points", "80",
        "--dz-nm", "0.2",
        "--kt-points", "5",
        "--electron-states", "1",
        "--hole-states", "2",
        "--out-json", str(json_path),
        "--out-csv", str(csv_path),
        "--plot", str(png_path),
    )

    data = json.loads(json_path.read_text(encoding="utf-8"))
    peak_nm = data["gain"]["peak_te_wavelength_nm"]
    # Coarse grid shifts peak; just verify it's well beyond O-band
    assert 1350.0 < peak_nm < 1650.0
