"""Tests for load_design_json and --design-json CLI integration."""

import csv
import json
import subprocess
import sys
from pathlib import Path

import pytest

from BasicMQWDesign import design_default, load_design_json, write_design_json

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


# --- Unit tests for load_design_json ---


def test_round_trip_write_load(tmp_path: Path) -> None:
    design = design_default(family="ingaasp")
    json_path = tmp_path / "design.json"
    write_design_json(design, json_path)
    loaded = load_design_json(json_path)
    assert loaded["family"] == "ingaasp"
    assert loaded["wells"] == design["wells"]
    assert loaded["well_nm"] == design["well_nm"]
    assert loaded["barrier_nm"] == design["barrier_nm"]
    assert loaded["qc"] == design["qc"]


def test_load_fixture_oband() -> None:
    design = load_design_json(FIXTURES / "ingaasp_oband_design.json")
    assert design["family"] == "ingaasp"
    wavelength_nm = design["transition"]["lambda_transition_hh_um"] * 1000.0
    assert 1250.0 < wavelength_nm < 1370.0


def test_load_fixture_cband() -> None:
    design = load_design_json(FIXTURES / "ingaasp_cband_design.json")
    assert design["family"] == "ingaasp"
    wavelength_nm = design["transition"]["lambda_transition_hh_um"] * 1000.0
    assert 1450.0 < wavelength_nm < 1600.0


def test_load_rejects_missing_keys(tmp_path: Path) -> None:
    bad = {"family": "ingaasp"}
    path = tmp_path / "bad.json"
    path.write_text(json.dumps(bad), encoding="utf-8")
    with pytest.raises(ValueError, match="missing required keys"):
        load_design_json(path)


def test_load_rejects_non_object(tmp_path: Path) -> None:
    path = tmp_path / "bad.json"
    path.write_text("[1,2,3]", encoding="utf-8")
    with pytest.raises(ValueError, match="JSON object"):
        load_design_json(path)


# --- CLI integration: MQWGainDesign.py --design-json ---


def test_gain_cli_with_design_json(tmp_path: Path) -> None:
    json_path = tmp_path / "gain.json"
    csv_path = tmp_path / "gain.csv"
    png_path = tmp_path / "gain.png"

    _run_script(
        "MQWGainDesign.py",
        "--design-json",
        str(FIXTURES / "ingaasp_oband_design.json"),
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
    assert png_path.stat().st_size > 0


def test_gain_cli_cband_design_json(tmp_path: Path) -> None:
    json_path = tmp_path / "gain.json"
    csv_path = tmp_path / "gain.csv"
    png_path = tmp_path / "gain.png"

    _run_script(
        "MQWGainDesign.py",
        "--design-json",
        str(FIXTURES / "ingaasp_cband_design.json"),
        "--dz-nm", "0.2",
        "--kt-points", "5",
        "--energy-points", "80",
        "--energy-min-eV", "0.75",
        "--energy-max-eV", "0.95",
        "--electron-states", "1",
        "--hole-states", "2",
        "--out-json", str(json_path),
        "--out-csv", str(csv_path),
        "--plot", str(png_path),
    )

    data = json.loads(json_path.read_text(encoding="utf-8"))
    assert data["design"]["family"] == "ingaasp"
    peak_nm = data["gain"]["peak_te_wavelength_nm"]
    # Coarse grid; just verify peak is redshifted vs O-band (~1310 nm)
    assert 1300.0 < peak_nm < 1650.0


# --- CLI integration: MQWGainSweep.py --design-json ---


def test_sweep_cli_with_design_json(tmp_path: Path) -> None:
    json_path = tmp_path / "sweep.json"
    csv_path = tmp_path / "sweep.csv"
    plot_path = tmp_path / "sweep.png"
    spectra_csv = tmp_path / "spectra.csv"
    spectra_plot = tmp_path / "spectra.png"

    _run_script(
        "MQWGainSweep.py",
        "--design-json",
        str(FIXTURES / "ingaasp_oband_design.json"),
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
    assert plot_path.stat().st_size > 0
