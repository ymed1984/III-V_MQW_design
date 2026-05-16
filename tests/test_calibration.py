import argparse
import json
import subprocess
import sys
from pathlib import Path

import pytest

from BasicMQWDesign import design_default
from calibration import load_calibration, resolve_calibration

REPO_ROOT = Path(__file__).resolve().parents[1]


def _args(**overrides):
    defaults = {
        "qc": None,
        "eg_offset_well_eV": None,
        "eg_offset_barrier_eV": None,
        "broadening_eV": None,
        "line_shape": None,
        "gain_scale_cm": None,
    }
    defaults.update(overrides)
    return argparse.Namespace(**defaults)


def test_load_calibration_json(tmp_path: Path) -> None:
    path = tmp_path / "calibration.json"
    path.write_text(
        json.dumps(
            {
                "name": "unit_cal",
                "band": {
                    "qc": 0.43,
                    "Eg_offset_well_eV": -0.01,
                    "Eg_offset_barrier_eV": 0.02,
                },
                "gain": {
                    "broadening_eV": 0.04,
                    "line_shape": "gaussian",
                    "gain_scale_cm": 3100.0,
                },
            }
        ),
        encoding="utf-8",
    )

    calibration = load_calibration(path)
    resolved = resolve_calibration(_args(), calibration)

    assert calibration.name == "unit_cal"
    assert resolved.q_c == pytest.approx(0.43)
    assert resolved.Eg_offset_well_eV == pytest.approx(-0.01)
    assert resolved.Eg_offset_barrier_eV == pytest.approx(0.02)
    assert resolved.broadening_eV == pytest.approx(0.04)
    assert resolved.line_shape == "gaussian"
    assert resolved.gain_scale_cm == pytest.approx(3100.0)


def test_cli_override_takes_precedence(tmp_path: Path) -> None:
    path = tmp_path / "calibration.json"
    path.write_text(
        json.dumps(
            {
                "band": {"qc": 0.35},
                "gain": {"broadening_eV": 0.04, "gain_scale_cm": 2000.0},
            }
        ),
        encoding="utf-8",
    )

    calibration = load_calibration(path)
    resolved = resolve_calibration(
        _args(qc=0.5, broadening_eV=0.025, gain_scale_cm=3000.0),
        calibration,
    )

    assert resolved.q_c == pytest.approx(0.5)
    assert resolved.broadening_eV == pytest.approx(0.025)
    assert resolved.gain_scale_cm == pytest.approx(3000.0)
    assert resolved.overrides["qc"] is True
    assert resolved.overrides["broadening_eV"] is True
    assert resolved.overrides["gain_scale_cm"] is True


def test_eg_offset_changes_transition_wavelength() -> None:
    nominal = design_default(family="ingaasp")
    shifted = design_default(family="ingaasp", eg_offset_well_eV=-0.02)

    assert shifted["transition"]["lambda_transition_hh_um"] > nominal["transition"][
        "lambda_transition_hh_um"
    ]


@pytest.mark.parametrize(
    "payload",
    [
        {"band": {"qc": -0.1}},
        {"band": {"qc": 1.1}},
        {"gain": {"broadening_eV": 0.0}},
        {"gain": {"gain_scale_cm": -1.0}},
        {"gain": {"line_shape": "bad"}},
    ],
)
def test_invalid_calibration_values_raise(tmp_path: Path, payload: dict) -> None:
    path = tmp_path / "bad.json"
    path.write_text(json.dumps(payload), encoding="utf-8")

    with pytest.raises(ValueError):
        load_calibration(path)


def test_gain_cli_writes_calibration_summary(tmp_path: Path) -> None:
    calibration_path = tmp_path / "calibration.json"
    json_path = tmp_path / "gain.json"
    csv_path = tmp_path / "gain.csv"
    png_path = tmp_path / "gain.png"
    calibration_path.write_text(
        json.dumps(
            {
                "name": "cli_cal",
                "band": {
                    "qc": 0.37,
                    "Eg_offset_well_eV": -0.01,
                    "Eg_offset_barrier_eV": 0.0,
                },
                "gain": {
                    "broadening_eV": 0.04,
                    "line_shape": "gaussian",
                    "gain_scale_cm": 2100.0,
                },
            }
        ),
        encoding="utf-8",
    )

    subprocess.run(
        [
            sys.executable,
            "-B",
            str(REPO_ROOT / "src" / "MQWGainDesign.py"),
            "--calibration",
            str(calibration_path),
            "--gain-scale-cm",
            "3000",
            "--dz-nm",
            "0.2",
            "--kt-points",
            "5",
            "--energy-points",
            "60",
            "--electron-states",
            "1",
            "--hole-states",
            "2",
            "--out-json",
            str(json_path),
            "--out-csv",
            str(csv_path),
            "--plot",
            str(png_path),
        ],
        check=True,
        cwd=REPO_ROOT,
        text=True,
        capture_output=True,
    )

    data = json.loads(json_path.read_text(encoding="utf-8"))
    applied = data["calibration"]["applied"]
    overrides = data["calibration"]["overrides"]

    assert data["calibration"]["name"] == "cli_cal"
    assert applied["qc"] == pytest.approx(0.37)
    assert applied["Eg_offset_well_eV"] == pytest.approx(-0.01)
    assert applied["broadening_eV"] == pytest.approx(0.04)
    assert applied["line_shape"] == "gaussian"
    assert applied["gain_scale_cm"] == pytest.approx(3000.0)
    assert overrides["gain_scale_cm"] is True
