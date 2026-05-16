import csv
import json
import subprocess
import sys
from pathlib import Path

import pytest

from FitCalibration import build_arg_parser, fit_calibration, resolve_fit_targets


REPO_ROOT = Path(__file__).resolve().parents[1]


def _base_args(tmp_path: Path):
    return build_arg_parser().parse_args(
        [
            "--calibration-in",
            str(REPO_ROOT / "calibrations" / "ingaasp_oband_example.json"),
            "--target-peak-wavelength-nm",
            "1260",
            "--target-te-peak-gain-cm",
            "30",
            "--target-fwhm-meV",
            "80",
            "--dz-nm",
            "0.3",
            "--kt-points",
            "3",
            "--energy-points",
            "50",
            "--electron-states",
            "1",
            "--hole-states",
            "1",
            "--eg-offset-min-eV",
            "-0.04",
            "--eg-offset-max-eV",
            "0.04",
            "--fit-maxiter",
            "5",
            "--out",
            str(tmp_path / "fit.json"),
        ]
    )


def test_fit_calibration_returns_positive_gain_scale(tmp_path: Path) -> None:
    output, state = fit_calibration(_base_args(tmp_path))

    assert output["band"]["Eg_offset_well_eV"] == state.Eg_offset_well_eV
    assert output["gain"]["gain_scale_cm"] > 0.0
    assert output["fit_result"]["metrics"]["peak_gain_cm"] > 0.0


def test_fit_targets_can_be_derived_from_reference_csv(tmp_path: Path) -> None:
    reference_path = tmp_path / "reference.csv"
    rows = [
        {"energy_eV": 1.0, "wavelength_nm": 1239.8, "gain_TE_cm-1": 0.0, "gain_TM_cm-1": 0.0},
        {"energy_eV": 1.1, "wavelength_nm": 1127.1, "gain_TE_cm-1": 20.0, "gain_TM_cm-1": 8.0},
        {"energy_eV": 1.2, "wavelength_nm": 1033.2, "gain_TE_cm-1": 0.0, "gain_TM_cm-1": 0.0},
    ]
    with reference_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    args = build_arg_parser().parse_args(
        [
            "--reference-csv",
            str(reference_path),
            "--reference-polarization",
            "TM",
            "--out",
            str(tmp_path / "fit.json"),
        ]
    )

    targets = resolve_fit_targets(args)

    assert targets.polarization == "TM"
    assert targets.peak_wavelength_nm == pytest.approx(1127.1)
    assert targets.peak_gain_cm == pytest.approx(8.0)
    assert targets.fwhm_meV == pytest.approx(100.0)


def test_fit_calibration_cli_writes_json(tmp_path: Path) -> None:
    out_path = tmp_path / "fit.json"

    subprocess.run(
        [
            sys.executable,
            "-B",
            str(REPO_ROOT / "src" / "FitCalibration.py"),
            "--calibration-in",
            str(REPO_ROOT / "calibrations" / "ingaasp_oband_example.json"),
            "--target-peak-wavelength-nm",
            "1260",
            "--target-te-peak-gain-cm",
            "30",
            "--target-fwhm-meV",
            "80",
            "--dz-nm",
            "0.3",
            "--kt-points",
            "3",
            "--energy-points",
            "50",
            "--electron-states",
            "1",
            "--hole-states",
            "1",
            "--eg-offset-min-eV",
            "-0.04",
            "--eg-offset-max-eV",
            "0.04",
            "--fit-maxiter",
            "5",
            "--out",
            str(out_path),
        ],
        check=True,
        cwd=REPO_ROOT,
        text=True,
        capture_output=True,
    )

    data = json.loads(out_path.read_text(encoding="utf-8"))
    assert data["gain"]["gain_scale_cm"] > 0.0
    assert "fit_result" in data
