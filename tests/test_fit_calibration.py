import json
import subprocess
import sys
from pathlib import Path

from FitCalibration import build_arg_parser, fit_calibration


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
