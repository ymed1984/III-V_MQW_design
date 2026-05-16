import csv
import json
import subprocess
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]


def _run_script(script: str, *args: str) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [sys.executable, "-B", str(REPO_ROOT / "src" / script), *args],
        check=True,
        cwd=REPO_ROOT,
        text=True,
        capture_output=True,
    )


def test_basic_design_cli_writes_outputs(tmp_path: Path) -> None:
    json_path = tmp_path / "design.json"
    lsf_path = tmp_path / "design.lsf"

    _run_script(
        "BasicMQWDesign.py",
        "--json",
        str(json_path),
        "--lsf",
        str(lsf_path),
    )

    data = json.loads(json_path.read_text(encoding="utf-8"))
    assert data["family"] == "ingaasp"
    assert json_path.stat().st_size > 0
    assert lsf_path.stat().st_size > 0


def test_critical_film_stress_cli_writes_json(tmp_path: Path) -> None:
    json_path = tmp_path / "stress.json"

    _run_script(
        "CriticalFilmStress.py",
        "--family",
        "ingaasp",
        "--strain",
        "-0.006",
        "--as-frac",
        "0.567",
        "--thickness-nm",
        "7.0",
        "--json",
        str(json_path),
    )

    data = json.loads(json_path.read_text(encoding="utf-8"))
    assert data["family"] == "ingaasp"
    assert data["film_stress_GPa"] < 0.0


def test_gain_cli_writes_outputs(tmp_path: Path) -> None:
    json_path = tmp_path / "gain.json"
    csv_path = tmp_path / "gain.csv"
    png_path = tmp_path / "gain.png"
    band_path = tmp_path / "band.png"
    wavefunction_path = tmp_path / "wavefunctions.png"
    dispersion_path = tmp_path / "dispersion.png"

    _run_script(
        "MQWGainDesign.py",
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
        "--band-plot",
        str(band_path),
        "--wavefunction-plot",
        str(wavefunction_path),
        "--dispersion-plot",
        str(dispersion_path),
    )

    data = json.loads(json_path.read_text(encoding="utf-8"))
    assert data["design"]["family"] == "ingaasp"
    assert data["gain"]["transition_count"] > 0
    with csv_path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    assert len(rows) == 60
    assert png_path.stat().st_size > 0
    assert band_path.stat().st_size > 0
    assert wavefunction_path.stat().st_size > 0
    assert dispersion_path.stat().st_size > 0


def test_gain_sweep_cli_accepts_calibration(tmp_path: Path) -> None:
    json_path = tmp_path / "sweep.json"
    csv_path = tmp_path / "sweep.csv"
    plot_path = tmp_path / "sweep.png"
    spectra_csv_path = tmp_path / "spectra.csv"
    spectra_plot_path = tmp_path / "spectra.png"
    reference_csv_path = tmp_path / "reference.csv"
    comparison_csv_path = tmp_path / "comparison.csv"
    with reference_csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["energy_eV", "wavelength_nm", "gain_TE_cm-1", "gain_TM_cm-1"],
        )
        writer.writeheader()
        writer.writerows(
            [
                {
                    "energy_eV": 0.95,
                    "wavelength_nm": 1305.0,
                    "gain_TE_cm-1": 1.0,
                    "gain_TM_cm-1": 0.5,
                },
                {
                    "energy_eV": 1.00,
                    "wavelength_nm": 1240.0,
                    "gain_TE_cm-1": 2.0,
                    "gain_TM_cm-1": 1.0,
                },
            ]
        )

    _run_script(
        "MQWGainSweep.py",
        "--calibration",
        str(REPO_ROOT / "calibrations" / "ingaasp_oband_example.json"),
        "--sweep",
        "qc",
        "--values",
        "0.35,0.40",
        "--dz-nm",
        "0.2",
        "--kt-points",
        "3",
        "--energy-points",
        "40",
        "--electron-states",
        "1",
        "--hole-states",
        "1",
        "--out-json",
        str(json_path),
        "--out-csv",
        str(csv_path),
        "--plot",
        str(plot_path),
        "--spectra-csv",
        str(spectra_csv_path),
        "--spectra-plot",
        str(spectra_plot_path),
        "--reference-csv",
        str(reference_csv_path),
        "--comparison-polarization",
        "TE",
        "--comparison-csv",
        str(comparison_csv_path),
    )

    data = json.loads(json_path.read_text(encoding="utf-8"))
    assert data["calibration"]["name"] == "ingaasp_oband_example"
    assert data["sweep"]["name"] == "qc"
    assert data["sweep"]["overrides_calibration_field"] == "qc"
    assert [row["qc"] for row in data["rows"]] == [0.35, 0.4]
    assert data["reference_comparison"]["best"]["score_rmse_gain_cm"] >= 0.0
    assert csv_path.stat().st_size > 0
    assert plot_path.stat().st_size > 0
    assert spectra_csv_path.stat().st_size > 0
    assert spectra_plot_path.stat().st_size > 0
    assert comparison_csv_path.stat().st_size > 0


def test_gain_sweep_broadening_overrides_calibration(tmp_path: Path) -> None:
    json_path = tmp_path / "sweep.json"
    csv_path = tmp_path / "sweep.csv"
    plot_path = tmp_path / "sweep.png"
    spectra_csv_path = tmp_path / "spectra.csv"
    spectra_plot_path = tmp_path / "spectra.png"

    _run_script(
        "MQWGainSweep.py",
        "--calibration",
        str(REPO_ROOT / "calibrations" / "ingaasp_oband_example.json"),
        "--sweep",
        "broadening-eV",
        "--values",
        "0.02,0.04",
        "--dz-nm",
        "0.2",
        "--kt-points",
        "3",
        "--energy-points",
        "40",
        "--electron-states",
        "1",
        "--hole-states",
        "1",
        "--out-json",
        str(json_path),
        "--out-csv",
        str(csv_path),
        "--plot",
        str(plot_path),
        "--spectra-csv",
        str(spectra_csv_path),
        "--spectra-plot",
        str(spectra_plot_path),
    )

    data = json.loads(json_path.read_text(encoding="utf-8"))
    assert data["sweep"]["name"] == "broadening-eV"
    assert data["sweep"]["overrides_calibration_field"] == "broadening_eV"
    assert [row["broadening_eV"] for row in data["rows"]] == [0.02, 0.04]
