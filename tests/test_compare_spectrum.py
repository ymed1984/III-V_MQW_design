import csv
import json
import subprocess
import sys
from pathlib import Path

import pytest

from CompareSpectrum import compare_spectra, filter_wavelength_range, read_spectrum_csv


REPO_ROOT = Path(__file__).resolve().parents[1]


def _rows(shift_nm: float = 0.0, scale: float = 1.0) -> list[dict[str, float]]:
    return [
        {
            "energy_eV": 0.9,
            "wavelength_nm": 1377.6 + shift_nm,
            "gain_TE_cm-1": 0.0 * scale,
            "gain_TM_cm-1": 0.0 * scale,
        },
        {
            "energy_eV": 1.0,
            "wavelength_nm": 1239.8 + shift_nm,
            "gain_TE_cm-1": 10.0 * scale,
            "gain_TM_cm-1": 4.0 * scale,
        },
        {
            "energy_eV": 1.1,
            "wavelength_nm": 1127.1 + shift_nm,
            "gain_TE_cm-1": 20.0 * scale,
            "gain_TM_cm-1": 8.0 * scale,
        },
        {
            "energy_eV": 1.2,
            "wavelength_nm": 1033.2 + shift_nm,
            "gain_TE_cm-1": 10.0 * scale,
            "gain_TM_cm-1": 4.0 * scale,
        },
        {
            "energy_eV": 1.3,
            "wavelength_nm": 953.7 + shift_nm,
            "gain_TE_cm-1": 0.0 * scale,
            "gain_TM_cm-1": 0.0 * scale,
        },
    ]


def _write_csv(path: Path, rows: list[dict[str, float]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def test_compare_spectra_reports_peak_deltas_and_rmse() -> None:
    result = compare_spectra(_rows(), _rows(shift_nm=2.0, scale=0.5), ["TE"])

    te = result["TE"]
    assert te["delta"]["peak_wavelength_nm"] == pytest.approx(-2.0)
    assert te["delta"]["peak_gain_cm"] == pytest.approx(10.0)
    assert te["rmse_gain_cm"] > 0.0


def test_read_spectrum_csv_validates_required_columns(tmp_path: Path) -> None:
    path = tmp_path / "bad.csv"
    path.write_text("energy_eV,wavelength_nm,gain_TE_cm-1\n1.0,1239.8,1.0\n")

    with pytest.raises(ValueError, match="missing required columns"):
        read_spectrum_csv(path)


def test_filter_wavelength_range_rejects_empty_range() -> None:
    with pytest.raises(ValueError, match="removed all rows"):
        filter_wavelength_range(_rows(), 1400.0, 1500.0)


def test_compare_spectrum_cli_writes_json_and_plot(tmp_path: Path) -> None:
    predicted = tmp_path / "predicted.csv"
    reference = tmp_path / "reference.csv"
    out_json = tmp_path / "comparison.json"
    out_plot = tmp_path / "comparison.png"
    _write_csv(predicted, _rows())
    _write_csv(reference, _rows(scale=0.8))

    completed = subprocess.run(
        [
            sys.executable,
            "-B",
            str(REPO_ROOT / "src" / "CompareSpectrum.py"),
            "--predicted",
            str(predicted),
            "--reference",
            str(reference),
            "--polarization",
            "TE",
            "--out-json",
            str(out_json),
            "--out-plot",
            str(out_plot),
        ],
        check=True,
        cwd=REPO_ROOT,
        text=True,
        capture_output=True,
    )

    data = json.loads(out_json.read_text(encoding="utf-8"))
    assert "TE peak delta" in completed.stdout
    assert data["comparisons"]["TE"]["delta"]["peak_gain_cm"] == pytest.approx(4.0)
    assert out_plot.stat().st_size > 0
