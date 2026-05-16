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
    )

    data = json.loads(json_path.read_text(encoding="utf-8"))
    assert data["design"]["family"] == "ingaasp"
    assert data["gain"]["transition_count"] > 0
    with csv_path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    assert len(rows) == 60
    assert png_path.stat().st_size > 0
