"""CSV helpers for wavelength-domain MQW gain spectra."""

from __future__ import annotations

import csv
from pathlib import Path

REQUIRED_SPECTRUM_COLUMNS = {
    "energy_eV",
    "wavelength_nm",
    "gain_TE_cm-1",
    "gain_TM_cm-1",
}


def read_spectrum_csv(path: Path) -> list[dict[str, float]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        rows = [{key: float(value) for key, value in row.items()} for row in reader]
    if not rows:
        raise ValueError(f"{path} contains no rows")
    missing = REQUIRED_SPECTRUM_COLUMNS.difference(rows[0])
    if missing:
        raise ValueError(f"{path} is missing required columns: {sorted(missing)}")
    return rows


def filter_wavelength_range(
    rows: list[dict[str, float]],
    min_nm: float | None,
    max_nm: float | None,
) -> list[dict[str, float]]:
    result = []
    for row in rows:
        wavelength = row["wavelength_nm"]
        if min_nm is not None and wavelength < min_nm:
            continue
        if max_nm is not None and wavelength > max_nm:
            continue
        result.append(row)
    if not result:
        raise ValueError("wavelength filter removed all rows")
    return result
