#!/usr/bin/env python3
"""Calibration file handling for compact MQW gain calculations."""

from __future__ import annotations

import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

DEFAULT_BROADENING_EV = 0.030
DEFAULT_LINE_SHAPE = "lorentzian"
DEFAULT_GAIN_SCALE_CM = 2400.0
LINE_SHAPES = {"lorentzian", "gaussian"}


@dataclass(frozen=True)
class CalibrationBand:
    qc: float | None = None
    Eg_offset_well_eV: float = 0.0
    Eg_offset_barrier_eV: float = 0.0


@dataclass(frozen=True)
class CalibrationGain:
    broadening_eV: float | None = None
    line_shape: str | None = None
    gain_scale_cm: float | None = None


@dataclass(frozen=True)
class Calibration:
    name: str | None
    description: str | None
    path: str | None
    band: CalibrationBand
    gain: CalibrationGain
    raw: dict[str, Any]


@dataclass(frozen=True)
class ResolvedCalibration:
    q_c: float | None
    Eg_offset_well_eV: float
    Eg_offset_barrier_eV: float
    broadening_eV: float
    line_shape: str
    gain_scale_cm: float
    overrides: dict[str, bool]


def _optional_float(value: Any, field_name: str) -> float | None:
    if value is None:
        return None
    if isinstance(value, bool):
        raise ValueError(f"{field_name} must be numeric")
    result = float(value)
    if not math.isfinite(result):
        raise ValueError(f"{field_name} must be finite")
    return result


def _float_with_default(value: Any, default: float, field_name: str) -> float:
    result = _optional_float(value, field_name)
    return default if result is None else result


def _validate_qc(qc: float | None) -> None:
    if qc is not None and not (0.0 <= qc <= 1.0):
        raise ValueError("band.qc must be in [0, 1]")


def _validate_gain(gain: CalibrationGain) -> None:
    if gain.broadening_eV is not None and gain.broadening_eV <= 0.0:
        raise ValueError("gain.broadening_eV must be positive")
    if gain.gain_scale_cm is not None and gain.gain_scale_cm <= 0.0:
        raise ValueError("gain.gain_scale_cm must be positive")
    if gain.line_shape is not None and gain.line_shape not in LINE_SHAPES:
        raise ValueError("gain.line_shape must be 'lorentzian' or 'gaussian'")


def load_calibration(path: str | Path | None) -> Calibration:
    if path is None:
        return Calibration(
            name=None,
            description=None,
            path=None,
            band=CalibrationBand(),
            gain=CalibrationGain(),
            raw={},
        )

    resolved_path = Path(path)
    raw = json.loads(resolved_path.read_text(encoding="utf-8"))
    if not isinstance(raw, dict):
        raise ValueError("calibration file root must be a JSON object")

    band_raw = raw.get("band", {})
    gain_raw = raw.get("gain", {})
    if not isinstance(band_raw, dict):
        raise ValueError("calibration.band must be an object")
    if not isinstance(gain_raw, dict):
        raise ValueError("calibration.gain must be an object")

    band = CalibrationBand(
        qc=_optional_float(band_raw.get("qc"), "band.qc"),
        Eg_offset_well_eV=_float_with_default(
            band_raw.get("Eg_offset_well_eV"),
            0.0,
            "band.Eg_offset_well_eV",
        ),
        Eg_offset_barrier_eV=_float_with_default(
            band_raw.get("Eg_offset_barrier_eV"),
            0.0,
            "band.Eg_offset_barrier_eV",
        ),
    )
    gain = CalibrationGain(
        broadening_eV=_optional_float(gain_raw.get("broadening_eV"), "gain.broadening_eV"),
        line_shape=gain_raw.get("line_shape"),
        gain_scale_cm=_optional_float(gain_raw.get("gain_scale_cm"), "gain.gain_scale_cm"),
    )
    if gain.line_shape is not None and not isinstance(gain.line_shape, str):
        raise ValueError("gain.line_shape must be a string")
    _validate_qc(band.qc)
    _validate_gain(gain)

    return Calibration(
        name=raw.get("name"),
        description=raw.get("description"),
        path=str(resolved_path),
        band=band,
        gain=gain,
        raw=raw,
    )


def _resolve(
    cli_value: Any,
    calibration_value: Any,
    default_value: Any,
) -> tuple[Any, bool]:
    if cli_value is not None:
        return cli_value, True
    if calibration_value is not None:
        return calibration_value, False
    return default_value, False


def resolve_calibration(args: Any, calibration: Calibration) -> ResolvedCalibration:
    q_c, qc_override = _resolve(args.qc, calibration.band.qc, None)
    eg_well, eg_well_override = _resolve(
        args.eg_offset_well_eV,
        calibration.band.Eg_offset_well_eV,
        0.0,
    )
    eg_barrier, eg_barrier_override = _resolve(
        args.eg_offset_barrier_eV,
        calibration.band.Eg_offset_barrier_eV,
        0.0,
    )
    broadening, broadening_override = _resolve(
        args.broadening_eV,
        calibration.gain.broadening_eV,
        DEFAULT_BROADENING_EV,
    )
    line_shape, line_shape_override = _resolve(
        args.line_shape,
        calibration.gain.line_shape,
        DEFAULT_LINE_SHAPE,
    )
    gain_scale, gain_scale_override = _resolve(
        args.gain_scale_cm,
        calibration.gain.gain_scale_cm,
        DEFAULT_GAIN_SCALE_CM,
    )

    resolved = ResolvedCalibration(
        q_c=None if q_c is None else float(q_c),
        Eg_offset_well_eV=float(eg_well),
        Eg_offset_barrier_eV=float(eg_barrier),
        broadening_eV=float(broadening),
        line_shape=str(line_shape),
        gain_scale_cm=float(gain_scale),
        overrides={
            "qc": qc_override,
            "Eg_offset_well_eV": eg_well_override,
            "Eg_offset_barrier_eV": eg_barrier_override,
            "broadening_eV": broadening_override,
            "line_shape": line_shape_override,
            "gain_scale_cm": gain_scale_override,
        },
    )
    _validate_qc(resolved.q_c)
    _validate_gain(
        CalibrationGain(
            broadening_eV=resolved.broadening_eV,
            line_shape=resolved.line_shape,
            gain_scale_cm=resolved.gain_scale_cm,
        )
    )
    return resolved


def calibration_summary(
    calibration: Calibration,
    resolved: ResolvedCalibration,
    applied_qc: float,
) -> dict[str, Any]:
    applied = asdict(resolved)
    applied.pop("overrides")
    applied["qc"] = applied.pop("q_c")
    if applied["qc"] is None:
        applied["qc"] = applied_qc
    return {
        "path": calibration.path,
        "name": calibration.name,
        "description": calibration.description,
        "applied": applied,
        "overrides": resolved.overrides,
    }
