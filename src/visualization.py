#!/usr/bin/env python3
"""Plot helpers for compact MQW k.p gain calculations."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from kp_solver import KpSubbands, MQWProfile


def _save(fig, path: str | Path) -> Path:
    resolved = Path(path)
    resolved.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(resolved, dpi=180)
    plt.close(fig)
    return resolved


def plot_gain_spectrum(
    rows: list[dict[str, float]],
    path: str | Path,
    title: str | None = None,
) -> Path:
    wavelength = np.array([row["wavelength_nm"] for row in rows], dtype=float)
    gain_te = np.array([row["gain_TE_cm-1"] for row in rows], dtype=float)
    gain_tm = np.array([row["gain_TM_cm-1"] for row in rows], dtype=float)
    order = np.argsort(wavelength)
    te_peak = int(np.argmax(gain_te))
    tm_peak = int(np.argmax(gain_tm))

    fig, ax = plt.subplots(figsize=(7.0, 4.2), constrained_layout=True)
    ax.plot(wavelength[order], gain_te[order], label="TE")
    ax.plot(wavelength[order], gain_tm[order], label="TM")
    ax.scatter(
        [wavelength[te_peak], wavelength[tm_peak]],
        [gain_te[te_peak], gain_tm[tm_peak]],
        marker="o",
        s=28,
        zorder=3,
        color=["tab:blue", "tab:orange"],
    )
    ax.axhline(0.0, color="0.35", linewidth=0.8)
    ax.set_xlabel("Wavelength [nm]")
    ax.set_ylabel("Material gain [cm$^{-1}$]")
    if title:
        ax.set_title(title)
    ax.legend()
    ax.grid(True, alpha=0.25)
    return _save(fig, path)


def plot_band_diagram(
    profile: MQWProfile,
    subbands: KpSubbands,
    path: str | Path,
    title: str | None = None,
) -> Path:
    z = profile.z_nm
    ec = profile.V_e_eV
    hh_edge = -profile.V_hh_eV
    lh_edge = -profile.V_lh_eV

    fig, ax = plt.subplots(figsize=(8.0, 4.6), constrained_layout=True)
    ax.plot(z, ec, label="Ec profile", color="tab:blue")
    ax.plot(z, hh_edge, label="HH edge (-Vhh)", color="tab:red")
    ax.plot(z, lh_edge, label="LH edge (-Vlh)", color="tab:purple")
    if subbands.electron and subbands.electron[0]:
        ax.axhline(subbands.electron[0][0].energy_eV, color="tab:blue", linestyle="--", label="e1")
    if subbands.valence and subbands.valence[0]:
        ax.axhline(-subbands.valence[0][0].energy_eV, color="tab:red", linestyle="--", label="h1")
    ax.set_xlabel("z [nm]")
    ax.set_ylabel("Relative energy [eV]")
    if title:
        ax.set_title(title)
    ax.grid(True, alpha=0.25)
    ax.legend(fontsize=8)
    return _save(fig, path)


def plot_wavefunctions(
    profile: MQWProfile,
    subbands: KpSubbands,
    path: str | Path,
    title: str | None = None,
) -> Path:
    z = profile.z_nm
    fig, axes = plt.subplots(2, 1, figsize=(8.0, 5.8), sharex=True, constrained_layout=True)

    axes[0].plot(z, profile.V_e_eV, color="0.45", linewidth=1.0, label="Ec profile")
    if subbands.electron and subbands.electron[0]:
        for index, state in enumerate(subbands.electron[0][:3], start=1):
            density = np.abs(state.psi) ** 2
            if np.max(density) > 0:
                density = density / np.max(density) * 0.05
            axes[0].plot(z, state.energy_eV + density, label=f"e{index}")
            axes[0].axhline(state.energy_eV, color="0.7", linewidth=0.7)
    axes[0].set_ylabel("Electron energy [eV]")
    axes[0].legend(fontsize=8)
    axes[0].grid(True, alpha=0.25)

    axes[1].plot(z, -profile.V_hh_eV, color="0.45", linewidth=1.0, label="HH edge")
    if subbands.valence and subbands.valence[0]:
        for index, state in enumerate(subbands.valence[0][:3], start=1):
            density = np.abs(state.hh) ** 2 + np.abs(state.lh) ** 2
            if np.max(density) > 0:
                density = density / np.max(density) * 0.05
            energy = -state.energy_eV
            axes[1].plot(z, energy - density, label=f"h{index}")
            axes[1].axhline(energy, color="0.7", linewidth=0.7)
    axes[1].set_xlabel("z [nm]")
    axes[1].set_ylabel("Hole energy [eV]")
    axes[1].legend(fontsize=8)
    axes[1].grid(True, alpha=0.25)
    if title:
        fig.suptitle(title)
    return _save(fig, path)


def plot_subband_dispersion(
    subbands: KpSubbands,
    path: str | Path,
    title: str | None = None,
) -> Path:
    kt = subbands.kt_nm
    fig, (ax_e, ax_h) = plt.subplots(
        2,
        1,
        figsize=(7.0, 6.0),
        sharex=True,
        constrained_layout=True,
    )
    n_e = min(len(row) for row in subbands.electron) if subbands.electron else 0
    n_h = min(len(row) for row in subbands.valence) if subbands.valence else 0
    for index in range(n_e):
        ax_e.plot(kt, [row[index].energy_eV for row in subbands.electron], marker="o", label=f"e{index + 1}")
    for index in range(n_h):
        ax_h.plot(kt, [row[index].energy_eV for row in subbands.valence], marker="o", label=f"h{index + 1}")
    ax_e.set_ylabel("Electron confinement [eV]")
    ax_h.set_ylabel("Hole confinement [eV]")
    ax_h.set_xlabel("$k_t$ [nm$^{-1}$]")
    ax_e.grid(True, alpha=0.25)
    ax_h.grid(True, alpha=0.25)
    ax_e.legend(fontsize=8)
    ax_h.legend(fontsize=8)
    if title:
        fig.suptitle(title)
    return _save(fig, path)


def plot_sweep_summary(
    rows: list[dict[str, float]],
    sweep_label: str,
    path: str | Path,
    title: str | None = None,
) -> Path:
    x = np.array([row["sweep_value"] for row in rows], dtype=float)
    te_wavelength = np.array([row["peak_TE_wavelength_nm"] for row in rows], dtype=float)
    tm_wavelength = np.array([row["peak_TM_wavelength_nm"] for row in rows], dtype=float)
    te_gain = np.array([row["peak_TE_gain_cm-1"] for row in rows], dtype=float)
    tm_gain = np.array([row["peak_TM_gain_cm-1"] for row in rows], dtype=float)
    order = np.argsort(x)

    fig, (ax_wave, ax_gain) = plt.subplots(
        2,
        1,
        figsize=(7.0, 6.4),
        sharex=True,
        constrained_layout=True,
    )
    ax_wave.plot(x[order], te_wavelength[order], marker="o", label="TE peak")
    ax_wave.plot(x[order], tm_wavelength[order], marker="s", label="TM peak")
    ax_wave.set_ylabel("Peak wavelength [nm]")
    ax_wave.grid(True, alpha=0.25)
    ax_wave.legend()

    ax_gain.plot(x[order], te_gain[order], marker="o", label="TE peak")
    ax_gain.plot(x[order], tm_gain[order], marker="s", label="TM peak")
    ax_gain.axhline(0.0, color="0.35", linewidth=0.8)
    ax_gain.set_xlabel(sweep_label)
    ax_gain.set_ylabel("Peak material gain [cm$^{-1}$]")
    ax_gain.grid(True, alpha=0.25)
    ax_gain.legend()
    if title:
        fig.suptitle(title)
    return _save(fig, path)
