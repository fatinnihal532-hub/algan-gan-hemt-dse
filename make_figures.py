"""Regenerates every figure and CSV in results/ from the model in src/.
Run: python3 make_figures.py   (about 5 s)
"""
import csv
import os

import matplotlib.pyplot as plt
import numpy as np

from src.materials import AlGaN
from src.charge_control import BarrierStack
from src.current_model import HEMTDevice, HEMTGeometry
from src.extraction import extract_metrics
from src.design_space import run_sweep, ron_bv_pareto_front
from src.plotstyle import save_light_dark

REFERENCE = dict(gate_length_m=0.25e-6, gate_drain_spacing_m=2e-6,
                  thickness_m=20e-9, x_al=0.25)


def reference_device():
    stack = BarrierStack(AlGaN(REFERENCE["x_al"]), REFERENCE["thickness_m"])
    geom = HEMTGeometry(gate_length_m=REFERENCE["gate_length_m"],
                         gate_drain_spacing_m=REFERENCE["gate_drain_spacing_m"])
    return HEMTDevice(stack, geom)


def fig_output_iv(theme):
    dev = reference_device()
    vth = dev.threshold_voltage(1.0)
    vds = np.linspace(0, 10, 60)
    fig, ax = plt.subplots(figsize=(5.6, 3.8))
    for i, vov in enumerate([1, 2, 3, 4]):
        vgs = vth + vov
        ids = [dev.drain_current(vgs, v) / (dev.geom.width_m * 1e3) for v in vds]  # A/mm
        ax.plot(vds, ids, color=theme["series"][i % 4], lw=1.8,
                label=f"$V_{{GS}}-V_{{th}}$={vov:.0f} V")
    ax.set_xlabel("$V_{DS}$ (V)")
    ax.set_ylabel("$I_D$ (A/mm)")
    ax.set_title("Reference device output characteristics")
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    return fig


def fig_transfer(theme):
    dev = reference_device()
    vth = dev.threshold_voltage(5.0)
    vgs = np.linspace(vth - 2, vth + 5, 80)
    ids = np.array([dev.drain_current(v, 5.0) / (dev.geom.width_m * 1e3) for v in vgs]) * 1e3  # mA/mm
    gm = np.gradient(ids, vgs)  # mA/mm per V = mS/mm
    fig, ax1 = plt.subplots(figsize=(5.6, 3.8))
    ax1.plot(vgs, ids, color=theme["series"][0], lw=1.8, label="$I_D$")
    ax1.set_xlabel("$V_{GS}$ (V)")
    ax1.set_ylabel("$I_D$ (mA/mm)", color=theme["series"][0])
    ax2 = ax1.twinx()
    ax2.plot(vgs, gm, color=theme["series"][1], lw=1.4, ls="--", label="$g_m$")
    ax2.set_ylabel("$g_m$ (mS/mm)", color=theme["series"][1])
    ax2.grid(False)
    ax1.set_title(f"Transfer characteristic at $V_{{DS}}$=5 V  ($V_{{th}}$={vth:.2f} V)")
    fig.tight_layout()
    return fig


def fig_ron_bv_pareto(theme, rows, front):
    fig, ax = plt.subplots(figsize=(5.6, 4.2))
    ax.scatter([r["ron_sp_ohm_mm"] for r in rows],
               [r["breakdown_V"] for r in rows],
               s=10, color=theme["axis"], alpha=0.55, label="all swept configurations")
    fr = sorted(front, key=lambda r: r["ron_sp_ohm_mm"])
    ax.plot([r["ron_sp_ohm_mm"] for r in fr], [r["breakdown_V"] for r in fr],
            color=theme["series"][1], marker="o", lw=1.8, ms=5,
            label="Pareto front (min $R_{on,sp}$, max BV)")
    ax.set_xlabel(r"$R_{on,sp}$ (Ohm.mm)")
    ax.set_ylabel("Breakdown voltage estimate (V)")
    ax.set_title("On-resistance vs. breakdown-voltage trade-off")
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    return fig


def fig_ion_vs_length_by_alfrac(theme, rows):
    fig, ax = plt.subplots(figsize=(5.6, 3.8))
    al_values = sorted(set(r["al_fraction"] for r in rows))
    lgd_fixed = min(r["gate_drain_spacing_um"] for r in rows)
    d_fixed = sorted(set(r["barrier_thickness_nm"] for r in rows))[len(set(r["barrier_thickness_nm"] for r in rows)) // 2]
    for i, x in enumerate(al_values):
        subset = sorted(
            [r for r in rows if r["al_fraction"] == x
             and r["gate_drain_spacing_um"] == lgd_fixed
             and r["barrier_thickness_nm"] == d_fixed],
            key=lambda r: r["gate_length_nm"])
        if not subset:
            continue
        ax.plot([r["gate_length_nm"] for r in subset],
                [r["ion_A_per_mm"] for r in subset],
                marker="o", ms=4, lw=1.6, color=theme["series"][i % 4],
                label=f"x$_{{Al}}$={x:.2f}")
    ax.set_xlabel("Gate length (nm)")
    ax.set_ylabel("$I_{on}$ (A/mm)")
    ax.set_title(f"On-current vs. gate length ($L_{{GD}}$={lgd_fixed:.0f} um, d={d_fixed:.0f} nm)")
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    return fig


def fig_ss_dibl_vs_length(theme, rows):
    fig, ax1 = plt.subplots(figsize=(5.6, 3.8))
    d_fixed = sorted(set(r["barrier_thickness_nm"] for r in rows))[0]
    x_fixed = sorted(set(r["al_fraction"] for r in rows))[0]
    lgd_fixed = min(r["gate_drain_spacing_um"] for r in rows)
    subset = sorted(
        [r for r in rows if r["barrier_thickness_nm"] == d_fixed
         and r["al_fraction"] == x_fixed
         and r["gate_drain_spacing_um"] == lgd_fixed],
        key=lambda r: r["gate_length_nm"])
    lg = [r["gate_length_nm"] for r in subset]
    ax1.plot(lg, [r["ss_mV_per_dec"] for r in subset], color=theme["series"][0],
              marker="o", ms=4, lw=1.8, label="SS")
    ax1.axhline(60, color=theme["secondary"], lw=1, ls=":", label="60 mV/dec (Boltzmann limit)")
    ax1.set_xlabel("Gate length (nm)")
    ax1.set_ylabel("Subthreshold swing (mV/dec)", color=theme["series"][0])
    ax2 = ax1.twinx()
    ax2.plot(lg, [r["dibl_mV_per_V"] for r in subset], color=theme["series"][1],
              marker="s", ms=4, lw=1.4, ls="--", label="DIBL")
    ax2.set_ylabel("DIBL (mV/V)", color=theme["series"][1])
    ax2.grid(False)
    ax1.set_title(f"Short-channel degradation (d={d_fixed:.0f} nm, x$_{{Al}}$={x_fixed:.2f})")
    fig.tight_layout()
    return fig


def write_csv(path, rows, fieldnames):
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for r in rows:
            w.writerow(r)


def main():
    os.makedirs("results", exist_ok=True)  # absent on a fresh CI checkout
    print("Running design-space sweep...")
    rows = run_sweep(
        gate_lengths_m=[x * 1e-9 for x in (100, 150, 200, 250, 350, 500, 750, 1000)],
        gate_drain_spacings_m=[x * 1e-6 for x in (1, 2, 3, 4, 6, 8)],
        thicknesses_m=[x * 1e-9 for x in (15, 20, 25)],
        al_fractions=[0.20, 0.25, 0.30],
    )
    print(f"  {len(rows)} configurations evaluated")

    front = ron_bv_pareto_front(rows)
    print(f"  {len(front)} Pareto-optimal configurations (min Ron_sp, max BV)")

    fieldnames = list(rows[0].keys())
    write_csv("results/design_space_sweep.csv", rows, fieldnames)
    write_csv("results/pareto_front.csv", front, fieldnames)

    ref = extract_metrics(reference_device())
    with open("results/reference_device_metrics.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["parameter", "value", "unit"])
        for field, unit in [
            ("vth_V", "V"), ("gm_max_S_per_mm", "S/mm"),
            ("ss_mV_per_dec", "mV/dec"), ("dibl_mV_per_V", "mV/V"),
            ("ion_A_per_mm", "A/mm"), ("ioff_A_per_mm", "A/mm"),
            ("ion_ioff_ratio", "-"), ("ron_sp_ohm_mm", "Ohm.mm"),
            ("breakdown_V", "V"), ("baliga_fom_V2_per_ohm_mm", "V^2/(Ohm.mm)"),
        ]:
            w.writerow([field, getattr(ref, field), unit])

    save_light_dark(fig_output_iv, "output_characteristics")
    save_light_dark(fig_transfer, "transfer_characteristic")
    save_light_dark(lambda th: fig_ron_bv_pareto(th, rows, front), "ron_bv_pareto")
    save_light_dark(lambda th: fig_ion_vs_length_by_alfrac(th, rows), "ion_vs_length")
    save_light_dark(lambda th: fig_ss_dibl_vs_length(th, rows), "ss_dibl_vs_length")
    print("Wrote figures and CSVs to results/")

    print("\nReference device (Lg=250nm, Lgd=2um, d=20nm, x_Al=0.25):")
    for field in ["vth_V", "gm_max_S_per_mm", "ss_mV_per_dec", "dibl_mV_per_V",
                  "ion_A_per_mm", "ron_sp_ohm_mm", "breakdown_V"]:
        print(f"  {field:28s} {getattr(ref, field):.4g}")


if __name__ == "__main__":
    main()
