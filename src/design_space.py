"""
Design-space sweep across gate length, gate-drain spacing, barrier
thickness and Al composition, plus Pareto-front extraction.

Everything here is downstream of the analytical/self-consistent model in
materials.py, charge_control.py and current_model.py: this module only
enumerates device configurations, evaluates extract_metrics() on each,
and post-processes the resulting table. No new physics is introduced.
"""
from dataclasses import asdict
from itertools import product

from .materials import AlGaN
from .charge_control import BarrierStack
from .current_model import HEMTDevice, HEMTGeometry
from .extraction import extract_metrics


def build_device(gate_length_m, gate_drain_spacing_m, thickness_m, x_al):
    stack = BarrierStack(AlGaN(x_al), thickness_m)
    geom = HEMTGeometry(gate_length_m=gate_length_m,
                         gate_drain_spacing_m=gate_drain_spacing_m)
    return HEMTDevice(stack, geom)


def run_sweep(gate_lengths_m, gate_drain_spacings_m, thicknesses_m, al_fractions):
    """Full factorial sweep. Returns a list of flat dicts, one per
    configuration, combining the swept parameters with extract_metrics()."""
    rows = []
    for lg, lgd, d, x in product(gate_lengths_m, gate_drain_spacings_m,
                                  thicknesses_m, al_fractions):
        dev = build_device(lg, lgd, d, x)
        metrics = extract_metrics(dev)
        row = {
            # rounded only to strip float noise from the unit conversion
            # (e.g. 100.00000000000001 nm); the model itself uses the
            # unrounded SI values above
            "gate_length_nm": round(lg * 1e9, 6),
            "gate_drain_spacing_um": round(lgd * 1e6, 6),
            "barrier_thickness_nm": round(d * 1e9, 6),
            "al_fraction": round(x, 6),
        }
        row.update(asdict(metrics))
        rows.append(row)
    return rows


def pareto_front(rows, maximize_key, also_maximize_key):
    """Non-dominated set for two objectives, both maximized.

    A configuration is Pareto-optimal if no other configuration in the
    sweep is at least as good on both objectives and strictly better on
    at least one -- the standard definition of Pareto dominance for a
    two-objective trade-off (see e.g. any multi-objective optimization
    text, e.g. Deb, "Multi-Objective Optimization Using Evolutionary
    Algorithms," Wiley, 2001, Ch. 2).
    """
    front = []
    for i, row in enumerate(rows):
        a, b = row[maximize_key], row[also_maximize_key]
        dominated = False
        for j, other in enumerate(rows):
            if i == j:
                continue
            oa, ob = other[maximize_key], other[also_maximize_key]
            if oa >= a and ob >= b and (oa > a or ob > b):
                dominated = True
                break
        if not dominated:
            front.append(row)
    front.sort(key=lambda r: r[maximize_key])
    return front


def ron_bv_pareto_front(rows):
    """Pareto front for the two objectives that actually trade off against
    each other in this model: minimizing Ron_sp and maximizing breakdown
    voltage (the Baliga trade-off). Implemented via pareto_front() on
    -Ron_sp (so "maximize" is the right sense for both objectives)."""
    augmented = [dict(r, _neg_ron_sp=-r["ron_sp_ohm_mm"]) for r in rows]
    front = pareto_front(augmented, "_neg_ron_sp", "breakdown_V")
    for row in front:
        row.pop("_neg_ron_sp", None)
    return front
