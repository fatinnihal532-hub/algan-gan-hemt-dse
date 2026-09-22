"""Checks the AlGaN/GaN HEMT design-space model against physical sanity
and closed-form theory. Exits non-zero if any check fails.

Run: python3 verify.py   (under a second)
"""
import sys

from src.materials import AlGaN, Q
from src.charge_control import BarrierStack
from src.current_model import HEMTDevice, HEMTGeometry
from src.extraction import extract_metrics
from src.design_space import run_sweep, ron_bv_pareto_front

PASS, FAIL = [], []


def check(name, condition, detail=""):
    if condition:
        PASS.append(name)
        print(f"  [pass] {name}")
    else:
        FAIL.append(name)
        print(f"  [FAIL] {name}")
    if detail:
        print(f"         {detail}")


def main():
    print("AlGaN/GaN HEMT design-space model - checks against theory")
    print("=" * 62)

    stack = BarrierStack(AlGaN(0.25), 20e-9)
    vth = stack.threshold_voltage_V()

    # 1. ns(Vth) should be ~0 by construction of the threshold formula.
    ns_at_vth = stack.sheet_density_m2(vth)
    check("ns(Vth) is ~0 by construction of the Vth formula",
          ns_at_vth < 1e13,  # << typical operating ns of ~1e17 m^-2
          f"ns(Vth)={ns_at_vth:.3e} m^-2")

    # 2. ns(Vgs) is monotonically non-decreasing in Vgs.
    vgs_grid = [-8 + 0.25 * i for i in range(60)]
    ns_grid = [stack.sheet_density_m2(v) for v in vgs_grid]
    check("ns(Vgs) is monotonically non-decreasing",
          all(b >= a - 1e10 for a, b in zip(ns_grid, ns_grid[1:])))

    # 3. Higher Al fraction increases the 2DEG density at a fixed gate bias
    #    (more polarization charge from a larger lattice/piezo mismatch).
    ns_low_x = BarrierStack(AlGaN(0.15), 20e-9).sheet_density_m2(0.0)
    ns_high_x = BarrierStack(AlGaN(0.35), 20e-9).sheet_density_m2(0.0)
    check("higher Al fraction increases ns at fixed Vgs=0",
          ns_high_x > ns_low_x,
          f"x=0.15: {ns_low_x:.3e} m^-2, x=0.35: {ns_high_x:.3e} m^-2")

    # 4. Thinner barrier increases the gate-to-channel capacitive coupling
    #    (Cg = q dns/dVgs), consistent with a parallel-plate capacitor.
    def cg(d):
        s = BarrierStack(AlGaN(0.25), d)
        dv = 1e-3
        return Q * (s.sheet_density_m2(dv) - s.sheet_density_m2(-dv)) / (2 * dv)
    check("thinner barrier gives higher gate capacitance Cg",
          cg(10e-9) > cg(30e-9),
          f"d=10nm: Cg={cg(10e-9):.3e} F/m^2, d=30nm: Cg={cg(30e-9):.3e} F/m^2")

    # 5. Velocity saturation reduces Id below the ideal (non-saturating)
    #    square-law prediction at the same bias.
    geom = HEMTGeometry(gate_length_m=0.25e-6, gate_drain_spacing_m=2e-6)
    dev = HEMTDevice(stack, geom)
    vgs_test, vds_test = vth + 3.0, 3.0
    id_actual = dev._intrinsic_current(vgs_test, vds_test)
    ns_test = stack.sheet_density_m2(vgs_test)
    cg_test = dev._gate_capacitance(vgs_test)
    id_ideal_square_law = ((geom.width_m * dev.proc.mobility_m2Vs / geom.gate_length_m)
                            * (Q * ns_test * vds_test - 0.5 * cg_test * vds_test ** 2))
    check("velocity saturation reduces Id below the ideal square law",
          0 < id_actual < id_ideal_square_law,
          f"Id={id_actual:.4g} A vs ideal square law {id_ideal_square_law:.4g} A")

    # 6. Current continuity of the INTRINSIC (access-resistance-free)
    #    triode/saturation piecewise model: Id at the triode/saturation
    #    boundary matches the current well into saturation. This is
    #    checked on the intrinsic model specifically because the
    #    terminal characteristic (with access resistance) has its own,
    #    physically expected "knee stretch" -- part of the terminal Vds
    #    is dropped resistively before the intrinsic device even reaches
    #    Vdsat, which is a real effect, not a discontinuity.
    vdsat = dev._vdsat(vgs_test)
    id_at_vdsat = dev._intrinsic_current(vgs_test, vdsat)
    id_past_vdsat = dev._intrinsic_current(vgs_test, vdsat + 3.0)
    check("intrinsic Id is continuous across the triode/saturation boundary",
          abs(id_past_vdsat - id_at_vdsat) < 1e-6 * max(id_at_vdsat, 1e-6),
          f"Id(Vdsat)={id_at_vdsat:.4g} A, Id(Vdsat+3V)={id_past_vdsat:.4g} A")

    # 7. Subthreshold swing approaches the room-temperature Boltzmann
    #    limit (~60 mV/dec) as the gate length becomes long (aspect
    #    ratio d_eff/Lg -> 0, i.e. good long-channel electrostatics).
    dev_long = HEMTDevice(stack, HEMTGeometry(gate_length_m=5e-6, gate_drain_spacing_m=2e-6))
    ss_long = extract_metrics(dev_long).ss_mV_per_dec
    check("SS approaches the 60 mV/dec limit for a long gate",
          59.5 < ss_long < 62.0,
          f"SS(Lg=5um)={ss_long:.3f} mV/dec")

    # 8. DIBL vanishes as the gate length becomes long.
    dibl_long = extract_metrics(dev_long).dibl_mV_per_V
    check("DIBL is small for a long gate",
          abs(dibl_long) < 5.0,
          f"DIBL(Lg=5um)={dibl_long:.3f} mV/V")

    # 9. Longer gate-drain spacing increases both the breakdown-voltage
    #    estimate and the specific on-resistance (the access-resistance
    #    penalty is the mechanism that makes this a real trade-off rather
    #    than a free parameter -- see docs/methodology.md).
    dev_short_lgd = HEMTDevice(stack, HEMTGeometry(gate_length_m=0.25e-6, gate_drain_spacing_m=1e-6))
    dev_long_lgd = HEMTDevice(stack, HEMTGeometry(gate_length_m=0.25e-6, gate_drain_spacing_m=6e-6))
    m_short, m_long = extract_metrics(dev_short_lgd), extract_metrics(dev_long_lgd)
    check("longer Lgd increases BV and Ron_sp together",
          m_long.breakdown_V > m_short.breakdown_V and m_long.ron_sp_ohm_mm > m_short.ron_sp_ohm_mm,
          f"Lgd=1um: BV={m_short.breakdown_V:.0f}V Ron={m_short.ron_sp_ohm_mm:.3f}; "
          f"Lgd=6um: BV={m_long.breakdown_V:.0f}V Ron={m_long.ron_sp_ohm_mm:.3f}")

    # 10. The Ron_sp/BV Pareto front is genuinely non-dominated: for every
    #     pair of front points, neither dominates the other.
    rows = run_sweep(
        gate_lengths_m=[x * 1e-9 for x in (150, 250, 500)],
        gate_drain_spacings_m=[x * 1e-6 for x in (1, 2, 4, 8)],
        thicknesses_m=[20e-9], al_fractions=[0.25],
    )
    front = ron_bv_pareto_front(rows)
    non_dominated = True
    for i, a in enumerate(front):
        for j, b in enumerate(front):
            if i == j:
                continue
            if (b["ron_sp_ohm_mm"] <= a["ron_sp_ohm_mm"] and b["breakdown_V"] >= a["breakdown_V"]
                    and (b["ron_sp_ohm_mm"] < a["ron_sp_ohm_mm"] or b["breakdown_V"] > a["breakdown_V"])):
                non_dominated = False
    check("Ron_sp/BV Pareto front is genuinely non-dominated",
          non_dominated and len(front) >= 2,
          f"front size={len(front)}")

    # 11. Determinism / reproducibility: identical inputs give identical
    #     outputs (no hidden randomness anywhere in the model).
    m1 = extract_metrics(HEMTDevice(BarrierStack(AlGaN(0.25), 20e-9), geom))
    m2 = extract_metrics(HEMTDevice(BarrierStack(AlGaN(0.25), 20e-9), geom))
    check("identical inputs reproduce identical extracted metrics",
          m1 == m2)

    print("=" * 62)
    print(f"{len(PASS)} checks passed, {len(FAIL)} failed")
    return 1 if FAIL else 0


if __name__ == "__main__":
    sys.exit(main())
