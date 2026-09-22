"""
DC drain-current model for a AlGaN/GaN HEMT, built directly on top of the
polarization charge-control model in charge_control.py.

Core (long-channel) physics: a velocity-saturated "charge-sheet" model,
the same structure as the standard velocity-saturated square-law MOSFET
model (see e.g. Sze & Ng, "Physics of Semiconductor Devices," 3rd ed.,
Ch. 6; Muller & Kamins, "Device Electronics for Integrated Circuits"),
except that the channel charge ns(Vgs) and its incremental capacitance
Cg(Vgs) = q*dns/dVgs come directly from the self-consistent polarization
charge-control relation, not from an assumed constant oxide/MOS
capacitance.

Two effects that a 1D charge-control model cannot capture without a full
2D (TCAD) solve -- subthreshold-swing degradation and drain-induced
barrier lowering (DIBL) as the gate shrinks relative to the barrier
thickness -- are added as clearly-labelled empirical proxies keyed to the
aspect ratio d_eff/Lg, which is the same "scale length" idea used in
rigorous short-channel scaling theory (Yan, Ourmazd and Lee, "Scaling the
Si MOSFET: from bulk to SOI to bulk," IEEE Trans. Electron Devices 39,
1704 (1992)) without repeating its full derivation. These proxies are the explicitly-approximate part of the
model; everything upstream of them (polarization charge, Vth0, velocity
saturation) is closed-form device physics.

Performance note: ns(Vgs) and Cg(Vgs) only depend on Vgs, not Vds, but the
terminal current solve (drain_current) holds Vgs fixed while it searches
over an *intrinsic* Vds by bisection. Recomputing the self-consistent
charge-control iteration inside every bisection step would multiply two
already-iterative solves together for no reason, so _at_vgs() computes
ns/Cg/Vdsat once per (Vgs) and the bisection loop reuses that cached
result -- the physics is unchanged, only the redundant recomputation is
removed.
"""
import math
from dataclasses import dataclass

from .charge_control import BarrierStack
from .materials import Q

KB = 1.380649e-23          # Boltzmann constant, J/K
T_DEFAULT = 300.0          # K


@dataclass(frozen=True)
class HEMTGeometry:
    gate_length_m: float          # Lg
    gate_drain_spacing_m: float   # Lgd (access-region / field-plate proxy)
    source_gate_spacing_m: float = 1.0e-6  # Lsg, held fixed across the sweep
    width_m: float = 100e-6       # gate width W (fixed; not swept)


@dataclass(frozen=True)
class ShortChannelParams:
    """Empirical proxies for 2D short-channel electrostatics, keyed to
    the barrier-thickness / gate-length aspect ratio d_eff/Lg. See module
    docstring: this is the explicitly-approximate half of the model."""
    ss_ideality_n0: float = 0.8     # subthreshold ideality-factor slope
    dibl_coeff_V_per_V: float = 0.6  # DIBL scale coefficient


@dataclass(frozen=True)
class ProcessParams:
    mobility_m2Vs: float = 0.15    # low-field 2DEG mobility, m^2/(V.s)
                                    # (1500 cm^2/V.s, representative AlGaN/GaN)
    v_sat_mps: float = 1.3e5       # saturation velocity, m/s (1.3e7 cm/s)
    e_crit_Vpm: float = 3.3e8      # GaN critical (avalanche) field, V/m (3.3 MV/cm)
    temperature_K: float = T_DEFAULT


@dataclass(frozen=True)
class _VgsPoint:
    """Everything about the intrinsic device that depends on Vgs alone,
    computed once and reused across an entire bisection solve."""
    ns: float
    cg: float
    vdsat: float


class HEMTDevice:
    """A complete DC device model: barrier stack + geometry + process."""

    def __init__(self, stack: BarrierStack, geom: HEMTGeometry,
                 proc: ProcessParams = ProcessParams(),
                 sc: ShortChannelParams = ShortChannelParams()):
        self.stack = stack
        self.geom = geom
        self.proc = proc
        self.sc = sc
        self._vth0 = stack.threshold_voltage_V()
        # Ungated 2DEG sheet density under the source/drain access regions,
        # taken at Vgs=0 as a representative (fixed) value -- the access
        # regions are not gated, so their charge does not track Vgs the
        # way the intrinsic channel does. Sheet resistance follows
        # directly: Rsh = 1/(q*ns0*mu) (Ohms per square).
        self._ns_access = stack.sheet_density_m2(0.0)
        self._ec_l = (proc.v_sat_mps / proc.mobility_m2Vs) * geom.gate_length_m
        self._vgs_cache: dict = {}

    # -- aspect ratio driving the short-channel proxies -------------------
    @property
    def _aspect_ratio(self) -> float:
        d_eff = self.stack.thickness_m + 0.5e-9
        return d_eff / self.geom.gate_length_m

    def threshold_voltage(self, vds: float) -> float:
        """Vth(Vds), including the DIBL proxy roll-off with drain bias."""
        dibl = self.sc.dibl_coeff_V_per_V * self._aspect_ratio * vds
        return self._vth0 - dibl

    def _gate_capacitance(self, vgs: float, dv: float = 1e-3) -> float:
        """Incremental gate capacitance Cg = q * dns/dVgs (F/m^2),
        evaluated by central finite difference of the self-consistent
        charge-control relation."""
        ns_p = self.stack.sheet_density_m2(vgs + dv)
        ns_m = self.stack.sheet_density_m2(vgs - dv)
        return Q * (ns_p - ns_m) / (2 * dv)

    def _at_vgs(self, vgs: float) -> _VgsPoint:
        """ns(Vgs), Cg(Vgs) and Vdsat(Vgs), computed once per distinct Vgs
        value and cached for the lifetime of this device instance."""
        cached = self._vgs_cache.get(vgs)
        if cached is not None:
            return cached
        ns = self.stack.sheet_density_m2(vgs)
        cg = self._gate_capacitance(vgs)
        if ns <= 0 or cg <= 0:
            vdsat = 0.0
        else:
            vov = Q * ns / cg  # q*ns/Cg: effective overdrive voltage
            vdsat = self._ec_l * (math.sqrt(1 + 2 * vov / self._ec_l) - 1)
        point = _VgsPoint(ns=ns, cg=cg, vdsat=vdsat)
        self._vgs_cache[vgs] = point
        return point

    def _vdsat(self, vgs: float) -> float:
        return self._at_vgs(vgs).vdsat

    def _id_triode_from_point(self, point: "_VgsPoint", vds: float) -> float:
        if point.ns <= 0:
            return 0.0
        w, L = self.geom.width_m, self.geom.gate_length_m
        numer = Q * point.ns * vds - 0.5 * point.cg * vds ** 2
        denom = 1.0 + vds / self._ec_l
        return (w * self.proc.mobility_m2Vs / L) * numer / denom

    def _id_triode(self, vgs: float, vds: float) -> float:
        """Velocity-saturated triode-region current (A), for vds <= vdsat."""
        return self._id_triode_from_point(self._at_vgs(vgs), vds)

    def _id_subthreshold(self, vgs: float, vds: float) -> float:
        """Exponential subthreshold model, anchored at Vth to the
        (small, non-zero) triode current evaluated just above threshold,
        with an ideality factor that grows as the gate shrinks relative
        to the barrier (see module docstring)."""
        vth = self.threshold_voltage(vds)
        n = 1.0 + self.sc.ss_ideality_n0 * self._aspect_ratio
        vt = KB * self.proc.temperature_K / Q
        anchor_vgs = vth + 0.05  # 50 mV above threshold, still in the
                                  # charge-control model's valid range
        anchor_point = self._at_vgs(anchor_vgs)
        i_anchor = max(self._id_triode_from_point(anchor_point, min(vds, anchor_point.vdsat)), 1e-18)
        return i_anchor * math.exp((vgs - anchor_vgs) / (n * vt))

    def _access_resistance_ohm(self) -> float:
        """Total ungated source+drain access resistance (source-gate plus
        gate-drain spacing), from the sheet resistance of the ungated 2DEG.
        Both access regions are simple 2DEG resistors in series with the
        intrinsic (gated) channel."""
        if self._ns_access <= 0:
            return 0.0
        r_sheet = 1.0 / (Q * self._ns_access * self.proc.mobility_m2Vs)  # Ohm/sq
        length_total = self.geom.source_gate_spacing_m + self.geom.gate_drain_spacing_m
        return r_sheet * length_total / self.geom.width_m

    def _intrinsic_current(self, vgs: float, vds_intrinsic: float) -> float:
        """Id(Vgs, Vds) across the intrinsic (gated) device only, ignoring
        access resistance -- this is the model described above."""
        vth = self.threshold_voltage(vds_intrinsic)
        if vgs < vth:
            return self._id_subthreshold(vgs, vds_intrinsic)
        point = self._at_vgs(vgs)
        if vds_intrinsic <= point.vdsat:
            return self._id_triode_from_point(point, vds_intrinsic)
        return self._id_triode_from_point(point, point.vdsat)  # saturation:
        # current continuity, channel-length modulation neglected

    def drain_current(self, vgs: float, vds: float, n_iter: int = 40) -> float:
        """Terminal Id(Vgs, Vds) in amps, including the voltage dropped
        across the source+drain access resistance.

        Id must satisfy Id = f(Vgs, Vds - Id*Racc), i.e. it is the root of
        g(Id) = f(Vgs, Vds - Id*Racc) - Id = 0. Naive fixed-point
        iteration on this relation can oscillate (a large Racc overshoots
        past Vds_intrinsic=0 and back every step), so this is solved by
        bisection instead: f is non-increasing in Id (more current means
        more drop across Racc means less intrinsic Vds means less or
        equal intrinsic current), so g is strictly decreasing and has a
        single root bracketed by Id=0 (g>=0) and Id=Vds/Racc, the current
        that would drop all of Vds across Racc (g<=0 there). Vgs is fixed
        throughout this search, so _at_vgs(vgs) is computed once (via the
        cache) rather than once per bisection step.
        """
        if vds <= 0:
            return 0.0
        racc = self._access_resistance_ohm()
        if racc <= 0:
            return self._intrinsic_current(vgs, vds)

        def g(id_val):
            vds_intrinsic = max(vds - id_val * racc, 0.0)
            return self._intrinsic_current(vgs, vds_intrinsic) - id_val

        lo, hi = 0.0, vds / racc
        g_lo = g(lo)
        if g_lo <= 0:
            return 0.0
        for _ in range(n_iter):
            mid = 0.5 * (lo + hi)
            if g(mid) > 0:
                lo = mid
            else:
                hi = mid
        return 0.5 * (lo + hi)

    def breakdown_voltage_estimate(self) -> float:
        """A first-order, parallel-plate estimate: BV ~ Ecrit * Lgd.

        This neglects field crowding at the gate edge, the field plate
        (if any), and vertical buffer leakage -- all of which set the
        breakdown voltage of a real device and require a 2D field solve
        (TCAD) to capture. It is used here purely as a monotonic proxy
        for how the gate-drain access-region length trades off against
        on-state performance in the design-space sweep, not as a
        prediction of a specific device's rating.
        """
        return self.proc.e_crit_Vpm * self.geom.gate_drain_spacing_m
