"""
Polarization-induced 2DEG charge-control model.

Implements the linearized charge-control relation of Ambacher et al.
(J. Appl. Phys. 85, 3222, 1999):

    ns(Vg) = sigma(x)/q
             - eps0*eps_r/(q^2*(d+dd)) * (q*phib + EF - dEc)
             + eps0*eps_r/(q*(d+dd)) * Vg

closed self-consistently against a single-subband, triangular-well
approximation for the Fermi level of a degenerate 2D electron gas:

    EF(ns) = (pi * hbar^2 / m*) * ns

This is the standard simplification used in analytical/compact
AlGaN/GaN HEMT models when a full self-consistent Schrodinger-Poisson
solve (i.e. what TCAD performs) is not being run. It is solved here by
straightforward fixed-point iteration, which converges in a handful of
steps because the EF correction is a small perturbation on ns.
"""
import math
from dataclasses import dataclass

from .materials import AlGaN, EPS0, Q, HBAR, M0, M_STAR_GAN


DELTA_D_M = 0.5e-9   # effective-thickness offset (m): the 2DEG charge
                     # centroid sits a short distance below the AlGaN/GaN
                     # interface, not exactly at it, so the gate-to-channel
                     # capacitance is eps/(d + dd) rather than eps/d. A fixed
                     # 0.5 nm is a representative value, kept small relative
                     # to d; it is a modelling choice, not a fitted one.


@dataclass(frozen=True)
class BarrierStack:
    """A gate-controlled AlGaN barrier over a GaN buffer."""
    alloy: AlGaN
    thickness_m: float          # AlGaN barrier thickness, d
    m_star: float = M_STAR_GAN  # effective mass, units of m0

    @property
    def eps(self) -> float:
        return EPS0 * self.alloy.eps_r

    def _fermi_level_eV(self, ns_m2: float) -> float:
        """EF above the GaN conduction-band edge for a degenerate 2DEG
        of sheet density ns (single-subband triangular-well DOS)."""
        if ns_m2 <= 0:
            return 0.0
        m_eff = self.m_star * M0
        ef_J = (math.pi * HBAR ** 2 / m_eff) * ns_m2
        return ef_J / Q

    def sheet_density_m2(self, vgs: float, n_iter: int = 10) -> float:
        """Self-consistent 2DEG sheet density ns (m^-2) at gate bias vgs.

        Returns 0 if the gate bias depletes the channel (below threshold),
        rather than a spurious negative density.
        """
        d_eff = self.thickness_m + DELTA_D_M
        sigma = self.alloy.total_polarization_charge_Cm2
        phib = self.alloy.schottky_barrier_eV()
        dEc = self.alloy.delta_ec_eV
        cap_term = self.eps / (Q * d_eff)  # eps/(q*d), units 1/(V.m^2) * ... see below

        ef = 0.0
        ns = 0.0
        for _ in range(n_iter):
            ns_new = (sigma / Q
                      - (self.eps / (Q ** 2 * d_eff)) * (Q * phib + Q * ef - Q * dEc)
                      + cap_term * vgs)
            ns_new = max(ns_new, 0.0)
            ef_new = self._fermi_level_eV(ns_new)
            if abs(ns_new - ns) < 1e-6 * max(ns_new, 1.0):
                ns, ef = ns_new, ef_new
                break
            ns, ef = ns_new, ef_new
        return ns

    def threshold_voltage_V(self) -> float:
        """Analytical Vth: the gate bias at which ns -> 0 and EF -> 0
        simultaneously (the ns = 0 limit of the relation above)."""
        d_eff = self.thickness_m + DELTA_D_M
        sigma = self.alloy.total_polarization_charge_Cm2
        phib = self.alloy.schottky_barrier_eV()
        dEc = self.alloy.delta_ec_eV
        return phib - dEc - sigma * d_eff / self.eps
