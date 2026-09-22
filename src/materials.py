"""
AlxGa1-xN / GaN material parameters and polarization-induced sheet charge.

All formulas and reference parameters follow the standard analytical
treatment of polarization-induced 2DEG formation in AlGaN/GaN heterostructures:

  Ambacher, O. et al., "Two-dimensional electron gases induced by
  spontaneous and piezoelectric polarization in undoped and doped
  AlGaN/GaN heterostructures," J. Appl. Phys. 85, 3222 (1999).

  Ambacher, O. et al., "Pyroelectric properties of Al(In)GaN/GaN
  hetero- and quantum well structures," J. Phys.: Condens. Matter
  14, 3399 (2002).  [elastic constants, piezoelectric coefficients]

  Vurgaftman, I. & Meyer, J. R., "Band parameters for nitrogen-
  containing semiconductors," J. Appl. Phys. 89, 5815 (2001).
  [bandgap bowing]

This is a deliberately simplified, closed-form / self-consistent-iteration
model, NOT a TCAD (Schrodinger-Poisson) solve. It is meant to expose clean,
physically-motivated trends across a device design space, not to reproduce
a specific measured wafer to high accuracy. Every constant below is a
representative literature value, not a fit to a specific fabricated device.
"""
from dataclasses import dataclass

Q = 1.602176634e-19          # elementary charge, C
EPS0 = 8.8541878128e-12      # vacuum permittivity, F/m
HBAR = 1.054571817e-34       # reduced Planck constant, J.s
M0 = 9.1093837015e-31        # free electron mass, kg

# --- Binary endpoint parameters (GaN, AlN) -----------------------------
# Relative dielectric constant (static)
EPS_R_GAN, EPS_R_ALN = 8.9, 8.5
# Bandgap (eV) and bowing parameter (Vurgaftman & Meyer 2001)
EG_GAN, EG_ALN, EG_BOW = 3.42, 6.20, 1.00
# In-plane lattice constant (m)  (Ambacher 2002)
A0_GAN, A0_ALN = 3.189e-10, 3.112e-10
# Elastic constants (Pa) (Ambacher 2002)
C13_GAN, C13_ALN = 103e9, 108e9
C33_GAN, C33_ALN = 405e9, 373e9
# Piezoelectric coefficients (C/m^2) (Ambacher 2002)
E31_GAN, E31_ALN = -0.49, -0.60
E33_GAN, E33_ALN = 0.73, 1.46
# Spontaneous polarization (C/m^2) and bowing (Ambacher 1999/2002)
PSP_GAN, PSP_ALN, PSP_BOW = -0.029, -0.081, -0.021
# Conduction-band effective mass (in-plane, GaN 2DEG), in units of m0
M_STAR_GAN = 0.22
# Conduction-band offset fraction of the bandgap difference that appears
# as a barrier for electrons -- the commonly used 70/30 split of the
# AlGaN/GaN band offset: DeltaEc = 0.7 * DeltaEg
DELTA_EC_FRACTION = 0.70


def _lerp(x, a, b):
    """Linear (Vegard-law style) interpolation between binary endpoints."""
    return x * b + (1 - x) * a


@dataclass(frozen=True)
class AlGaN:
    """A single AlxGa1-xN barrier composition, with all derived parameters
    computed once at construction time."""
    x: float  # Al mole fraction, 0 <= x <= 1

    def __post_init__(self):
        if not (0.0 <= self.x <= 1.0):
            raise ValueError(f"Al mole fraction x={self.x} out of [0, 1]")

    @property
    def eps_r(self) -> float:
        return _lerp(self.x, EPS_R_GAN, EPS_R_ALN)

    @property
    def bandgap_eV(self) -> float:
        # Quadratic bowing: Eg(x) = x*Eg_AlN + (1-x)*Eg_GaN - b*x*(1-x)
        return _lerp(self.x, EG_GAN, EG_ALN) - EG_BOW * self.x * (1 - self.x)

    @property
    def delta_ec_eV(self) -> float:
        """Conduction-band offset (eV) presented to the 2DEG at the
        AlGaN/GaN interface."""
        delta_eg = self.bandgap_eV - EG_GAN
        return DELTA_EC_FRACTION * delta_eg

    @property
    def lattice_a_m(self) -> float:
        return _lerp(self.x, A0_GAN, A0_ALN)

    @property
    def piezo_polarization_Cm2(self) -> float:
        """Piezoelectric polarization (C/m^2) of AlGaN pseudomorphically
        strained in-plane to relaxed GaN (Ambacher 1999)."""
        a = self.lattice_a_m
        c13 = _lerp(self.x, C13_GAN, C13_ALN)
        c33 = _lerp(self.x, C33_GAN, C33_ALN)
        e31 = _lerp(self.x, E31_GAN, E31_ALN)
        e33 = _lerp(self.x, E33_GAN, E33_ALN)
        strain_inplane = (A0_GAN - a) / a
        return 2 * strain_inplane * (e31 - e33 * c13 / c33)

    @property
    def spontaneous_polarization_Cm2(self) -> float:
        return (_lerp(self.x, PSP_GAN, PSP_ALN)
                + PSP_BOW * self.x * (1 - self.x))

    @property
    def total_polarization_charge_Cm2(self) -> float:
        """Net bound polarization sheet charge (C/m^2) at the AlGaN/GaN
        interface, relative to an unstrained GaN buffer (whose
        polarization difference is zero by construction)."""
        delta_psp = self.spontaneous_polarization_Cm2 - PSP_GAN
        return abs(delta_psp + self.piezo_polarization_Cm2)

    def schottky_barrier_eV(self) -> float:
        """Metal/AlGaN Schottky barrier height.

        Modelled as increasing linearly with Al content, consistent with
        the qualitative trend reported for Ni/Au contacts on AlGaN
        (higher Al content -> larger barrier). This is a representative
        modelling choice, not a fit to one specific reported contact, and
        is documented here explicitly as such.
        """
        return 0.84 + 1.3 * self.x
