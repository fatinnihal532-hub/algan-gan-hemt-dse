# Methodology: full derivations

This document gives every equation used in `src/`, in the order the code
evaluates them, and states explicitly which parts are closed-form device
physics and which are documented simplifications standing in for a full
2D/TCAD solve. See the README's Limitations section for the short version.

## 1. Alloy interpolation (`materials.py`)

For AlxGa1-xN, every binary-endpoint property `P` (bandgap, lattice
constant, elastic constants, piezoelectric coefficients, dielectric
constant) is interpolated linearly (Vegard's law) between the GaN (x=0)
and AlN (x=1) values:

```
P(x) = x * P_AlN + (1 - x) * P_GaN
```

Two properties additionally include a quadratic bowing correction, per the
cited literature:

**Bandgap** (Vurgaftman & Meyer, 2001):
```
Eg(x) = x*Eg_AlN + (1-x)*Eg_GaN - b_Eg * x * (1-x),   b_Eg = 1.00 eV
```

**Spontaneous polarization** (Ambacher, 1999/2002):
```
P_sp(x) = x*P_sp,AlN + (1-x)*P_sp,GaN + b_Psp * x * (1-x),   b_Psp = -0.021 C/m^2
```

## 2. Polarization sheet charge (`materials.py`)

**Piezoelectric polarization** of AlGaN pseudomorphically strained in-plane
to relaxed GaN (Ambacher 1999):

```
P_pe(x) = 2 * [a_GaN - a(x)] / a(x) * [e31(x) - e33(x) * C13(x)/C33(x)]
```

where `a(x)` is the AlGaN lattice constant and `C13`, `C33`, `e31`, `e33`
are the interpolated elastic and piezoelectric constants.

**Total interface polarization charge** (magnitude, relative to an
unstrained GaN buffer whose own polarization difference is zero):

```
sigma(x) = | [P_sp(x) - P_sp,GaN] + P_pe(x) |
```

## 3. Conduction-band offset and Schottky barrier

```
dEc(x) = 0.70 * [Eg(x) - Eg_GaN]           (the commonly used 70/30 band-offset split)
phi_b(x) = 0.84 + 1.3 * x    [eV]           (representative linear fit — see
                                              materials.py docstring; not a
                                              measured value for a specific
                                              reported contact)
```

## 4. Self-consistent 2DEG sheet density (`charge_control.py`)

Ambacher's linearized charge-control relation (Ambacher et al., 1999):

```
ns(Vg) = sigma(x)/q
         - eps0*eps_r(x) / [q^2 (d + Δd)] * [q*phi_b(x) + EF - dEc(x)]
         + eps0*eps_r(x) / [q (d + Δd)] * Vg
```

`Δd = 0.5 nm` is a small, fixed effective-thickness offset: the 2DEG's charge
centroid sits a short distance below the interface rather than exactly at
it, so the gate sees a slightly thicker dielectric than the barrier alone.
0.5 nm is a representative modelling choice, not a fitted value.

`EF` (the electron gas's Fermi level above the GaN conduction-band edge)
is *itself* a function of `ns`, via the single-subband, triangular-well
density-of-states relation for a degenerate 2D electron gas:

```
EF(ns) = (pi * hbar^2 / m*) * ns,    m* = 0.22 m0 (GaN conduction band)
```

Because `ns` and `EF` depend on each other, the code solves them together by
fixed-point iteration (`BarrierStack.sheet_density_m2`): guess `EF=0`,
compute `ns`, recompute `EF(ns)`, repeat. This converges in a handful of
iterations because `EF` is a small (tens of meV) correction to a charge
relation whose dominant term is linear in `Vg`. This iterative EF closure —
rather than a full self-consistent Schrödinger-Poisson solve — is the
standard simplification made in analytical/compact HEMT charge-control
models when a full quantum solve (TCAD) is not being run.

**Threshold voltage** falls out in closed form by setting `ns=0` and,
self-consistently, `EF(0)=0`:

```
Vth(x, d) = phi_b(x) - dEc(x) - sigma(x) * (d+Δd) / [eps0*eps_r(x)]
```

## 5. Drain current (`current_model.py`)

### 5.1 Incremental gate capacitance

Rather than assume a constant MOS-style oxide capacitance, the incremental
gate capacitance is extracted directly and numerically from the
self-consistent charge-control relation:

```
Cg(Vgs) = q * dns/dVgs   (central finite difference, step 1 mV)
```

This ties the current model's charge term to the same polarization physics
used everywhere else in the framework, rather than introducing a second,
inconsistent notion of "gate capacitance."

### 5.2 Velocity-saturated triode current

Using the standard velocity-saturated charge-sheet model (structurally
identical to the classic velocity-saturated square-law MOSFET model, e.g.
Sze & Ng, *Physics of Semiconductor Devices*, 3rd ed., Ch. 6), but with
`ns(Vgs)` and `Cg(Vgs)` from the polarization model above instead of an
assumed constant capacitance:

```
Id,triode(Vgs,Vds) = (W*mu/L) * [q*ns(Vgs)*Vds - 0.5*Cg(Vgs)*Vds^2] / (1 + Vds/(Ec*L))

Ec*L = (v_sat / mu) * L    [a voltage; "critical field times channel length"]
```

### 5.3 Saturation drain voltage

Setting the overdrive in terms of the actual charge-control charge (rather
than an idealized `Vgs - Vth`):

```
V_ov = q*ns(Vgs) / Cg(Vgs)     [an effective overdrive voltage]
Vdsat = Ec*L * [ sqrt(1 + 2*V_ov/(Ec*L)) - 1 ]
```

This is the standard closed-form result for the velocity-saturated
long-channel MOSFET model (e.g. Muller & Kamins, *Device Electronics for
Integrated Circuits*), evaluated here with the polarization-model overdrive
instead of the ideal square-law one.

### 5.4 Saturation-region current

```
Id(Vgs,Vds>=Vdsat) = Id,triode(Vgs, Vdsat)
```

i.e. current continuity at the triode/saturation boundary, with
channel-length modulation deliberately neglected (documented limitation).

### 5.5 Subthreshold conduction

The charge-control model gives `ns=0` (and hence `Id=0`) for `Vgs<Vth`,
which is not physically realistic — real devices show finite thermionic
subthreshold conduction. This is modelled as an exponential, anchored to the
triode current 50 mV above threshold:

```
Id,sub(Vgs,Vds) = Id(Vth+0.05, Vds) * exp[(Vgs - Vth - 0.05) / (n * Vt)]

Vt = kT/q  (thermal voltage, ~25.85 mV at 300 K)
n(Lg) = 1 + n0 * (d_eff/Lg),   n0 = 0.8    [explicitly-approximate proxy]
```

`n(Lg)` is an empirical stand-in for genuine 2D short-channel electrostatic
loss of gate control — the same physical idea captured rigorously by a
"scale length" analysis (Yan, Ourmazd and Lee, IEEE Trans. Electron
Devices, 1992), applied here as a simple monotonic proxy rather than
its full derivation for this specific HEMT cross-section.

### 5.6 Drain-induced barrier lowering (DIBL)

Modelled as a Vds-dependent threshold roll-off, using the same aspect-ratio
proxy as the subthreshold-swing ideality factor:

```
Vth(Lg, Vds) = Vth0(x,d) - k_DIBL * (d_eff/Lg) * Vds,    k_DIBL = 0.6
```

Again, this is the explicitly-approximate half of the model — see
Limitations in the README.

### 5.7 Source/drain access resistance

The source-gate and gate-drain regions are ungated 2DEG, with sheet density
fixed at its `Vgs=0` value (a common simplification: the access-region
charge is not gate-modulated):

```
R_sheet = 1 / (q * ns(Vgs=0) * mu)                [Ohm/square]
R_access = R_sheet * (Lsg + Lgd) / W              [Lsg fixed at 1 um]
```

The terminal current must satisfy `Id = f(Vgs, Vds - Id*R_access)`, where
`f` is the intrinsic (access-resistance-free) current model above. Naive
fixed-point iteration on this relation can oscillate when `R_access` is
large (the correction overshoots past `Vds_intrinsic=0` and back every
step), so it is solved instead by bisection on `g(Id) = f(...) - Id`, which
is guaranteed to have a unique root because `f` is non-increasing in the
implied `Id` (more current means more voltage dropped across `R_access`
means less intrinsic drain bias means less or equal intrinsic current).

### 5.8 Breakdown-voltage estimate

A first-order, parallel-plate estimate used purely as a monotonic proxy for
how the gate-drain access length trades off against on-state performance,
**not** a prediction of a specific device's rating (see Limitations):

```
BV ≈ E_crit * Lgd,    E_crit(GaN) = 3.3 MV/cm
```

## 6. Figure-of-merit extraction (`extraction.py`)

- **Vth**: `threshold_voltage(Vds=5V)` — the closed-form Vth above, including
  the DIBL roll-off at that bias.
- **gm,max**: maximum of the finite-difference derivative of the simulated
  `Id(Vgs, 5V)` transfer characteristic.
- **Subthreshold swing**: `dVgs / d(log10 Id)`, fit over the window
  `[Vth-0.9, Vth-0.3]` using the *intrinsic* current (access-resistance
  drop is many orders of magnitude below Vds at these current levels, so
  including it there only invites floating-point underflow without
  changing the physics).
- **DIBL**: `[Vth(Vds=0.1V) - Vth(Vds=5V)] / (5V - 0.1V)`.
- **Ion / Ioff**: read at `Vgs = Vth+4V` / `Vgs = Vth-1V`, both at
  `Vds=5V`. Ioff is floored at 1 µA/mm, a representative order-of-magnitude
  gate/buffer leakage floor for real AlGaN/GaN HEMTs — the model's
  thermionic-only subthreshold branch has no leakage path and would
  otherwise predict an unphysically small Ioff, making Ion/Ioff meaningless.
- **Ron,sp**: `Vds / Id` at `Vgs=Vth+4V`, `Vds=0.1V` (deep linear region),
  normalized per unit gate width (Ohm·mm — the HEMT convention, as opposed
  to the Ohm·cm² area-normalized convention used for vertical Si/SiC power
  devices, since HEMT layouts trade off periphery, not die area).
- **Baliga figure of merit**: `BV^2 / Ron,sp`, reported in V²/(Ohm·mm) — see
  the docstring in `extraction.py` for why this is not numerically
  comparable to a textbook area-normalized Baliga FOM without knowing an
  actual device pitch.

## 7. Design-space sweep and Pareto front (`design_space.py`)

A full-factorial sweep over `(Lg, Lgd, d, x)` evaluates every metric above
for each configuration. The Pareto front for two objectives (both
maximized) is the standard non-dominated set: configuration A is removed if
some other configuration B is at least as good on both objectives and
strictly better on at least one. For the Ron,sp/BV trade-off, this is
applied to `(-Ron,sp, BV)` so that both objectives are "maximize."
