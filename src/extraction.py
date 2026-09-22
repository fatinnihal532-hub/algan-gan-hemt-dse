"""
Figure-of-merit extraction from a HEMTDevice's DC characteristics.

All extraction conventions are documented explicitly here because, unlike
a measured or TCAD I-V curve, an analytical model's "Vth" or "Ioff" is
only meaningful once you say exactly how it was read off the curve.
"""
import math
from dataclasses import dataclass

from .current_model import HEMTDevice

VDS_HIGH = 5.0          # V, bias used for the nominal transfer characteristic
VDS_LOW = 0.1           # V, low-Vds bias used for the DIBL comparison AND for Ron_sp
OVERDRIVE_ON_V = 4.0    # Vgs_on  = Vth(Vds_high) + OVERDRIVE_ON_V
UNDERDRIVE_OFF_V = 1.0  # Vgs_off = Vth(Vds_high) - UNDERDRIVE_OFF_V
IOFF_FLOOR_A_PER_MM = 1e-6  # representative gate/buffer leakage floor
                            # (typical order of magnitude reported for
                            # AlGaN/GaN HEMTs); the thermionic-only
                            # subthreshold model has no leakage path and
                            # would otherwise predict an unphysically
                            # small Ioff, making Ion/Ioff meaningless.


@dataclass
class DeviceMetrics:
    vth_V: float
    gm_max_S_per_mm: float
    ss_mV_per_dec: float
    dibl_mV_per_V: float
    ion_A_per_mm: float
    ioff_A_per_mm: float
    ion_ioff_ratio: float
    ron_sp_ohm_mm: float
    breakdown_V: float
    baliga_fom_V2_per_ohm_mm: float


def _gm_max(dev: HEMTDevice, vds: float, vgs_grid) -> float:
    ids = [dev.drain_current(v, vds) for v in vgs_grid]
    gm = [(ids[i + 1] - ids[i]) / (vgs_grid[i + 1] - vgs_grid[i])
          for i in range(len(vgs_grid) - 1)]
    return max(gm) if gm else 0.0


def _subthreshold_swing(dev: HEMTDevice, vds: float, vth: float) -> float:
    """SS (mV/decade), read from the slope of log10(Id) vs Vgs over a
    window strictly below threshold where the model's exponential
    subthreshold branch applies.

    Uses the intrinsic current directly (bypassing the access-resistance
    solve): at subthreshold current levels the I*Racc drop is many
    orders of magnitude below Vds, so it is physically negligible here,
    and skipping it avoids gratuitous floating-point underflow in the
    bisection solve at currents far below any physically meaningful
    floor.
    """
    v_lo, v_hi = vth - 0.9, vth - 0.3
    i_lo = dev._intrinsic_current(v_lo, vds)
    i_hi = dev._intrinsic_current(v_hi, vds)
    if i_lo <= 0 or i_hi <= 0 or i_hi <= i_lo:
        return float("nan")
    decades = math.log10(i_hi / i_lo)
    return (v_hi - v_lo) / decades * 1000.0  # V -> mV


def extract_metrics(dev: HEMTDevice) -> DeviceMetrics:
    w_mm = dev.geom.width_m * 1e3

    vth_high = dev.threshold_voltage(VDS_HIGH)
    vth_low = dev.threshold_voltage(VDS_LOW)
    dibl = (vth_low - vth_high) / (VDS_HIGH - VDS_LOW) * 1000.0  # mV/V

    vgs_grid = [vth_high + 0.05 * i for i in range(-4, 100)]
    gm_max = _gm_max(dev, VDS_HIGH, vgs_grid) / w_mm

    ss = _subthreshold_swing(dev, VDS_HIGH, vth_high)

    vgs_on = vth_high + OVERDRIVE_ON_V
    vgs_off = vth_high - UNDERDRIVE_OFF_V
    ion = dev.drain_current(vgs_on, VDS_HIGH) / w_mm
    ioff = max(dev.drain_current(vgs_off, VDS_HIGH) / w_mm, IOFF_FLOOR_A_PER_MM)
    ratio = ion / ioff if ioff > 0 else float("inf")

    bv = dev.breakdown_voltage_estimate()

    # Specific on-resistance: Vds/Id in the linear region at a fixed
    # small Vds, gate driven on -- the standard way Ron is reported for
    # power devices (as opposed to gm, which is read from the
    # saturation-region transfer characteristic above).
    id_lin = dev.drain_current(vgs_on, VDS_LOW)
    ron_sp = (VDS_LOW / id_lin) * w_mm if id_lin > 0 else float("inf")

    # Width-normalized analogue of the Baliga figure of merit, BV^2/Ron_sp
    # (Baliga, B. J., "Power semiconductor device figure of merit for
    # high-frequency applications," IEEE Electron Device Lett. 10, 455
    # (1989)). The textbook Baliga FOM normalizes Ron by chip AREA
    # (Ohm.cm^2); HEMTs are conventionally specified per unit gate WIDTH
    # (Ohm.mm) instead, because periphery, not die area, is what a HEMT
    # layout trades off. Reported here in V^2/(Ohm.mm) -- consistent
    # for comparing configurations within this sweep, but not
    # numerically comparable to an area-normalized Si/SiC Baliga FOM
    # quoted in the literature without knowing the actual device pitch.
    baliga_fom = (bv ** 2) / ron_sp if ron_sp > 0 else float("inf")

    return DeviceMetrics(
        vth_V=vth_high,
        gm_max_S_per_mm=gm_max,
        ss_mV_per_dec=ss,
        dibl_mV_per_V=dibl,
        ion_A_per_mm=ion,
        ioff_A_per_mm=ioff,
        ion_ioff_ratio=ratio,
        ron_sp_ohm_mm=ron_sp,
        breakdown_V=bv,
        baliga_fom_V2_per_ohm_mm=baliga_fom,
    )
