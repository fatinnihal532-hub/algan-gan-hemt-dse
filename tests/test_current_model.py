import unittest
from src.materials import AlGaN
from src.charge_control import BarrierStack
from src.current_model import HEMTDevice, HEMTGeometry


class TestCurrentModel(unittest.TestCase):
    def setUp(self):
        self.stack = BarrierStack(AlGaN(0.25), 20e-9)
        self.geom = HEMTGeometry(gate_length_m=0.25e-6, gate_drain_spacing_m=2e-6)
        self.dev = HEMTDevice(self.stack, self.geom)

    def test_zero_current_at_zero_drain_bias(self):
        self.assertEqual(self.dev.drain_current(2.0, 0.0), 0.0)

    def test_zero_current_below_threshold_far_out(self):
        # Far below threshold the exponential subthreshold model predicts
        # a current many orders of magnitude below anything physically
        # meaningful (< 1 pA/mm here); the exact floor is sensitive to the
        # charge-control solver's convergence tolerance at these extreme
        # tails, so this checks "negligible relative to Ion", not an
        # arbitrary absolute floor.
        vth = self.dev.threshold_voltage(5.0)
        ion = self.dev.drain_current(vth + 4.0, 5.0)
        ioff_far = self.dev.drain_current(vth - 4.0, 5.0)
        self.assertLess(ioff_far, 1e-9 * ion)

    def test_current_increases_with_gate_overdrive(self):
        vth = self.dev.threshold_voltage(5.0)
        i1 = self.dev.drain_current(vth + 1.0, 5.0)
        i2 = self.dev.drain_current(vth + 3.0, 5.0)
        self.assertGreater(i2, i1)

    def test_current_saturates_with_drain_bias(self):
        vgs = self.dev.threshold_voltage(5.0) + 2.0
        i_mid = self.dev.drain_current(vgs, 5.0)
        i_high = self.dev.drain_current(vgs, 15.0)
        # Saturation (no channel-length modulation): current should be
        # essentially flat, not still rising linearly.
        self.assertLess(abs(i_high - i_mid), 0.05 * i_mid)

    def test_longer_access_region_increases_ron(self):
        short = HEMTDevice(self.stack, HEMTGeometry(gate_length_m=0.25e-6, gate_drain_spacing_m=1e-6))
        long_ = HEMTDevice(self.stack, HEMTGeometry(gate_length_m=0.25e-6, gate_drain_spacing_m=6e-6))
        self.assertGreater(long_._access_resistance_ohm(), short._access_resistance_ohm())

    def test_breakdown_estimate_scales_with_gate_drain_spacing(self):
        short = HEMTDevice(self.stack, HEMTGeometry(gate_length_m=0.25e-6, gate_drain_spacing_m=1e-6))
        long_ = HEMTDevice(self.stack, HEMTGeometry(gate_length_m=0.25e-6, gate_drain_spacing_m=4e-6))
        self.assertAlmostEqual(long_.breakdown_voltage_estimate(),
                                4 * short.breakdown_voltage_estimate(), places=6)


if __name__ == "__main__":
    unittest.main()
