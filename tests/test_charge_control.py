import unittest
from src.materials import AlGaN
from src.charge_control import BarrierStack


class TestChargeControl(unittest.TestCase):
    def setUp(self):
        self.stack = BarrierStack(AlGaN(0.25), 20e-9)

    def test_ns_is_zero_below_threshold(self):
        vth = self.stack.threshold_voltage_V()
        self.assertEqual(self.stack.sheet_density_m2(vth - 2.0), 0.0)

    def test_ns_is_positive_and_reasonable_above_threshold(self):
        ns = self.stack.sheet_density_m2(0.0)
        # A representative undoped AlGaN/GaN 2DEG sits in the low-to-mid
        # 1e12-1e13 cm^-2 range (1e16-1e17 m^-2).
        self.assertGreater(ns, 1e15)
        self.assertLess(ns, 1e18)

    def test_thicker_barrier_gives_more_negative_threshold(self):
        # Vth = phib - dEc/q - sigma*d/eps: a thicker barrier makes Vth
        # more negative (larger polarization contribution to charge
        # neutrality has to be compensated by a more negative gate bias).
        vth_thin = BarrierStack(AlGaN(0.25), 10e-9).threshold_voltage_V()
        vth_thick = BarrierStack(AlGaN(0.25), 30e-9).threshold_voltage_V()
        self.assertLess(vth_thick, vth_thin)

    def test_deterministic(self):
        a = self.stack.sheet_density_m2(1.0)
        b = self.stack.sheet_density_m2(1.0)
        self.assertEqual(a, b)


if __name__ == "__main__":
    unittest.main()
